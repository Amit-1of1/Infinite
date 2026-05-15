#!/usr/bin/env python3

import argparse
import sys
import time
import requests
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Sequence & Ortholog Search using NCBI BLAST API")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--email", required=True, help="Email address for NCBI Entrez")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file containing homologs")
    parser.add_argument("--e_value", type=float, default=0.001, help="E-value threshold for BLAST")
    parser.add_argument("--max_hits", type=int, default=10, help="Maximum number of homolog hits to retrieve")
    args = parser.parse_args()

    Entrez.email = args.email

    try:
        record = SeqIO.read(args.input, format="fasta")
    except Exception as e:
        print(f"Error reading input FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Running BLAST search for {record.id} via NCBI WWW (this may take a few minutes)...")

    # In a real scenario, running NCBIWWW.qblast takes a while. We wrap it for safety.
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=args.max_hits, expect=args.e_value)
    except Exception as e:
        print(f"BLAST search failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("Parsing BLAST results...")
    blast_record = NCBIXML.read(result_handle)

    homolog_ids = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < args.e_value:
                # The title typically contains the accession/ID
                accession = alignment.accession
                homolog_ids.append(accession)
                break # Only need one per alignment

    if not homolog_ids:
        print("No homologs found within the E-value threshold.")
        # Write the original sequence as the only output just to continue the pipeline
        with open(args.output_fasta, "w") as f_out:
            SeqIO.write(record, f_out, "fasta")
        sys.exit(0)

    print(f"Found {len(homolog_ids)} homologs. Fetching sequences from NCBI...")

    # Deduplicate
    homolog_ids = list(set(homolog_ids))

    fetched_records = []
    # Always include the original sequence
    fetched_records.append(record)

    for homolog_id in homolog_ids:
        try:
            print(f"Fetching {homolog_id}...")
            # Use Entrez efetch to get the sequence
            handle = Entrez.efetch(db="protein", id=homolog_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            fetched_records.append(seq_record)
            handle.close()
            time.sleep(0.34) # Respect NCBI API limits (3 requests per second)
        except Exception as e:
            print(f"Warning: Failed to fetch sequence for {homolog_id}: {e}")

    print(f"Writing {len(fetched_records)} sequences to {args.output_fasta}...")
    with open(args.output_fasta, "w") as f_out:
        SeqIO.write(fetched_records, f_out, "fasta")

    print("Stage 1 complete.")

if __name__ == "__main__":
    main()
