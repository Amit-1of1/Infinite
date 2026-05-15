#!/usr/bin/env python3

import argparse
import sys
import os
import requests
from Bio import SeqIO
import time

def fetch_esmfold_structure(sequence, output_pdb):
    """Fetches structure prediction from ESMFold API."""
    # Using ESMFold API (which is very fast for sequences < 400 aa)
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(url, data=str(sequence), headers={"Content-Type": "text/plain"})
        if response.status_code == 200:
            with open(output_pdb, "w") as f:
                f.write(response.text)
            return True
        else:
            print(f"ESMFold API Error: {response.status_code} - {response.text}", file=sys.stderr)
            return False
    except Exception as e:
        print(f"Exception during ESMFold API call: {e}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="Stage 2: Structure Prediction using ESMFold API")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file with monomer sequences")
    parser.add_argument("--outdir", required=True, help="Directory to save output PDB files")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    try:
        records = list(SeqIO.parse(args.input_fasta, "fasta"))
    except Exception as e:
        print(f"Error reading input FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(records)} sequences. Starting structure prediction...")

    successful = 0
    for record in records:
        pdb_filename = os.path.join(args.outdir, f"{record.id.replace('|', '_')}.pdb")
        print(f"Predicting structure for {record.id} (Length: {len(record.seq)})...")

        # ESMFold API might fail for very long sequences
        if len(record.seq) > 400:
            print(f"Warning: Sequence {record.id} is long ({len(record.seq)}). API might timeout or fail.")

        success = fetch_esmfold_structure(record.seq, pdb_filename)
        if success:
            print(f"Saved structure to {pdb_filename}")
            successful += 1
        else:
            print(f"Failed to predict structure for {record.id}.")

        # Add a delay to be polite to the public API
        time.sleep(1)

    print(f"Stage 2 complete. Successfully predicted {successful}/{len(records)} structures.")

if __name__ == "__main__":
    main()
