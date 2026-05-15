#!/usr/bin/env python3

import argparse
import sys
import os
import glob
from Bio.PDB import PDBParser, PDBIO, Select

def calculate_interface_contacts(structure, distance_cutoff=5.0):
    """
    Identifies residues at the interface of a multimer.
    Returns a set of residue identifiers (chain_id, res_id).
    """
    interface_residues = set()
    chains = list(structure[0])

    if len(chains) < 2:
        return interface_residues

    for i in range(len(chains)):
        for j in range(i+1, len(chains)):
            chain1 = chains[i]
            chain2 = chains[j]

            for res1 in chain1:
                if 'CA' not in res1: continue
                for res2 in chain2:
                    if 'CA' not in res2: continue

                    dist = res1['CA'] - res2['CA']
                    if dist < distance_cutoff:
                        interface_residues.add((chain1.id, res1.id[1]))
                        interface_residues.add((chain2.id, res2.id[1]))

    return interface_residues

class InterfaceSelect(Select):
    def __init__(self, interface_residues):
        self.interface_residues = interface_residues

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        res_id = residue.id[1]
        if (chain_id, res_id) in self.interface_residues:
            return 1
        return 0

def main():
    parser = argparse.ArgumentParser(description="Stage 5: Interface Evolution Analysis")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing multimer PDBs")
    parser.add_argument("--outdir", required=True, help="Output directory for interface PDBs")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
    parser_pdb = PDBParser(QUIET=True)
    io = PDBIO()

    if not pdb_files:
        print("No PDB files found for interface analysis.", file=sys.stderr)
        sys.exit(0)

    for f in pdb_files:
        print(f"Analyzing interfaces for {f}...")
        try:
            structure = parser_pdb.get_structure('struct', f)
            interface_res = calculate_interface_contacts(structure)

            print(f"Found {len(interface_res)} interface residues across chains.")

            # Save a PDB containing ONLY the interface residues
            base_name = os.path.basename(f)
            out_pdb = os.path.join(args.outdir, base_name.replace(".pdb", "_interface.pdb"))

            io.set_structure(structure)
            io.save(out_pdb, InterfaceSelect(interface_res))
            print(f"Saved interface structure to {out_pdb}")

        except Exception as e:
            print(f"Error processing {f}: {e}", file=sys.stderr)

    print("Stage 5 complete.")

if __name__ == "__main__":
    main()
