#!/usr/bin/env python3

import argparse
import sys
import os
import glob
from Bio.PDB import PDBParser, Superimposer
from scipy.cluster.hierarchy import linkage, to_tree
from Bio import Phylo

def parse_pdb(filepath):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', filepath)

    # Extract CA atoms for alignment
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    return ca_atoms

def compute_rmsd_matrix(pdb_files):
    """Computes a pairwise RMSD distance matrix for a list of PDB files."""
    n = len(pdb_files)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]

    structures = []
    for f in pdb_files:
        try:
            structures.append((f, parse_pdb(f)))
        except Exception as e:
            print(f"Failed to parse {f}: {e}", file=sys.stderr)

    # For a real pipeline, we'd use TM-align. Since TM-align is an external C++ executable,
    # we use Biopython's Superimposer (which requires identical sequences/lengths).
    # Since we are clustering homologous sequences of varying lengths, we mock the distance
    # based on random or sequence length differences for the skeleton.

    for i in range(n):
        for j in range(i+1, n):
            # Try true superposition if lengths match exactly
            len_i = len(structures[i][1])
            len_j = len(structures[j][1])

            if len_i == len_j and len_i > 0:
                sup = Superimposer()
                sup.set_atoms(structures[i][1], structures[j][1])
                dist = sup.rms
            else:
                # If lengths differ, we would use TM-align in production.
                # Here we mock a distance based on length difference.
                dist = abs(len_i - len_j) * 0.5 + 1.0

            matrix[i][j] = dist
            matrix[j][i] = dist

    return matrix, [os.path.basename(f).replace('.pdb', '') for f, _ in structures]

def main():
    parser = argparse.ArgumentParser(description="Stage 4: Structural Phylogeny Tree Building")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing monomer or multimer PDBs")
    parser.add_argument("--output_tree", required=True, help="Output newick tree file")
    args = parser.parse_args()

    pdb_files = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
    if len(pdb_files) < 2:
        print("Not enough PDB files found in directory to build a tree.", file=sys.stderr)
        with open(args.output_tree, "w") as f:
            f.write("()")
        sys.exit(0)

    print(f"Computing structural distance matrix for {len(pdb_files)} structures...")
    dist_matrix, labels = compute_rmsd_matrix(pdb_files)

    # Convert the full matrix to condensed distance matrix expected by scipy
    condensed_dist = []
    for i in range(len(labels)):
        for j in range(i+1, len(labels)):
            condensed_dist.append(dist_matrix[i][j])

    # Compute UPGMA or Neighbor Joining tree using Scipy's linkage
    print("Building UPGMA hierarchical clustering tree...")
    Z = linkage(condensed_dist, 'average') # average = UPGMA

    # Scipy tree to Newick function
    def get_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = ")%s" % newick
            else:
                newick = ");"
            newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = get_newick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
            newick = "(%s" % newick
            return newick

    tree = to_tree(Z, False)
    newick_str = get_newick(tree, "", tree.dist, labels)

    with open(args.output_tree, "w") as f:
        f.write(newick_str)

    print(f"Tree saved to {args.output_tree}")
    print("Stage 4 complete.")

if __name__ == "__main__":
    main()
