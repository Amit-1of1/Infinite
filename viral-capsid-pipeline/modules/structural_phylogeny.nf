process STRUCTURAL_PHYLOGENY {
    publishDir "${params.outdir}/04_structural_phylogeny", mode: 'copy'

    input:
    path pdb_dir

    output:
    path "structural_tree.nwk", emit: tree

    script:
    """
    04_structural_phylogeny.py \\
        --pdb_dir ${pdb_dir} \\
        --output_tree structural_tree.nwk
    """
}
