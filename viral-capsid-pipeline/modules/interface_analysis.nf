process INTERFACE_ANALYSIS {
    publishDir "${params.outdir}/05_interface_analysis", mode: 'copy'

    input:
    path pdb_dir
    path msa_file

    output:
    path "interfaces", emit: interface_pdbs

    script:
    """
    05_interface_analysis.py \\
        --pdb_dir ${pdb_dir} \\
        --outdir interfaces
    """
}
