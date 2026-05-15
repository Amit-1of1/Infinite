process STRUCTURE_PREDICTION {
    publishDir "${params.outdir}/02_structure_prediction", mode: 'copy'

    input:
    path fasta_file

    output:
    path "pdbs", emit: monomer_pdbs

    script:
    """
    02_structure_prediction.py \\
        --input_fasta ${fasta_file} \\
        --outdir pdbs
    """
}
