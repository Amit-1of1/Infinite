process SEQUENCE_SEARCH {
    publishDir "${params.outdir}/01_sequence_search", mode: 'copy'

    input:
    path input_fasta
    val email

    output:
    path "homologs.fasta", emit: homologs_fasta
    // A real pipeline would align them into an MSA here
    path "homologs.fasta", emit: msa

    script:
    """
    01_sequence_search.py \\
        --input ${input_fasta} \\
        --email ${email} \\
        --output_fasta homologs.fasta \\
        --e_value 0.001 \\
        --max_hits 5
    """
}
