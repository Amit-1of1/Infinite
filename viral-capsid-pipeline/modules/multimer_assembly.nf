process MULTIMER_ASSEMBLY {
    publishDir "${params.outdir}/03_multimer_assembly", mode: 'copy'

    input:
    // In a real scenario we might take the FASTA instead of monomer PDBs for AlphaFold multimer,
    // but the input channel comes from the previous steps. Let's just pass the original sequences from stage 1
    // For this skeleton, we pass the directory just as a trigger, but we actually need the fasta.
    // We will update the main.nf to pass the fasta to multimer assembly, not the monomer pdbs.
    path input_fasta

    output:
    path "multimer_pdbs", emit: multimer_pdbs

    script:
    def gpu_flag = params.use_gpu ? "--use_gpu" : ""
    """
    03_multimer_assembly.py \\
        --input_fasta ${input_fasta} \\
        --outdir multimer_pdbs \\
        --copies 3 \\
        ${gpu_flag}
    """
}
