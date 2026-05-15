nextflow.enable.dsl=2

/*
 * Viral Capsid / Large Complex Evolution & Assembly Analysis Pipeline
 *
 * This pipeline combines phylogenetics with 3D structural analysis for multimeric complexes.
 */

// Define input parameters
params.input_fasta = "data/input.fasta"
params.email = "example@example.com" // Required for NCBI Entrez
params.outdir = "results"
params.threads = 4

// Include modules
include { SEQUENCE_SEARCH } from './modules/sequence_search.nf'
include { STRUCTURE_PREDICTION } from './modules/structure_prediction.nf'
include { MULTIMER_ASSEMBLY } from './modules/multimer_assembly.nf'
include { STRUCTURAL_PHYLOGENY } from './modules/structural_phylogeny.nf'
include { INTERFACE_ANALYSIS } from './modules/interface_analysis.nf'

workflow {
    // Stage 1: Sequence & Ortholog Search
    fasta_ch = Channel.fromPath(params.input_fasta, checkIfExists: true)
    SEQUENCE_SEARCH(fasta_ch, params.email)

    // Stage 2: Structure Prediction / Homology Modeling (Monomers)
    STRUCTURE_PREDICTION(SEQUENCE_SEARCH.out.homologs_fasta)

    // Stage 3: Multimer Assembly Prediction (needs FASTA sequence, not monomer PDBs)
    MULTIMER_ASSEMBLY(SEQUENCE_SEARCH.out.homologs_fasta)

    // Stage 4: Structural Phylogeny (Using monomer PDBs)
    STRUCTURAL_PHYLOGENY(STRUCTURE_PREDICTION.out.monomer_pdbs)

    // Stage 5: Interface Evolution Analysis (Using Multimer PDBs and MSA)
    INTERFACE_ANALYSIS(MULTIMER_ASSEMBLY.out.multimer_pdbs, SEQUENCE_SEARCH.out.msa)
}
