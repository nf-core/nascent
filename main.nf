#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/groseq
========================================================================================
    Github : https://github.com/nf-core/groseq
    Website: https://nf-co.re/groseq
    Slack  : https://nfcore.slack.com/channels/groseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf       = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff       = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed  = WorkflowMain.getGenomeAttribute(params, 'bed12')
// FIXME Issue with $INDEX in bwa and igenomes hosting multiple versions
// So we have to generate our own index for now
// params.bwa_index = Checks.get_genome_attribute(params, 'bwa')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { GROSEQ } from './workflows/groseq'

//
// WORKFLOW: Run main nf-core/groseq analysis pipeline
//
workflow NFCORE_GROSEQ {
    GROSEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_GROSEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
