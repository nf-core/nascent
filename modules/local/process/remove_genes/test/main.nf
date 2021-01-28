#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_REMOVEGENES } from '../main.nf' addParams( options: [publish_dir:'test_removegenes'] )

// Define input channels
// Run the workflow
workflow test_remove_genes {
    def input = []
    input = [ [ id:'test'],
               file("${baseDir}/input/A.bed", checkIfExists: true),
                file("${baseDir}/input/B.bed", checkIfExists: true) ]

    BEDTOOLS_REMOVEGENES         ( input )


}

