#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GROHMM } from '../../../../subworkflows/local/grohmm' addParams( options: [ publish_dir:'test_removegenes' ], slop_options: [args: '-l 15 -r 30'] )

workflow test_grohmm {
    def input = [[ id: 'test' ],
            [ file(params.test_data['grohmm']['s0mR1'], checkIfExists: true),
            file(params.test_data['grohmm']['s40mR1'], checkIfExists: true) ]]
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    GROHMM (
        input,
        gtf
    )
}
