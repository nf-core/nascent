#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_FORWARD} from '../process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm', args:'--strand=+'] )
include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_REVERSE } from '../process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm', args:'--strand=-'] )
include { GROHMM_TRANSCRIPTCALLING } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_PARAMETERTUNING } from '../process/grohmm/parametertuning/main.nf' addParams( options: [publish_dir:'test_grohmm'] )


/*
 * Test with forward strand
 */
workflow GROHMM {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://raw.githack.com/Kraus-Lab/groHMM/master/inst/extdata/S0mR1.bam', checkIfExists: true),
               file("${launchDir}/modules/local/process/grohmm/test/tune.csv", checkIfExists: true) ]]
    GROHMM_PARAMETERTUNING ( input )
}

/*
 * Test with reverse strand

workflow test_grohmm_makeucscfile_reverse {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_MAKEUCSCFILE_REVERSE ( input )
}




workflow test_grohmm_transcriptcalling {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_TRANSCRIPTCALLING ( input )
}

workflow test_grohmm_parametertuning {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_PARAMETERTUNING ( input )
}
*/
