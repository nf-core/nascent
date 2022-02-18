/*
 * TODO
 */

workflow COVERAGE_GRAPHS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    emit:

    versions = ch_versions                      // channel: [ versions.yml ]
}
