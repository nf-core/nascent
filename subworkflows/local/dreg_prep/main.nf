include { BEDTOOLS_BAMTOBED } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { GAWK } from '../../../modules/nf-core/gawk/main'

workflow DREG_PREP {
    take:
    bam
    bai
    sizes

    main:
    BEDTOOLS_BAMTOBED (
        bam
    )

    GAWK (
        BEDTOOLS_BAMTOBED.out.bed,
        file("${projectDir}/assets/dreg.awk")
    )

}
