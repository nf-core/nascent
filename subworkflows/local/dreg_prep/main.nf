include { BEDTOOLS_BAMTOBED } from '../modules/nf-core/bedtools/bamtobed/main'

workflow DREG_PREP {
    take:
    bam
    bai
    sizes

    main:
    BEDTOOLS_BAMTOBED(
        bam
    )

}
