include { BEDTOOLS_BAMTOBED } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { GAWK } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SORT } from '../../../modules/nf-core/bedtools/sort/main'
include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PLUS
} from '../../../modules/nf-core/bedtools/genomecov/main'

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

    BEDTOOLS_SORT (
        GAWK.out.output,
        sizes
    )

    BEDTOOLS_GENOMECOV_PLUS (
        BEDTOOLS_SORT.out.sorted.combine(Channel.from(0)),
        sizes,
        'bedGraph'
    )

}
