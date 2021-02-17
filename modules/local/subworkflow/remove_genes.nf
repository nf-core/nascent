<<<<<<< HEAD
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_REMOVEGENES } from '../process/remove_genes/main.nf' addParams( options: [publish_dir:'test_removegenes'] )

// Define input channels
// Run the workflow
workflow test_remove_genes {
    def input = []
    input = [ [ id:'test'],
               file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/A.bed", checkIfExists: true),
                file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/B.bed", checkIfExists: true) ]

    BEDTOOLS_REMOVEGENES         ( input )


=======
include { BEDTOOLS_SORT } from '../../nf-core/software/bedtools/sort/main' addParams( options: [:] )
include { BEDTOOLS_SLOP } from '../../nf-core/software/bedtools/slop/main' addParams( options: [:] )
include { BEDTOOLS_INTERSECT } from '../../nf-core/software/bedtools/intersect/main' addParams( options: [:] )

workflow REMOVE_GENES {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    refseq // file: /path/to/refseq
    chromsizes

    main:
    BEDTOOLS_SORT ( refseq )
    BEDTOOLS_SLOP ( BEDTOOLS_SORT.out.bed, chromsizes )
    BEDTOOLS_INTERSECT ( BEDTOOLS_SLOP.out.bed, reads )


    emit:
    transcrips = BEDTOOLS_INERSECT.out.bed
>>>>>>> 5e8a4e3a3fe9c7ec5dd1f74e7d49c5016eb82a95
}

