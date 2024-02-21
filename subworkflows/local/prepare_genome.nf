//
// Uncompress and prepare reference genome files
//

include { GTF2BED } from '../../modules/local/gtf2bed'

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED } from '../../modules/nf-core/gunzip/main'
include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_DRAGMAP } from '../../modules/nf-core/untar/main'
include { GFFREAD } from '../../modules/nf-core/gffread/main'
include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE } from '../../modules/nf-core/dragmap/hashtable/main'
include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main'

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices
    fasta
    gtf
    gff

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff))
            }
            ch_gtf      = GFFREAD ( ch_gff ).gtf
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    ch_dragmap = Channel.empty()
    ch_bowtie2_index = Channel.empty()
    // TODO Turn this into a switch
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
                ch_versions = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                // TODO Give the meta from basename or genome?
                ch_bwa_index = [ [meta: "Genome"], file(params.bwa_index) ]
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    } else if ('bwamem2' in prepare_tool_indices) {
        if (params.bwamem2_index) {
            if (params.bwamem2_index.endsWith('.tar.gz') || params.bwamem2_index.endsWith('.tgz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwamem2_index ] ).untar
                ch_versions = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                // TODO Give the meta from basename or genome?
                ch_bwa_index = [ [meta: "Genome"], file(params.bwamem2_index) ]
            }
        } else {
            ch_bwa_index = BWAMEM2_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }
    } else if ('dragmap' in prepare_tool_indices) {
        if (params.dragmap) {
            if (params.dragmap.endsWith('.tar.gz')) {
                ch_dragmap = UNTAR_DRAGMAP_INDEX ( params.dragmap ).untar
                ch_versions = ch_versions.mix(UNTAR_DRAGMAP_INDEX.out.versions)
            } else {
                // TODO Give the meta from basename or genome?
                ch_dragmap = [ [meta: "Genome"], file(params.dragmap) ]
            }
        } else {
            ch_dragmap = DRAGMAP_HASHTABLE( ch_fasta.map { [ [:], it ] } ).hashmap
            ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
        }
    } else if ('bowtie2' in prepare_tool_indices) {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.bowtie2_index ).untar
                ch_versions = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                // TODO Give the meta from basename or genome?
                ch_bowtie2_index = [ [meta: "Genome"], file(params.bowtie2_index) ]
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    emit:
    fasta = ch_fasta
    fai = ch_fai
    gtf = ch_gtf
    gene_bed = ch_gene_bed
    chrom_sizes = ch_chrom_sizes
    bwa_index = ch_bwa_index
    dragmap = ch_dragmap
    bowtie2_index = ch_bowtie2_index

    versions = ch_versions.ifEmpty(null)
}
