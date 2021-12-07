//
// Uncompress and prepare reference genome files
//

include { GTF2BED                     } from '../../modules/local/gtf2bed'
include { CAT_ADDITIONAL_FASTA        } from '../../modules/local/cat_additional_fasta'
include { GTF_GENE_FILTER             } from '../../modules/local/gtf_gene_filter'

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/modules/gunzip/main'
include { UNTAR as UNTAR_BWA_INDEX    } from '../../modules/nf-core/modules/untar/main'
include { GFFREAD                     } from '../../modules/nf-core/modules/gffread/main'
include { BWA_INDEX                   } from '../../modules/nf-core/modules/bwa/index/main'
include { BWAMEM2_INDEX               } from '../../modules/nf-core/modules/bwamem2/index/main'
include { CUSTOM_GETCHROMSIZES        } from '../../modules/nf-core/modules/custom/getchromsizes/main'

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( params.gene_bed ).gunzip
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf )
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES ( ch_fasta ).sizes

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index   = Channel.empty()
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index   = BWA_INDEX ( ch_fasta ).index
            ch_versions    = ch_versions.mix(BWA_INDEX.out.versions)
        }
    } else if ('bwamem2' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index = BWAMEM2_INDEX ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        }
    }

    emit:
    fasta       = ch_fasta       // path: genome.fasta
    gtf         = ch_gtf         // path: genome.gtf
    gene_bed    = ch_gene_bed    // path: gene.bed
    chrom_sizes = ch_chrom_sizes // path: genome.sizes
    bwa_index   = ch_bwa_index   // path: star/index/

    versions    = ch_versions    // channel: [ versions.yml ]
}
