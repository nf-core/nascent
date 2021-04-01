/*
 * Uncompress and prepare reference genome files
*/

params.genome_options       = [:]
params.index_options        = [:]
params.gffread_options      = [:]
params.bwa_index_options    = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/local/gunzip'               addParams( options: params.genome_options       )
include { GTF2BED                     } from '../../modules/local/gtf2bed'              addParams( options: params.genome_options       )
include { CAT_ADDITIONAL_FASTA        } from '../../modules/local/cat_additional_fasta' addParams( options: params.genome_options       )
include { GTF_GENE_FILTER             } from '../../modules/local/gtf_gene_filter'      addParams( options: params.genome_options       )
include { GET_CHROM_SIZES             } from '../../modules/local/get_chrom_sizes'      addParams( options: params.genome_options       )
include { UNTAR as UNTAR_BWA_INDEX    } from '../../modules/local/untar'                addParams( options: params.bwa_index_options   )

include { GFFREAD                     } from '../../modules/nf-core/software/gffread/main'   addParams( options: params.gffread_options      )
include { BWA_INDEX                   } from '../../modules/nf-core/software/bwa/index/main' addParams( options: params.bwa_index_options )
include { BWAMEM2_INDEX               } from '../../modules/nf-core/software/bwamem2/index/main' addParams( options: params.bwa_index_options )

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:
    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * Uncompress GTF annotation file or create from GFF3 if required
     */
    ch_gffread_version = Channel.empty()
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
        ch_gffread_version = GFFREAD.out.version
    }

    /*
     * Uncompress additional fasta file and concatenate with reference fasta and gtf files
     */
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( params.additional_fasta ).gunzip
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta )
        ch_fasta = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf   = CAT_ADDITIONAL_FASTA.out.gtf
    }

    /*
     * Uncompress gene BED annotation file or create from GTF if required
     */
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( params.gene_bed ).gunzip
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf )
    }

    /*
     * Create chromosome sizes file
     */
    ch_chrom_sizes = GET_CHROM_SIZES ( ch_fasta ).sizes

    /*
     * Uncompress BWA index or generate from scratch if required
     */
    ch_bwa_index   = Channel.empty()
    ch_bwa_version = Channel.empty()
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index   = BWA_INDEX ( ch_fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }
    } else if ('bwamem2' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {
            ch_bwa_index   = BWAMEM2_INDEX ( ch_fasta ).index
            ch_bwa_version = BWAMEM2_INDEX.out.version
        }
    }

    emit:
    fasta            = ch_fasta            // path: genome.fasta
    gtf              = ch_gtf              // path: genome.gtf
    gene_bed         = ch_gene_bed         // path: gene.bed
    chrom_sizes      = ch_chrom_sizes      // path: genome.sizes
    bwa_index        = ch_bwa_index       // path: star/index/
    gffread_version  = ch_gffread_version  // path: *.version.txt
}
