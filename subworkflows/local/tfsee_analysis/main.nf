/*
 * TFSee analysis workflow for transcription factor-enhancer prediction
 */

include { TFSEE_ANALYSIS                      } from '../../../modules/local/tfsee/main'
include { BEDTOOLS_GETFASTA                   } from '../../../modules/nf-core/bedtools/getfasta/main'
include { SAMTOOLS_FAIDX                      } from '../../../modules/nf-core/samtools/faidx/main'

workflow TFSEE_ANALYSIS_WORKFLOW {
    take:
    enhancer_beds      // channel: [ meta, bed ]
    coverage_tracks    // channel: [ meta, [forward.bw, reverse.bw] ]
    fasta              // channel: [ fasta ]
    motif_database     // channel: [ path ]
    tf_expression      // channel: [ path ] (optional)
    tf_chip_peaks      // channel: [ path ] (optional)

    main:
    ch_versions = Channel.empty()

    // Prepare genome fasta index if needed
    ch_fasta_fai = fasta.map { fasta_file -> 
        def fai_file = file("${fasta_file}.fai")
        if (fai_file.exists()) {
            [[:], fasta_file, fai_file]
        } else {
            [[:], fasta_file, []]
        }
    }

    // Index fasta if not already indexed
    ch_fasta_to_index = ch_fasta_fai.filter { meta, fasta_file, fai_file -> 
        fai_file.size() == 0 
    }.map { meta, fasta_file, fai_file -> 
        [meta, fasta_file] 
    }

    if (!ch_fasta_to_index.isEmpty()) {
        SAMTOOLS_FAIDX(ch_fasta_to_index, [[], []])
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        
        ch_fasta_indexed = SAMTOOLS_FAIDX.out.fai.map { meta, fasta_file, fai_file ->
            [meta, fasta_file, fai_file]
        }
    } else {
        ch_fasta_indexed = ch_fasta_fai.filter { meta, fasta_file, fai_file -> 
            fai_file.size() > 0 
        }
    }

    // Combine enhancer beds with coverage tracks
    ch_enhancers_coverage = enhancer_beds.join(coverage_tracks, by: [0])

    // Prepare inputs for TFSee analysis
    ch_tfsee_input = ch_enhancers_coverage.map { meta, bed, coverage ->
        def forward_bw = coverage.find { it.name.contains('_pl.bw') || it.name.contains('_forward.bw') || it.name.contains('_plus.bw') }
        def reverse_bw = coverage.find { it.name.contains('_mn.bw') || it.name.contains('_reverse.bw') || it.name.contains('_minus.bw') }
        
        if (!forward_bw || !reverse_bw) {
            // If we can't identify forward/reverse, assume first is forward, second is reverse
            forward_bw = coverage[0]
            reverse_bw = coverage.size() > 1 ? coverage[1] : coverage[0]
        }
        
        [meta, forward_bw, reverse_bw, bed]
    }

    // Set up optional inputs
    ch_motif_db = motif_database ?: Channel.empty()
    ch_tf_expr = tf_expression ?: Channel.empty()
    ch_tf_chip = tf_chip_peaks ?: Channel.empty()
    
    // Default TFSee configuration
    def tfsee_config = [
        window_size: params.tfsee_window_size ?: 2000,
        min_peak_height: params.tfsee_min_peak_height ?: 0.1,
        smoothing_sigma: params.tfsee_smoothing_sigma ?: 2.0,
        n_clusters: params.tfsee_n_clusters ?: 10,
        motif_threshold: params.tfsee_motif_threshold ?: 0.7,
        distance_threshold: params.tfsee_distance_threshold ?: 1000000,
        enable_clustering: params.tfsee_enable_clustering ?: true,
        calculate_statistics: params.tfsee_calculate_statistics ?: true
    ]

    // Run TFSee analysis
    TFSEE_ANALYSIS(
        ch_tfsee_input.map { meta, forward_bw, reverse_bw, bed -> 
            [meta, forward_bw, reverse_bw] 
        },
        ch_tfsee_input.map { meta, forward_bw, reverse_bw, bed -> bed }.first(),
        ch_motif_db.ifEmpty([]),
        ch_tf_expr.ifEmpty([]),
        ch_tf_chip.ifEmpty([]),
        tfsee_config
    )
    ch_versions = ch_versions.mix(TFSEE_ANALYSIS.out.versions.first())

    emit:
    results                = TFSEE_ANALYSIS.out.results
    motif_enrichment      = TFSEE_ANALYSIS.out.motif_enrichment
    tf_enhancer_scores    = TFSEE_ANALYSIS.out.tf_enhancer_scores
    clustering_results    = TFSEE_ANALYSIS.out.clustering_results
    feature_matrix        = TFSEE_ANALYSIS.out.features
    statistics            = TFSEE_ANALYSIS.out.statistics
    versions              = ch_versions
}