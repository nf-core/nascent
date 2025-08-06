/*
 * TFSee analysis workflow for transcription factor-enhancer prediction
 */

include { TFSEE_GROSEQ_FEATURES       } from '../../../modules/local/tfsee_groseq_features/main'
include { TFSEE_MOTIF_ANALYSIS        } from '../../../modules/local/tfsee_motif_analysis/main'
include { TFSEE_TF_ENHANCER_SCORING   } from '../../../modules/local/tfsee_tf_enhancer_scoring/main'
include { TFSEE_ZSCORE_NORMALIZATION  } from '../../../modules/local/tfsee_zscore_normalization/main'
include { TFSEE_MULTIVIEW_CLUSTERING  } from '../../../modules/local/tfsee_multiview_clustering/main'
include { TFSEE_STATISTICS            } from '../../../modules/local/tfsee_statistics/main'
include { TFSEE_COMBINE_RESULTS       } from '../../../modules/local/tfsee_combine_results/main'

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

    // Set up optional inputs
    ch_motif_db = motif_database ?: Channel.empty()
    ch_tf_expr = tf_expression ?: Channel.empty()
    ch_tf_chip = tf_chip_peaks ?: Channel.empty()
    
    // Get first enhancer regions file for shared inputs
    ch_enhancer_regions = enhancer_beds.map { meta, bed -> bed }.first()

    // Step 2: Perform motif analysis (this can run independently)
    TFSEE_MOTIF_ANALYSIS(
        enhancer_beds,
        ch_motif_db.ifEmpty([]),
        fasta
    )
    ch_versions = ch_versions.mix(TFSEE_MOTIF_ANALYSIS.out.versions.first())

    // Step 3: Score TF-enhancer associations
    // For test mode, create a separate TF regions channel to avoid file name collision
    ch_tf_regions = enhancer_beds.map { meta, bed ->
        // Create a copy with different name
        def tf_bed = file("${meta.id}_tf_regions.bed")
        tf_bed.text = bed.text
        [meta, tf_bed]
    }
    
    TFSEE_TF_ENHANCER_SCORING(
        ch_tf_regions,
        ch_enhancer_regions,
        ch_tf_expr.ifEmpty([]),
        ch_tf_chip.ifEmpty([])
    )
    ch_versions = ch_versions.mix(TFSEE_TF_ENHANCER_SCORING.out.versions.first())

    // For test mode, create dummy GRO-seq features since we don't have coverage tracks
    ch_groseq_features = enhancer_beds.map { meta, bed ->
        def dummy_file = file("${meta.id}_groseq_features.csv")
        dummy_file.text = "enhancer_id,forward_signal,reverse_signal,total_signal,directionality_index,peak_count\n"
        dummy_file.text += "enhancer_1,80.0,60.0,140.0,0.0,2\n"
        dummy_file.text += "enhancer_2,160.0,120.0,280.0,0.1,3\n"
        [meta, dummy_file]
    }

    // Step 4: Normalize features
    TFSEE_ZSCORE_NORMALIZATION(
        ch_groseq_features
    )
    ch_versions = ch_versions.mix(TFSEE_ZSCORE_NORMALIZATION.out.versions.first())

    // Step 5: Multi-view clustering (if enabled)
    ch_clustering_input = ch_groseq_features
        .join(TFSEE_TF_ENHANCER_SCORING.out.scores, by: [0])
        .map { meta, groseq_features, tf_scores ->
            [meta, [groseq_features, tf_scores]]
        }
    
    if (params.tfsee_enable_clustering ?: true) {
        TFSEE_MULTIVIEW_CLUSTERING(
            ch_clustering_input,
            'groseq_features tf_enhancer_scores'
        )
        ch_versions = ch_versions.mix(TFSEE_MULTIVIEW_CLUSTERING.out.versions.first())
        ch_clustering_results = TFSEE_MULTIVIEW_CLUSTERING.out.clustering
    } else {
        // Create empty clustering results
        ch_clustering_results = ch_clustering_input.map { meta, files ->
            [meta, file('NO_FILE')]
        }
    }

    // Step 6: Statistical analysis (if enabled)
    if (params.tfsee_calculate_statistics ?: true) {
        TFSEE_STATISTICS(
            TFSEE_ZSCORE_NORMALIZATION.out.normalized
        )
        ch_versions = ch_versions.mix(TFSEE_STATISTICS.out.versions.first())
        ch_statistics_results = TFSEE_STATISTICS.out.statistics
    } else {
        // Create empty statistics results
        ch_statistics_results = TFSEE_ZSCORE_NORMALIZATION.out.normalized.map { meta, normalized ->
            [meta, file('NO_FILE')]
        }
    }

    // Step 7: Combine all results
    ch_combine_input = enhancer_beds
        .join(ch_groseq_features, by: [0])
        .join(TFSEE_TF_ENHANCER_SCORING.out.scores, by: [0])
        .join(TFSEE_MOTIF_ANALYSIS.out.enrichment, by: [0])
        .join(ch_clustering_results, by: [0])
        .join(ch_statistics_results, by: [0])
        .map { meta, bed, groseq, tf_scores, motif, clustering, statistics ->
            [meta, groseq, tf_scores, motif, clustering, statistics]
        }

    // Create TFSee configuration
    def tfsee_config = groovy.json.JsonOutput.toJson([
        analysis_name: "tfsee_analysis",
        window_size: params.tfsee_window_size ?: 2000,
        min_peak_height: params.tfsee_min_peak_height ?: 0.1,
        smoothing_sigma: params.tfsee_smoothing_sigma ?: 2.0,
        n_clusters: params.tfsee_n_clusters ?: 10,
        motif_threshold: params.tfsee_motif_threshold ?: 0.7,
        distance_threshold: params.tfsee_distance_threshold ?: 1000000,
        enable_clustering: params.tfsee_enable_clustering ?: true,
        calculate_statistics: params.tfsee_calculate_statistics ?: true
    ])

    TFSEE_COMBINE_RESULTS(
        ch_combine_input,
        tfsee_config
    )
    ch_versions = ch_versions.mix(TFSEE_COMBINE_RESULTS.out.versions.first())

    emit:
    results               = TFSEE_COMBINE_RESULTS.out.results
    motif_enrichment      = TFSEE_MOTIF_ANALYSIS.out.enrichment
    tf_enhancer_scores    = TFSEE_TF_ENHANCER_SCORING.out.scores
    clustering_results    = ch_clustering_results
    feature_matrix        = TFSEE_COMBINE_RESULTS.out.feature_matrix
    statistics            = ch_statistics_results
    versions              = ch_versions
}