process TFSEE_ANALYSIS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tfsee:1.0.0--py_0' :
        'biocontainers/tfsee:1.0.0--py_0' }"

    input:
    tuple val(meta), path(forward_bigwig), path(reverse_bigwig)
    path enhancer_regions
    path motif_database
    path expression_data
    path tf_chip_peaks
    val tfsee_config

    output:
    tuple val(meta), path("*_tfsee_results.csv")      , emit: results
    tuple val(meta), path("*_motif_enrichment.csv")   , emit: motif_enrichment
    tuple val(meta), path("*_tf_enhancer_scores.csv") , emit: tf_enhancer_scores
    tuple val(meta), path("*_clustering_results.csv") , emit: clustering_results
    tuple val(meta), path("*_feature_matrix.csv")     , emit: features
    tuple val(meta), path("*_statistics.csv")         , emit: statistics
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config_args = tfsee_config ? "--config ${tfsee_config}" : ""
    def expression_args = expression_data ? "--expression-data ${expression_data}" : ""
    def chip_args = tf_chip_peaks ? "--chip-peaks ${tf_chip_peaks}" : ""
    
    """
    # Create TFSee configuration
    cat > tfsee_config.json << EOF
    {
        "analysis_name": "${prefix}",
        "window_size": ${task.ext.window_size ?: 2000},
        "min_peak_height": ${task.ext.min_peak_height ?: 0.1},
        "smoothing_sigma": ${task.ext.smoothing_sigma ?: 2.0},
        "n_clusters": ${task.ext.n_clusters ?: 10},
        "chunk_size": ${task.ext.chunk_size ?: 1000},
        "n_workers": ${task.cpus},
        "memory_limit_mb": ${task.memory ? task.memory.toMega() * 0.8 : 8000},
        "motif_score_threshold": ${task.ext.motif_threshold ?: 0.7},
        "distance_threshold": ${task.ext.distance_threshold ?: 1000000},
        "use_sparse": ${task.ext.use_sparse ?: true},
        "enable_clustering": ${task.ext.enable_clustering ?: true},
        "calculate_statistics": ${task.ext.calculate_statistics ?: true}
    }
    EOF

    # Step 1: Extract GRO-seq features from enhancer regions
    echo "Extracting GRO-seq features..."
    tfsee_groseq_features.py \\
        --enhancers ${enhancer_regions} \\
        --forward-bigwig ${forward_bigwig} \\
        --reverse-bigwig ${reverse_bigwig} \\
        --output ${prefix}_groseq_features.csv \\
        --window-size \$(jq -r '.window_size' tfsee_config.json) \\
        --min-peak-height \$(jq -r '.min_peak_height' tfsee_config.json) \\
        --smoothing-sigma \$(jq -r '.smoothing_sigma' tfsee_config.json) \\
        --calculate-activity \\
        ${args}

    # Step 2: Perform motif analysis and enrichment
    echo "Performing motif analysis..."
    tfsee_motif_analysis.py \\
        --foreground ${enhancer_regions} \\
        --background ${enhancer_regions} \\
        --motifs ${motif_database} \\
        --output ${prefix}_motif_enrichment.csv \\
        --min-score \$(jq -r '.motif_score_threshold' tfsee_config.json) \\
        ${args}

    # Step 3: Extract enhancer sequences for motif scoring
    bedtools getfasta -fi \${GENOME_FASTA:-/dev/null} -bed ${enhancer_regions} -fo ${prefix}_enhancer_sequences.fa || \\
        echo "Warning: GENOME_FASTA not set, skipping sequence extraction"

    # Step 4: Score TF-enhancer associations
    echo "Scoring TF-enhancer associations..."
    tfsee_tf_enhancer_scoring.py \\
        --tf-regions ${enhancer_regions} \\
        --enhancer-regions ${enhancer_regions} \\
        --output ${prefix}_tf_enhancer_scores.csv \\
        --max-distance \$(jq -r '.distance_threshold' tfsee_config.json) \\
        --top-k 10000 \\
        ${expression_args} \\
        ${chip_args} \\
        ${args}

    # Step 5: Normalize features using Z-score normalization
    echo "Normalizing features..."
    tfsee_zscore_normalization.py \\
        --input ${prefix}_groseq_features.csv \\
        --data-type expression \\
        --output ${prefix}_normalized_features.csv \\
        --log-transform \\
        ${args}

    # Step 6: Perform multi-view clustering (if enabled)
    if \$(jq -r '.enable_clustering' tfsee_config.json); then
        echo "Performing multi-view clustering..."
        tfsee_multiview_clustering.py \\
            --view-data ${prefix}_groseq_features.csv ${prefix}_tf_enhancer_scores.csv \\
            --view-names groseq_features tf_enhancer_scores \\
            --n-clusters \$(jq -r '.n_clusters' tfsee_config.json) \\
            --output-prefix ${prefix}_clustering \\
            --method multiview \\
            --n-runs 5 \\
            ${args}
        
        # Copy clustering results to expected output name
        cp ${prefix}_clustering_clustering_results.csv ${prefix}_clustering_results.csv || \\
            touch ${prefix}_clustering_results.csv
    else
        echo "Clustering disabled, creating empty results file"
        echo "region_id,cluster_label" > ${prefix}_clustering_results.csv
    fi

    # Step 7: Statistical analysis (if enabled)
    if \$(jq -r '.calculate_statistics' tfsee_config.json); then
        echo "Performing statistical analysis..."
        tfsee_statistics.py \\
            --data ${prefix}_normalized_features.csv \\
            --output ${prefix}_statistics.csv \\
            --alpha 0.05 \\
            --correction-method fdr_bh \\
            --n-permutations 1000 \\
            ${args}
    else
        echo "Statistical analysis disabled, creating empty results file"
        echo "test_name,method,statistic,p_value,effect_size" > ${prefix}_statistics.csv
    fi

    # Step 8: Combine results into final output
    echo "Combining results..."
    python3 << 'EOF'
import pandas as pd
import json
import numpy as np

# Load configuration
with open('tfsee_config.json', 'r') as f:
    config = json.load(f)

prefix = config['analysis_name']

try:
    # Load all results
    groseq_features = pd.read_csv(f'{prefix}_groseq_features.csv')
    tf_enhancer_scores = pd.read_csv(f'{prefix}_tf_enhancer_scores.csv')
    motif_enrichment = pd.read_csv(f'{prefix}_motif_enrichment.csv')
    
    # Try to load clustering results
    try:
        clustering_results = pd.read_csv(f'{prefix}_clustering_results.csv')
    except:
        clustering_results = pd.DataFrame()
    
    # Try to load statistics
    try:
        statistics = pd.read_csv(f'{prefix}_statistics.csv')
    except:
        statistics = pd.DataFrame()
    
    # Create comprehensive results summary
    summary_data = {
        'analysis_id': [prefix],
        'n_enhancers': [len(groseq_features)],
        'n_tf_enhancer_associations': [len(tf_enhancer_scores)],
        'n_significant_motifs': [len(motif_enrichment[motif_enrichment['p_value_corrected'] < 0.05]) if 'p_value_corrected' in motif_enrichment.columns else 0],
        'mean_enhancer_activity': [groseq_features['activity_score'].mean() if 'activity_score' in groseq_features.columns else np.nan],
        'n_clusters': [len(clustering_results['cluster_label'].unique()) if len(clustering_results) > 0 else 0],
        'n_statistical_tests': [len(statistics) if len(statistics) > 0 else 0]
    }
    
    # Add top associations
    if len(tf_enhancer_scores) > 0:
        top_association = tf_enhancer_scores.iloc[0]
        summary_data['top_tf'] = [top_association.get('tf_name', 'N/A')]
        summary_data['top_enhancer'] = [top_association.get('enhancer_id', 'N/A')]
        summary_data['top_association_score'] = [top_association.get('association_score', np.nan)]
    else:
        summary_data['top_tf'] = ['N/A']
        summary_data['top_enhancer'] = ['N/A']
        summary_data['top_association_score'] = [np.nan]
    
    # Add configuration parameters
    for key, value in config.items():
        if key not in summary_data:
            summary_data[f'config_{key}'] = [value]
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f'{prefix}_tfsee_results.csv', index=False)
    
    print(f"TFSee analysis completed successfully for {prefix}")
    print(f"Processed {len(groseq_features)} enhancers")
    print(f"Generated {len(tf_enhancer_scores)} TF-enhancer associations")
    
except Exception as e:
    print(f"Error in results combination: {e}")
    # Create minimal results file
    minimal_results = pd.DataFrame({
        'analysis_id': [prefix],
        'status': ['error'],
        'error_message': [str(e)]
    })
    minimal_results.to_csv(f'{prefix}_tfsee_results.csv', index=False)

EOF

    # Create feature matrix combining all features
    echo "Creating feature matrix..."
    python3 << 'EOF'
import pandas as pd
import numpy as np

prefix = "${prefix}"

try:
    # Load GRO-seq features
    groseq_features = pd.read_csv(f'{prefix}_groseq_features.csv', index_col='region_id')
    
    # Initialize feature matrix with GRO-seq features
    feature_matrix = groseq_features.copy()
    
    # Add TF-enhancer association scores as features
    try:
        tf_scores = pd.read_csv(f'{prefix}_tf_enhancer_scores.csv')
        if len(tf_scores) > 0:
            # Pivot to get TF scores as columns for each enhancer
            tf_pivot = tf_scores.pivot_table(
                index='enhancer_id', 
                columns='tf_name', 
                values='association_score',
                fill_value=0
            )
            
            # Add TF_ prefix to column names
            tf_pivot.columns = [f'TF_{col}' for col in tf_pivot.columns]
            
            # Merge with feature matrix
            feature_matrix = feature_matrix.join(tf_pivot, how='left')
            feature_matrix = feature_matrix.fillna(0)
    except Exception as e:
        print(f"Warning: Could not add TF association features: {e}")
    
    # Add clustering labels if available
    try:
        clustering = pd.read_csv(f'{prefix}_clustering_results.csv')
        if len(clustering) > 0 and 'region_id' in clustering.columns:
            clustering_indexed = clustering.set_index('region_id')
            feature_matrix = feature_matrix.join(clustering_indexed, how='left')
    except Exception as e:
        print(f"Warning: Could not add clustering features: {e}")
    
    # Save feature matrix
    feature_matrix.to_csv(f'{prefix}_feature_matrix.csv')
    print(f"Feature matrix created with shape: {feature_matrix.shape}")
    
except Exception as e:
    print(f"Error creating feature matrix: {e}")
    # Create minimal feature matrix
    pd.DataFrame({'error': [str(e)]}).to_csv(f'{prefix}_feature_matrix.csv')

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
        tfsee: "1.0.0"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tfsee_results.csv
    touch ${prefix}_motif_enrichment.csv
    touch ${prefix}_tf_enhancer_scores.csv
    touch ${prefix}_clustering_results.csv
    touch ${prefix}_feature_matrix.csv
    touch ${prefix}_statistics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g' || echo "3.8.0")
        tfsee: "1.0.0"
    END_VERSIONS
    """
}