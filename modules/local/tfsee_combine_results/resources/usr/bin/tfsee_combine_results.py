#!/usr/bin/env python3
"""
Combine TFSee analysis results into a comprehensive summary.
"""

import argparse
import pandas as pd
import numpy as np
import json
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Combine TFSee analysis results')
    parser.add_argument('--groseq-features', required=True, help='GRO-seq features CSV file')
    parser.add_argument('--tf-enhancer-scores', required=True, help='TF-enhancer scores CSV file')
    parser.add_argument('--motif-enrichment', required=True, help='Motif enrichment CSV file')
    parser.add_argument('--clustering-results', help='Clustering results CSV file')
    parser.add_argument('--statistics', help='Statistics CSV file')
    parser.add_argument('--config', required=True, help='TFSee configuration JSON file')
    parser.add_argument('--output-summary', required=True, help='Output summary CSV file')
    parser.add_argument('--output-features', required=True, help='Output feature matrix CSV file')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    prefix = config.get('analysis_name', 'tfsee')
    
    try:
        # Load all results
        groseq_features = pd.read_csv(args.groseq_features)
        tf_enhancer_scores = pd.read_csv(args.tf_enhancer_scores)
        motif_enrichment = pd.read_csv(args.motif_enrichment)
        
        # Try to load optional files
        clustering_results = pd.DataFrame()
        if args.clustering_results and Path(args.clustering_results).exists():
            try:
                clustering_results = pd.read_csv(args.clustering_results)
            except Exception as e:
                print(f"Warning: Could not load clustering results: {e}")
        
        statistics = pd.DataFrame()
        if args.statistics and Path(args.statistics).exists():
            try:
                statistics = pd.read_csv(args.statistics)
            except Exception as e:
                print(f"Warning: Could not load statistics: {e}")
        
        # Create comprehensive results summary
        summary_data = {
            'analysis_id': [prefix],
            'n_enhancers': [len(groseq_features)],
            'n_tf_enhancer_associations': [len(tf_enhancer_scores)],
            'n_significant_motifs': [len(motif_enrichment[motif_enrichment['p_value_corrected'] < 0.05]) if 'p_value_corrected' in motif_enrichment.columns else 0],
            'mean_enhancer_activity': [groseq_features['activity_score'].mean() if 'activity_score' in groseq_features.columns else np.nan],
            'n_clusters': [len(clustering_results['cluster_label'].unique()) if len(clustering_results) > 0 and 'cluster_label' in clustering_results.columns else 0],
            'n_statistical_tests': [len(statistics) if len(statistics) > 0 else 0]
        }
        
        # Add top associations
        if len(tf_enhancer_scores) > 0 and 'association_score' in tf_enhancer_scores.columns:
            tf_enhancer_scores_sorted = tf_enhancer_scores.sort_values('association_score', ascending=False)
            top_association = tf_enhancer_scores_sorted.iloc[0]
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
        summary_df.to_csv(args.output_summary, index=False)
        
        print(f"TFSee analysis completed successfully for {prefix}")
        print(f"Processed {len(groseq_features)} enhancers")
        print(f"Generated {len(tf_enhancer_scores)} TF-enhancer associations")
        
        # Create feature matrix combining all features
        try:
            # Set index on GRO-seq features
            if 'region_id' in groseq_features.columns:
                groseq_features = groseq_features.set_index('region_id')
            else:
                # Create a region_id column if it doesn't exist
                groseq_features['region_id'] = [f'region_{i}' for i in range(len(groseq_features))]
                groseq_features = groseq_features.set_index('region_id')
            
            # Initialize feature matrix with GRO-seq features
            feature_matrix = groseq_features.copy()
            
            # Add TF-enhancer association scores as features
            if len(tf_enhancer_scores) > 0 and 'enhancer_id' in tf_enhancer_scores.columns and 'tf_name' in tf_enhancer_scores.columns:
                # Pivot to get TF scores as columns for each enhancer
                tf_pivot = tf_enhancer_scores.pivot_table(
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
            
            # Add clustering labels if available
            if len(clustering_results) > 0 and 'region_id' in clustering_results.columns:
                clustering_indexed = clustering_results.set_index('region_id')
                feature_matrix = feature_matrix.join(clustering_indexed, how='left')
            
            # Save feature matrix
            feature_matrix.to_csv(args.output_features)
            print(f"Feature matrix created with shape: {feature_matrix.shape}")
            
        except Exception as e:
            print(f"Error creating feature matrix: {e}")
            # Create minimal feature matrix
            pd.DataFrame({'error': [str(e)]}).to_csv(args.output_features)
        
    except Exception as e:
        print(f"Error in results combination: {e}")
        # Create minimal results file
        minimal_results = pd.DataFrame({
            'analysis_id': [prefix],
            'status': ['error'],
            'error_message': [str(e)]
        })
        minimal_results.to_csv(args.output_summary, index=False)
        
        # Create minimal feature matrix
        pd.DataFrame({'error': [str(e)]}).to_csv(args.output_features)
        
        sys.exit(1)

if __name__ == "__main__":
    main()