#!/usr/bin/env python
"""
TFSee Main Analysis Pipeline

This is the main orchestrator script for TFSee analysis, integrating all
machine learning and algorithmic components for comprehensive transcription
factor-enhancer relationship analysis from GRO-seq data.

Usage:
    tfsee_main.py --config config.json --output-dir results/

Features integrated:
- GRO-seq feature engineering
- Motif probability calculation with Stouffer method
- TF-enhancer association scoring
- Multi-view clustering
- Z-score normalization
- Statistical significance testing
- Performance optimization
"""

import argparse
import logging
import json
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
import time
import warnings

# Import TFSee modules
try:
    from tfsee_groseq_features import GROSeqSignalProcessor, EnhancerActivityQuantifier
    from tfsee_motif_analysis import PWMScorer, MotifEnrichmentAnalyzer, StoufferCombination
    from tfsee_tf_enhancer_scoring import IntegratedTFEnhancerScorer, DistanceScorer
    from tfsee_multiview_clustering import MultiViewSpectralClustering, ConsensusClusteringEvaluator
    from tfsee_zscore_normalization import DataTypeSpecificNormalizer
    from tfsee_statistics import MultipleTestingCorrector, PermutationTester
    from tfsee_optimization import optimize_tfsee_pipeline, OptimizationConfig
except ImportError as e:
    logging.error(f"Could not import TFSee modules: {e}")
    logging.error("Make sure all TFSee Python modules are in your PATH")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    level=logging.INFO
)
logger = logging.getLogger("TFSee")


class TFSeeConfig:
    """Configuration management for TFSee analysis."""
    
    def __init__(self, config_file: Optional[str] = None):
        """Initialize with default configuration."""
        self.config = {
            # Input/Output
            "analysis_name": "tfsee_analysis",
            "output_directory": "./tfsee_results",
            
            # GRO-seq processing
            "window_size": 2000,
            "min_peak_height": 0.1,
            "smoothing_sigma": 2.0,
            
            # Motif analysis
            "motif_score_threshold": 0.7,
            "min_motif_score": 5.0,
            "pseudocount": 0.001,
            
            # TF-enhancer association
            "max_distance": 1000000,
            "distance_decay": "exponential",
            "association_weights": {
                "distance": 0.3,
                "motif": 0.4,
                "expression": 0.2,
                "chip": 0.1
            },
            
            # Clustering
            "n_clusters": 10,
            "clustering_method": "multiview",
            "n_clustering_runs": 10,
            "optimize_clusters": True,
            
            # Normalization
            "normalization_method": "robust_zscore",
            "outlier_method": "iqr",
            "outlier_threshold": 3.0,
            
            # Statistics
            "alpha": 0.05,
            "correction_method": "fdr_bh",
            "n_permutations": 1000,
            
            # Optimization
            "chunk_size": 1000,
            "n_workers": None,
            "memory_limit_mb": None,
            "use_sparse": True,
            "enable_caching": True,
            
            # Analysis options
            "calculate_activity_scores": True,
            "perform_clustering": True,
            "calculate_statistics": True,
            "generate_plots": True
        }
        
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
    
    def load_config(self, config_file: str):
        """Load configuration from JSON file."""
        try:
            with open(config_file, 'r') as f:
                user_config = json.load(f)
            
            # Update default config with user settings
            self._update_config(self.config, user_config)
            logger.info(f"Loaded configuration from {config_file}")
            
        except Exception as e:
            logger.error(f"Error loading config file {config_file}: {e}")
            raise
    
    def _update_config(self, base_config: dict, user_config: dict):
        """Recursively update configuration."""
        for key, value in user_config.items():
            if key in base_config and isinstance(base_config[key], dict) and isinstance(value, dict):
                self._update_config(base_config[key], value)
            else:
                base_config[key] = value
    
    def save_config(self, output_file: str):
        """Save current configuration to file."""
        with open(output_file, 'w') as f:
            json.dump(self.config, f, indent=2)
    
    def get(self, key: str, default=None):
        """Get configuration value."""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value


class TFSeeAnalysisPipeline:
    """Main TFSee analysis pipeline."""
    
    def __init__(self, config: TFSeeConfig):
        """Initialize pipeline with configuration."""
        self.config = config
        self.results = {}
        self.optimization_components = None
        
        # Create output directory
        output_dir = Path(self.config.get('output_directory'))
        output_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = output_dir
        
        # Initialize optimization components
        self._setup_optimization()
        
        logger.info(f"TFSee pipeline initialized")
        logger.info(f"Output directory: {self.output_dir}")
    
    def _setup_optimization(self):
        """Setup optimization components."""
        opt_config = OptimizationConfig(
            chunk_size=self.config.get('chunk_size'),
            n_workers=self.config.get('n_workers'),
            use_sparse=self.config.get('use_sparse'),
            memory_limit_mb=self.config.get('memory_limit_mb'),
            cache_size=128 if self.config.get('enable_caching') else 0
        )
        
        self.optimization_components = optimize_tfsee_pipeline(opt_config)
    
    def run_analysis(self, 
                    forward_bigwig: str,
                    reverse_bigwig: str,
                    enhancer_regions: str,
                    motif_database: str,
                    expression_data: Optional[str] = None,
                    chip_peaks: Optional[str] = None) -> Dict:
        """
        Run complete TFSee analysis pipeline.
        
        Args:
            forward_bigwig: Path to forward strand BigWig file
            reverse_bigwig: Path to reverse strand BigWig file
            enhancer_regions: Path to enhancer regions BED file
            motif_database: Path to motif database file
            expression_data: Optional path to expression data CSV
            chip_peaks: Optional path to ChIP-seq peaks file
            
        Returns:
            Dictionary with analysis results
        """
        logger.info("Starting TFSee analysis pipeline")
        start_time = time.time()
        
        try:
            # Step 1: Extract GRO-seq features
            logger.info("Step 1: Extracting GRO-seq features")
            groseq_features = self._extract_groseq_features(
                forward_bigwig, reverse_bigwig, enhancer_regions
            )
            
            # Step 2: Motif analysis
            logger.info("Step 2: Performing motif analysis")
            motif_results = self._analyze_motifs(enhancer_regions, motif_database)
            
            # Step 3: TF-enhancer association scoring
            logger.info("Step 3: Scoring TF-enhancer associations")
            association_scores = self._score_tf_enhancer_associations(
                enhancer_regions, motif_results, expression_data, chip_peaks
            )
            
            # Step 4: Data normalization
            logger.info("Step 4: Normalizing features")
            normalized_features = self._normalize_features(groseq_features)
            
            # Step 5: Multi-view clustering (optional)
            clustering_results = None
            if self.config.get('perform_clustering'):
                logger.info("Step 5: Performing multi-view clustering")
                clustering_results = self._perform_clustering(
                    normalized_features, association_scores
                )
            
            # Step 6: Statistical analysis
            statistics_results = None
            if self.config.get('calculate_statistics'):
                logger.info("Step 6: Performing statistical analysis")
                statistics_results = self._calculate_statistics(
                    normalized_features, association_scores, clustering_results
                )
            
            # Step 7: Generate comprehensive results
            logger.info("Step 7: Generating final results")
            final_results = self._generate_final_results(
                groseq_features, motif_results, association_scores,
                normalized_features, clustering_results, statistics_results
            )
            
            # Step 8: Save results
            logger.info("Step 8: Saving results")
            self._save_results(final_results)
            
            end_time = time.time()
            logger.info(f"TFSee analysis completed in {end_time - start_time:.2f} seconds")
            
            return final_results
            
        except Exception as e:
            logger.error(f"Error in TFSee analysis: {e}")
            raise
    
    def _extract_groseq_features(self, 
                                forward_bigwig: str, 
                                reverse_bigwig: str, 
                                enhancer_regions: str) -> pd.DataFrame:
        """Extract GRO-seq features from enhancer regions."""
        # This would typically call the GRO-seq feature extraction module
        # For now, return placeholder
        features_file = self.output_dir / "groseq_features.csv"
        
        # Call GRO-seq feature extraction
        cmd = [
            "tfsee_groseq_features.py",
            "--enhancers", enhancer_regions,
            "--forward-bigwig", forward_bigwig,
            "--reverse-bigwig", reverse_bigwig,
            "--output", str(features_file),
            "--window-size", str(self.config.get('window_size')),
            "--min-peak-height", str(self.config.get('min_peak_height')),
            "--smoothing-sigma", str(self.config.get('smoothing_sigma'))
        ]
        
        if self.config.get('calculate_activity_scores'):
            cmd.append("--calculate-activity")
        
        self._run_command(cmd)
        
        return pd.read_csv(features_file)
    
    def _analyze_motifs(self, enhancer_regions: str, motif_database: str) -> pd.DataFrame:
        """Perform motif analysis on enhancer regions."""
        motif_file = self.output_dir / "motif_enrichment.csv"
        
        cmd = [
            "tfsee_motif_analysis.py",
            "--foreground", enhancer_regions,
            "--background", enhancer_regions,  # Would use appropriate background
            "--motifs", motif_database,
            "--output", str(motif_file),
            "--min-score", str(self.config.get('min_motif_score')),
            "--pseudocount", str(self.config.get('pseudocount'))
        ]
        
        self._run_command(cmd)
        
        return pd.read_csv(motif_file)
    
    def _score_tf_enhancer_associations(self, 
                                      enhancer_regions: str,
                                      motif_results: pd.DataFrame,
                                      expression_data: Optional[str],
                                      chip_peaks: Optional[str]) -> pd.DataFrame:
        """Score TF-enhancer associations."""
        scores_file = self.output_dir / "tf_enhancer_scores.csv"
        
        cmd = [
            "tfsee_tf_enhancer_scoring.py",
            "--tf-regions", enhancer_regions,  # Would use actual TF regions
            "--enhancer-regions", enhancer_regions,
            "--output", str(scores_file),
            "--max-distance", str(self.config.get('max_distance')),
            "--top-k", "10000"
        ]
        
        if expression_data:
            cmd.extend(["--expression-data", expression_data])
        
        if chip_peaks:
            cmd.extend(["--chip-peaks", chip_peaks])
        
        self._run_command(cmd)
        
        return pd.read_csv(scores_file)
    
    def _normalize_features(self, features: pd.DataFrame) -> pd.DataFrame:
        """Normalize feature data."""
        features_file = self.output_dir / "normalized_features.csv"
        input_file = self.output_dir / "temp_features.csv"
        
        # Save features temporarily
        features.to_csv(input_file)
        
        cmd = [
            "tfsee_zscore_normalization.py",
            "--input", str(input_file),
            "--data-type", "expression",
            "--output", str(features_file),
            "--log-transform"
        ]
        
        self._run_command(cmd)
        
        # Clean up temporary file
        input_file.unlink()
        
        return pd.read_csv(features_file, index_col=0)
    
    def _perform_clustering(self, 
                          features: pd.DataFrame, 
                          associations: pd.DataFrame) -> pd.DataFrame:
        """Perform multi-view clustering."""
        clustering_file = self.output_dir / "clustering_results.csv"
        features_file = self.output_dir / "temp_clustering_features.csv"
        associations_file = self.output_dir / "temp_clustering_associations.csv"
        
        # Save data temporarily
        features.to_csv(features_file)
        associations.to_csv(associations_file)
        
        cmd = [
            "tfsee_multiview_clustering.py",
            "--view-data", str(features_file), str(associations_file),
            "--view-names", "features", "associations",
            "--n-clusters", str(self.config.get('n_clusters')),
            "--output-prefix", str(self.output_dir / "clustering"),
            "--method", self.config.get('clustering_method'),
            "--n-runs", str(self.config.get('n_clustering_runs'))
        ]
        
        if self.config.get('optimize_clusters'):
            cmd.append("--optimize-k")
        
        self._run_command(cmd)
        
        # Clean up temporary files
        features_file.unlink()
        associations_file.unlink()
        
        return pd.read_csv(clustering_file)
    
    def _calculate_statistics(self, 
                            features: pd.DataFrame,
                            associations: pd.DataFrame,
                            clustering: Optional[pd.DataFrame]) -> pd.DataFrame:
        """Perform statistical analysis."""
        stats_file = self.output_dir / "statistics.csv"
        data_file = self.output_dir / "temp_stats_data.csv"
        
        # Combine data for analysis
        combined_data = features.copy()
        
        # Save temporarily
        combined_data.to_csv(data_file)
        
        cmd = [
            "tfsee_statistics.py",
            "--data", str(data_file),
            "--output", str(stats_file),
            "--alpha", str(self.config.get('alpha')),
            "--correction-method", self.config.get('correction_method'),
            "--n-permutations", str(self.config.get('n_permutations'))
        ]
        
        self._run_command(cmd)
        
        # Clean up temporary file
        data_file.unlink()
        
        return pd.read_csv(stats_file)
    
    def _generate_final_results(self, 
                              groseq_features: pd.DataFrame,
                              motif_results: pd.DataFrame,
                              association_scores: pd.DataFrame,
                              normalized_features: pd.DataFrame,
                              clustering_results: Optional[pd.DataFrame],
                              statistics_results: Optional[pd.DataFrame]) -> Dict:
        """Generate comprehensive final results."""
        
        results = {
            'analysis_config': self.config.config,
            'summary': {
                'n_enhancers': len(groseq_features),
                'n_motifs_tested': len(motif_results),
                'n_tf_enhancer_associations': len(association_scores),
                'n_features': normalized_features.shape[1],
                'analysis_timestamp': time.time()
            },
            'groseq_features': groseq_features,
            'motif_enrichment': motif_results,
            'tf_enhancer_associations': association_scores,
            'normalized_features': normalized_features,
            'clustering_results': clustering_results,
            'statistics': statistics_results
        }
        
        # Add summary statistics
        if 'activity_score' in groseq_features.columns:
            results['summary']['mean_activity_score'] = groseq_features['activity_score'].mean()
            results['summary']['high_activity_enhancers'] = (groseq_features['activity_score'] > 0.8).sum()
        
        if len(motif_results) > 0 and 'p_value_corrected' in motif_results.columns:
            results['summary']['significant_motifs'] = (motif_results['p_value_corrected'] < 0.05).sum()
        
        if clustering_results is not None and len(clustering_results) > 0:
            results['summary']['n_clusters'] = clustering_results['cluster_label'].nunique()
        
        return results
    
    def _save_results(self, results: Dict):
        """Save all results to files."""
        # Save main results summary
        summary_file = self.output_dir / f"{self.config.get('analysis_name')}_results.json"
        
        # Convert DataFrames to dictionaries for JSON serialization
        json_results = {
            'analysis_config': results['analysis_config'],
            'summary': results['summary']
        }
        
        with open(summary_file, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        # Save individual DataFrames
        for key, data in results.items():
            if isinstance(data, pd.DataFrame):
                filename = self.output_dir / f"{self.config.get('analysis_name')}_{key}.csv"
                data.to_csv(filename, index=False)
        
        # Save configuration
        config_file = self.output_dir / f"{self.config.get('analysis_name')}_config.json"
        self.config.save_config(str(config_file))
        
        logger.info(f"Results saved to {self.output_dir}")
    
    def _run_command(self, cmd: List[str]):
        """Run external command."""
        import subprocess
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stdout:
                logger.debug(f"Command output: {result.stdout}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed: {' '.join(cmd)}")
            logger.error(f"Error: {e.stderr}")
            raise


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee: Comprehensive Transcription Factor-Enhancer Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required inputs
    parser.add_argument(
        '--forward-bigwig', required=True,
        help="Forward strand GRO-seq BigWig file"
    )
    parser.add_argument(
        '--reverse-bigwig', required=True,
        help="Reverse strand GRO-seq BigWig file"
    )
    parser.add_argument(
        '--enhancer-regions', required=True,
        help="BED file with enhancer regions"
    )
    parser.add_argument(
        '--motif-database', required=True,
        help="MEME format motif database file"
    )
    
    # Optional inputs
    parser.add_argument(
        '--expression-data',
        help="CSV file with gene expression data"
    )
    parser.add_argument(
        '--chip-peaks',
        help="Pickle file with ChIP-seq peak data"
    )
    parser.add_argument(
        '--config', '-c',
        help="JSON configuration file"
    )
    
    # Output options
    parser.add_argument(
        '--output-dir', '-o', default='./tfsee_results',
        help="Output directory for results"
    )
    parser.add_argument(
        '--analysis-name', default='tfsee_analysis',
        help="Name for this analysis"
    )
    
    # Analysis options
    parser.add_argument(
        '--n-clusters', type=int, default=10,
        help="Number of clusters for clustering analysis"
    )
    parser.add_argument(
        '--n-workers', type=int,
        help="Number of parallel workers (default: auto)"
    )
    parser.add_argument(
        '--memory-limit', type=float,
        help="Memory limit in MB (default: auto)"
    )
    parser.add_argument(
        '--skip-clustering', action='store_true',
        help="Skip clustering analysis"
    )
    parser.add_argument(
        '--skip-statistics', action='store_true',
        help="Skip statistical analysis"
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize configuration
    config = TFSeeConfig(args.config)
    
    # Update config with command line arguments
    config.config['output_directory'] = args.output_dir
    config.config['analysis_name'] = args.analysis_name
    config.config['n_clusters'] = args.n_clusters
    
    if args.n_workers:
        config.config['n_workers'] = args.n_workers
    if args.memory_limit:
        config.config['memory_limit_mb'] = args.memory_limit
    if args.skip_clustering:
        config.config['perform_clustering'] = False
    if args.skip_statistics:
        config.config['calculate_statistics'] = False
    
    # Initialize and run pipeline
    try:
        pipeline = TFSeeAnalysisPipeline(config)
        
        results = pipeline.run_analysis(
            forward_bigwig=args.forward_bigwig,
            reverse_bigwig=args.reverse_bigwig,
            enhancer_regions=args.enhancer_regions,
            motif_database=args.motif_database,
            expression_data=args.expression_data,
            chip_peaks=args.chip_peaks
        )
        
        # Print summary
        print("\nTFSee Analysis Summary:")
        print("=" * 50)
        for key, value in results['summary'].items():
            print(f"{key}: {value}")
        
        print(f"\nResults saved to: {args.output_dir}")
        print("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"TFSee analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()