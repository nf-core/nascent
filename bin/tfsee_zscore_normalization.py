#!/usr/bin/env python
"""
TFSee Z-Score Normalization Module

This module implements comprehensive Z-score normalization strategies for
different genomic data types used in TFSee analysis.

Key Features:
- Data type-specific normalization strategies
- Robust Z-score calculation with outlier handling
- Batch effect correction
- Missing value imputation
- Distribution transformation methods
- Quality control metrics

Supported data types:
- RNA-seq expression data
- ChIP-seq peak intensities
- Chromatin accessibility (ATAC-seq/DNase-seq)
- Motif scores
- Distance-based features
- Multi-condition comparisons
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import zscore, normaltest, boxcox, yeojohnson
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer
from sklearn.impute import KNNImputer, IterativeImputer
from sklearn.decomposition import PCA
import warnings
from typing import List, Dict, Tuple, Optional, Union, Literal
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class NormalizationResult:
    """Container for normalization results and metadata."""
    normalized_data: np.ndarray
    original_data: np.ndarray
    method: str
    parameters: Dict
    quality_metrics: Dict
    outliers_detected: np.ndarray
    missing_imputed: np.ndarray


class RobustZScoreNormalizer:
    """Robust Z-score normalization with outlier handling."""
    
    def __init__(self, 
                 outlier_method: str = "iqr",
                 outlier_threshold: float = 3.0,
                 missing_strategy: str = "median"):
        """
        Initialize robust Z-score normalizer.
        
        Args:
            outlier_method: Method for outlier detection ('iqr', 'zscore', 'isolation')
            outlier_threshold: Threshold for outlier detection
            missing_strategy: Strategy for handling missing values ('median', 'mean', 'knn')
        """
        self.outlier_method = outlier_method
        self.outlier_threshold = outlier_threshold
        self.missing_strategy = missing_strategy
    
    def detect_outliers(self, data: np.ndarray) -> np.ndarray:
        """
        Detect outliers in data.
        
        Args:
            data: Input data array
            
        Returns:
            Boolean array indicating outliers
        """
        if self.outlier_method == "iqr":
            q1 = np.percentile(data, 25)
            q3 = np.percentile(data, 75)
            iqr = q3 - q1
            lower_bound = q1 - self.outlier_threshold * iqr
            upper_bound = q3 + self.outlier_threshold * iqr
            outliers = (data < lower_bound) | (data > upper_bound)
        
        elif self.outlier_method == "zscore":
            z_scores = np.abs(zscore(data, nan_policy='omit'))
            outliers = z_scores > self.outlier_threshold
        
        elif self.outlier_method == "modified_zscore":
            # Modified Z-score using median
            median = np.median(data)
            mad = np.median(np.abs(data - median))
            modified_z_scores = 0.6745 * (data - median) / mad
            outliers = np.abs(modified_z_scores) > self.outlier_threshold
        
        else:
            raise ValueError(f"Unknown outlier method: {self.outlier_method}")
        
        return outliers
    
    def handle_missing_values(self, data: np.ndarray) -> np.ndarray:
        """
        Handle missing values in data.
        
        Args:
            data: Input data with potential missing values
            
        Returns:
            Data with missing values imputed
        """
        if not np.any(np.isnan(data)):
            return data
        
        if self.missing_strategy == "median":
            imputer = lambda x: np.nanmedian(x)
        elif self.missing_strategy == "mean":
            imputer = lambda x: np.nanmean(x)
        elif self.missing_strategy == "knn":
            if data.ndim > 1:
                knn_imputer = KNNImputer(n_neighbors=5)
                return knn_imputer.fit_transform(data)
            else:
                imputer = lambda x: np.nanmedian(x)
        else:
            raise ValueError(f"Unknown missing strategy: {self.missing_strategy}")
        
        if data.ndim == 1:
            imputed_value = imputer(data)
            data = np.where(np.isnan(data), imputed_value, data)
        else:
            for i in range(data.shape[1]):
                column = data[:, i]
                if np.any(np.isnan(column)):
                    imputed_value = imputer(column)
                    data[:, i] = np.where(np.isnan(column), imputed_value, column)
        
        return data
    
    def robust_zscore(self, data: np.ndarray) -> Tuple[np.ndarray, Dict]:
        """
        Calculate robust Z-scores.
        
        Args:
            data: Input data array
            
        Returns:
            Tuple of (normalized_data, parameters)
        """
        # Use median and MAD for robust normalization
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        
        # Avoid division by zero
        if mad == 0:
            mad = np.std(data)
            if mad == 0:
                mad = 1.0
        
        # Calculate robust Z-scores
        z_scores = (data - median) / (1.4826 * mad)  # 1.4826 makes MAD consistent with std
        
        parameters = {
            'median': median,
            'mad': mad,
            'scale_factor': 1.4826
        }
        
        return z_scores, parameters


class DataTypeSpecificNormalizer:
    """Data type-specific normalization strategies."""
    
    def __init__(self):
        self.robust_normalizer = RobustZScoreNormalizer()
    
    def normalize_expression_data(self, 
                                expression_data: np.ndarray,
                                log_transform: bool = True,
                                quantile_normalize: bool = False) -> NormalizationResult:
        """
        Normalize RNA-seq expression data.
        
        Args:
            expression_data: Expression matrix (genes x samples)
            log_transform: Apply log2 transformation
            quantile_normalize: Apply quantile normalization
            
        Returns:
            NormalizationResult object
        """
        original_data = expression_data.copy()
        data = expression_data.copy()
        
        # Handle zeros for log transformation
        if log_transform:
            # Add pseudocount to avoid log(0)
            pseudocount = 1.0
            data = np.log2(data + pseudocount)
        
        # Quantile normalization
        if quantile_normalize:
            qt = QuantileTransformer(output_distribution='normal')
            if data.ndim > 1:
                data = qt.fit_transform(data.T).T
            else:
                data = qt.fit_transform(data.reshape(-1, 1)).flatten()
        
        # Handle missing values
        missing_mask = np.isnan(data)
        data = self.robust_normalizer.handle_missing_values(data)
        
        # Detect outliers
        outliers = np.zeros_like(data, dtype=bool)
        if data.ndim > 1:
            for i in range(data.shape[0]):
                outliers[i, :] = self.robust_normalizer.detect_outliers(data[i, :])
        else:
            outliers = self.robust_normalizer.detect_outliers(data)
        
        # Z-score normalization (per gene)
        normalized_data = np.zeros_like(data)
        if data.ndim > 1:
            for i in range(data.shape[0]):
                normalized_data[i, :], _ = self.robust_normalizer.robust_zscore(data[i, :])
        else:
            normalized_data, _ = self.robust_normalizer.robust_zscore(data)
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(original_data, normalized_data)
        
        return NormalizationResult(
            normalized_data=normalized_data,
            original_data=original_data,
            method="expression_specific",
            parameters={'log_transform': log_transform, 'quantile_normalize': quantile_normalize},
            quality_metrics=quality_metrics,
            outliers_detected=outliers,
            missing_imputed=missing_mask
        )
    
    def normalize_chipseq_data(self, 
                             chipseq_data: np.ndarray,
                             background_correction: bool = True) -> NormalizationResult:
        """
        Normalize ChIP-seq peak intensity data.
        
        Args:
            chipseq_data: ChIP-seq intensity matrix
            background_correction: Apply background correction
            
        Returns:
            NormalizationResult object
        """
        original_data = chipseq_data.copy()
        data = chipseq_data.copy()
        
        # Background correction (subtract minimum)
        if background_correction:
            min_val = np.min(data[data > 0]) if np.any(data > 0) else 0
            data = np.maximum(data - min_val, 0)
        
        # Log transformation for ChIP-seq data
        data = np.log1p(data)  # log(1 + x) to handle zeros
        
        # Handle missing values
        missing_mask = np.isnan(data)
        data = self.robust_normalizer.handle_missing_values(data)
        
        # Detect outliers
        outliers = self.robust_normalizer.detect_outliers(data.flatten()).reshape(data.shape)
        
        # Robust normalization
        if data.ndim > 1:
            normalized_data = np.zeros_like(data)
            for i in range(data.shape[1]):
                normalized_data[:, i], _ = self.robust_normalizer.robust_zscore(data[:, i])
        else:
            normalized_data, _ = self.robust_normalizer.robust_zscore(data)
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(original_data, normalized_data)
        
        return NormalizationResult(
            normalized_data=normalized_data,
            original_data=original_data,
            method="chipseq_specific",
            parameters={'background_correction': background_correction},
            quality_metrics=quality_metrics,
            outliers_detected=outliers,
            missing_imputed=missing_mask
        )
    
    def normalize_accessibility_data(self, 
                                   accessibility_data: np.ndarray,
                                   peak_calling_threshold: Optional[float] = None) -> NormalizationResult:
        """
        Normalize chromatin accessibility data (ATAC-seq/DNase-seq).
        
        Args:
            accessibility_data: Accessibility signal matrix
            peak_calling_threshold: Threshold for binarizing peaks
            
        Returns:
            NormalizationResult object
        """
        original_data = accessibility_data.copy()
        data = accessibility_data.copy()
        
        # Apply peak calling threshold if provided
        if peak_calling_threshold is not None:
            data = (data > peak_calling_threshold).astype(float)
        else:
            # Quantile transformation for continuous accessibility data
            qt = QuantileTransformer(output_distribution='normal')
            if data.ndim > 1:
                data = qt.fit_transform(data)
            else:
                data = qt.fit_transform(data.reshape(-1, 1)).flatten()
        
        # Handle missing values
        missing_mask = np.isnan(data)
        data = self.robust_normalizer.handle_missing_values(data)
        
        # For binary data, no outlier detection needed
        if peak_calling_threshold is not None:
            outliers = np.zeros_like(data, dtype=bool)
            normalized_data = data  # Already binary
        else:
            # Detect outliers for continuous data
            outliers = self.robust_normalizer.detect_outliers(data.flatten()).reshape(data.shape)
            
            # Standard normalization
            if data.ndim > 1:
                scaler = StandardScaler()
                normalized_data = scaler.fit_transform(data)
            else:
                normalized_data, _ = self.robust_normalizer.robust_zscore(data)
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(original_data, normalized_data)
        
        return NormalizationResult(
            normalized_data=normalized_data,
            original_data=original_data,
            method="accessibility_specific",
            parameters={'peak_calling_threshold': peak_calling_threshold},
            quality_metrics=quality_metrics,
            outliers_detected=outliers,
            missing_imputed=missing_mask
        )
    
    def normalize_motif_scores(self, 
                             motif_scores: np.ndarray,
                             score_threshold: Optional[float] = None) -> NormalizationResult:
        """
        Normalize motif scores.
        
        Args:
            motif_scores: Motif score matrix
            score_threshold: Threshold for significant motifs
            
        Returns:
            NormalizationResult object
        """
        original_data = motif_scores.copy()
        data = motif_scores.copy()
        
        # Apply threshold if provided
        if score_threshold is not None:
            data = np.where(data >= score_threshold, data, 0)
        
        # Handle missing values
        missing_mask = np.isnan(data)
        data = self.robust_normalizer.handle_missing_values(data)
        
        # Detect outliers
        outliers = self.robust_normalizer.detect_outliers(data.flatten()).reshape(data.shape)
        
        # Min-max normalization for motif scores (preserve relative ranking)
        data_min = np.min(data)
        data_max = np.max(data)
        
        if data_max > data_min:
            normalized_data = (data - data_min) / (data_max - data_min)
        else:
            normalized_data = np.ones_like(data) * 0.5
        
        # Optional: Apply Z-score on top of min-max
        if normalized_data.std() > 0:
            normalized_data = (normalized_data - normalized_data.mean()) / normalized_data.std()
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(original_data, normalized_data)
        
        return NormalizationResult(
            normalized_data=normalized_data,
            original_data=original_data,
            method="motif_specific",
            parameters={'score_threshold': score_threshold},
            quality_metrics=quality_metrics,
            outliers_detected=outliers,
            missing_imputed=missing_mask
        )
    
    def normalize_distance_features(self, 
                                  distance_data: np.ndarray,
                                  max_distance: float = 1000000) -> NormalizationResult:
        """
        Normalize distance-based features.
        
        Args:
            distance_data: Distance values
            max_distance: Maximum distance for normalization
            
        Returns:
            NormalizationResult object
        """
        original_data = distance_data.copy()
        data = distance_data.copy()
        
        # Cap distances at maximum
        data = np.minimum(data, max_distance)
        
        # Log transformation (distance + 1 to handle zeros)
        data = np.log1p(data)
        
        # Handle missing values
        missing_mask = np.isnan(data)
        data = self.robust_normalizer.handle_missing_values(data)
        
        # Outlier detection
        outliers = self.robust_normalizer.detect_outliers(data.flatten()).reshape(data.shape)
        
        # Inverse transformation (closer = higher score)
        max_log_dist = np.log1p(max_distance)
        normalized_data = 1 - (data / max_log_dist)
        normalized_data = np.clip(normalized_data, 0, 1)
        
        # Apply Z-score normalization
        if normalized_data.std() > 0:
            normalized_data = (normalized_data - normalized_data.mean()) / normalized_data.std()
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(original_data, normalized_data)
        
        return NormalizationResult(
            normalized_data=normalized_data,
            original_data=original_data,
            method="distance_specific",
            parameters={'max_distance': max_distance},
            quality_metrics=quality_metrics,
            outliers_detected=outliers,
            missing_imputed=missing_mask
        )
    
    def _calculate_quality_metrics(self, 
                                 original_data: np.ndarray, 
                                 normalized_data: np.ndarray) -> Dict:
        """Calculate quality metrics for normalization."""
        metrics = {}
        
        # Basic statistics
        metrics['original_mean'] = np.mean(original_data)
        metrics['original_std'] = np.std(original_data)
        metrics['normalized_mean'] = np.mean(normalized_data)
        metrics['normalized_std'] = np.std(normalized_data)
        
        # Normality tests (on subset if data is large)
        sample_size = min(5000, len(normalized_data.flatten()))
        sample_indices = np.random.choice(
            len(normalized_data.flatten()), 
            sample_size, 
            replace=False
        )
        sample_data = normalized_data.flatten()[sample_indices]
        
        # Shapiro-Wilk test for small samples, Anderson-Darling for larger
        if sample_size <= 5000:
            try:
                stat, p_value = normaltest(sample_data)
                metrics['normality_test_stat'] = stat
                metrics['normality_test_pvalue'] = p_value
                metrics['is_normal'] = p_value > 0.05
            except:
                metrics['normality_test_stat'] = np.nan
                metrics['normality_test_pvalue'] = np.nan
                metrics['is_normal'] = False
        
        # Skewness and kurtosis
        metrics['skewness'] = stats.skew(normalized_data.flatten())
        metrics['kurtosis'] = stats.kurtosis(normalized_data.flatten())
        
        # Outlier percentage
        if normalized_data.ndim > 1:
            z_scores = np.abs(zscore(normalized_data, axis=None, nan_policy='omit'))
        else:
            z_scores = np.abs(zscore(normalized_data, nan_policy='omit'))
        
        outlier_mask = z_scores > 3
        metrics['outlier_percentage'] = np.sum(outlier_mask) / len(normalized_data.flatten()) * 100
        
        return metrics


class BatchEffectCorrector:
    """Correct for batch effects in genomic data."""
    
    def __init__(self, method: str = "combat"):
        """
        Initialize batch effect corrector.
        
        Args:
            method: Method for batch correction ('combat', 'pca', 'linear')
        """
        self.method = method
    
    def correct_batch_effects(self, 
                            data: np.ndarray, 
                            batch_labels: np.ndarray) -> np.ndarray:
        """
        Correct batch effects in data.
        
        Args:
            data: Data matrix (features x samples)
            batch_labels: Batch labels for each sample
            
        Returns:
            Batch-corrected data
        """
        if self.method == "linear":
            return self._linear_batch_correction(data, batch_labels)
        elif self.method == "pca":
            return self._pca_batch_correction(data, batch_labels)
        else:
            logger.warning(f"Batch correction method '{self.method}' not implemented. Returning original data.")
            return data
    
    def _linear_batch_correction(self, 
                               data: np.ndarray, 
                               batch_labels: np.ndarray) -> np.ndarray:
        """Linear model-based batch correction."""
        corrected_data = data.copy()
        
        # For each feature, fit linear model and remove batch effect
        unique_batches = np.unique(batch_labels)
        
        if len(unique_batches) <= 1:
            return data
        
        for i in range(data.shape[0]):
            feature_data = data[i, :]
            
            # Calculate batch means
            batch_means = {}
            overall_mean = np.mean(feature_data)
            
            for batch in unique_batches:
                batch_mask = batch_labels == batch
                batch_means[batch] = np.mean(feature_data[batch_mask])
            
            # Correct each sample
            for j, batch in enumerate(batch_labels):
                corrected_data[i, j] = feature_data[j] - (batch_means[batch] - overall_mean)
        
        return corrected_data
    
    def _pca_batch_correction(self, 
                            data: np.ndarray, 
                            batch_labels: np.ndarray) -> np.ndarray:
        """PCA-based batch correction."""
        # Simple approach: remove first PC if it correlates with batch
        pca = PCA()
        transformed = pca.fit_transform(data.T)
        
        # Check correlation between first PC and batch
        unique_batches = np.unique(batch_labels)
        batch_numeric = np.array([np.where(unique_batches == b)[0][0] for b in batch_labels])
        
        correlation = np.corrcoef(transformed[:, 0], batch_numeric)[0, 1]
        
        if abs(correlation) > 0.5:  # If first PC is batch-related
            # Remove first PC
            transformed[:, 0] = 0
            corrected_data = pca.inverse_transform(transformed).T
        else:
            corrected_data = data
        
        return corrected_data


def visualize_normalization_results(result: NormalizationResult, 
                                  output_prefix: str = "normalization"):
    """
    Visualize normalization results.
    
    Args:
        result: NormalizationResult object
        output_prefix: Prefix for output files
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Original data distribution
    axes[0, 0].hist(result.original_data.flatten(), bins=50, alpha=0.7, density=True)
    axes[0, 0].set_title('Original Data Distribution')
    axes[0, 0].set_xlabel('Value')
    axes[0, 0].set_ylabel('Density')
    
    # Normalized data distribution
    axes[0, 1].hist(result.normalized_data.flatten(), bins=50, alpha=0.7, density=True)
    axes[0, 1].set_title('Normalized Data Distribution')
    axes[0, 1].set_xlabel('Value')
    axes[0, 1].set_ylabel('Density')
    
    # Q-Q plot for normality check
    stats.probplot(result.normalized_data.flatten()[:1000], dist="norm", plot=axes[1, 0])
    axes[1, 0].set_title('Q-Q Plot (Normality Check)')
    
    # Outlier visualization
    if np.any(result.outliers_detected):
        outlier_data = result.normalized_data.copy()
        outlier_data[~result.outliers_detected] = np.nan
        
        axes[1, 1].scatter(range(len(result.normalized_data.flatten())), 
                          result.normalized_data.flatten(), 
                          alpha=0.5, s=1, label='Normal')
        axes[1, 1].scatter(range(len(outlier_data.flatten())), 
                          outlier_data.flatten(), 
                          color='red', s=1, label='Outliers')
        axes[1, 1].set_title('Outlier Detection')
        axes[1, 1].set_xlabel('Index')
        axes[1, 1].set_ylabel('Normalized Value')
        axes[1, 1].legend()
    else:
        axes[1, 1].text(0.5, 0.5, 'No outliers detected', 
                       transform=axes[1, 1].transAxes, ha='center', va='center')
        axes[1, 1].set_title('Outlier Detection')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_visualization.png", dpi=300, bbox_inches='tight')
    plt.close()


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee Z-Score Normalization for Genomic Data"
    )
    parser.add_argument(
        '--input', '-i', required=True,
        help="Input data file (CSV format)"
    )
    parser.add_argument(
        '--data-type', required=True,
        choices=['expression', 'chipseq', 'accessibility', 'motif', 'distance'],
        help="Type of genomic data"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output file for normalized data"
    )
    parser.add_argument(
        '--batch-file',
        help="File with batch labels for batch effect correction"
    )
    parser.add_argument(
        '--log-transform', action='store_true',
        help="Apply log transformation (for expression data)"
    )
    parser.add_argument(
        '--quantile-normalize', action='store_true',
        help="Apply quantile normalization (for expression data)"
    )
    parser.add_argument(
        '--threshold', type=float,
        help="Threshold for data filtering"
    )
    parser.add_argument(
        '--visualize', action='store_true',
        help="Generate visualization plots"
    )
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading {args.data_type} data from {args.input}")
    data = pd.read_csv(args.input, index_col=0).values
    
    # Initialize normalizer
    normalizer = DataTypeSpecificNormalizer()
    
    # Perform normalization based on data type
    if args.data_type == 'expression':
        result = normalizer.normalize_expression_data(
            data, 
            log_transform=args.log_transform,
            quantile_normalize=args.quantile_normalize
        )
    
    elif args.data_type == 'chipseq':
        result = normalizer.normalize_chipseq_data(data)
    
    elif args.data_type == 'accessibility':
        result = normalizer.normalize_accessibility_data(
            data, 
            peak_calling_threshold=args.threshold
        )
    
    elif args.data_type == 'motif':
        result = normalizer.normalize_motif_scores(
            data, 
            score_threshold=args.threshold
        )
    
    elif args.data_type == 'distance':
        max_dist = args.threshold if args.threshold else 1000000
        result = normalizer.normalize_distance_features(data, max_distance=max_dist)
    
    # Batch effect correction if requested
    if args.batch_file:
        logger.info("Applying batch effect correction...")
        batch_labels = pd.read_csv(args.batch_file, header=None).values.flatten()
        
        batch_corrector = BatchEffectCorrector()
        result.normalized_data = batch_corrector.correct_batch_effects(
            result.normalized_data, batch_labels
        )
    
    # Save results
    logger.info(f"Saving normalized data to {args.output}")
    
    # Save normalized data
    normalized_df = pd.DataFrame(result.normalized_data)
    normalized_df.to_csv(args.output)
    
    # Save quality metrics
    metrics_file = args.output.replace('.csv', '_metrics.json')
    import json
    with open(metrics_file, 'w') as f:
        json.dump(result.quality_metrics, f, indent=2)
    
    # Save outlier information
    if np.any(result.outliers_detected):
        outliers_file = args.output.replace('.csv', '_outliers.csv')
        outliers_df = pd.DataFrame(result.outliers_detected)
        outliers_df.to_csv(outliers_file)
    
    # Generate visualizations if requested
    if args.visualize:
        output_prefix = args.output.replace('.csv', '')
        visualize_normalization_results(result, output_prefix)
    
    # Print summary
    logger.info("Normalization completed successfully!")
    logger.info(f"Method: {result.method}")
    logger.info(f"Original data - Mean: {result.quality_metrics['original_mean']:.3f}, "
               f"Std: {result.quality_metrics['original_std']:.3f}")
    logger.info(f"Normalized data - Mean: {result.quality_metrics['normalized_mean']:.3f}, "
               f"Std: {result.quality_metrics['normalized_std']:.3f}")
    logger.info(f"Outliers detected: {result.quality_metrics['outlier_percentage']:.1f}%")
    
    if 'is_normal' in result.quality_metrics:
        logger.info(f"Data appears normal: {result.quality_metrics['is_normal']}")


if __name__ == "__main__":
    main()