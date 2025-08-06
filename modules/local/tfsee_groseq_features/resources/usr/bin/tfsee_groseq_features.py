#!/usr/bin/env python
"""
TFSee GRO-seq Feature Engineering Module

This module implements feature engineering for enhancer activity prediction
from GRO-seq (Global Run-On sequencing) signals and related nascent transcription data.

Key Features:
- Signal intensity features from GRO-seq data
- Bidirectional transcription detection
- Temporal dynamics analysis
- Peak shape characterization
- Enhancer activity quantification
- Integration with chromatin accessibility
- Multi-condition comparison features

GRO-seq specific features:
- Transcriptional burst characteristics
- Pause site identification
- Elongation rate estimates
- Directional transcription bias
- Signal-to-noise ratios
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats, signal
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.ndimage import gaussian_filter1d
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import warnings
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
import pysam
import pyBigWig
from collections import defaultdict

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class EnhancerRegion:
    """Represents an enhancer region with coordinates."""
    chromosome: str
    start: int
    end: int
    name: str = ""
    strand: str = "."
    
    @property
    def length(self) -> int:
        return self.end - self.start
    
    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


@dataclass
class GROSeqFeatures:
    """Container for GRO-seq derived features."""
    region_id: str
    
    # Basic intensity features
    mean_signal: float
    max_signal: float
    total_signal: float
    signal_variance: float
    
    # Bidirectional features
    forward_signal: float
    reverse_signal: float
    bidirectional_ratio: float
    
    # Peak characteristics
    n_peaks: int
    peak_heights: List[float]
    peak_widths: List[float]
    peak_prominences: List[float]
    
    # Shape features
    skewness: float
    kurtosis: float
    signal_entropy: float
    
    # Spatial features
    signal_spread: float
    center_bias: float
    
    # Quality metrics
    snr: float  # Signal-to-noise ratio
    coverage_fraction: float


class GROSeqSignalProcessor:
    """Process GRO-seq signals and extract features."""
    
    def __init__(self, 
                 window_size: int = 2000,
                 min_peak_height: float = 0.1,
                 smoothing_sigma: float = 2.0):
        """
        Initialize GRO-seq signal processor.
        
        Args:
            window_size: Size of window around enhancer center
            min_peak_height: Minimum height for peak detection
            smoothing_sigma: Sigma for Gaussian smoothing
        """
        self.window_size = window_size
        self.min_peak_height = min_peak_height
        self.smoothing_sigma = smoothing_sigma
    
    def extract_signal_from_bigwig(self, 
                                 bigwig_file: str, 
                                 region: EnhancerRegion,
                                 strand: str = "+") -> np.ndarray:
        """
        Extract signal from BigWig file for a genomic region.
        
        Args:
            bigwig_file: Path to BigWig file
            region: Enhancer region
            strand: Strand information ("+", "-", or "both")
            
        Returns:
            Signal array
        """
        try:
            bw = pyBigWig.open(bigwig_file)
            
            # Define extraction window
            center = region.center
            start = max(0, center - self.window_size // 2)
            end = center + self.window_size // 2
            
            # Ensure coordinates are within chromosome bounds
            chrom_size = bw.chroms().get(region.chromosome, 0)
            if chrom_size > 0:
                end = min(end, chrom_size)
            
            # Extract signal
            signal_values = bw.values(region.chromosome, start, end)
            
            # Handle missing values
            signal_values = np.array(signal_values)
            signal_values = np.nan_to_num(signal_values, nan=0.0)
            
            bw.close()
            
            return signal_values
            
        except Exception as e:
            logger.warning(f"Error extracting signal from {bigwig_file}: {e}")
            return np.zeros(self.window_size)
    
    def extract_bidirectional_signals(self, 
                                    forward_bigwig: str, 
                                    reverse_bigwig: str, 
                                    region: EnhancerRegion) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract forward and reverse strand signals.
        
        Args:
            forward_bigwig: Forward strand BigWig file
            reverse_bigwig: Reverse strand BigWig file
            region: Enhancer region
            
        Returns:
            Tuple of (forward_signal, reverse_signal)
        """
        forward_signal = self.extract_signal_from_bigwig(forward_bigwig, region, "+")
        reverse_signal = self.extract_signal_from_bigwig(reverse_bigwig, region, "-")
        
        # Reverse the reverse signal for consistent orientation
        reverse_signal = np.flip(reverse_signal)
        
        return forward_signal, reverse_signal
    
    def smooth_signal(self, signal: np.ndarray) -> np.ndarray:
        """Apply Gaussian smoothing to signal."""
        if self.smoothing_sigma > 0:
            return gaussian_filter1d(signal, sigma=self.smoothing_sigma)
        return signal
    
    def detect_peaks(self, signal: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Detect peaks in GRO-seq signal.
        
        Args:
            signal: Input signal array
            
        Returns:
            Dictionary with peak information
        """
        # Smooth signal
        smoothed_signal = self.smooth_signal(signal)
        
        # Find peaks
        peaks, properties = find_peaks(
            smoothed_signal, 
            height=self.min_peak_height,
            distance=10  # Minimum distance between peaks
        )
        
        # Calculate peak properties
        if len(peaks) > 0:
            # Peak widths at half maximum
            widths = peak_widths(smoothed_signal, peaks, rel_height=0.5)[0]
            
            # Peak prominences
            prominences = peak_prominences(smoothed_signal, peaks)[0]
            
            # Peak heights
            heights = smoothed_signal[peaks]
        else:
            widths = np.array([])
            prominences = np.array([])
            heights = np.array([])
        
        return {
            'peak_positions': peaks,
            'peak_heights': heights,
            'peak_widths': widths,
            'peak_prominences': prominences
        }
    
    def calculate_signal_entropy(self, signal: np.ndarray) -> float:
        """Calculate signal entropy as a measure of complexity."""
        # Normalize signal to probability distribution
        signal_positive = signal - np.min(signal) + 1e-10
        prob_dist = signal_positive / np.sum(signal_positive)
        
        # Calculate entropy
        entropy = -np.sum(prob_dist * np.log2(prob_dist + 1e-10))
        
        return entropy
    
    def calculate_center_bias(self, signal: np.ndarray) -> float:
        """Calculate bias towards center of the region."""
        center_idx = len(signal) // 2
        window_quarter = len(signal) // 4
        
        # Signal in center quarter vs. edges
        center_signal = np.mean(signal[center_idx - window_quarter:center_idx + window_quarter])
        edge_signal = np.mean(np.concatenate([
            signal[:window_quarter], 
            signal[-window_quarter:]
        ]))
        
        if edge_signal > 0:
            return center_signal / edge_signal
        else:
            return 1.0
    
    def calculate_snr(self, signal: np.ndarray) -> float:
        """Calculate signal-to-noise ratio."""
        if len(signal) == 0:
            return 0.0
        
        # Use robust estimates
        signal_power = np.percentile(signal, 90)
        noise_power = np.std(signal[signal < np.percentile(signal, 25)])
        
        if noise_power > 0:
            return signal_power / noise_power
        else:
            return float('inf') if signal_power > 0 else 0.0
    
    def extract_features(self, 
                        forward_signal: np.ndarray, 
                        reverse_signal: np.ndarray, 
                        region_id: str) -> GROSeqFeatures:
        """
        Extract comprehensive features from GRO-seq signals.
        
        Args:
            forward_signal: Forward strand signal
            reverse_signal: Reverse strand signal
            region_id: Identifier for the region
            
        Returns:
            GROSeqFeatures object
        """
        # Combine signals for overall statistics
        combined_signal = forward_signal + reverse_signal
        
        # Basic intensity features
        mean_signal = np.mean(combined_signal)
        max_signal = np.max(combined_signal)
        total_signal = np.sum(combined_signal)
        signal_variance = np.var(combined_signal)
        
        # Bidirectional features
        forward_total = np.sum(forward_signal)
        reverse_total = np.sum(reverse_signal)
        total_both = forward_total + reverse_total
        
        if total_both > 0:
            bidirectional_ratio = min(forward_total, reverse_total) / total_both
        else:
            bidirectional_ratio = 0.0
        
        # Peak detection on combined signal
        peak_info = self.detect_peaks(combined_signal)
        n_peaks = len(peak_info['peak_positions'])
        
        # Shape features
        if len(combined_signal) > 0 and np.std(combined_signal) > 0:
            skewness = stats.skew(combined_signal)
            kurtosis = stats.kurtosis(combined_signal)
        else:
            skewness = 0.0
            kurtosis = 0.0
        
        # Information content
        signal_entropy = self.calculate_signal_entropy(combined_signal)
        
        # Spatial features
        signal_spread = np.std(np.arange(len(combined_signal)) * combined_signal) if total_signal > 0 else 0
        center_bias = self.calculate_center_bias(combined_signal)
        
        # Quality metrics
        snr = self.calculate_snr(combined_signal)
        coverage_fraction = np.sum(combined_signal > 0) / len(combined_signal)
        
        return GROSeqFeatures(
            region_id=region_id,
            mean_signal=mean_signal,
            max_signal=max_signal,
            total_signal=total_signal,
            signal_variance=signal_variance,
            forward_signal=forward_total,
            reverse_signal=reverse_total,
            bidirectional_ratio=bidirectional_ratio,
            n_peaks=n_peaks,
            peak_heights=peak_info['peak_heights'].tolist(),
            peak_widths=peak_info['peak_widths'].tolist(),
            peak_prominences=peak_info['peak_prominences'].tolist(),
            skewness=skewness,
            kurtosis=kurtosis,
            signal_entropy=signal_entropy,
            signal_spread=signal_spread,
            center_bias=center_bias,
            snr=snr,
            coverage_fraction=coverage_fraction
        )


class TemporalFeatureExtractor:
    """Extract temporal dynamics features from time-series GRO-seq data."""
    
    def __init__(self):
        self.time_points = None
    
    def extract_temporal_features(self, 
                                signals_time_series: List[np.ndarray], 
                                time_points: List[float],
                                region_id: str) -> Dict[str, float]:
        """
        Extract temporal dynamics features.
        
        Args:
            signals_time_series: List of signal arrays at different time points
            time_points: Time points corresponding to signals
            region_id: Region identifier
            
        Returns:
            Dictionary of temporal features
        """
        if len(signals_time_series) < 2:
            return {'region_id': region_id}
        
        # Calculate signal intensity over time
        intensities = [np.sum(signal) for signal in signals_time_series]
        
        # Temporal features
        features = {
            'region_id': region_id,
            'max_intensity': max(intensities),
            'min_intensity': min(intensities),
            'intensity_range': max(intensities) - min(intensities),
            'peak_time': time_points[np.argmax(intensities)],
            'intensity_trend': self._calculate_trend(time_points, intensities),
            'intensity_variability': np.std(intensities),
            'response_delay': self._calculate_response_delay(time_points, intensities),
            'decay_rate': self._calculate_decay_rate(time_points, intensities)
        }
        
        return features
    
    def _calculate_trend(self, time_points: List[float], intensities: List[float]) -> float:
        """Calculate overall trend in intensity over time."""
        if len(time_points) < 2:
            return 0.0
        
        slope, _, _, _, _ = stats.linregress(time_points, intensities)
        return slope
    
    def _calculate_response_delay(self, time_points: List[float], intensities: List[float]) -> float:
        """Calculate delay to reach maximum response."""
        baseline = intensities[0]
        max_intensity = max(intensities)
        
        if max_intensity <= baseline:
            return float('inf')
        
        # Find first time point where signal exceeds 50% of max response
        threshold = baseline + 0.5 * (max_intensity - baseline)
        
        for i, intensity in enumerate(intensities):
            if intensity >= threshold:
                return time_points[i] - time_points[0]
        
        return float('inf')
    
    def _calculate_decay_rate(self, time_points: List[float], intensities: List[float]) -> float:
        """Calculate signal decay rate after peak."""
        peak_idx = np.argmax(intensities)
        
        if peak_idx >= len(intensities) - 1:
            return 0.0
        
        # Fit exponential decay to post-peak data
        post_peak_times = time_points[peak_idx:]
        post_peak_intensities = intensities[peak_idx:]
        
        if len(post_peak_times) < 2:
            return 0.0
        
        # Log-linear fit for exponential decay
        log_intensities = np.log(np.maximum(post_peak_intensities, 1e-10))
        slope, _, _, _, _ = stats.linregress(post_peak_times, log_intensities)
        
        return -slope  # Return positive decay rate


class EnhancerActivityQuantifier:
    """Quantify enhancer activity from GRO-seq features."""
    
    def __init__(self, 
                 activity_weights: Optional[Dict[str, float]] = None):
        """
        Initialize enhancer activity quantifier.
        
        Args:
            activity_weights: Weights for different features in activity calculation
        """
        if activity_weights is None:
            self.activity_weights = {
                'total_signal': 0.3,
                'bidirectional_ratio': 0.25,
                'n_peaks': 0.15,
                'snr': 0.15,
                'signal_entropy': 0.15
            }
        else:
            self.activity_weights = activity_weights
    
    def calculate_activity_score(self, features: GROSeqFeatures) -> float:
        """
        Calculate enhancer activity score from features.
        
        Args:
            features: GROSeqFeatures object
            
        Returns:
            Activity score (0-1 range)
        """
        # Normalize individual components
        normalized_features = {}
        
        # Total signal (log-scaled)
        normalized_features['total_signal'] = np.tanh(np.log1p(features.total_signal) / 10)
        
        # Bidirectional ratio (already 0-1)
        normalized_features['bidirectional_ratio'] = features.bidirectional_ratio
        
        # Number of peaks (scaled)
        normalized_features['n_peaks'] = np.tanh(features.n_peaks / 5.0)
        
        # Signal-to-noise ratio (log-scaled)
        normalized_features['snr'] = np.tanh(np.log1p(features.snr) / 5)
        
        # Signal entropy (scaled)
        normalized_features['signal_entropy'] = np.tanh(features.signal_entropy / 10)
        
        # Calculate weighted sum
        activity_score = 0.0
        for feature_name, weight in self.activity_weights.items():
            if feature_name in normalized_features:
                activity_score += weight * normalized_features[feature_name]
        
        return np.clip(activity_score, 0.0, 1.0)
    
    def classify_enhancer_activity(self, activity_score: float) -> str:
        """
        Classify enhancer activity level.
        
        Args:
            activity_score: Activity score from calculate_activity_score
            
        Returns:
            Activity classification string
        """
        if activity_score >= 0.8:
            return "high"
        elif activity_score >= 0.5:
            return "medium"
        elif activity_score >= 0.2:
            return "low"
        else:
            return "inactive"


class MultiConditionComparator:
    """Compare GRO-seq features across multiple conditions."""
    
    def __init__(self):
        pass
    
    def compare_conditions(self, 
                         features_by_condition: Dict[str, List[GROSeqFeatures]]) -> pd.DataFrame:
        """
        Compare GRO-seq features across conditions.
        
        Args:
            features_by_condition: Dictionary mapping condition names to feature lists
            
        Returns:
            DataFrame with comparison statistics
        """
        conditions = list(features_by_condition.keys())
        
        if len(conditions) < 2:
            raise ValueError("Need at least 2 conditions for comparison")
        
        # Collect all region IDs
        all_regions = set()
        for features_list in features_by_condition.values():
            all_regions.update([f.region_id for f in features_list])
        
        comparison_results = []
        
        for region_id in all_regions:
            region_data = {}
            
            # Extract features for this region across conditions
            for condition in conditions:
                features_list = features_by_condition[condition]
                region_features = next((f for f in features_list if f.region_id == region_id), None)
                
                if region_features:
                    region_data[condition] = {
                        'total_signal': region_features.total_signal,
                        'bidirectional_ratio': region_features.bidirectional_ratio,
                        'n_peaks': region_features.n_peaks,
                        'snr': region_features.snr
                    }
            
            # Skip if not present in all conditions
            if len(region_data) != len(conditions):
                continue
            
            # Calculate differences and statistics
            result = {'region_id': region_id}
            
            # Pairwise comparisons
            for i, cond1 in enumerate(conditions):
                for j, cond2 in enumerate(conditions[i+1:], i+1):
                    for feature in ['total_signal', 'bidirectional_ratio', 'n_peaks', 'snr']:
                        val1 = region_data[cond1][feature]
                        val2 = region_data[cond2][feature]
                        
                        # Fold change
                        if val2 > 0:
                            fold_change = val1 / val2
                        else:
                            fold_change = float('inf') if val1 > 0 else 1.0
                        
                        # Log2 fold change
                        log2_fc = np.log2(fold_change) if fold_change > 0 and fold_change != float('inf') else 0
                        
                        result[f"{feature}_{cond1}_vs_{cond2}_log2fc"] = log2_fc
            
            comparison_results.append(result)
        
        return pd.DataFrame(comparison_results)


def load_enhancer_regions(bed_file: str) -> List[EnhancerRegion]:
    """Load enhancer regions from BED file."""
    regions = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
                strand = fields[5] if len(fields) > 5 else "."
                
                regions.append(EnhancerRegion(chrom, start, end, name, strand))
    
    return regions


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee GRO-seq Feature Engineering for Enhancer Activity"
    )
    parser.add_argument(
        '--enhancers', '-e', required=True,
        help="BED file with enhancer regions"
    )
    parser.add_argument(
        '--forward-bigwig', required=True,
        help="BigWig file with forward strand GRO-seq signal"
    )
    parser.add_argument(
        '--reverse-bigwig', required=True,
        help="BigWig file with reverse strand GRO-seq signal"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output CSV file for features"
    )
    parser.add_argument(
        '--window-size', type=int, default=2000,
        help="Window size around enhancer center (default: 2000)"
    )
    parser.add_argument(
        '--min-peak-height', type=float, default=0.1,
        help="Minimum peak height for detection (default: 0.1)"
    )
    parser.add_argument(
        '--smoothing-sigma', type=float, default=2.0,
        help="Gaussian smoothing sigma (default: 2.0)"
    )
    parser.add_argument(
        '--calculate-activity', action='store_true',
        help="Calculate enhancer activity scores"
    )
    
    args = parser.parse_args()
    
    # Load enhancer regions
    logger.info(f"Loading enhancer regions from {args.enhancers}")
    enhancer_regions = load_enhancer_regions(args.enhancers)
    logger.info(f"Loaded {len(enhancer_regions)} enhancer regions")
    
    # Initialize signal processor
    processor = GROSeqSignalProcessor(
        window_size=args.window_size,
        min_peak_height=args.min_peak_height,
        smoothing_sigma=args.smoothing_sigma
    )
    
    # Initialize activity quantifier if requested
    if args.calculate_activity:
        activity_quantifier = EnhancerActivityQuantifier()
    
    # Extract features for each enhancer
    logger.info("Extracting GRO-seq features...")
    all_features = []
    
    for i, region in enumerate(enhancer_regions):
        if i % 100 == 0:
            logger.info(f"Processing region {i+1}/{len(enhancer_regions)}")
        
        try:
            # Extract bidirectional signals
            forward_signal, reverse_signal = processor.extract_bidirectional_signals(
                args.forward_bigwig, args.reverse_bigwig, region
            )
            
            # Extract features
            features = processor.extract_features(forward_signal, reverse_signal, region.name)
            all_features.append(features)
            
        except Exception as e:
            logger.warning(f"Error processing region {region.name}: {e}")
            continue
    
    # Convert to DataFrame
    logger.info("Converting features to DataFrame...")
    feature_data = []
    
    for features in all_features:
        row = {
            'region_id': features.region_id,
            'mean_signal': features.mean_signal,
            'max_signal': features.max_signal,
            'total_signal': features.total_signal,
            'signal_variance': features.signal_variance,
            'forward_signal': features.forward_signal,
            'reverse_signal': features.reverse_signal,
            'bidirectional_ratio': features.bidirectional_ratio,
            'n_peaks': features.n_peaks,
            'mean_peak_height': np.mean(features.peak_heights) if features.peak_heights else 0,
            'mean_peak_width': np.mean(features.peak_widths) if features.peak_widths else 0,
            'mean_peak_prominence': np.mean(features.peak_prominences) if features.peak_prominences else 0,
            'skewness': features.skewness,
            'kurtosis': features.kurtosis,
            'signal_entropy': features.signal_entropy,
            'signal_spread': features.signal_spread,
            'center_bias': features.center_bias,
            'snr': features.snr,
            'coverage_fraction': features.coverage_fraction
        }
        
        # Add activity score if calculated
        if args.calculate_activity:
            activity_score = activity_quantifier.calculate_activity_score(features)
            activity_class = activity_quantifier.classify_enhancer_activity(activity_score)
            row['activity_score'] = activity_score
            row['activity_class'] = activity_class
        
        feature_data.append(row)
    
    # Create DataFrame and save
    features_df = pd.DataFrame(feature_data)
    features_df.to_csv(args.output, index=False)
    
    logger.info(f"Features saved to {args.output}")
    logger.info(f"Extracted features for {len(features_df)} regions")
    
    # Print summary statistics
    if args.calculate_activity:
        activity_counts = features_df['activity_class'].value_counts()
        logger.info("Activity classification summary:")
        for activity_class, count in activity_counts.items():
            logger.info(f"  {activity_class}: {count} regions")
    
    # Print feature statistics
    numeric_cols = features_df.select_dtypes(include=[np.number]).columns
    logger.info("Feature statistics (mean ± std):")
    for col in numeric_cols[:10]:  # Show first 10 features
        mean_val = features_df[col].mean()
        std_val = features_df[col].std()
        logger.info(f"  {col}: {mean_val:.3f} ± {std_val:.3f}")


if __name__ == "__main__":
    main()