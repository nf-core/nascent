#!/usr/bin/env python
"""
TFSee TF-Enhancer Association Scoring Module

This module implements algorithms for scoring associations between transcription
factors and enhancer regions based on multiple evidence sources including:
- Motif presence and strength
- Chromatin accessibility
- Gene expression correlation
- Distance-based proximity scoring
- ChIP-seq peak overlap (when available)

Key Features:
- Multi-evidence integration scoring
- Distance-weighted association scoring
- Correlation-based TF-target gene relationships
- Robust statistical frameworks for association testing
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr, spearmanr
import warnings
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
import pickle

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class GenomicRegion:
    """Represents a genomic region with coordinates and metadata."""
    chromosome: str
    start: int
    end: int
    name: str = ""
    score: float = 0.0
    strand: str = "."
    
    def __post_init__(self):
        if self.end <= self.start:
            raise ValueError(f"End position {self.end} must be greater than start {self.start}")
    
    @property
    def length(self) -> int:
        return self.end - self.start
    
    @property
    def center(self) -> int:
        return (self.start + self.end) // 2
    
    def distance_to(self, other: 'GenomicRegion') -> int:
        """Calculate distance to another genomic region."""
        if self.chromosome != other.chromosome:
            return float('inf')
        
        # If regions overlap, distance is 0
        if self.end > other.start and self.start < other.end:
            return 0
        
        # Calculate minimum distance between regions
        return min(abs(self.start - other.end), abs(self.end - other.start))
    
    def overlaps_with(self, other: 'GenomicRegion', min_overlap: int = 1) -> bool:
        """Check if this region overlaps with another."""
        if self.chromosome != other.chromosome:
            return False
        
        overlap = min(self.end, other.end) - max(self.start, other.start)
        return overlap >= min_overlap


class DistanceScorer:
    """Distance-based scoring for TF-enhancer associations."""
    
    def __init__(self, max_distance: int = 1000000, decay_function: str = "exponential"):
        """
        Initialize distance scorer.
        
        Args:
            max_distance: Maximum distance to consider for associations
            decay_function: Type of distance decay ('exponential', 'linear', 'power')
        """
        self.max_distance = max_distance
        self.decay_function = decay_function
    
    def calculate_distance_weight(self, distance: float) -> float:
        """
        Calculate distance weight based on decay function.
        
        Args:
            distance: Distance between TF and enhancer
            
        Returns:
            Weight between 0 and 1
        """
        if distance > self.max_distance:
            return 0.0
        
        normalized_distance = distance / self.max_distance
        
        if self.decay_function == "exponential":
            # Exponential decay: weight = exp(-λ * distance)
            lambda_param = 3.0  # Adjust for steepness
            return np.exp(-lambda_param * normalized_distance)
        
        elif self.decay_function == "linear":
            # Linear decay: weight = 1 - distance/max_distance
            return 1.0 - normalized_distance
        
        elif self.decay_function == "power":
            # Power decay: weight = (1 - distance/max_distance)^α
            alpha = 2.0  # Power parameter
            return (1.0 - normalized_distance) ** alpha
        
        else:
            raise ValueError(f"Unknown decay function: {self.decay_function}")
    
    def score_tf_enhancer_distance(self, 
                                 tf_regions: List[GenomicRegion], 
                                 enhancer_regions: List[GenomicRegion]) -> np.ndarray:
        """
        Score TF-enhancer associations based on distance.
        
        Args:
            tf_regions: List of TF binding regions
            enhancer_regions: List of enhancer regions
            
        Returns:
            Matrix of shape (n_tfs, n_enhancers) with distance scores
        """
        n_tfs = len(tf_regions)
        n_enhancers = len(enhancer_regions)
        scores = np.zeros((n_tfs, n_enhancers))
        
        for i, tf_region in enumerate(tf_regions):
            for j, enhancer_region in enumerate(enhancer_regions):
                distance = tf_region.distance_to(enhancer_region)
                scores[i, j] = self.calculate_distance_weight(distance)
        
        return scores


class MotifScorer:
    """Motif-based scoring for TF-enhancer associations."""
    
    def __init__(self, motif_scores: Dict[str, Dict[str, float]], 
                 score_threshold: float = 0.7):
        """
        Initialize motif scorer.
        
        Args:
            motif_scores: Nested dict {tf_name: {enhancer_id: score}}
            score_threshold: Minimum score to consider significant
        """
        self.motif_scores = motif_scores
        self.score_threshold = score_threshold
    
    def normalize_scores(self, scores: np.ndarray) -> np.ndarray:
        """Normalize motif scores to 0-1 range."""
        min_score = np.min(scores)
        max_score = np.max(scores)
        
        if max_score == min_score:
            return np.ones_like(scores) * 0.5
        
        return (scores - min_score) / (max_score - min_score)
    
    def score_tf_enhancer_motifs(self, 
                               tf_names: List[str], 
                               enhancer_ids: List[str]) -> np.ndarray:
        """
        Score TF-enhancer associations based on motif presence.
        
        Args:
            tf_names: List of TF names
            enhancer_ids: List of enhancer IDs
            
        Returns:
            Matrix of shape (n_tfs, n_enhancers) with motif scores
        """
        n_tfs = len(tf_names)
        n_enhancers = len(enhancer_ids)
        scores = np.zeros((n_tfs, n_enhancers))
        
        for i, tf_name in enumerate(tf_names):
            if tf_name in self.motif_scores:
                for j, enhancer_id in enumerate(enhancer_ids):
                    if enhancer_id in self.motif_scores[tf_name]:
                        raw_score = self.motif_scores[tf_name][enhancer_id]
                        # Apply threshold
                        scores[i, j] = raw_score if raw_score >= self.score_threshold else 0.0
        
        # Normalize scores for each TF
        for i in range(n_tfs):
            if np.max(scores[i, :]) > 0:
                scores[i, :] = self.normalize_scores(scores[i, :])
        
        return scores


class ExpressionCorrelationScorer:
    """Expression correlation-based scoring for TF-target relationships."""
    
    def __init__(self, expression_data: pd.DataFrame, correlation_method: str = "pearson"):
        """
        Initialize expression correlation scorer.
        
        Args:
            expression_data: DataFrame with genes as rows, samples as columns
            correlation_method: 'pearson' or 'spearman'
        """
        self.expression_data = expression_data
        self.correlation_method = correlation_method
    
    def calculate_tf_target_correlations(self, 
                                       tf_genes: List[str], 
                                       target_genes: List[str]) -> np.ndarray:
        """
        Calculate correlations between TF and target gene expression.
        
        Args:
            tf_genes: List of TF gene names
            target_genes: List of target gene names
            
        Returns:
            Correlation matrix of shape (n_tfs, n_targets)
        """
        n_tfs = len(tf_genes)
        n_targets = len(target_genes)
        correlations = np.zeros((n_tfs, n_targets))
        
        for i, tf_gene in enumerate(tf_genes):
            if tf_gene not in self.expression_data.index:
                logger.warning(f"TF gene {tf_gene} not found in expression data")
                continue
                
            tf_expression = self.expression_data.loc[tf_gene].values
            
            for j, target_gene in enumerate(target_genes):
                if target_gene not in self.expression_data.index:
                    continue
                
                target_expression = self.expression_data.loc[target_gene].values
                
                # Calculate correlation
                if self.correlation_method == "pearson":
                    corr, _ = pearsonr(tf_expression, target_expression)
                elif self.correlation_method == "spearman":
                    corr, _ = spearmanr(tf_expression, target_expression)
                else:
                    raise ValueError(f"Unknown correlation method: {self.correlation_method}")
                
                correlations[i, j] = corr if not np.isnan(corr) else 0.0
        
        return correlations
    
    def score_tf_enhancer_via_targets(self, 
                                    tf_genes: List[str], 
                                    enhancer_target_map: Dict[str, List[str]]) -> Dict[str, np.ndarray]:
        """
        Score TF-enhancer associations via target gene correlations.
        
        Args:
            tf_genes: List of TF gene names
            enhancer_target_map: Map from enhancer IDs to target gene lists
            
        Returns:
            Dictionary mapping enhancer IDs to correlation scores
        """
        scores = {}
        
        for enhancer_id, target_genes in enhancer_target_map.items():
            if not target_genes:
                continue
            
            correlations = self.calculate_tf_target_correlations(tf_genes, target_genes)
            
            # Aggregate correlations across targets (mean absolute correlation)
            enhancer_scores = np.mean(np.abs(correlations), axis=1)
            scores[enhancer_id] = enhancer_scores
        
        return scores


class ChIPSeqOverlapScorer:
    """ChIP-seq peak overlap scoring for TF-enhancer associations."""
    
    def __init__(self, chip_peaks: Dict[str, List[GenomicRegion]], 
                 min_overlap: int = 1):
        """
        Initialize ChIP-seq overlap scorer.
        
        Args:
            chip_peaks: Dictionary mapping TF names to peak regions
            min_overlap: Minimum overlap in base pairs
        """
        self.chip_peaks = chip_peaks
        self.min_overlap = min_overlap
    
    def score_tf_enhancer_overlap(self, 
                                tf_names: List[str], 
                                enhancer_regions: List[GenomicRegion]) -> np.ndarray:
        """
        Score TF-enhancer associations based on ChIP-seq peak overlap.
        
        Args:
            tf_names: List of TF names
            enhancer_regions: List of enhancer regions
            
        Returns:
            Binary overlap matrix of shape (n_tfs, n_enhancers)
        """
        n_tfs = len(tf_names)
        n_enhancers = len(enhancer_regions)
        overlaps = np.zeros((n_tfs, n_enhancers))
        
        for i, tf_name in enumerate(tf_names):
            if tf_name not in self.chip_peaks:
                continue
            
            peaks = self.chip_peaks[tf_name]
            
            for j, enhancer in enumerate(enhancer_regions):
                # Check if any peak overlaps with enhancer
                for peak in peaks:
                    if enhancer.overlaps_with(peak, self.min_overlap):
                        overlaps[i, j] = 1.0
                        break
        
        return overlaps


class IntegratedTFEnhancerScorer:
    """Integrated scoring combining multiple evidence sources."""
    
    def __init__(self, 
                 distance_weight: float = 0.3,
                 motif_weight: float = 0.4, 
                 expression_weight: float = 0.2,
                 chip_weight: float = 0.1):
        """
        Initialize integrated scorer with evidence weights.
        
        Args:
            distance_weight: Weight for distance evidence
            motif_weight: Weight for motif evidence
            expression_weight: Weight for expression correlation evidence
            chip_weight: Weight for ChIP-seq evidence
        """
        self.weights = {
            'distance': distance_weight,
            'motif': motif_weight,
            'expression': expression_weight,
            'chip': chip_weight
        }
        
        # Normalize weights
        total_weight = sum(self.weights.values())
        self.weights = {k: v/total_weight for k, v in self.weights.items()}
    
    def integrate_scores(self, 
                        score_matrices: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Integrate multiple score matrices using weighted combination.
        
        Args:
            score_matrices: Dictionary mapping evidence types to score matrices
            
        Returns:
            Integrated score matrix
        """
        # Check matrix dimensions
        shapes = [matrix.shape for matrix in score_matrices.values()]
        if not all(shape == shapes[0] for shape in shapes):
            raise ValueError("All score matrices must have the same shape")
        
        integrated_scores = np.zeros(shapes[0])
        total_weight = 0.0
        
        for evidence_type, matrix in score_matrices.items():
            if evidence_type in self.weights:
                weight = self.weights[evidence_type]
                integrated_scores += weight * matrix
                total_weight += weight
        
        # Normalize by total weight used
        if total_weight > 0:
            integrated_scores /= total_weight
        
        return integrated_scores
    
    def calculate_confidence_scores(self, 
                                  score_matrices: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Calculate confidence scores based on evidence agreement.
        
        Args:
            score_matrices: Dictionary mapping evidence types to score matrices
            
        Returns:
            Confidence score matrix
        """
        if len(score_matrices) < 2:
            return np.ones(list(score_matrices.values())[0].shape)
        
        # Calculate pairwise correlations between evidence types
        evidence_stack = np.stack(list(score_matrices.values()), axis=2)
        n_evidence = evidence_stack.shape[2]
        
        confidence_scores = np.zeros(evidence_stack.shape[:2])
        
        for i in range(evidence_stack.shape[0]):
            for j in range(evidence_stack.shape[1]):
                values = evidence_stack[i, j, :]
                
                # Calculate mean pairwise correlation
                correlations = []
                for k in range(n_evidence):
                    for l in range(k+1, n_evidence):
                        if np.std(values[[k, l]]) > 0:
                            corr = np.corrcoef(values[[k, l]])[0, 1]
                            if not np.isnan(corr):
                                correlations.append(abs(corr))
                
                confidence_scores[i, j] = np.mean(correlations) if correlations else 0.0
        
        return confidence_scores


def rank_tf_enhancer_associations(score_matrix: np.ndarray, 
                                tf_names: List[str], 
                                enhancer_ids: List[str],
                                top_k: int = 100) -> pd.DataFrame:
    """
    Rank TF-enhancer associations by score.
    
    Args:
        score_matrix: Matrix of association scores
        tf_names: List of TF names
        enhancer_ids: List of enhancer IDs
        top_k: Number of top associations to return
        
    Returns:
        DataFrame with ranked associations
    """
    # Flatten matrix and get indices
    flat_scores = score_matrix.flatten()
    tf_indices, enhancer_indices = np.unravel_index(
        np.argsort(flat_scores)[::-1], score_matrix.shape
    )
    
    # Create results dataframe
    results = []
    for i in range(min(top_k, len(flat_scores))):
        tf_idx = tf_indices[i]
        enh_idx = enhancer_indices[i]
        score = score_matrix[tf_idx, enh_idx]
        
        if score > 0:  # Only include non-zero scores
            results.append({
                'tf_name': tf_names[tf_idx],
                'enhancer_id': enhancer_ids[enh_idx],
                'association_score': score,
                'rank': i + 1
            })
    
    return pd.DataFrame(results)


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee TF-Enhancer Association Scoring"
    )
    parser.add_argument(
        '--tf-regions', required=True,
        help="BED file with TF binding regions"
    )
    parser.add_argument(
        '--enhancer-regions', required=True,
        help="BED file with enhancer regions"
    )
    parser.add_argument(
        '--motif-scores',
        help="Pickle file with motif scores dictionary"
    )
    parser.add_argument(
        '--expression-data',
        help="CSV file with gene expression data"
    )
    parser.add_argument(
        '--chip-peaks',
        help="Pickle file with ChIP-seq peaks dictionary"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output file for association scores"
    )
    parser.add_argument(
        '--max-distance', type=int, default=1000000,
        help="Maximum distance for TF-enhancer associations (default: 1000000)"
    )
    parser.add_argument(
        '--top-k', type=int, default=1000,
        help="Number of top associations to output (default: 1000)"
    )
    
    args = parser.parse_args()
    
    # Read genomic regions
    def read_bed_file(filename):
        regions = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
                    score = float(fields[4]) if len(fields) > 4 else 0.0
                    strand = fields[5] if len(fields) > 5 else "."
                    
                    regions.append(GenomicRegion(chrom, start, end, name, score, strand))
        return regions
    
    logger.info("Reading TF and enhancer regions...")
    tf_regions = read_bed_file(args.tf_regions)
    enhancer_regions = read_bed_file(args.enhancer_regions)
    
    tf_names = [region.name for region in tf_regions]
    enhancer_ids = [region.name for region in enhancer_regions]
    
    logger.info(f"Loaded {len(tf_regions)} TF regions and {len(enhancer_regions)} enhancer regions")
    
    # Initialize scorers and calculate scores
    score_matrices = {}
    
    # Distance scoring
    logger.info("Calculating distance scores...")
    distance_scorer = DistanceScorer(max_distance=args.max_distance)
    score_matrices['distance'] = distance_scorer.score_tf_enhancer_distance(
        tf_regions, enhancer_regions
    )
    
    # Motif scoring
    if args.motif_scores:
        logger.info("Loading motif scores...")
        with open(args.motif_scores, 'rb') as f:
            motif_scores = pickle.load(f)
        
        motif_scorer = MotifScorer(motif_scores)
        score_matrices['motif'] = motif_scorer.score_tf_enhancer_motifs(
            tf_names, enhancer_ids
        )
    
    # Expression correlation scoring
    if args.expression_data:
        logger.info("Loading expression data...")
        expression_data = pd.read_csv(args.expression_data, index_col=0)
        
        # For simplicity, assume enhancer targets are nearby genes
        # In practice, this would come from enhancer-target mapping
        enhancer_target_map = {enh_id: [enh_id.split('_')[0]] for enh_id in enhancer_ids}
        
        expr_scorer = ExpressionCorrelationScorer(expression_data)
        expr_scores = expr_scorer.score_tf_enhancer_via_targets(tf_names, enhancer_target_map)
        
        # Convert to matrix format
        expr_matrix = np.zeros((len(tf_names), len(enhancer_ids)))
        for j, enh_id in enumerate(enhancer_ids):
            if enh_id in expr_scores:
                expr_matrix[:, j] = expr_scores[enh_id]
        
        score_matrices['expression'] = expr_matrix
    
    # ChIP-seq overlap scoring
    if args.chip_peaks:
        logger.info("Loading ChIP-seq peaks...")
        with open(args.chip_peaks, 'rb') as f:
            chip_peaks = pickle.load(f)
        
        chip_scorer = ChIPSeqOverlapScorer(chip_peaks)
        score_matrices['chip'] = chip_scorer.score_tf_enhancer_overlap(
            tf_names, enhancer_regions
        )
    
    # Integrate scores
    logger.info("Integrating evidence sources...")
    integrated_scorer = IntegratedTFEnhancerScorer()
    integrated_scores = integrated_scorer.integrate_scores(score_matrices)
    confidence_scores = integrated_scorer.calculate_confidence_scores(score_matrices)
    
    # Rank associations
    logger.info("Ranking associations...")
    results_df = rank_tf_enhancer_associations(
        integrated_scores, tf_names, enhancer_ids, top_k=args.top_k
    )
    
    # Add confidence scores
    if len(results_df) > 0:
        confidence_values = []
        for _, row in results_df.iterrows():
            tf_idx = tf_names.index(row['tf_name'])
            enh_idx = enhancer_ids.index(row['enhancer_id'])
            confidence_values.append(confidence_scores[tf_idx, enh_idx])
        
        results_df['confidence_score'] = confidence_values
    
    # Save results
    results_df.to_csv(args.output, index=False)
    logger.info(f"Results saved to {args.output}")
    
    # Print summary
    logger.info(f"Generated {len(results_df)} TF-enhancer associations")
    if len(results_df) > 0:
        logger.info(f"Top association: {results_df.iloc[0]['tf_name']} -> {results_df.iloc[0]['enhancer_id']} "
                   f"(score: {results_df.iloc[0]['association_score']:.3f})")


if __name__ == "__main__":
    main()