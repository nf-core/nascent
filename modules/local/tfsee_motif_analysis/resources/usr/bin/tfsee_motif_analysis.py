#!/usr/bin/env python
"""
TFSee Motif Analysis Module

This module implements motif probability calculation and integration methods
for transcription factor binding site analysis, including the Stouffer method
for combining multiple p-values.

Key Features:
- Position Weight Matrix (PWM) scoring
- Stouffer's method for p-value combination
- Motif enrichment analysis
- Background sequence modeling
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm
import warnings
from typing import List, Dict, Tuple, Optional, Union
import re
from collections import defaultdict

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


class PWMScorer:
    """Position Weight Matrix scorer for motif analysis."""
    
    def __init__(self, pseudocount: float = 0.001):
        """
        Initialize PWM scorer.
        
        Args:
            pseudocount: Small value to avoid log(0) in PWM calculations
        """
        self.pseudocount = pseudocount
        self.nucleotides = ['A', 'C', 'G', 'T']
        
    def build_pwm_from_counts(self, count_matrix: np.ndarray) -> np.ndarray:
        """
        Build Position Weight Matrix from count matrix.
        
        Args:
            count_matrix: Matrix of shape (4, motif_length) with nucleotide counts
            
        Returns:
            PWM matrix of shape (4, motif_length)
        """
        # Add pseudocounts
        count_matrix = count_matrix + self.pseudocount
        
        # Convert to frequencies
        freq_matrix = count_matrix / np.sum(count_matrix, axis=0)
        
        # Background frequencies (uniform for simplicity, can be genome-specific)
        background = np.array([0.25, 0.25, 0.25, 0.25])
        
        # Calculate log-odds ratios
        pwm = np.log2(freq_matrix / background.reshape(-1, 1))
        
        return pwm
    
    def score_sequence(self, sequence: str, pwm: np.ndarray) -> float:
        """
        Score a sequence using PWM.
        
        Args:
            sequence: DNA sequence string
            pwm: Position Weight Matrix
            
        Returns:
            PWM score for the sequence
        """
        if len(sequence) != pwm.shape[1]:
            raise ValueError(f"Sequence length {len(sequence)} doesn't match PWM length {pwm.shape[1]}")
        
        score = 0.0
        nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        for i, nucleotide in enumerate(sequence.upper()):
            if nucleotide in nuc_to_idx:
                score += pwm[nuc_to_idx[nucleotide], i]
            else:
                # Handle ambiguous nucleotides by using average score
                score += np.mean(pwm[:, i])
        
        return score
    
    def scan_sequence(self, sequence: str, pwm: np.ndarray) -> List[Tuple[int, float]]:
        """
        Scan sequence for motif matches using sliding window.
        
        Args:
            sequence: DNA sequence to scan
            pwm: Position Weight Matrix
            
        Returns:
            List of (position, score) tuples
        """
        motif_length = pwm.shape[1]
        matches = []
        
        for i in range(len(sequence) - motif_length + 1):
            subseq = sequence[i:i + motif_length]
            score = self.score_sequence(subseq, pwm)
            matches.append((i, score))
        
        return matches


class StoufferCombination:
    """Implementation of Stouffer's method for combining p-values."""
    
    @staticmethod
    def combine_pvalues(pvalues: np.ndarray, weights: Optional[np.ndarray] = None) -> Tuple[float, float]:
        """
        Combine p-values using Stouffer's method.
        
        Args:
            pvalues: Array of p-values to combine
            weights: Optional weights for each p-value
            
        Returns:
            Tuple of (combined_z_score, combined_p_value)
        """
        # Remove NaN values
        valid_mask = ~np.isnan(pvalues)
        pvalues = pvalues[valid_mask]
        
        if len(pvalues) == 0:
            return np.nan, np.nan
        
        if weights is not None:
            weights = weights[valid_mask]
            if len(weights) != len(pvalues):
                raise ValueError("Weights and p-values must have same length after removing NaNs")
        
        # Convert p-values to z-scores
        z_scores = norm.ppf(1 - pvalues)
        
        # Handle extreme p-values
        z_scores = np.clip(z_scores, -10, 10)
        
        if weights is None:
            # Unweighted Stouffer's method
            combined_z = np.sum(z_scores) / np.sqrt(len(z_scores))
        else:
            # Weighted Stouffer's method
            weights = weights / np.sum(weights)  # Normalize weights
            combined_z = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))
        
        # Convert back to p-value
        combined_p = 1 - norm.cdf(combined_z)
        
        return combined_z, combined_p


class MotifEnrichmentAnalyzer:
    """Motif enrichment analysis using PWM scanning and statistical testing."""
    
    def __init__(self, pwm_scorer: PWMScorer, min_score_threshold: float = 5.0):
        """
        Initialize motif enrichment analyzer.
        
        Args:
            pwm_scorer: PWMScorer instance
            min_score_threshold: Minimum PWM score to consider a hit
        """
        self.pwm_scorer = pwm_scorer
        self.min_score_threshold = min_score_threshold
        self.stouffer = StoufferCombination()
    
    def find_motif_hits(self, sequences: List[str], pwm: np.ndarray) -> List[List[Tuple[int, float]]]:
        """
        Find motif hits in a list of sequences.
        
        Args:
            sequences: List of DNA sequences
            pwm: Position Weight Matrix
            
        Returns:
            List of hits for each sequence
        """
        all_hits = []
        
        for seq in sequences:
            hits = self.pwm_scorer.scan_sequence(seq, pwm)
            # Filter by threshold
            filtered_hits = [(pos, score) for pos, score in hits if score >= self.min_score_threshold]
            all_hits.append(filtered_hits)
        
        return all_hits
    
    def calculate_enrichment(self, 
                           foreground_sequences: List[str], 
                           background_sequences: List[str], 
                           pwm: np.ndarray) -> Dict[str, float]:
        """
        Calculate motif enrichment between foreground and background sequences.
        
        Args:
            foreground_sequences: Sequences of interest (e.g., near enhancers)
            background_sequences: Background sequences
            pwm: Position Weight Matrix
            
        Returns:
            Dictionary with enrichment statistics
        """
        # Find hits in both sets
        fg_hits = self.find_motif_hits(foreground_sequences, pwm)
        bg_hits = self.find_motif_hits(background_sequences, pwm)
        
        # Count sequences with hits
        fg_with_hits = sum(1 for hits in fg_hits if hits)
        bg_with_hits = sum(1 for hits in bg_hits if hits)
        
        fg_total = len(foreground_sequences)
        bg_total = len(background_sequences)
        
        # Calculate rates
        fg_rate = fg_with_hits / fg_total if fg_total > 0 else 0
        bg_rate = bg_with_hits / bg_total if bg_total > 0 else 0
        
        # Fisher's exact test
        contingency_table = [
            [fg_with_hits, fg_total - fg_with_hits],
            [bg_with_hits, bg_total - bg_with_hits]
        ]
        
        _, p_value = stats.fisher_exact(contingency_table, alternative='greater')
        
        # Calculate fold enrichment
        fold_enrichment = (fg_rate / bg_rate) if bg_rate > 0 else float('inf')
        
        return {
            'foreground_rate': fg_rate,
            'background_rate': bg_rate,
            'fold_enrichment': fold_enrichment,
            'p_value': p_value,
            'fg_hits': fg_with_hits,
            'fg_total': fg_total,
            'bg_hits': bg_with_hits,
            'bg_total': bg_total
        }
    
    def multi_motif_analysis(self, 
                           foreground_sequences: List[str], 
                           background_sequences: List[str], 
                           pwm_dict: Dict[str, np.ndarray]) -> pd.DataFrame:
        """
        Perform enrichment analysis for multiple motifs.
        
        Args:
            foreground_sequences: Sequences of interest
            background_sequences: Background sequences
            pwm_dict: Dictionary mapping motif names to PWMs
            
        Returns:
            DataFrame with enrichment results for each motif
        """
        results = []
        
        for motif_name, pwm in pwm_dict.items():
            try:
                enrichment = self.calculate_enrichment(
                    foreground_sequences, background_sequences, pwm
                )
                enrichment['motif_name'] = motif_name
                results.append(enrichment)
            except Exception as e:
                logger.warning(f"Error analyzing motif {motif_name}: {e}")
                continue
        
        df = pd.DataFrame(results)
        
        # Multiple testing correction
        if len(df) > 1:
            _, corrected_p, _, _ = stats.multipletests(df['p_value'], method='fdr_bh')
            df['p_value_corrected'] = corrected_p
        else:
            df['p_value_corrected'] = df['p_value']
        
        # Sort by significance
        df = df.sort_values('p_value')
        
        return df
    
    def combine_motif_evidence(self, 
                             motif_scores: Dict[str, List[float]], 
                             motif_weights: Optional[Dict[str, float]] = None) -> Tuple[float, float]:
        """
        Combine evidence from multiple motifs using Stouffer's method.
        
        Args:
            motif_scores: Dictionary mapping motif names to lists of p-values
            motif_weights: Optional weights for each motif
            
        Returns:
            Tuple of (combined_z_score, combined_p_value)
        """
        all_pvalues = []
        all_weights = []
        
        for motif_name, pvalues in motif_scores.items():
            weight = motif_weights.get(motif_name, 1.0) if motif_weights else 1.0
            
            all_pvalues.extend(pvalues)
            all_weights.extend([weight] * len(pvalues))
        
        if not all_pvalues:
            return np.nan, np.nan
        
        return self.stouffer.combine_pvalues(
            np.array(all_pvalues), 
            np.array(all_weights)
        )


def parse_meme_format(meme_file: str) -> Dict[str, np.ndarray]:
    """
    Parse PWMs from MEME format file.
    
    Args:
        meme_file: Path to MEME format file
        
    Returns:
        Dictionary mapping motif names to PWM matrices
    """
    pwms = {}
    
    with open(meme_file, 'r') as f:
        content = f.read()
    
    # Split into motif blocks
    motif_blocks = re.split(r'MOTIF\s+', content)[1:]  # Skip header
    
    for block in motif_blocks:
        lines = block.strip().split('\n')
        
        # Extract motif name
        motif_name = lines[0].split()[0]
        
        # Find letter-probability matrix
        matrix_start = None
        for i, line in enumerate(lines):
            if 'letter-probability matrix' in line:
                matrix_start = i + 1
                break
        
        if matrix_start is None:
            continue
        
        # Parse matrix
        matrix_data = []
        for line in lines[matrix_start:]:
            if line.strip() and not line.startswith('//'):
                try:
                    values = [float(x) for x in line.strip().split()]
                    if len(values) == 4:  # A, C, G, T
                        matrix_data.append(values)
                    else:
                        break
                except ValueError:
                    break
        
        if matrix_data:
            pwm_matrix = np.array(matrix_data).T  # Transpose to get (4, length)
            # Convert to log-odds
            background = np.array([0.25, 0.25, 0.25, 0.25])
            pwm_matrix = np.log2(pwm_matrix / background.reshape(-1, 1))
            pwms[motif_name] = pwm_matrix
    
    return pwms


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee Motif Analysis - PWM scoring and enrichment analysis"
    )
    parser.add_argument(
        '--foreground', '-f', required=True,
        help="FASTA file with foreground sequences"
    )
    parser.add_argument(
        '--background', '-b', required=True,
        help="FASTA file with background sequences"
    )
    parser.add_argument(
        '--motifs', '-m', required=True,
        help="MEME format file with PWMs"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output file for enrichment results"
    )
    parser.add_argument(
        '--min-score', type=float, default=5.0,
        help="Minimum PWM score threshold (default: 5.0)"
    )
    parser.add_argument(
        '--pseudocount', type=float, default=0.001,
        help="Pseudocount for PWM calculation (default: 0.001)"
    )
    
    args = parser.parse_args()
    
    # Read sequences
    def read_fasta(filename):
        sequences = []
        with open(filename, 'r') as f:
            sequence = ""
            for line in f:
                if line.startswith('>'):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                sequences.append(sequence)
        return sequences
    
    logger.info("Reading sequences...")
    fg_sequences = read_fasta(args.foreground)
    bg_sequences = read_fasta(args.background)
    
    logger.info(f"Loaded {len(fg_sequences)} foreground and {len(bg_sequences)} background sequences")
    
    # Parse motifs
    logger.info("Parsing motifs...")
    pwm_dict = parse_meme_format(args.motifs)
    logger.info(f"Loaded {len(pwm_dict)} motifs")
    
    # Initialize analyzer
    pwm_scorer = PWMScorer(pseudocount=args.pseudocount)
    analyzer = MotifEnrichmentAnalyzer(pwm_scorer, min_score_threshold=args.min_score)
    
    # Perform analysis
    logger.info("Performing enrichment analysis...")
    results_df = analyzer.multi_motif_analysis(fg_sequences, bg_sequences, pwm_dict)
    
    # Save results
    results_df.to_csv(args.output, index=False)
    logger.info(f"Results saved to {args.output}")
    
    # Print summary
    significant = results_df[results_df['p_value_corrected'] < 0.05]
    logger.info(f"Found {len(significant)} significantly enriched motifs")


if __name__ == "__main__":
    main()