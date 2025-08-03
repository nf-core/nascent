#!/usr/bin/env python
"""
TFSee Statistical Methods Module

This module implements comprehensive statistical methods for significance testing
in TFSee analysis, including multiple testing correction, permutation tests,
and specialized genomic statistics.

Key Features:
- Multiple testing correction (FDR, Bonferroni, etc.)
- Permutation and bootstrap testing
- Genomic-specific statistical tests
- Effect size calculations
- Confidence interval estimation
- Power analysis
- Model selection and validation
- Robust statistical methods

Statistical tests implemented:
- Enrichment tests (Fisher's exact, hypergeometric)
- Association tests (correlation, regression)
- Comparative tests (t-test, Mann-Whitney, Kolmogorov-Smirnov)
- Clustering validation (silhouette, gap statistic)
- Gene set enrichment analysis
- Spatial statistics for genomic data
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import (
    ttest_ind, mannwhitneyu, ks_2samp, fisher_exact, chi2_contingency,
    pearsonr, spearmanr, kendalltau, multivariate_normal, norm,
    false_discovery_control, combine_pvalues
)
from scipy.special import comb
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.power import ttest_power
from sklearn.model_selection import PermutationTestScore, cross_val_score
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import warnings
from typing import List, Dict, Tuple, Optional, Union, Callable
from dataclasses import dataclass
from collections import defaultdict
import itertools

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class StatisticalResult:
    """Container for statistical test results."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: Optional[float] = None
    confidence_interval: Optional[Tuple[float, float]] = None
    sample_size: Optional[int] = None
    power: Optional[float] = None
    corrected_p_value: Optional[float] = None
    metadata: Optional[Dict] = None


@dataclass
class EnrichmentResult:
    """Container for enrichment analysis results."""
    category: str
    observed: int
    expected: float
    fold_enrichment: float
    p_value: float
    odds_ratio: float
    confidence_interval: Tuple[float, float]
    background_size: int
    category_size: int


class MultipleTestingCorrector:
    """Handle multiple testing correction with various methods."""
    
    def __init__(self, alpha: float = 0.05):
        """
        Initialize multiple testing corrector.
        
        Args:
            alpha: Family-wise error rate
        """
        self.alpha = alpha
    
    def correct_pvalues(self, 
                       pvalues: np.ndarray, 
                       method: str = 'fdr_bh') -> Tuple[np.ndarray, np.ndarray]:
        """
        Correct p-values for multiple testing.
        
        Args:
            pvalues: Array of p-values to correct
            method: Correction method ('fdr_bh', 'fdr_by', 'bonferroni', 'holm', etc.)
            
        Returns:
            Tuple of (rejected_hypotheses, corrected_pvalues)
        """
        # Remove NaN values
        valid_mask = ~np.isnan(pvalues)
        valid_pvalues = pvalues[valid_mask]
        
        if len(valid_pvalues) == 0:
            return np.zeros_like(pvalues, dtype=bool), np.full_like(pvalues, np.nan)
        
        # Perform correction
        rejected, corrected_p, _, _ = multipletests(
            valid_pvalues, 
            alpha=self.alpha, 
            method=method
        )
        
        # Reconstruct full arrays
        full_rejected = np.zeros(len(pvalues), dtype=bool)
        full_corrected = np.full(len(pvalues), np.nan)
        
        full_rejected[valid_mask] = rejected
        full_corrected[valid_mask] = corrected_p
        
        return full_rejected, full_corrected
    
    def stouffer_combination(self, 
                           pvalues: np.ndarray, 
                           weights: Optional[np.ndarray] = None) -> Tuple[float, float]:
        """
        Combine p-values using Stouffer's method.
        
        Args:
            pvalues: Array of p-values to combine
            weights: Optional weights for each p-value
            
        Returns:
            Tuple of (combined_statistic, combined_p_value)
        """
        # Remove NaN values
        valid_mask = ~np.isnan(pvalues)
        valid_pvalues = pvalues[valid_mask]
        
        if len(valid_pvalues) == 0:
            return np.nan, np.nan
        
        if weights is not None:
            valid_weights = weights[valid_mask]
        else:
            valid_weights = None
        
        return combine_pvalues(valid_pvalues, method='stouffer', weights=valid_weights)
    
    def fisher_combination(self, pvalues: np.ndarray) -> Tuple[float, float]:
        """
        Combine p-values using Fisher's method.
        
        Args:
            pvalues: Array of p-values to combine
            
        Returns:
            Tuple of (combined_statistic, combined_p_value)
        """
        valid_pvalues = pvalues[~np.isnan(pvalues)]
        
        if len(valid_pvalues) == 0:
            return np.nan, np.nan
        
        return combine_pvalues(valid_pvalues, method='fisher')


class PermutationTester:
    """Perform permutation tests for various genomic analyses."""
    
    def __init__(self, n_permutations: int = 10000, random_state: int = 42):
        """
        Initialize permutation tester.
        
        Args:
            n_permutations: Number of permutations to perform
            random_state: Random seed for reproducibility
        """
        self.n_permutations = n_permutations
        self.random_state = random_state
        np.random.seed(random_state)
    
    def permutation_test(self, 
                        group1: np.ndarray, 
                        group2: np.ndarray, 
                        test_statistic: Callable = None) -> StatisticalResult:
        """
        Perform permutation test between two groups.
        
        Args:
            group1: First group data
            group2: Second group data
            test_statistic: Function to compute test statistic
            
        Returns:
            StatisticalResult object
        """
        if test_statistic is None:
            test_statistic = lambda x, y: np.mean(x) - np.mean(y)
        
        # Observed test statistic
        observed_stat = test_statistic(group1, group2)
        
        # Combine groups for permutation
        combined = np.concatenate([group1, group2])
        n1, n2 = len(group1), len(group2)
        
        # Permutation distribution
        perm_stats = []
        for _ in range(self.n_permutations):
            # Shuffle and split
            np.random.shuffle(combined)
            perm_group1 = combined[:n1]
            perm_group2 = combined[n1:]
            
            perm_stat = test_statistic(perm_group1, perm_group2)
            perm_stats.append(perm_stat)
        
        perm_stats = np.array(perm_stats)
        
        # Calculate p-value (two-tailed)
        p_value = np.mean(np.abs(perm_stats) >= np.abs(observed_stat))
        
        # Effect size (Cohen's d for mean difference)
        if test_statistic == lambda x, y: np.mean(x) - np.mean(y):
            pooled_std = np.sqrt(((n1-1)*np.var(group1) + (n2-1)*np.var(group2)) / (n1+n2-2))
            effect_size = observed_stat / pooled_std if pooled_std > 0 else 0
        else:
            effect_size = None
        
        return StatisticalResult(
            test_name="permutation_test",
            statistic=observed_stat,
            p_value=p_value,
            effect_size=effect_size,
            sample_size=n1 + n2,
            metadata={'n_permutations': self.n_permutations}
        )
    
    def permutation_correlation_test(self, 
                                   x: np.ndarray, 
                                   y: np.ndarray) -> StatisticalResult:
        """
        Permutation test for correlation significance.
        
        Args:
            x: First variable
            y: Second variable
            
        Returns:
            StatisticalResult object
        """
        # Observed correlation
        observed_corr, _ = pearsonr(x, y)
        
        # Permutation distribution
        perm_corrs = []
        for _ in range(self.n_permutations):
            y_permuted = np.random.permutation(y)
            perm_corr, _ = pearsonr(x, y_permuted)
            perm_corrs.append(perm_corr)
        
        perm_corrs = np.array(perm_corrs)
        
        # Calculate p-value
        p_value = np.mean(np.abs(perm_corrs) >= np.abs(observed_corr))
        
        return StatisticalResult(
            test_name="permutation_correlation",
            statistic=observed_corr,
            p_value=p_value,
            effect_size=observed_corr**2,  # R-squared
            sample_size=len(x)
        )
    
    def permutation_enrichment_test(self, 
                                  query_set: set, 
                                  background_set: set, 
                                  annotation_sets: Dict[str, set]) -> Dict[str, StatisticalResult]:
        """
        Permutation test for set enrichment.
        
        Args:
            query_set: Set of query items
            background_set: Background set of items
            annotation_sets: Dictionary of annotation sets
            
        Returns:
            Dictionary of StatisticalResult objects for each annotation
        """
        results = {}
        
        for annotation_name, annotation_set in annotation_sets.items():
            # Observed overlap
            observed_overlap = len(query_set.intersection(annotation_set))
            
            # Permutation distribution
            perm_overlaps = []
            query_size = len(query_set)
            
            for _ in range(self.n_permutations):
                # Sample random set of same size from background
                random_query = set(np.random.choice(
                    list(background_set), 
                    size=query_size, 
                    replace=False
                ))
                perm_overlap = len(random_query.intersection(annotation_set))
                perm_overlaps.append(perm_overlap)
            
            perm_overlaps = np.array(perm_overlaps)
            
            # Calculate p-value (one-tailed, testing for enrichment)
            p_value = np.mean(perm_overlaps >= observed_overlap)
            
            # Expected overlap under null
            expected_overlap = np.mean(perm_overlaps)
            
            # Fold enrichment
            fold_enrichment = observed_overlap / expected_overlap if expected_overlap > 0 else float('inf')
            
            results[annotation_name] = StatisticalResult(
                test_name="permutation_enrichment",
                statistic=observed_overlap,
                p_value=p_value,
                effect_size=fold_enrichment,
                sample_size=query_size,
                metadata={
                    'expected_overlap': expected_overlap,
                    'annotation_size': len(annotation_set),
                    'background_size': len(background_set)
                }
            )
        
        return results


class EnrichmentAnalyzer:
    """Perform enrichment analysis for genomic data."""
    
    def __init__(self):
        pass
    
    def hypergeometric_test(self, 
                          observed: int, 
                          query_size: int, 
                          category_size: int, 
                          background_size: int) -> EnrichmentResult:
        """
        Hypergeometric test for enrichment.
        
        Args:
            observed: Observed overlap
            query_size: Size of query set
            category_size: Size of category
            background_size: Total background size
            
        Returns:
            EnrichmentResult object
        """
        # Expected overlap under null hypothesis
        expected = (query_size * category_size) / background_size
        
        # Hypergeometric test
        p_value = stats.hypergeom.sf(observed - 1, background_size, category_size, query_size)
        
        # Odds ratio and confidence interval using Fisher's exact test
        # Construct 2x2 contingency table
        in_category_in_query = observed
        in_category_not_query = category_size - observed
        not_category_in_query = query_size - observed
        not_category_not_query = background_size - category_size - query_size + observed
        
        contingency_table = [
            [in_category_in_query, in_category_not_query],
            [not_category_in_query, not_category_not_query]
        ]
        
        odds_ratio, fisher_p = fisher_exact(contingency_table, alternative='greater')
        
        # Confidence interval for odds ratio (log scale)
        log_or = np.log(odds_ratio) if odds_ratio > 0 else -np.inf
        se_log_or = np.sqrt(sum(1/x for x in [in_category_in_query, in_category_not_query, 
                                            not_category_in_query, not_category_not_query] if x > 0))
        
        if se_log_or > 0 and np.isfinite(log_or):
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
        else:
            ci_lower, ci_upper = 0, np.inf
        
        # Fold enrichment
        fold_enrichment = observed / expected if expected > 0 else float('inf')
        
        return EnrichmentResult(
            category="test_category",
            observed=observed,
            expected=expected,
            fold_enrichment=fold_enrichment,
            p_value=p_value,
            odds_ratio=odds_ratio,
            confidence_interval=(ci_lower, ci_upper),
            background_size=background_size,
            category_size=category_size
        )
    
    def gsea_analysis(self, 
                     gene_scores: Dict[str, float], 
                     gene_sets: Dict[str, set], 
                     n_permutations: int = 1000) -> Dict[str, StatisticalResult]:
        """
        Gene Set Enrichment Analysis (GSEA).
        
        Args:
            gene_scores: Dictionary mapping gene IDs to scores
            gene_sets: Dictionary mapping set names to gene sets
            n_permutations: Number of permutations for p-value calculation
            
        Returns:
            Dictionary of StatisticalResult objects for each gene set
        """
        results = {}
        
        # Sort genes by score
        sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
        gene_list = [gene for gene, _ in sorted_genes]
        score_list = [score for _, score in sorted_genes]
        
        for set_name, gene_set in gene_sets.items():
            # Calculate enrichment score
            es = self._calculate_enrichment_score(gene_list, score_list, gene_set)
            
            # Permutation test for significance
            perm_scores = []
            for _ in range(n_permutations):
                # Permute gene scores
                perm_scores_dict = dict(zip(gene_scores.keys(), np.random.permutation(list(gene_scores.values()))))
                perm_sorted = sorted(perm_scores_dict.items(), key=lambda x: x[1], reverse=True)
                perm_gene_list = [gene for gene, _ in perm_sorted]
                perm_score_list = [score for _, score in perm_sorted]
                
                perm_es = self._calculate_enrichment_score(perm_gene_list, perm_score_list, gene_set)
                perm_scores.append(perm_es)
            
            # Calculate p-value
            p_value = np.mean(np.abs(perm_scores) >= np.abs(es))
            
            results[set_name] = StatisticalResult(
                test_name="gsea",
                statistic=es,
                p_value=p_value,
                effect_size=es,
                sample_size=len(gene_set),
                metadata={'gene_set_size': len(gene_set)}
            )
        
        return results
    
    def _calculate_enrichment_score(self, 
                                  gene_list: List[str], 
                                  score_list: List[float], 
                                  gene_set: set) -> float:
        """Calculate GSEA enrichment score."""
        n_genes = len(gene_list)
        n_set = len(gene_set)
        
        # Create running sum
        running_sum = 0
        max_deviation = 0
        
        # Sum of absolute scores for genes in set
        sum_scores_in_set = sum(abs(score_list[i]) for i, gene in enumerate(gene_list) if gene in gene_set)
        
        for i, gene in enumerate(gene_list):
            if gene in gene_set:
                # Gene is in set
                if sum_scores_in_set > 0:
                    running_sum += abs(score_list[i]) / sum_scores_in_set
            else:
                # Gene is not in set
                running_sum -= 1 / (n_genes - n_set)
            
            # Track maximum deviation
            if abs(running_sum) > abs(max_deviation):
                max_deviation = running_sum
        
        return max_deviation


class EffectSizeCalculator:
    """Calculate effect sizes for various statistical tests."""
    
    @staticmethod
    def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
        """
        Calculate Cohen's d for two groups.
        
        Args:
            group1: First group
            group2: Second group
            
        Returns:
            Cohen's d effect size
        """
        n1, n2 = len(group1), len(group2)
        
        # Pooled standard deviation
        pooled_std = np.sqrt(((n1-1)*np.var(group1, ddof=1) + (n2-1)*np.var(group2, ddof=1)) / (n1+n2-2))
        
        if pooled_std == 0:
            return 0
        
        return (np.mean(group1) - np.mean(group2)) / pooled_std
    
    @staticmethod
    def glass_delta(group1: np.ndarray, group2: np.ndarray) -> float:
        """
        Calculate Glass's delta (uses control group SD).
        
        Args:
            group1: Treatment group
            group2: Control group
            
        Returns:
            Glass's delta effect size
        """
        control_std = np.std(group2, ddof=1)
        
        if control_std == 0:
            return 0
        
        return (np.mean(group1) - np.mean(group2)) / control_std
    
    @staticmethod
    def correlation_effect_size(r: float) -> str:
        """
        Interpret correlation effect size.
        
        Args:
            r: Correlation coefficient
            
        Returns:
            Effect size interpretation
        """
        abs_r = abs(r)
        
        if abs_r < 0.1:
            return "negligible"
        elif abs_r < 0.3:
            return "small"
        elif abs_r < 0.5:
            return "medium"
        else:
            return "large"
    
    @staticmethod
    def odds_ratio_to_cohens_d(odds_ratio: float) -> float:
        """
        Convert odds ratio to Cohen's d.
        
        Args:
            odds_ratio: Odds ratio
            
        Returns:
            Approximate Cohen's d
        """
        if odds_ratio <= 0:
            return 0
        
        return np.log(odds_ratio) * np.sqrt(3) / np.pi


class PowerAnalyzer:
    """Perform power analysis for study design."""
    
    def __init__(self, alpha: float = 0.05):
        """
        Initialize power analyzer.
        
        Args:
            alpha: Significance level
        """
        self.alpha = alpha
    
    def t_test_power(self, 
                    effect_size: float, 
                    sample_size: int, 
                    alternative: str = 'two-sided') -> float:
        """
        Calculate power for t-test.
        
        Args:
            effect_size: Cohen's d effect size
            sample_size: Sample size per group
            alternative: Type of test ('two-sided', 'larger', 'smaller')
            
        Returns:
            Statistical power
        """
        return ttest_power(effect_size, sample_size, self.alpha, alternative=alternative)
    
    def sample_size_for_power(self, 
                            effect_size: float, 
                            power: float = 0.8, 
                            alternative: str = 'two-sided') -> int:
        """
        Calculate required sample size for desired power.
        
        Args:
            effect_size: Expected effect size
            power: Desired power
            alternative: Type of test
            
        Returns:
            Required sample size per group
        """
        from statsmodels.stats.power import ttest_power
        
        # Binary search for sample size
        min_n, max_n = 1, 10000
        
        while max_n - min_n > 1:
            mid_n = (min_n + max_n) // 2
            calculated_power = ttest_power(effect_size, mid_n, self.alpha, alternative=alternative)
            
            if calculated_power < power:
                min_n = mid_n
            else:
                max_n = mid_n
        
        return max_n
    
    def correlation_power(self, 
                         effect_size: float, 
                         sample_size: int) -> float:
        """
        Calculate power for correlation test.
        
        Args:
            effect_size: Expected correlation coefficient
            sample_size: Sample size
            
        Returns:
            Statistical power
        """
        # Fisher's z-transformation
        z_r = 0.5 * np.log((1 + effect_size) / (1 - effect_size))
        se = 1 / np.sqrt(sample_size - 3)
        
        # Critical value
        z_critical = norm.ppf(1 - self.alpha/2)
        
        # Power calculation
        power = 1 - norm.cdf(z_critical - abs(z_r) / se) + norm.cdf(-z_critical - abs(z_r) / se)
        
        return power


class ModelValidator:
    """Validate statistical models and predictions."""
    
    def __init__(self):
        pass
    
    def cross_validation_test(self, 
                            X: np.ndarray, 
                            y: np.ndarray, 
                            model=None, 
                            cv: int = 5) -> StatisticalResult:
        """
        Perform cross-validation test for model performance.
        
        Args:
            X: Feature matrix
            y: Target variable
            model: Model to test (default: RandomForest)
            cv: Number of cross-validation folds
            
        Returns:
            StatisticalResult with cross-validation scores
        """
        if model is None:
            if len(np.unique(y)) == 2:  # Binary classification
                model = LogisticRegression(random_state=42)
            else:  # Multiclass or regression
                model = RandomForestClassifier(random_state=42)
        
        # Perform cross-validation
        cv_scores = cross_val_score(model, X, y, cv=cv, scoring='roc_auc' if len(np.unique(y)) == 2 else 'accuracy')
        
        # Statistical test for significance
        t_stat, p_value = ttest_ind(cv_scores, [0.5] * len(cv_scores))  # Test against chance
        
        return StatisticalResult(
            test_name="cross_validation",
            statistic=np.mean(cv_scores),
            p_value=p_value,
            effect_size=np.mean(cv_scores) - 0.5,  # Improvement over chance
            sample_size=X.shape[0],
            metadata={
                'cv_scores': cv_scores.tolist(),
                'cv_std': np.std(cv_scores),
                'cv_folds': cv
            }
        )
    
    def permutation_importance_test(self, 
                                  X: np.ndarray, 
                                  y: np.ndarray, 
                                  model=None, 
                                  n_permutations: int = 100) -> Dict[int, StatisticalResult]:
        """
        Test feature importance using permutation.
        
        Args:
            X: Feature matrix
            y: Target variable
            model: Trained model
            n_permutations: Number of permutations
            
        Returns:
            Dictionary mapping feature indices to StatisticalResult objects
        """
        if model is None:
            model = RandomForestClassifier(random_state=42)
            model.fit(X, y)
        
        # Baseline score
        baseline_score = model.score(X, y)
        
        feature_importance = {}
        
        for feature_idx in range(X.shape[1]):
            # Permute feature and calculate score drop
            perm_scores = []
            
            for _ in range(n_permutations):
                X_perm = X.copy()
                X_perm[:, feature_idx] = np.random.permutation(X_perm[:, feature_idx])
                perm_score = model.score(X_perm, y)
                perm_scores.append(baseline_score - perm_score)
            
            # Statistical test
            importance = np.mean(perm_scores)
            t_stat, p_value = ttest_ind(perm_scores, [0] * len(perm_scores))
            
            feature_importance[feature_idx] = StatisticalResult(
                test_name="permutation_importance",
                statistic=importance,
                p_value=p_value,
                effect_size=importance / baseline_score if baseline_score > 0 else 0,
                sample_size=X.shape[0],
                metadata={'baseline_score': baseline_score}
            )
        
        return feature_importance


def comprehensive_statistical_analysis(data: pd.DataFrame, 
                                     target_column: str = None,
                                     group_column: str = None) -> Dict[str, StatisticalResult]:
    """
    Perform comprehensive statistical analysis on dataset.
    
    Args:
        data: Input dataframe
        target_column: Target variable for supervised analysis
        group_column: Grouping variable for comparative analysis
        
    Returns:
        Dictionary of statistical results
    """
    results = {}
    
    # Initialize components
    corrector = MultipleTestingCorrector()
    permutation_tester = PermutationTester()
    enrichment_analyzer = EnrichmentAnalyzer()
    effect_calculator = EffectSizeCalculator()
    power_analyzer = PowerAnalyzer()
    
    # Basic descriptive statistics
    numeric_columns = data.select_dtypes(include=[np.number]).columns
    
    # Normality tests
    for col in numeric_columns:
        if col not in [target_column, group_column]:
            values = data[col].dropna()
            if len(values) > 3:
                stat, p_val = stats.normaltest(values)
                results[f"{col}_normality"] = StatisticalResult(
                    test_name="normality_test",
                    statistic=stat,
                    p_value=p_val,
                    sample_size=len(values)
                )
    
    # Correlation analysis
    if len(numeric_columns) > 1:
        corr_matrix = data[numeric_columns].corr()
        
        # Test correlations for significance
        n_tests = 0
        corr_pvalues = []
        
        for i, col1 in enumerate(numeric_columns):
            for j, col2 in enumerate(numeric_columns[i+1:], i+1):
                if col1 != col2:
                    x = data[col1].dropna()
                    y = data[col2].dropna()
                    
                    # Find common indices
                    common_idx = data[col1].notna() & data[col2].notna()
                    x_common = data.loc[common_idx, col1]
                    y_common = data.loc[common_idx, col2]
                    
                    if len(x_common) > 3:
                        corr, p_val = pearsonr(x_common, y_common)
                        corr_pvalues.append(p_val)
                        
                        results[f"{col1}_vs_{col2}_correlation"] = StatisticalResult(
                            test_name="correlation",
                            statistic=corr,
                            p_value=p_val,
                            effect_size=corr**2,
                            sample_size=len(x_common)
                        )
                        n_tests += 1
        
        # Multiple testing correction for correlations
        if corr_pvalues:
            _, corrected_p = corrector.correct_pvalues(np.array(corr_pvalues))
            
            corr_test_idx = 0
            for i, col1 in enumerate(numeric_columns):
                for j, col2 in enumerate(numeric_columns[i+1:], i+1):
                    if col1 != col2:
                        key = f"{col1}_vs_{col2}_correlation"
                        if key in results:
                            results[key].corrected_p_value = corrected_p[corr_test_idx]
                            corr_test_idx += 1
    
    # Group comparisons
    if group_column and group_column in data.columns:
        groups = data[group_column].unique()
        
        if len(groups) == 2:
            for col in numeric_columns:
                if col != group_column:
                    group_data = {}
                    for group in groups:
                        group_data[group] = data[data[group_column] == group][col].dropna()
                    
                    if len(group_data[groups[0]]) > 0 and len(group_data[groups[1]]) > 0:
                        # T-test
                        t_stat, t_p = ttest_ind(group_data[groups[0]], group_data[groups[1]])
                        
                        # Mann-Whitney U test
                        u_stat, u_p = mannwhitneyu(group_data[groups[0]], group_data[groups[1]])
                        
                        # Effect size
                        cohens_d = effect_calculator.cohens_d(
                            group_data[groups[0]], 
                            group_data[groups[1]]
                        )
                        
                        results[f"{col}_ttest_{groups[0]}_vs_{groups[1]}"] = StatisticalResult(
                            test_name="t_test",
                            statistic=t_stat,
                            p_value=t_p,
                            effect_size=cohens_d,
                            sample_size=len(group_data[groups[0]]) + len(group_data[groups[1]])
                        )
                        
                        results[f"{col}_mannwhitney_{groups[0]}_vs_{groups[1]}"] = StatisticalResult(
                            test_name="mann_whitney",
                            statistic=u_stat,
                            p_value=u_p,
                            sample_size=len(group_data[groups[0]]) + len(group_data[groups[1]])
                        )
    
    return results


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee Statistical Analysis"
    )
    parser.add_argument(
        '--data', '-d', required=True,
        help="Input data file (CSV format)"
    )
    parser.add_argument(
        '--target-column',
        help="Target column for supervised analysis"
    )
    parser.add_argument(
        '--group-column',
        help="Grouping column for comparative analysis"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output file for statistical results"
    )
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help="Significance level (default: 0.05)"
    )
    parser.add_argument(
        '--correction-method', default='fdr_bh',
        choices=['fdr_bh', 'fdr_by', 'bonferroni', 'holm'],
        help="Multiple testing correction method"
    )
    parser.add_argument(
        '--n-permutations', type=int, default=1000,
        help="Number of permutations for permutation tests"
    )
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.data}")
    data = pd.read_csv(args.data, index_col=0)
    logger.info(f"Data shape: {data.shape}")
    
    # Perform comprehensive analysis
    logger.info("Performing statistical analysis...")
    results = comprehensive_statistical_analysis(
        data, 
        target_column=args.target_column,
        group_column=args.group_column
    )
    
    # Convert results to DataFrame
    results_data = []
    for test_name, result in results.items():
        row = {
            'test_name': test_name,
            'method': result.test_name,
            'statistic': result.statistic,
            'p_value': result.p_value,
            'effect_size': result.effect_size,
            'sample_size': result.sample_size,
            'corrected_p_value': result.corrected_p_value
        }
        results_data.append(row)
    
    results_df = pd.DataFrame(results_data)
    
    # Apply multiple testing correction if not already done
    if 'corrected_p_value' not in results_df.columns or results_df['corrected_p_value'].isna().all():
        corrector = MultipleTestingCorrector(alpha=args.alpha)
        _, corrected_p = corrector.correct_pvalues(
            results_df['p_value'].values, 
            method=args.correction_method
        )
        results_df['corrected_p_value'] = corrected_p
    
    # Save results
    results_df.to_csv(args.output, index=False)
    logger.info(f"Results saved to {args.output}")
    
    # Print summary
    significant_tests = results_df[results_df['corrected_p_value'] < args.alpha]
    logger.info(f"Performed {len(results_df)} statistical tests")
    logger.info(f"Found {len(significant_tests)} significant results after correction")
    
    if len(significant_tests) > 0:
        logger.info("Most significant results:")
        top_results = significant_tests.nsmallest(5, 'corrected_p_value')
        for _, row in top_results.iterrows():
            logger.info(f"  {row['test_name']}: p = {row['corrected_p_value']:.2e}, "
                       f"effect = {row['effect_size']:.3f}")


if __name__ == "__main__":
    main()