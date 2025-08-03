#!/usr/bin/env python
"""
TFSee Multi-View Clustering Module

This module implements multi-view clustering algorithms for discovering
TF-enhancer relationship patterns from heterogeneous data sources.

Key Features:
- Multi-view spectral clustering
- Consensus clustering across views
- View-specific and shared embedding discovery
- Co-clustering of TFs and enhancers
- Robust clustering evaluation metrics

Multi-view clustering integrates:
- Motif presence/absence patterns
- Distance-based proximity networks
- Expression correlation patterns
- Chromatin accessibility patterns
- ChIP-seq overlap networks
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA, NMF
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.neighbors import kneighbors_graph
import warnings
from typing import List, Dict, Tuple, Optional, Union
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
import pickle

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class ClusteringResult:
    """Container for clustering results."""
    labels: np.ndarray
    n_clusters: int
    silhouette_score: float
    view_weights: Optional[np.ndarray] = None
    embedding: Optional[np.ndarray] = None
    consensus_matrix: Optional[np.ndarray] = None


class MultiViewSpectralClustering:
    """Multi-view spectral clustering for TF-enhancer relationship discovery."""
    
    def __init__(self, 
                 n_clusters: int = 10,
                 view_weights: Optional[np.ndarray] = None,
                 gamma: float = 1.0,
                 normalize_views: bool = True):
        """
        Initialize multi-view spectral clustering.
        
        Args:
            n_clusters: Number of clusters to find
            view_weights: Weights for each view (if None, equal weights)
            gamma: Parameter for view weight optimization
            normalize_views: Whether to normalize each view
        """
        self.n_clusters = n_clusters
        self.view_weights = view_weights
        self.gamma = gamma
        self.normalize_views = normalize_views
    
    def _normalize_affinity_matrix(self, affinity: np.ndarray) -> np.ndarray:
        """Normalize affinity matrix to [0, 1] range."""
        if np.max(affinity) == np.min(affinity):
            return np.ones_like(affinity) * 0.5
        
        normalized = (affinity - np.min(affinity)) / (np.max(affinity) - np.min(affinity))
        
        # Ensure symmetry
        normalized = (normalized + normalized.T) / 2
        
        # Set diagonal to 1
        np.fill_diagonal(normalized, 1.0)
        
        return normalized
    
    def _build_knn_graph(self, features: np.ndarray, k: int = 10) -> np.ndarray:
        """Build k-nearest neighbor graph from feature matrix."""
        if features.shape[0] < k:
            k = features.shape[0] - 1
        
        if k <= 0:
            return np.eye(features.shape[0])
        
        # Build KNN graph
        knn_graph = kneighbors_graph(features, n_neighbors=k, mode='connectivity')
        
        # Make symmetric
        knn_graph = (knn_graph + knn_graph.T) / 2
        
        return knn_graph.toarray()
    
    def _compute_laplacian(self, affinity: np.ndarray) -> np.ndarray:
        """Compute normalized graph Laplacian."""
        # Degree matrix
        degree = np.sum(affinity, axis=1)
        degree[degree == 0] = 1  # Avoid division by zero
        
        # Normalized Laplacian: L = I - D^(-1/2) * A * D^(-1/2)
        d_sqrt_inv = np.diag(1.0 / np.sqrt(degree))
        normalized_affinity = d_sqrt_inv @ affinity @ d_sqrt_inv
        
        laplacian = np.eye(affinity.shape[0]) - normalized_affinity
        
        return laplacian
    
    def _optimize_view_weights(self, view_embeddings: List[np.ndarray]) -> np.ndarray:
        """Optimize view weights based on embedding quality."""
        n_views = len(view_embeddings)
        
        if self.view_weights is not None:
            return self.view_weights
        
        # Initialize equal weights
        weights = np.ones(n_views) / n_views
        
        # Simple heuristic: weight by explained variance in embedding
        view_qualities = []
        for embedding in view_embeddings:
            if embedding.shape[1] > 1:
                # Use first PC explained variance as quality measure
                pca = PCA(n_components=1)
                pca.fit(embedding)
                quality = pca.explained_variance_ratio_[0]
            else:
                quality = np.var(embedding[:, 0])
            view_qualities.append(quality)
        
        # Normalize qualities to get weights
        view_qualities = np.array(view_qualities)
        if np.sum(view_qualities) > 0:
            weights = view_qualities / np.sum(view_qualities)
        
        return weights
    
    def fit_predict(self, view_data: List[np.ndarray]) -> ClusteringResult:
        """
        Perform multi-view spectral clustering.
        
        Args:
            view_data: List of data matrices, one per view
            
        Returns:
            ClusteringResult object
        """
        n_views = len(view_data)
        n_samples = view_data[0].shape[0]
        
        # Verify all views have same number of samples
        for i, data in enumerate(view_data):
            if data.shape[0] != n_samples:
                raise ValueError(f"View {i} has {data.shape[0]} samples, expected {n_samples}")
        
        # Normalize views if requested
        if self.normalize_views:
            view_data = [StandardScaler().fit_transform(data) for data in view_data]
        
        # Build affinity matrices for each view
        affinity_matrices = []
        view_embeddings = []
        
        for i, data in enumerate(view_data):
            logger.info(f"Processing view {i+1}/{n_views}")
            
            # Build affinity matrix (using correlation or RBF kernel)
            if data.shape[1] > 1:
                # Use correlation for multi-dimensional data
                corr_matrix = np.corrcoef(data)
                # Convert to positive similarity
                affinity = (corr_matrix + 1) / 2
            else:
                # Use RBF kernel for 1D data
                distances = squareform(pdist(data.reshape(-1, 1)))
                sigma = np.median(distances)
                affinity = np.exp(-distances**2 / (2 * sigma**2))
            
            # Normalize affinity matrix
            affinity = self._normalize_affinity_matrix(affinity)
            affinity_matrices.append(affinity)
            
            # Compute Laplacian and eigendecomposition
            laplacian = self._compute_laplacian(affinity)
            
            # Compute eigenvectors (smallest eigenvalues)
            eigenvals, eigenvecs = np.linalg.eigh(laplacian)
            
            # Select embedding dimensions (k smallest eigenvalues, skip first if ~0)
            start_idx = 1 if eigenvals[0] < 1e-8 else 0
            end_idx = min(start_idx + self.n_clusters, len(eigenvals))
            
            embedding = eigenvecs[:, start_idx:end_idx]
            view_embeddings.append(embedding)
        
        # Optimize view weights
        weights = self._optimize_view_weights(view_embeddings)
        logger.info(f"View weights: {weights}")
        
        # Combine embeddings using weighted average
        combined_embedding = np.zeros((n_samples, self.n_clusters))
        
        for i, (embedding, weight) in enumerate(zip(view_embeddings, weights)):
            # Pad or truncate embedding to match n_clusters
            if embedding.shape[1] < self.n_clusters:
                padded_embedding = np.zeros((n_samples, self.n_clusters))
                padded_embedding[:, :embedding.shape[1]] = embedding
                embedding = padded_embedding
            elif embedding.shape[1] > self.n_clusters:
                embedding = embedding[:, :self.n_clusters]
            
            combined_embedding += weight * embedding
        
        # Normalize rows of combined embedding
        combined_embedding = normalize(combined_embedding, norm='l2', axis=1)
        
        # Perform K-means clustering on combined embedding
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(combined_embedding)
        
        # Calculate silhouette score
        if self.n_clusters > 1:
            silhouette = silhouette_score(combined_embedding, labels)
        else:
            silhouette = 0.0
        
        # Build consensus matrix
        consensus_matrix = self._build_consensus_matrix(affinity_matrices, weights)
        
        return ClusteringResult(
            labels=labels,
            n_clusters=self.n_clusters,
            silhouette_score=silhouette,
            view_weights=weights,
            embedding=combined_embedding,
            consensus_matrix=consensus_matrix
        )
    
    def _build_consensus_matrix(self, 
                              affinity_matrices: List[np.ndarray], 
                              weights: np.ndarray) -> np.ndarray:
        """Build consensus matrix from multiple views."""
        n_samples = affinity_matrices[0].shape[0]
        consensus = np.zeros((n_samples, n_samples))
        
        for affinity, weight in zip(affinity_matrices, weights):
            consensus += weight * affinity
        
        return consensus


class CoClusteringTFEnhancer:
    """Co-clustering algorithm for TF-enhancer relationship discovery."""
    
    def __init__(self, 
                 n_tf_clusters: int = 10,
                 n_enhancer_clusters: int = 15,
                 max_iter: int = 100,
                 tol: float = 1e-4):
        """
        Initialize co-clustering algorithm.
        
        Args:
            n_tf_clusters: Number of TF clusters
            n_enhancer_clusters: Number of enhancer clusters
            max_iter: Maximum iterations for convergence
            tol: Convergence tolerance
        """
        self.n_tf_clusters = n_tf_clusters
        self.n_enhancer_clusters = n_enhancer_clusters
        self.max_iter = max_iter
        self.tol = tol
    
    def fit_predict(self, tf_enhancer_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Perform co-clustering on TF-enhancer association matrix.
        
        Args:
            tf_enhancer_matrix: Matrix of shape (n_tfs, n_enhancers)
            
        Returns:
            Tuple of (tf_labels, enhancer_labels)
        """
        n_tfs, n_enhancers = tf_enhancer_matrix.shape
        
        # Initialize cluster assignments randomly
        tf_labels = np.random.randint(0, self.n_tf_clusters, n_tfs)
        enhancer_labels = np.random.randint(0, self.n_enhancer_clusters, n_enhancers)
        
        prev_tf_labels = tf_labels.copy()
        prev_enhancer_labels = enhancer_labels.copy()
        
        for iteration in range(self.max_iter):
            # Update TF clusters based on enhancer clusters
            for tf_idx in range(n_tfs):
                best_cluster = 0
                best_score = -np.inf
                
                for cluster in range(self.n_tf_clusters):
                    # Calculate score for assigning TF to this cluster
                    score = 0.0
                    for enh_cluster in range(self.n_enhancer_clusters):
                        enhancer_mask = enhancer_labels == enh_cluster
                        if np.sum(enhancer_mask) > 0:
                            cluster_score = np.mean(tf_enhancer_matrix[tf_idx, enhancer_mask])
                            score += cluster_score
                    
                    if score > best_score:
                        best_score = score
                        best_cluster = cluster
                
                tf_labels[tf_idx] = best_cluster
            
            # Update enhancer clusters based on TF clusters
            for enh_idx in range(n_enhancers):
                best_cluster = 0
                best_score = -np.inf
                
                for cluster in range(self.n_enhancer_clusters):
                    # Calculate score for assigning enhancer to this cluster
                    score = 0.0
                    for tf_cluster in range(self.n_tf_clusters):
                        tf_mask = tf_labels == tf_cluster
                        if np.sum(tf_mask) > 0:
                            cluster_score = np.mean(tf_enhancer_matrix[tf_mask, enh_idx])
                            score += cluster_score
                    
                    if score > best_score:
                        best_score = score
                        best_cluster = cluster
                
                enhancer_labels[enh_idx] = best_cluster
            
            # Check convergence
            tf_changed = np.sum(tf_labels != prev_tf_labels)
            enh_changed = np.sum(enhancer_labels != prev_enhancer_labels)
            
            change_ratio = (tf_changed + enh_changed) / (n_tfs + n_enhancers)
            
            if change_ratio < self.tol:
                logger.info(f"Co-clustering converged after {iteration + 1} iterations")
                break
            
            prev_tf_labels = tf_labels.copy()
            prev_enhancer_labels = enhancer_labels.copy()
        
        return tf_labels, enhancer_labels


class ConsensusClusteringEvaluator:
    """Evaluate and compare clustering results across multiple views."""
    
    @staticmethod
    def calculate_consensus_matrix(cluster_results: List[np.ndarray]) -> np.ndarray:
        """
        Calculate consensus matrix from multiple clustering results.
        
        Args:
            cluster_results: List of cluster label arrays
            
        Returns:
            Consensus matrix indicating co-clustering frequency
        """
        n_samples = len(cluster_results[0])
        n_clusterings = len(cluster_results)
        
        consensus = np.zeros((n_samples, n_samples))
        
        for labels in cluster_results:
            for i in range(n_samples):
                for j in range(n_samples):
                    if labels[i] == labels[j]:
                        consensus[i, j] += 1
        
        # Normalize by number of clusterings
        consensus /= n_clusterings
        
        return consensus
    
    @staticmethod
    def evaluate_clustering_stability(cluster_results: List[np.ndarray]) -> Dict[str, float]:
        """
        Evaluate stability of clustering across multiple runs.
        
        Args:
            cluster_results: List of cluster label arrays from different runs
            
        Returns:
            Dictionary with stability metrics
        """
        n_clusterings = len(cluster_results)
        
        if n_clusterings < 2:
            return {'mean_ari': 1.0, 'std_ari': 0.0, 'mean_nmi': 1.0, 'std_nmi': 0.0}
        
        # Calculate pairwise ARI and NMI scores
        ari_scores = []
        nmi_scores = []
        
        for i in range(n_clusterings):
            for j in range(i + 1, n_clusterings):
                ari = adjusted_rand_score(cluster_results[i], cluster_results[j])
                nmi = normalized_mutual_info_score(cluster_results[i], cluster_results[j])
                
                ari_scores.append(ari)
                nmi_scores.append(nmi)
        
        return {
            'mean_ari': np.mean(ari_scores),
            'std_ari': np.std(ari_scores),
            'mean_nmi': np.mean(nmi_scores),
            'std_nmi': np.std(nmi_scores)
        }
    
    @staticmethod
    def optimal_cluster_number(data: np.ndarray, 
                             max_clusters: int = 20,
                             method: str = 'silhouette') -> Tuple[int, List[float]]:
        """
        Find optimal number of clusters using silhouette or elbow method.
        
        Args:
            data: Data matrix for clustering
            max_clusters: Maximum number of clusters to test
            method: 'silhouette' or 'elbow'
            
        Returns:
            Tuple of (optimal_k, scores)
        """
        scores = []
        k_range = range(2, min(max_clusters + 1, data.shape[0]))
        
        for k in k_range:
            if method == 'silhouette':
                kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                labels = kmeans.fit_predict(data)
                score = silhouette_score(data, labels)
                scores.append(score)
            
            elif method == 'elbow':
                kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                kmeans.fit(data)
                score = -kmeans.inertia_  # Negative for maximization
                scores.append(score)
        
        if method == 'silhouette':
            optimal_idx = np.argmax(scores)
        else:  # elbow method
            # Find elbow point using second derivative
            if len(scores) >= 3:
                second_deriv = np.diff(scores, 2)
                optimal_idx = np.argmax(second_deriv) + 1
            else:
                optimal_idx = 0
        
        optimal_k = list(k_range)[optimal_idx]
        
        return optimal_k, scores


def visualize_clustering_results(embedding: np.ndarray, 
                               labels: np.ndarray, 
                               output_file: str = None):
    """
    Visualize clustering results in 2D.
    
    Args:
        embedding: 2D embedding of data points
        labels: Cluster labels
        output_file: Optional output file path
    """
    # Reduce to 2D if needed
    if embedding.shape[1] > 2:
        pca = PCA(n_components=2)
        embedding_2d = pca.fit_transform(embedding)
    else:
        embedding_2d = embedding
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    unique_labels = np.unique(labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
    
    for i, label in enumerate(unique_labels):
        mask = labels == label
        plt.scatter(embedding_2d[mask, 0], embedding_2d[mask, 1], 
                   c=[colors[i]], label=f'Cluster {label}', alpha=0.7)
    
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.title('TF-Enhancer Clustering Results')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    plt.close()


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="TFSee Multi-View Clustering for TF-Enhancer Relationships"
    )
    parser.add_argument(
        '--view-data', nargs='+', required=True,
        help="List of data files for different views (CSV format)"
    )
    parser.add_argument(
        '--view-names', nargs='+',
        help="Names for each view (default: view1, view2, ...)"
    )
    parser.add_argument(
        '--n-clusters', type=int, default=10,
        help="Number of clusters (default: 10)"
    )
    parser.add_argument(
        '--output-prefix', '-o', required=True,
        help="Output prefix for results files"
    )
    parser.add_argument(
        '--method', choices=['multiview', 'coclustering'], default='multiview',
        help="Clustering method (default: multiview)"
    )
    parser.add_argument(
        '--optimize-k', action='store_true',
        help="Optimize number of clusters using silhouette score"
    )
    parser.add_argument(
        '--n-runs', type=int, default=10,
        help="Number of runs for stability evaluation (default: 10)"
    )
    
    args = parser.parse_args()
    
    # Load view data
    view_data = []
    view_names = args.view_names if args.view_names else [f"view{i+1}" for i in range(len(args.view_data))]
    
    logger.info(f"Loading {len(args.view_data)} views...")
    for i, data_file in enumerate(args.view_data):
        data = pd.read_csv(data_file, index_col=0).values
        view_data.append(data)
        logger.info(f"View {view_names[i]}: {data.shape}")
    
    # Optimize number of clusters if requested
    if args.optimize_k:
        logger.info("Optimizing number of clusters...")
        # Use first view for optimization
        optimal_k, silhouette_scores = ConsensusClusteringEvaluator.optimal_cluster_number(
            view_data[0], max_clusters=20
        )
        logger.info(f"Optimal number of clusters: {optimal_k}")
        n_clusters = optimal_k
    else:
        n_clusters = args.n_clusters
    
    # Perform clustering
    if args.method == 'multiview':
        logger.info("Performing multi-view spectral clustering...")
        
        # Multiple runs for stability evaluation
        cluster_results = []
        best_result = None
        best_silhouette = -1
        
        for run in range(args.n_runs):
            clusterer = MultiViewSpectralClustering(n_clusters=n_clusters)
            result = clusterer.fit_predict(view_data)
            cluster_results.append(result.labels)
            
            if result.silhouette_score > best_silhouette:
                best_silhouette = result.silhouette_score
                best_result = result
        
        # Evaluate stability
        stability_metrics = ConsensusClusteringEvaluator.evaluate_clustering_stability(cluster_results)
        
        logger.info(f"Best silhouette score: {best_silhouette:.3f}")
        logger.info(f"Clustering stability - ARI: {stability_metrics['mean_ari']:.3f} ± {stability_metrics['std_ari']:.3f}")
        
        # Save results
        results_df = pd.DataFrame({
            'cluster_label': best_result.labels,
            'silhouette_score': [best_result.silhouette_score] * len(best_result.labels)
        })
        
        # Add view weights
        view_weights_df = pd.DataFrame({
            'view_name': view_names,
            'weight': best_result.view_weights
        })
        
        # Save files
        results_df.to_csv(f"{args.output_prefix}_clustering_results.csv")
        view_weights_df.to_csv(f"{args.output_prefix}_view_weights.csv", index=False)
        
        # Save embedding and consensus matrix
        np.save(f"{args.output_prefix}_embedding.npy", best_result.embedding)
        np.save(f"{args.output_prefix}_consensus_matrix.npy", best_result.consensus_matrix)
        
        # Create visualization
        visualize_clustering_results(
            best_result.embedding, 
            best_result.labels,
            f"{args.output_prefix}_clustering_plot.png"
        )
        
    elif args.method == 'coclustering':
        logger.info("Performing co-clustering...")
        
        if len(view_data) != 1:
            raise ValueError("Co-clustering requires exactly one data matrix (TF x enhancer)")
        
        tf_enhancer_matrix = view_data[0]
        
        coclustering = CoClusteringTFEnhancer(
            n_tf_clusters=n_clusters,
            n_enhancer_clusters=n_clusters
        )
        
        tf_labels, enhancer_labels = coclustering.fit_predict(tf_enhancer_matrix)
        
        # Save results
        tf_results_df = pd.DataFrame({
            'tf_index': range(len(tf_labels)),
            'cluster_label': tf_labels
        })
        
        enhancer_results_df = pd.DataFrame({
            'enhancer_index': range(len(enhancer_labels)),
            'cluster_label': enhancer_labels
        })
        
        tf_results_df.to_csv(f"{args.output_prefix}_tf_clusters.csv", index=False)
        enhancer_results_df.to_csv(f"{args.output_prefix}_enhancer_clusters.csv", index=False)
    
    # Save stability metrics
    with open(f"{args.output_prefix}_stability_metrics.txt", 'w') as f:
        if args.method == 'multiview':
            f.write(f"Clustering Stability Metrics:\n")
            f.write(f"Mean ARI: {stability_metrics['mean_ari']:.3f} ± {stability_metrics['std_ari']:.3f}\n")
            f.write(f"Mean NMI: {stability_metrics['mean_nmi']:.3f} ± {stability_metrics['std_nmi']:.3f}\n")
            f.write(f"Best Silhouette Score: {best_silhouette:.3f}\n")
            f.write(f"Number of Clusters: {n_clusters}\n")
            f.write(f"Number of Runs: {args.n_runs}\n")
    
    logger.info(f"Results saved with prefix: {args.output_prefix}")


if __name__ == "__main__":
    main()