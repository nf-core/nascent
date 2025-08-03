#!/usr/bin/env python
"""
TFSee Optimization Module

This module implements optimization techniques for handling large-scale genomic
data in TFSee analysis, focusing on memory efficiency, computational speed,
and scalability.

Key Features:
- Memory-efficient data processing
- Parallel computation strategies
- Chunked data processing for large datasets
- Sparse matrix optimizations
- Caching and memoization
- GPU acceleration support (when available)
- Streaming data processing
- Incremental learning approaches

Optimization strategies:
- Batch processing for large datasets
- Memory mapping for large files
- Distributed computing support
- Efficient sparse operations
- Algorithmic complexity reduction
"""

import argparse
import logging
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
import multiprocessing as mp
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import h5py
import pickle
import warnings
from typing import List, Dict, Tuple, Optional, Union, Iterator, Callable
from dataclasses import dataclass
from functools import lru_cache, partial
import psutil
import gc
import time
from pathlib import Path
import tempfile
import mmap
import json

# Set up logging
logging.basicConfig(
    format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
    level=logging.INFO
)
logger = logging.getLogger(__file__)


@dataclass
class MemoryProfile:
    """Memory usage profiling information."""
    peak_memory_mb: float
    current_memory_mb: float
    available_memory_mb: float
    memory_usage_percentage: float


@dataclass
class OptimizationConfig:
    """Configuration for optimization parameters."""
    chunk_size: int = 1000
    n_workers: int = None
    use_sparse: bool = True
    cache_size: int = 128
    memory_limit_mb: float = None
    temp_directory: str = None
    enable_gpu: bool = False
    compression_level: int = 6


class MemoryMonitor:
    """Monitor and manage memory usage during computation."""
    
    def __init__(self, memory_limit_mb: Optional[float] = None):
        """
        Initialize memory monitor.
        
        Args:
            memory_limit_mb: Memory limit in MB (default: 80% of available)
        """
        if memory_limit_mb is None:
            total_memory = psutil.virtual_memory().total / (1024**2)
            self.memory_limit_mb = total_memory * 0.8
        else:
            self.memory_limit_mb = memory_limit_mb
        
        self.peak_memory = 0.0
    
    def get_memory_usage(self) -> MemoryProfile:
        """Get current memory usage profile."""
        memory = psutil.virtual_memory()
        current_memory = memory.used / (1024**2)
        available_memory = memory.available / (1024**2)
        
        self.peak_memory = max(self.peak_memory, current_memory)
        
        return MemoryProfile(
            peak_memory_mb=self.peak_memory,
            current_memory_mb=current_memory,
            available_memory_mb=available_memory,
            memory_usage_percentage=memory.percent
        )
    
    def check_memory_limit(self) -> bool:
        """Check if memory usage is within limits."""
        profile = self.get_memory_usage()
        return profile.current_memory_mb < self.memory_limit_mb
    
    def cleanup_memory(self):
        """Force garbage collection to free memory."""
        gc.collect()
    
    def memory_efficient_operation(self, operation: Callable, *args, **kwargs):
        """Execute operation with memory monitoring."""
        initial_memory = self.get_memory_usage().current_memory_mb
        
        try:
            result = operation(*args, **kwargs)
            
            final_memory = self.get_memory_usage().current_memory_mb
            memory_delta = final_memory - initial_memory
            
            logger.debug(f"Memory delta: {memory_delta:.2f} MB")
            
            return result
            
        except MemoryError:
            logger.error("Memory limit exceeded during operation")
            self.cleanup_memory()
            raise


class ChunkedDataProcessor:
    """Process large datasets in chunks to manage memory usage."""
    
    def __init__(self, chunk_size: int = 1000, overlap: int = 0):
        """
        Initialize chunked data processor.
        
        Args:
            chunk_size: Size of each chunk
            overlap: Overlap between chunks (for sliding window operations)
        """
        self.chunk_size = chunk_size
        self.overlap = overlap
    
    def chunk_array(self, array: np.ndarray) -> Iterator[Tuple[int, np.ndarray]]:
        """
        Split array into chunks.
        
        Args:
            array: Input array to chunk
            
        Yields:
            Tuple of (start_index, chunk_array)
        """
        n_samples = array.shape[0]
        
        for start in range(0, n_samples, self.chunk_size - self.overlap):
            end = min(start + self.chunk_size, n_samples)
            chunk = array[start:end]
            yield start, chunk
    
    def chunk_matrix(self, matrix: np.ndarray, axis: int = 0) -> Iterator[Tuple[int, np.ndarray]]:
        """
        Split matrix into chunks along specified axis.
        
        Args:
            matrix: Input matrix to chunk
            axis: Axis along which to chunk (0 for rows, 1 for columns)
            
        Yields:
            Tuple of (start_index, chunk_matrix)
        """
        n_elements = matrix.shape[axis]
        
        for start in range(0, n_elements, self.chunk_size):
            end = min(start + self.chunk_size, n_elements)
            
            if axis == 0:
                chunk = matrix[start:end, :]
            elif axis == 1:
                chunk = matrix[:, start:end]
            else:
                raise ValueError("Axis must be 0 or 1 for 2D matrices")
            
            yield start, chunk
    
    def process_chunks_parallel(self, 
                              data: np.ndarray, 
                              process_func: Callable,
                              n_workers: int = None,
                              combine_func: Callable = None) -> Union[List, np.ndarray]:
        """
        Process chunks in parallel.
        
        Args:
            data: Input data to process
            process_func: Function to apply to each chunk
            n_workers: Number of parallel workers
            combine_func: Function to combine results (default: concatenate)
            
        Returns:
            Combined results from all chunks
        """
        if n_workers is None:
            n_workers = min(mp.cpu_count(), 8)
        
        # Generate chunks
        chunks = list(self.chunk_array(data))
        
        # Process chunks in parallel
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(process_func, chunk) for _, chunk in chunks]
            results = [future.result() for future in futures]
        
        # Combine results
        if combine_func is None:
            if isinstance(results[0], np.ndarray):
                return np.concatenate(results, axis=0)
            else:
                return results
        else:
            return combine_func(results)


class SparseMatrixOptimizer:
    """Optimize operations on sparse matrices for memory and speed."""
    
    def __init__(self, density_threshold: float = 0.1):
        """
        Initialize sparse matrix optimizer.
        
        Args:
            density_threshold: Threshold below which to use sparse representation
        """
        self.density_threshold = density_threshold
    
    def optimize_matrix_format(self, matrix: Union[np.ndarray, sparse.spmatrix]) -> sparse.spmatrix:
        """
        Optimize matrix format for sparsity.
        
        Args:
            matrix: Input matrix
            
        Returns:
            Optimized sparse matrix
        """
        if sparse.issparse(matrix):
            density = matrix.nnz / (matrix.shape[0] * matrix.shape[1])
        else:
            density = np.count_nonzero(matrix) / matrix.size
        
        if density < self.density_threshold:
            if sparse.issparse(matrix):
                # Convert to most efficient sparse format
                if matrix.format != 'csr':
                    matrix = matrix.tocsr()
                return matrix
            else:
                return csr_matrix(matrix)
        else:
            # Dense matrix is more efficient
            if sparse.issparse(matrix):
                return matrix.toarray()
            else:
                return matrix
    
    def sparse_dot_product(self, 
                          matrix1: Union[np.ndarray, sparse.spmatrix],
                          matrix2: Union[np.ndarray, sparse.spmatrix]) -> Union[np.ndarray, sparse.spmatrix]:
        """
        Optimized sparse matrix multiplication.
        
        Args:
            matrix1: First matrix
            matrix2: Second matrix
            
        Returns:
            Matrix product
        """
        # Convert to optimal sparse formats
        if not sparse.issparse(matrix1):
            matrix1 = csr_matrix(matrix1)
        if not sparse.issparse(matrix2):
            matrix2 = csc_matrix(matrix2)
        
        # Ensure compatible formats for multiplication
        if matrix1.format != 'csr':
            matrix1 = matrix1.tocsr()
        if matrix2.format != 'csc':
            matrix2 = matrix2.tocsc()
        
        result = matrix1.dot(matrix2)
        
        # Convert back to dense if result is dense enough
        if hasattr(result, 'nnz'):
            density = result.nnz / (result.shape[0] * result.shape[1])
            if density > self.density_threshold:
                result = result.toarray()
        
        return result
    
    def sparse_correlation(self, matrix: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
        """
        Compute correlation matrix efficiently for sparse input.
        
        Args:
            matrix: Input matrix (features x samples)
            
        Returns:
            Correlation matrix
        """
        if not sparse.issparse(matrix):
            matrix = csr_matrix(matrix)
        
        # Center the data
        means = np.array(matrix.mean(axis=1)).flatten()
        matrix_centered = matrix.copy()
        
        # Subtract means (sparse-efficient)
        for i in range(matrix.shape[0]):
            matrix_centered.data[matrix_centered.indptr[i]:matrix_centered.indptr[i+1]] -= means[i]
        
        # Compute correlation
        norms = np.sqrt(np.array((matrix_centered.multiply(matrix_centered)).sum(axis=1)).flatten())
        norms[norms == 0] = 1  # Avoid division by zero
        
        # Normalize
        normalized_matrix = matrix_centered.copy()
        for i in range(matrix.shape[0]):
            normalized_matrix.data[normalized_matrix.indptr[i]:normalized_matrix.indptr[i+1]] /= norms[i]
        
        correlation_matrix = normalized_matrix.dot(normalized_matrix.T)
        
        return correlation_matrix.toarray()


class CacheManager:
    """Manage caching for frequently accessed computations."""
    
    def __init__(self, cache_size: int = 128, cache_dir: str = None):
        """
        Initialize cache manager.
        
        Args:
            cache_size: Maximum number of items to cache in memory
            cache_dir: Directory for disk-based caching
        """
        self.cache_size = cache_size
        self.memory_cache = {}
        self.cache_dir = Path(cache_dir) if cache_dir else Path(tempfile.gettempdir()) / "tfsee_cache"
        self.cache_dir.mkdir(exist_ok=True)
    
    def _get_cache_key(self, *args, **kwargs) -> str:
        """Generate cache key from arguments."""
        key_data = str(args) + str(sorted(kwargs.items()))
        return str(hash(key_data))
    
    def get_from_cache(self, key: str):
        """Retrieve item from cache."""
        # Check memory cache first
        if key in self.memory_cache:
            return self.memory_cache[key]
        
        # Check disk cache
        cache_file = self.cache_dir / f"{key}.pkl"
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    result = pickle.load(f)
                
                # Add to memory cache if there's space
                if len(self.memory_cache) < self.cache_size:
                    self.memory_cache[key] = result
                
                return result
            except:
                logger.warning(f"Failed to load cache file {cache_file}")
        
        return None
    
    def store_in_cache(self, key: str, value):
        """Store item in cache."""
        # Store in memory cache
        if len(self.memory_cache) >= self.cache_size:
            # Remove oldest item (simple FIFO)
            oldest_key = next(iter(self.memory_cache))
            del self.memory_cache[oldest_key]
        
        self.memory_cache[key] = value
        
        # Store in disk cache
        cache_file = self.cache_dir / f"{key}.pkl"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(value, f)
        except:
            logger.warning(f"Failed to save cache file {cache_file}")
    
    def cached_computation(self, func: Callable):
        """Decorator for caching function results."""
        def wrapper(*args, **kwargs):
            key = self._get_cache_key(*args, **kwargs)
            
            # Try to get from cache
            result = self.get_from_cache(key)
            if result is not None:
                logger.debug(f"Cache hit for {func.__name__}")
                return result
            
            # Compute and cache result
            logger.debug(f"Cache miss for {func.__name__}, computing...")
            result = func(*args, **kwargs)
            self.store_in_cache(key, result)
            
            return result
        
        return wrapper
    
    def clear_cache(self):
        """Clear all cache."""
        self.memory_cache.clear()
        
        # Clear disk cache
        for cache_file in self.cache_dir.glob("*.pkl"):
            try:
                cache_file.unlink()
            except:
                logger.warning(f"Failed to delete cache file {cache_file}")


class StreamingDataProcessor:
    """Process data in streaming fashion for very large datasets."""
    
    def __init__(self, buffer_size: int = 10000):
        """
        Initialize streaming data processor.
        
        Args:
            buffer_size: Size of internal buffer for streaming
        """
        self.buffer_size = buffer_size
    
    def stream_from_hdf5(self, file_path: str, dataset_name: str) -> Iterator[np.ndarray]:
        """
        Stream data from HDF5 file.
        
        Args:
            file_path: Path to HDF5 file
            dataset_name: Name of dataset in HDF5 file
            
        Yields:
            Chunks of data
        """
        with h5py.File(file_path, 'r') as f:
            dataset = f[dataset_name]
            n_samples = dataset.shape[0]
            
            for start in range(0, n_samples, self.buffer_size):
                end = min(start + self.buffer_size, n_samples)
                chunk = dataset[start:end]
                yield chunk
    
    def stream_from_csv(self, file_path: str, chunk_size: int = None) -> Iterator[pd.DataFrame]:
        """
        Stream data from CSV file.
        
        Args:
            file_path: Path to CSV file
            chunk_size: Size of chunks (default: use buffer_size)
            
        Yields:
            DataFrame chunks
        """
        if chunk_size is None:
            chunk_size = self.buffer_size
        
        for chunk in pd.read_csv(file_path, chunksize=chunk_size):
            yield chunk
    
    def streaming_correlation(self, data_stream: Iterator[np.ndarray]) -> np.ndarray:
        """
        Compute correlation matrix from streaming data.
        
        Args:
            data_stream: Iterator yielding data chunks
            
        Returns:
            Correlation matrix
        """
        # Online correlation computation using Welford's algorithm
        n_features = None
        n_samples = 0
        mean = None
        covariance = None
        
        for chunk in data_stream:
            if n_features is None:
                n_features = chunk.shape[1]
                mean = np.zeros(n_features)
                covariance = np.zeros((n_features, n_features))
            
            chunk_size = chunk.shape[0]
            
            # Update running statistics
            for sample in chunk:
                n_samples += 1
                delta = sample - mean
                mean += delta / n_samples
                
                # Update covariance matrix
                delta2 = sample - mean
                covariance += np.outer(delta, delta2)
        
        # Convert to correlation
        if n_samples > 1:
            covariance /= (n_samples - 1)
            std_devs = np.sqrt(np.diag(covariance))
            std_devs[std_devs == 0] = 1  # Avoid division by zero
            
            correlation = covariance / np.outer(std_devs, std_devs)
        else:
            correlation = np.eye(n_features)
        
        return correlation


class IncrementalLearning:
    """Implement incremental learning for large-scale data."""
    
    def __init__(self, learning_rate: float = 0.01):
        """
        Initialize incremental learning.
        
        Args:
            learning_rate: Learning rate for updates
        """
        self.learning_rate = learning_rate
        self.model_state = None
    
    def incremental_pca(self, data_chunk: np.ndarray, n_components: int = 10) -> np.ndarray:
        """
        Incremental Principal Component Analysis.
        
        Args:
            data_chunk: New data chunk
            n_components: Number of components to keep
            
        Returns:
            Updated principal components
        """
        if self.model_state is None:
            # Initialize with first chunk
            U, s, Vt = np.linalg.svd(data_chunk, full_matrices=False)
            self.model_state = {
                'components': Vt[:n_components],
                'explained_variance': s[:n_components]**2,
                'n_samples': data_chunk.shape[0]
            }
        else:
            # Update with new chunk
            # Project new data onto existing components
            projected = data_chunk @ self.model_state['components'].T
            
            # Update components using incremental update rule
            for i, component in enumerate(self.model_state['components']):
                gradient = np.mean((projected[:, i:i+1] * data_chunk), axis=0)
                component += self.learning_rate * gradient
                # Re-normalize
                component /= np.linalg.norm(component)
            
            self.model_state['n_samples'] += data_chunk.shape[0]
        
        return self.model_state['components']
    
    def incremental_clustering(self, data_chunk: np.ndarray, n_clusters: int = 10) -> np.ndarray:
        """
        Incremental clustering using online K-means.
        
        Args:
            data_chunk: New data chunk
            n_clusters: Number of clusters
            
        Returns:
            Updated cluster centers
        """
        if self.model_state is None:
            # Initialize cluster centers randomly
            n_features = data_chunk.shape[1]
            self.model_state = {
                'centers': np.random.randn(n_clusters, n_features),
                'counts': np.zeros(n_clusters)
            }
        
        centers = self.model_state['centers']
        counts = self.model_state['counts']
        
        # Assign points to clusters and update centers
        for point in data_chunk:
            # Find closest cluster
            distances = np.linalg.norm(centers - point, axis=1)
            closest_cluster = np.argmin(distances)
            
            # Update cluster center
            counts[closest_cluster] += 1
            learning_rate = 1.0 / counts[closest_cluster]
            centers[closest_cluster] += learning_rate * (point - centers[closest_cluster])
        
        return centers


class PerformanceProfiler:
    """Profile performance of TFSee operations."""
    
    def __init__(self):
        self.timing_data = {}
        self.memory_data = {}
    
    def time_operation(self, operation_name: str):
        """Decorator for timing operations."""
        def decorator(func):
            def wrapper(*args, **kwargs):
                start_time = time.time()
                memory_before = psutil.virtual_memory().used / (1024**2)
                
                result = func(*args, **kwargs)
                
                end_time = time.time()
                memory_after = psutil.virtual_memory().used / (1024**2)
                
                execution_time = end_time - start_time
                memory_delta = memory_after - memory_before
                
                # Store performance data
                if operation_name not in self.timing_data:
                    self.timing_data[operation_name] = []
                    self.memory_data[operation_name] = []
                
                self.timing_data[operation_name].append(execution_time)
                self.memory_data[operation_name].append(memory_delta)
                
                logger.info(f"{operation_name}: {execution_time:.3f}s, {memory_delta:.2f}MB")
                
                return result
            return wrapper
        return decorator
    
    def get_performance_summary(self) -> Dict[str, Dict[str, float]]:
        """Get summary of performance metrics."""
        summary = {}
        
        for operation in self.timing_data:
            times = self.timing_data[operation]
            memories = self.memory_data[operation]
            
            summary[operation] = {
                'mean_time': np.mean(times),
                'std_time': np.std(times),
                'total_time': np.sum(times),
                'mean_memory': np.mean(memories),
                'max_memory': np.max(memories),
                'n_calls': len(times)
            }
        
        return summary


def optimize_tfsee_pipeline(config: OptimizationConfig) -> Dict[str, object]:
    """
    Create optimized components for TFSee pipeline.
    
    Args:
        config: Optimization configuration
        
    Returns:
        Dictionary of optimized components
    """
    # Initialize components
    memory_monitor = MemoryMonitor(config.memory_limit_mb)
    chunk_processor = ChunkedDataProcessor(config.chunk_size)
    sparse_optimizer = SparseMatrixOptimizer()
    cache_manager = CacheManager(config.cache_size, config.temp_directory)
    streaming_processor = StreamingDataProcessor()
    incremental_learner = IncrementalLearning()
    profiler = PerformanceProfiler()
    
    # Set number of workers
    if config.n_workers is None:
        config.n_workers = min(mp.cpu_count(), 8)
    
    logger.info(f"Initialized optimization with {config.n_workers} workers")
    logger.info(f"Memory limit: {config.memory_limit_mb:.0f} MB")
    logger.info(f"Chunk size: {config.chunk_size}")
    
    return {
        'memory_monitor': memory_monitor,
        'chunk_processor': chunk_processor,
        'sparse_optimizer': sparse_optimizer,
        'cache_manager': cache_manager,
        'streaming_processor': streaming_processor,
        'incremental_learner': incremental_learner,
        'profiler': profiler,
        'config': config
    }


def main():
    """Main function for testing optimization components."""
    parser = argparse.ArgumentParser(
        description="TFSee Optimization Utilities"
    )
    parser.add_argument(
        '--test-data', required=True,
        help="Test data file (CSV or HDF5)"
    )
    parser.add_argument(
        '--chunk-size', type=int, default=1000,
        help="Chunk size for processing (default: 1000)"
    )
    parser.add_argument(
        '--n-workers', type=int,
        help="Number of workers (default: auto)"
    )
    parser.add_argument(
        '--memory-limit', type=float,
        help="Memory limit in MB (default: 80% of available)"
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help="Output file for performance results"
    )
    parser.add_argument(
        '--test-operation', 
        choices=['correlation', 'pca', 'clustering', 'all'],
        default='all',
        help="Test operation to perform"
    )
    
    args = parser.parse_args()
    
    # Initialize optimization config
    config = OptimizationConfig(
        chunk_size=args.chunk_size,
        n_workers=args.n_workers,
        memory_limit_mb=args.memory_limit
    )
    
    # Get optimized components
    components = optimize_tfsee_pipeline(config)
    
    # Load test data
    logger.info(f"Loading test data from {args.test_data}")
    if args.test_data.endswith('.h5') or args.test_data.endswith('.hdf5'):
        with h5py.File(args.test_data, 'r') as f:
            dataset_name = list(f.keys())[0]
            data = f[dataset_name][:]
    else:
        data = pd.read_csv(args.test_data, index_col=0).values
    
    logger.info(f"Data shape: {data.shape}")
    
    profiler = components['profiler']
    
    # Test operations
    results = {}
    
    if args.test_operation in ['correlation', 'all']:
        logger.info("Testing correlation computation...")
        
        @profiler.time_operation("standard_correlation")
        def standard_correlation():
            return np.corrcoef(data)
        
        @profiler.time_operation("sparse_correlation")
        def sparse_correlation():
            return components['sparse_optimizer'].sparse_correlation(data)
        
        @profiler.time_operation("chunked_correlation")
        def chunked_correlation():
            def compute_chunk_corr(chunk):
                return np.corrcoef(chunk)
            
            return components['chunk_processor'].process_chunks_parallel(
                data, compute_chunk_corr
            )
        
        results['standard_corr'] = standard_correlation()
        results['sparse_corr'] = sparse_correlation()
        # results['chunked_corr'] = chunked_correlation()
    
    if args.test_operation in ['pca', 'all']:
        logger.info("Testing PCA computation...")
        
        @profiler.time_operation("incremental_pca")
        def incremental_pca():
            components_list = []
            for _, chunk in components['chunk_processor'].chunk_array(data):
                components_arr = components['incremental_learner'].incremental_pca(chunk)
                components_list.append(components_arr)
            return components_list
        
        results['incremental_pca'] = incremental_pca()
    
    if args.test_operation in ['clustering', 'all']:
        logger.info("Testing clustering...")
        
        @profiler.time_operation("incremental_clustering")
        def incremental_clustering():
            centers_list = []
            for _, chunk in components['chunk_processor'].chunk_array(data):
                centers = components['incremental_learner'].incremental_clustering(chunk)
                centers_list.append(centers)
            return centers_list
        
        results['incremental_clustering'] = incremental_clustering()
    
    # Get performance summary
    performance_summary = profiler.get_performance_summary()
    
    # Save results
    logger.info(f"Saving results to {args.output}")
    with open(args.output, 'w') as f:
        json.dump(performance_summary, f, indent=2)
    
    # Print summary
    logger.info("Performance Summary:")
    for operation, metrics in performance_summary.items():
        logger.info(f"{operation}:")
        logger.info(f"  Mean time: {metrics['mean_time']:.3f}s")
        logger.info(f"  Peak memory: {metrics['max_memory']:.2f}MB")
        logger.info(f"  Calls: {metrics['n_calls']}")
    
    # Memory usage summary
    memory_profile = components['memory_monitor'].get_memory_usage()
    logger.info(f"Peak memory usage: {memory_profile.peak_memory_mb:.2f}MB")
    logger.info(f"Current memory usage: {memory_profile.current_memory_mb:.2f}MB")


if __name__ == "__main__":
    main()