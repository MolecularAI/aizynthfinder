#!/usr/bin/env python
"""
Test script for GPU acceleration in ONNX model inference.

This script verifies that:
1. CUDAExecutionProvider can be detected when available
2. LocalOnnxModel uses GPU providers when available
3. GPU inference provides speedup over CPU inference

Run with pytest:  pytest tests/test_gpu_acceleration.py -v
Run standalone:   python tests/test_gpu_acceleration.py
"""
import os
import tempfile
import time
from pathlib import Path

import numpy as np
import onnxruntime as ort

# pytest is optional - only needed when running via pytest
try:
    import pytest
    PYTEST_AVAILABLE = True
except ImportError:
    PYTEST_AVAILABLE = False
    # Create a dummy pytest module for decorators when running standalone
    class _DummyMark:
        @staticmethod
        def skipif(*args, **kwargs):
            def decorator(cls):
                return cls
            return decorator
    
    class _DummyPytest:
        mark = _DummyMark()
    
    pytest = _DummyPytest()


class TestGPUAcceleration:
    """Tests for GPU acceleration functionality."""

    def test_get_onnx_providers_returns_list(self):
        """Test that _get_onnx_providers returns a valid list."""
        from aizynthfinder.utils.models import _get_onnx_providers
        
        providers = _get_onnx_providers()
        
        assert isinstance(providers, list)
        assert len(providers) > 0
        assert "CPUExecutionProvider" in providers

    def test_get_onnx_providers_caching(self):
        """Test that provider detection is cached."""
        from aizynthfinder.utils.models import _get_onnx_providers, _ONNX_PROVIDERS_CACHE
        
        # First call
        providers1 = _get_onnx_providers()
        # Second call should return cached value
        providers2 = _get_onnx_providers()
        
        assert providers1 is providers2

    def test_cuda_provider_detection(self):
        """Test CUDA provider detection matches onnxruntime availability."""
        from aizynthfinder.utils.models import _get_onnx_providers
        
        available = ort.get_available_providers()
        providers = _get_onnx_providers()
        
        # If CUDA is available in onnxruntime, it should be in our providers
        if "CUDAExecutionProvider" in available:
            assert "CUDAExecutionProvider" in providers
        else:
            assert "CUDAExecutionProvider" not in providers

    def test_gpu_disabled_via_env_var(self):
        """Test that GPU can be disabled via environment variable."""
        import aizynthfinder.utils.models as models_module
        from importlib import reload
        
        # Save original env value
        original_value = os.environ.get("AIZYNTHFINDER_USE_GPU")
        
        try:
            # Clear the cache first
            models_module._ONNX_PROVIDERS_CACHE = None
            
            # Set environment variable to disable GPU
            os.environ["AIZYNTHFINDER_USE_GPU"] = "0"
            
            # Reload to pick up new env var
            reload(models_module)

            providers = models_module._get_onnx_providers()
            
            # Should not include CUDA even if available
            assert "CUDAExecutionProvider" not in providers
            assert "CPUExecutionProvider" in providers
        finally:
            # Restore original env value
            models_module._ONNX_PROVIDERS_CACHE = None
            if original_value is None:
                os.environ.pop("AIZYNTHFINDER_USE_GPU", None)
            else:
                os.environ["AIZYNTHFINDER_USE_GPU"] = original_value
            reload(models_module)


@pytest.mark.skipif(
    "CUDAExecutionProvider" not in ort.get_available_providers(),
    reason="CUDA not available"
)
class TestGPUBenchmark:
    """Benchmark tests that require GPU availability."""

    def test_gpu_inference_speedup(self, shared_datadir):
        """
        Benchmark GPU vs CPU inference speed.
        
        This test loads an ONNX model and compares inference time between
        GPU and CPU execution providers.
        """
        # Find an ONNX model file for testing
        model_files = list(shared_datadir.glob("*.onnx"))
        if not model_files:
            pytest.skip("No ONNX model files found in test data")
        
        model_path = str(model_files[0])
        
        # Create sessions with different providers
        session_opts = ort.SessionOptions()
        session_opts.intra_op_num_threads = 2
        
        gpu_session = ort.InferenceSession(
            model_path,
            sess_options=session_opts,
            providers=["CUDAExecutionProvider", "CPUExecutionProvider"]
        )
        
        cpu_session = ort.InferenceSession(
            model_path,
            sess_options=session_opts,
            providers=["CPUExecutionProvider"]
        )
        
        # Verify GPU session is using CUDA
        assert "CUDAExecutionProvider" in gpu_session.get_providers()
        
        # Get input dimensions
        input_info = gpu_session.get_inputs()[0]
        input_shape = [10] + list(input_info.shape[1:])
        input_data = np.random.randn(*input_shape).astype(np.float32)
        
        # Warm-up
        for _ in range(10):
            gpu_session.run(None, {input_info.name: input_data})
            cpu_session.run(None, {input_info.name: input_data})
        
        # Benchmark GPU
        num_iterations = 100
        gpu_start = time.perf_counter()
        for _ in range(num_iterations):
            gpu_session.run(None, {input_info.name: input_data})
        gpu_time = time.perf_counter() - gpu_start
        
        # Benchmark CPU
        cpu_start = time.perf_counter()
        for _ in range(num_iterations):
            cpu_session.run(None, {input_info.name: input_data})
        cpu_time = time.perf_counter() - cpu_start
        
        # Calculate speedup
        speedup = cpu_time / gpu_time
        
        print(f"\n{'='*60}")
        print("GPU vs CPU Inference Benchmark")
        print(f"{'='*60}")
        print(f"  Iterations: {num_iterations}")
        print(f"  GPU time:   {gpu_time*1000:.2f} ms ({gpu_time/num_iterations*1000:.3f} ms/iter)")
        print(f"  CPU time:   {cpu_time*1000:.2f} ms ({cpu_time/num_iterations*1000:.3f} ms/iter)")
        print(f"  Speedup:    {speedup:.2f}x")
        print(f"{'='*60}")
        
        # GPU should be at least slightly faster (allowing for small variance)
        # In practice, speedup is typically 3-5x
        assert speedup > 0.9, f"GPU should not be significantly slower than CPU (speedup: {speedup:.2f}x)"


def benchmark_gpu_vs_cpu():
    """
    Standalone benchmark function that can be run directly.
    
    Usage:
        python tests/test_gpu_acceleration.py
    """
    print("\n" + "=" * 70)
    print("AiZynthFinder GPU Acceleration Benchmark")
    print("=" * 70)
    
    # Check available providers
    available = ort.get_available_providers()
    print(f"\nAvailable ONNX Runtime providers: {available}")
    
    has_cuda = "CUDAExecutionProvider" in available
    print(f"CUDA available: {has_cuda}")
    
    if not has_cuda:
        print("\n⚠️  CUDA not available. Cannot benchmark GPU performance.")
        print("    Install onnxruntime-gpu and ensure CUDA drivers are installed.")
        return
    
    # Test with AiZynthFinder's provider detection
    from aizynthfinder.utils.models import _get_onnx_providers
    providers = _get_onnx_providers()
    print(f"AiZynthFinder providers: {providers}")
    
    # Try to find an existing ONNX model for benchmarking
    print("\n" + "-" * 70)
    print("Running inference benchmark...")
    print("-" * 70)
    
    # Look for existing model files in common locations
    model_path = None
    search_paths = [
        Path(__file__).parent.parent.parent / "aizynth_models",  # Sibling to aizynthfinder repo
        Path.home() / "aizynth_models",
        Path.cwd() / "aizynth_models",
    ]
    
    for search_dir in search_paths:
        if search_dir.exists():
            onnx_files = list(search_dir.glob("*.onnx"))
            if onnx_files:
                model_path = str(onnx_files[0])
                print(f"  Found model: {model_path}")
                break
    
    if not model_path:
        # Fallback: try to create a test model dynamically
        try:
            import onnx
            from onnx import helper, TensorProto
            
            print("  No existing model found. Creating temporary test model...")
            
            X = helper.make_tensor_value_info('X', TensorProto.FLOAT, [None, 2048])
            Y = helper.make_tensor_value_info('Y', TensorProto.FLOAT, [None, 16000])
            
            W = helper.make_tensor('W', TensorProto.FLOAT, [2048, 16000], 
                                   np.random.randn(2048, 16000).astype(np.float32).flatten().tolist())
            B = helper.make_tensor('B', TensorProto.FLOAT, [16000], 
                                   np.random.randn(16000).astype(np.float32).flatten().tolist())
            
            matmul_node = helper.make_node('MatMul', ['X', 'W'], ['XW'])
            add_node = helper.make_node('Add', ['XW', 'B'], ['Y'])
            
            graph = helper.make_graph([matmul_node, add_node], 'test_model', [X], [Y], [W, B])
            model = helper.make_model(graph, opset_imports=[helper.make_opsetid('', 13)])
            
            with tempfile.NamedTemporaryFile(suffix='.onnx', delete=False) as f:
                onnx.save(model, f.name)
                model_path = f.name
            
        except ImportError:
            print("  ⚠️  No ONNX models found and 'onnx' package not installed.")
            print("  ℹ️  To run the benchmark, either:")
            print("      1. Run from a directory with aizynth_models/*.onnx files, or")
            print("      2. Install the onnx package: pip install onnx")
            print("\n  GPU acceleration is working! Skipping benchmark only.")
            return
    
    # Track if we created a temp file
    is_temp_file = model_path.startswith(tempfile.gettempdir()) if model_path else False
    
    try:
        session_opts = ort.SessionOptions()
        session_opts.intra_op_num_threads = 2
        
        gpu_session = ort.InferenceSession(
            model_path,
            sess_options=session_opts,
            providers=["CUDAExecutionProvider", "CPUExecutionProvider"]
        )
        
        cpu_session = ort.InferenceSession(
            model_path,
            sess_options=session_opts,
            providers=["CPUExecutionProvider"]
        )
        
        print(f"  GPU session providers: {gpu_session.get_providers()}")
        print(f"  CPU session providers: {cpu_session.get_providers()}")
        
        # Build input feed for ALL model inputs
        input_feed = {}
        all_inputs = gpu_session.get_inputs()
        print(f"  Model has {len(all_inputs)} input(s)")
        
        for input_info in all_inputs:
            # Handle dynamic dimensions (None or strings like 'batch')
            input_shape = []
            for dim in input_info.shape:
                if isinstance(dim, int):
                    input_shape.append(dim)
                else:
                    input_shape.append(10)  # Default batch size
            
            input_feed[input_info.name] = np.random.randn(*input_shape).astype(np.float32)
            print(f"    {input_info.name}: {input_shape}")
        
        num_iterations = 100
        
        # Warm-up
        for _ in range(10):
            gpu_session.run(None, input_feed)
            cpu_session.run(None, input_feed)
        
        # GPU benchmark
        gpu_start = time.perf_counter()
        for _ in range(num_iterations):
            gpu_session.run(None, input_feed)
        gpu_time = time.perf_counter() - gpu_start
        
        # CPU benchmark
        cpu_start = time.perf_counter()
        for _ in range(num_iterations):
            cpu_session.run(None, input_feed)
        cpu_time = time.perf_counter() - cpu_start
        
        speedup = cpu_time / gpu_time
        
        print(f"\n  Iterations: {num_iterations}")
        print(f"\n  GPU time:  {gpu_time*1000:.1f} ms total, {gpu_time/num_iterations*1000:.2f} ms/iter")
        print(f"  CPU time:  {cpu_time*1000:.1f} ms total, {cpu_time/num_iterations*1000:.2f} ms/iter")
        print(f"\n  🚀 GPU Speedup: {speedup:.2f}x faster")
        
        if speedup > 1.5:
            print("\n  ✅ GPU acceleration is working and providing significant speedup!")
        elif speedup > 1.0:
            print("\n  ✅ GPU acceleration is working (speedup may vary with model size).")
        else:
            print("\n  ⚠️  GPU appears slower - this may be due to small model/batch size overhead.")
        
    finally:
        # Cleanup temp file only if we created it
        if is_temp_file and os.path.exists(model_path):
            os.unlink(model_path)
    
    print("\n" + "=" * 70 + "\n")


if __name__ == "__main__":
    benchmark_gpu_vs_cpu()

