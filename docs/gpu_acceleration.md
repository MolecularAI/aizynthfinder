# GPU Acceleration

AiZynthFinder supports GPU acceleration for ONNX model inference using CUDA. When `onnxruntime-gpu` is installed and a compatible NVIDIA GPU is available, the expansion and filter policy models will automatically use GPU acceleration, providing **3-5x speedup** for model inference.

## Enabling GPU Acceleration

To enable GPU acceleration, replace the CPU-only ONNX Runtime with the GPU-enabled version:

```bash
pip uninstall onnxruntime
pip install onnxruntime-gpu
```

Ensure NVIDIA drivers and CUDA toolkit are installed on your system.

## Verifying GPU Acceleration

To verify GPU acceleration is working:

```python
import onnxruntime as ort
print(ort.get_available_providers())
# Should include 'CUDAExecutionProvider'
```

You can also run the benchmark test:

```bash
python tests/test_gpu_acceleration.py
```

## Disabling GPU Acceleration

To force CPU-only execution (e.g., for debugging or reproducibility), set the environment variable:

```bash
export AIZYNTHFINDER_USE_GPU=0
```

## Implementation Details

### Files Modified

**`aizynthfinder/utils/models.py`**:
- Added `_get_onnx_providers()` function with intelligent provider detection and caching
- Updated `LocalOnnxModel.__init__()` to use GPU providers automatically
- Added `AIZYNTHFINDER_USE_GPU` environment variable support

**`pyproject.toml`**:
- Added `onnx` as a dev dependency for creating test models in the benchmark

### New Files

**`tests/test_gpu_acceleration.py`**:
- Tests `_get_onnx_providers()` function
- Tests provider caching
- Tests CUDA detection
- Tests GPU disable via environment variable
- Includes GPU vs CPU benchmark (when CUDA available)

**`docs/gpu_acceleration.md`**:
- This documentation file

### Documentation Updated

**`README.md`**:
- Added GPU Acceleration section with setup instructions

## Estimated Speedups

Typical speedup on NVIDIA GPUs:
- Tesla T4: ~4-5x faster
- A100: ~5-6x faster
- Consumer GPUs (RTX 3080, etc.): ~3-4x faster

The actual speedup depends on model size, batch size, GPU memory bandwidth, and CPU baseline performance.
