# Voroffset for Web Resin Slicer

*Based on the work of [geometryprocessing/voroffset](https://github.com/geometryprocessing/voroffset)*

**Forked and modified for WebAssembly integration on macOS (MacBook Air 2024, M3 chip).**

## Overview

This repository is a fork of the original Voroffset project. The goal of this fork is to adapt the 3D offset (hollowing) algorithm for use in a web-based light-curing resin slicer application. Our primary targets are:

- **WebAssembly Integration:** Produce a lightweight module for front-end usage.
- **Simplified External Dependencies:** Reduce or remove external libraries (e.g., Geogram) to lower the overhead when compiling to WebAssembly.
- **Clean API Interface:** Replace the original CLI interface with a structured API for better integration.

## Key Modifications

- **Build Configuration for macOS (M3 chip):**  
  The CMake files have been modified to disable external dependencies that fail to build on macOS—particularly Geogram. As a result, the dexelization functionality is currently stubbed out:
  - Functions such as `voroffset3d::CompressedVolume`, `voroffset3d::create_dexels`, and `voroffset3d::dexel_dump` (and related helpers in `dexelize.cpp`) are empty.
  - These will be replaced later with either alternative libraries or custom code.

- **API Refactoring:**  
  The original project used a CLI (defined in `main`) to allow user interaction. In this fork, the CLI is removed in favor of designing a dedicated API. This API will be the interface through which the front-end (via WebAssembly) interacts with the 3D offset functionality.

- **Minimized External Dependencies:**  
  We have intentionally reduced external dependencies to ease the future WebAssembly build process. Our focus is on maintaining the core 3D offset algorithm while providing our own implementations for parts that originally depended on Geogram.

## WebAssembly Build Instructions

This section documents how to set up and build the WebAssembly module (`MorphologyWASM`) on macOS.

### 1. Install Emscripten

Clone and install the Emscripten SDK:

```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

To persist the environment in every terminal session, add the following to your `~/.zshrc`:

```bash
source /add your path here.../emsdk/emsdk_env.sh
```

(Adjust the path based on where you cloned the SDK)

---

### 2. Build the WASM module

From the root of the `MorphologyWASM` project:

```bash
cd app/demos/MorphologyWASM
emcmake cmake -B build
cmake --build build
```

> The build output will be in: `build/dist/`  
> It includes:
> - `morphology.js` (JavaScript glue code)
> - `morphology.wasm` (WebAssembly binary)
> - `index.html` and `test.js` (copied automatically after build)

---

### 3. Run in a local server

You must serve the output via a web server (due to browser WASM restrictions):

```bash
cd build/dist
python3 -m http.server
```

Then open your browser to:

```
http://localhost:8000/index.html
```

If everything works correctly, you will see log output confirming that the WebAssembly module was loaded and executed.

## TODOs

1. [x] **Dexelization Replacement:**  
  Implement custom dexelization functions (or integrate another library) to replace the Geogram-dependent functionality.
  
2. [x] **API Development:**  
  Design and document the final API for users who will integrate this module with a web front-end.
  
3. [ ] **WebAssembly Compilation:**  
  Finalize the build process for compiling this project into a WebAssembly module to be used in the resin slicer application.
  
4. [ ] **Performance Optimization:**  
  Improve the performance and quality of the output meshes through multiple strategies:
  - **Memory Efficiency:** Refactor data structures and algorithms to reduce allocation overhead and minimize unnecessary copies.
  - **Parallel Processing:** Utilize multithreading or SIMD where feasible, especially during voxelization and morphological operations.
  - **Mesh Simplification:** Implement triangle decimation techniques to reduce polygon count while preserving shape fidelity.
  - **Surface Smoothing:** Apply algorithms such as Laplacian or Taubin smoothing to eliminate voxel-induced surface noise and generate more natural surfaces.
  - **GPU Acceleration (Experimental):** Explore offloading key stages to WebGL/WebGPU compute shaders to improve performance in browser environments.
  - **Web-Friendly Output:** Optimize resulting meshes for WebGL by balancing polygon count and visual accuracy for real-time rendering.

## Development Environment

- **Platform:** macOS (MacBook Air 2024, M3 chip)
- **Build Tools:** CMake, Visual Studio Code
- **Target:** WebAssembly for front-end integration in a light-curing web resin slicer project