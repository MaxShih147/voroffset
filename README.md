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

## Future Work

1. [ ] **Dexelization Replacement:**  
  Implement custom dexelization functions (or integrate another library) to replace the Geogram-dependent functionality.
  
2. [ ] **API Development:**  
  Design and document the final API for users who will integrate this module with a web front-end.
  
3. [ ] **WebAssembly Compilation:**  
  Finalize the build process for compiling this project into a WebAssembly module to be used in the resin slicer application.

## Development Environment

- **Platform:** macOS (MacBook Air 2024, M3 chip)
- **Build Tools:** CMake, Visual Studio Code
- **Target:** WebAssembly for front-end integration in a light-curing web resin slicer project