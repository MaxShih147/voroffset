cmake_minimum_required(VERSION 3.13)
project(MorphologyWASM)

set(CMAKE_CXX_STANDARD 17)

# Set output directory to dist/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist)

add_executable(morphology_wasm
    ../../../src/vor3d/MorphologyProcessor.cpp
    ../../../src/vor3d/AABB.cpp
    ../../../src/vor3d/Common.cpp
    ../../../src/vor3d/CompressedVolume.cpp
    ../../../src/vor3d/CompressedVolumeBase.cpp
    ../../../src/vor3d/CompressedVolumeWithRadii.cpp
    ../../../src/vor3d/Dexelize.cpp
    ../../../src/vor3d/HalfDilationOperator.cpp
    ../../../src/vor3d/Morpho2D.cpp
    ../../../src/vor3d/MorphologyOperators.cpp
    ../../../src/vor3d/MorphologyProcessor.cpp
    ../../../src/vor3d/SeparatePower2D.cpp
    ../../../src/vor3d/Timer.cpp
    ../../../src/vor3d/Voronoi.cpp
    ../../../src/vor3d/Voronoi2D.cpp
    ../../../src/vor3d/VoronoiBruteForce.cpp
    ../../../src/vor3d/VoronoiVorPower.cpp
    MorphologyCAPIImpl.cpp
    dummy_main.cpp
)

# Optional: include headers if needed
include_directories(../../../src)
include_directories(../../../external/eigen)

# Emscripten output: .wasm with .js wrapper
set_target_properties(morphology_wasm PROPERTIES
    OUTPUT_NAME "morphology"
    SUFFIX ".js"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist
)

# Allow outputting HTML if needed
set(CMAKE_EXECUTABLE_SUFFIX ".html")

# 清掉舊的 CXX_FLAGS linker 設定，避免 warning
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# 正確將 linker 專用 flags 放在這
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} \
  -s EXPORTED_FUNCTIONS=['_malloc','_free','_Morphology_Create','_Morphology_Run','_Morphology_Delete'] \
  -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap'] \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME='MorphologyModule' \
  -s STANDALONE_WASM=0 \
")

# Copy index.html and test.js into dist/ after build
add_custom_command(TARGET morphology_wasm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/index.html
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/index.html
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/test.js
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test.js
    COMMENT "Copying index.html and test.js into dist/"
)