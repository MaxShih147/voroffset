////////////////////////////////////////////////////////////////////////////////
// offset3d.cpp
//
// This file is the main entry point for mesh morphing using the vor3d module.
// It reads an input STL model, performs a morphological operation (e.g., dilation),
// and writes the output STL model. The STL read/write functions are integrated
// directly in this file for simplicity.
// 
// Author: Max Shih
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

// Include vor3d module headers (implementation in src/vor3d)
#include <vor3d/CompressedVolume.h>
#include <vor3d/Dexelize.h>
#include <vor3d/Logger.h>
#include <vor3d/VoronoiVorPower.h>
#include <vor3d/VoronoiBruteForce.h>
#include <vor3d/Timer.h>

#include <vor3d/MorphologyProcessor.h>

// -----------------------------------------------------------------------------
// STL IO Functions
// -----------------------------------------------------------------------------

/**
 * @brief Reads a binary STL file.
 * 
 * This function opens a binary STL file, reads the header, and then extracts
 * the triangle data. Each triangle is stored as 3 vertices (x, y, z interleaved)
 * and no vertex deduplication is performed (each triangle gets its own vertices).
 * The resulting data is stored in outVertices (as doubles) and outIndices
 * (3 indices per triangle, where each triangle's vertices are consecutive).
 *
 * @param filename The input STL file path.
 * @param outVertices [output] Flat vertex array (x, y, z interleaved).
 * @param outIndices [output] Triangle index array (3 indices per triangle).
 * @return true if reading was successful, false otherwise.
 */
bool ReadSTL(const std::string &filename, std::vector<float>& outVertices, std::vector<unsigned int>& outIndices) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Error: Cannot open STL file: " << filename << std::endl;
        return false;
    }

    // Read 80-byte header (ignored)
    char header[80] = {0};
    in.read(header, 80);
    if (in.gcount() != 80) {
        std::cerr << "Error: Failed to read STL header." << std::endl;
        return false;
    }

    // Read triangle count (4 bytes)
    uint32_t triangleCount = 0;
    in.read(reinterpret_cast<char*>(&triangleCount), sizeof(uint32_t));
    if (!in) {
        std::cerr << "Error: Failed to read triangle count." << std::endl;
        return false;
    }

    // Resize output arrays
    // Each triangle has 3 vertices, each vertex has 3 coordinates.
    outVertices.resize(triangleCount * 3 * 3);
    outIndices.resize(triangleCount * 3);

    // Each triangle is 50 bytes: 12 bytes for normal, 36 bytes for vertices, 2 bytes for attribute count.
    for (uint32_t i = 0; i < triangleCount; ++i) {
        // Read normal (ignored)
        float normal[3];
        in.read(reinterpret_cast<char*>(normal), 3 * sizeof(float));

        // Read triangle vertices (9 floats)
        float triangleVerts[9];
        in.read(reinterpret_cast<char*>(triangleVerts), 9 * sizeof(float));
        if (!in) {
            std::cerr << "Error: Failed to read triangle vertices." << std::endl;
            return false;
        }

        // Read attribute byte count (2 bytes, ignored)
        uint16_t attributeByteCount = 0;
        in.read(reinterpret_cast<char*>(&attributeByteCount), sizeof(uint16_t));

        // Copy vertices (convert to double)
        for (int j = 0; j < 9; ++j) {
            outVertices[i * 9 + j] = triangleVerts[j];
        }

        // Set triangle indices (each triangle gets new vertices)
        outIndices[i * 3 + 0] = i * 3 + 0;
        outIndices[i * 3 + 1] = i * 3 + 1;
        outIndices[i * 3 + 2] = i * 3 + 2;
    }

    return true;
}

/**
 * @brief Writes a binary STL file.
 * 
 * This function writes a binary STL file using the provided vertex and index arrays.
 * It assumes that the vertex array is flat (x, y, z interleaved) and that the index
 * array contains triangles (3 indices per triangle). The normal for each triangle is
 * computed using the cross product of two edges.
 *
 * @param filename The output STL file path.
 * @param vertices The flat vertex array (x, y, z interleaved).
 * @param indices The triangle index array (3 indices per triangle).
 * @return true if writing was successful, false otherwise.
 */
bool WriteSTL(const std::string &filename, const std::vector<float>& vertices, const std::vector<unsigned int>& indices) {
	// Check if the sizes are valid (each triangle should have 3 indices and each vertex has 3 coordinates)
    if (indices.size() % 3 != 0 || vertices.size() % 3 != 0) {
        std::cerr << "Error: Invalid vertices or indices size." << std::endl;
        return false;
    }
    
    std::ofstream out(filename, std::ios::binary | std::ios::out);
    if (!out) {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        return false;
    }
    
    // Write 80-byte header (fill with zeros)
    char header[80] = {0};
    out.write(header, 80);
    
    // Write triangle count (4 bytes)
    uint32_t triangleCount = static_cast<uint32_t>(indices.size() / 3);
    out.write(reinterpret_cast<const char*>(&triangleCount), sizeof(uint32_t));
    
    // For each triangle, compute the normal and write triangle data.
    for (uint32_t i = 0; i < triangleCount; ++i) {
        // Get vertex indices for triangle i
        unsigned int idx0 = indices[i * 3 + 0];
        unsigned int idx1 = indices[i * 3 + 1];
        unsigned int idx2 = indices[i * 3 + 2];
        
        // Extract vertices from the flat array and convert to float
        float v0[3], v1[3], v2[3];
        for (int j = 0; j < 3; ++j) {
            v0[j] = vertices[idx0 * 3 + j];
            v1[j] = vertices[idx1 * 3 + j];
            v2[j] = vertices[idx2 * 3 + j];
        }
        
        // Compute normal using cross product of (v1 - v0) and (v2 - v0)
        float e1[3] = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
        float e2[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };
        float normal[3] = { e1[1] * e2[2] - e1[2] * e2[1],
                            e1[2] * e2[0] - e1[0] * e2[2],
                            e1[0] * e2[1] - e1[1] * e2[0] };
        
        // Normalize the normal vector
        float norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        if (norm > 0) {
            normal[0] /= norm;
            normal[1] /= norm;
            normal[2] /= norm;
        } else {
            normal[0] = normal[1] = normal[2] = 0;
        }
        
        // Write normal (3 floats)
        out.write(reinterpret_cast<const char*>(normal), 3 * sizeof(float));
        // Write vertices (each triangle has 3 vertices, each vertex has 3 floats)
        out.write(reinterpret_cast<const char*>(v0), 3 * sizeof(float));
        out.write(reinterpret_cast<const char*>(v1), 3 * sizeof(float));
        out.write(reinterpret_cast<const char*>(v2), 3 * sizeof(float));
        // Write attribute byte count (2 bytes, usually 0)
        uint16_t attributeByteCount = 0;
        out.write(reinterpret_cast<const char*>(&attributeByteCount), sizeof(uint16_t));
    }
    
    return true;
}

// -----------------------------------------------------------------------------
// Main function
// -----------------------------------------------------------------------------

/**
 * @brief Main entry point for mesh morphing.
 * 
 * This function reads an input STL file, performs a morphological operation
 * (such as dilation) using the vor3d module, and writes the result to an output
 * STL file. It uses the integrated STL IO functions.
 */
int main() {

    LOG_TIME();

    std::string inputFilename = "input/01.stl";
    std::string outputFilename = "output/result.stl";

    // Load input mesh
    std::vector<float> inVertices;
    std::vector<unsigned int> inIndices;
    if (!ReadSTL(inputFilename, inVertices, inIndices)) {
        LOG_FAIL("Failed to read input STL...");
        return 1;
    }

    // Set parameters
    morpho3d::MorphologyParams params;
    // Type of morphological operation: options include "dilation", "erosion", "closing", "opening", "noop"
    params.operation = "erosion";
    // Size of each dexel (grid cell) in world units (e.g., mm)
    params.dexelSize = 1.0;    
    // Radius of morphological effect (in dexel units unless radiusInMM = true)
    params.radius = 8.0;   
    // Whether the radius is specified in millimeters (true) or dexel units (false)
    params.radiusInMM = false;   
    // Target resolution: number of dexels along the longest axis
    params.numDexels = 256;

    // Create processor
    morpho3d::MorphologyProcessor* processor = morpho3d::MorphologyProcessor::Create(params);

    // Output buffers
    float* outVertices = nullptr;
    unsigned int* outIndices = nullptr;
    size_t numOutVertices = 0, numOutIndices = 0;

    // Run morphology
    if (!processor->Run(
        inVertices.data(), inVertices.size(),
        inIndices.data(), inIndices.size(),
        &outVertices, &numOutVertices,
        &outIndices, &numOutIndices
    )) {
        LOG_FAIL("Morphological operation failed...");
        morpho3d::MorphologyProcessor::Delete(processor);
        return 1;
    }

    // Convert result to vectors for STL output
    std::vector<float> finalVerts(outVertices, outVertices + numOutVertices);
    std::vector<unsigned int> finalInds(outIndices, outIndices + numOutIndices);

    delete[] outVertices;
    delete[] outIndices;

    if (!WriteSTL(outputFilename, finalVerts, finalInds)) {
        LOG_FAIL("Failed to write output STL...");
        morpho3d::MorphologyProcessor::Delete(processor);
        return 1;
    }

    morpho3d::MorphologyProcessor::Delete(processor);
    LOG_PASS("Mesh morph operation completed.");
    return 0;
}