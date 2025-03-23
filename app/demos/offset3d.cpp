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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <memory>

// Include vor3d module headers (implementation in src/vor3d)
#include <vor3d/CompressedVolume.h>
#include <vor3d/VoronoiVorPower.h>
#include <vor3d/VoronoiBruteForce.h>
#include <vor3d/Dexelize.h>
#include <vor3d/Timer.h>

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

	std::string inputFilename = "input/01.stl";
	std::string outputFilename = "output/result.stl";

	// Read the input STL model into flat vertex array and triangle index array.
    std::vector<float> inVertices;
    std::vector<unsigned int> inIndices;
    if (!ReadSTL(inputFilename, inVertices, inIndices)) {
        std::cerr << "Error: Failed to read STL file: " << inputFilename << std::endl;
        return 1;
    }
    
    // -------------------------------------------------------------------------    
    // Convert the STL mesh into a CompressedVolume for dexel processing.
    // TODO: Replace this placeholder with your actual conversion function.
    vor3d::CompressedVolume inputVolume;
	double dexelsSize = 1.0;
	double padding = 0.0;
	int numDexels = 256;
    inputVolume = voroffset3d::CreateDexelsFromMeshBuffers(inVertices, inIndices, dexelsSize, padding, numDexels);
    
    // Set default parameters (adjust as necessary)
    double radius = 8.0;
    bool radiusInMM = false;
    
    // Convert radius from mm to dexel units if needed.
    if (radiusInMM) {
        radius /= inputVolume.spacing();
    }
    
	// Create the morphological (offset) operator based on the chosen method.
	std::unique_ptr<vor3d::VoronoiMorpho> offsetOp;
	std::string method = "ours";  // Options: "ours" or "brute_force"
	if (method == "ours") {
		offsetOp = std::make_unique<vor3d::VoronoiMorphoVorPower>();
	} else if (method == "brute_force") {
		offsetOp = std::make_unique<vor3d::VoronoiMorphoBruteForce>();
	} else {
		std::cerr << "Error: Invalid method: " << method << std::endl;
		return 1;
	}
	if (!offsetOp) {
		std::cerr << "Error: Failed to create offset operator." << std::endl;
		return 1;
	}

	// Apply the morphological operation (example: dilation)
	std::string operation = "dilation"; // Options: "noop", "erosion", "closing", "opening"
	double timeFirst = 0, timeSecond = 0;
	vor3d::CompressedVolume outputVolume;
	if (operation == "noop") {
		outputVolume = inputVolume;
	} else if (operation == "erosion") {
		offsetOp->erosion(inputVolume, outputVolume, radius, timeFirst, timeSecond);
	} else if (operation == "dilation") {
		offsetOp->dilation(inputVolume, outputVolume, radius, timeFirst, timeSecond);
	} else if (operation == "closing") {
		vor3d::CompressedVolume tmpVolume;
		offsetOp->dilation(inputVolume, tmpVolume, radius, timeFirst, timeSecond);
		offsetOp->erosion(tmpVolume, outputVolume, radius, timeFirst, timeSecond);
	} else if (operation == "opening") {
		vor3d::CompressedVolume tmpVolume;
		offsetOp->erosion(inputVolume, tmpVolume, radius, timeFirst, timeSecond);
		offsetOp->dilation(tmpVolume, outputVolume, radius, timeFirst, timeSecond);
	} else {
		std::cerr << "Error: Invalid operation: " << operation << std::endl;
		return 1;
	}

	// TODO: Convert the output CompressedVolume back to STL mesh arrays.
	// For now, we use placeholder empty vectors.
	std::vector<float> outVertices;         // Output vertex array (x, y, z interleaved)
	std::vector<unsigned int> outIndices;      // Output triangle index array (3 indices per triangle)

	voroffset3d::DumpDexelsIntoMeshBuffers(outputVolume, outVertices, outIndices);

	// Write the output STL model.
	if (!WriteSTL(outputFilename, outVertices, outIndices)) {
		std::cerr << "Error: Failed to write STL file: " << outputFilename << std::endl;
		return 1;
	}

	std::cout << "Mesh morph operation completed successfully!" << std::endl;
	return 0;
}