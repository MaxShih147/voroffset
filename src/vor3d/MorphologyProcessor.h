#pragma once

#include <string>
#include <vector>
#include <memory>
#include <vor3d/CompressedVolume.h>
// #include <vor3d/Dexelize.h>
#include <vor3d/VoronoiVorPower.h>
// #include <vor3d/VoronoiBruteForce.h>
// #include <vor3d/Timer.h>

namespace morpho3d {

/**
 * @brief Parameters for configuring the 3D morphological operation.
 *
 * This structure contains all configurable options for dexelization, 
 * offsetting, and mesh generation.
 */
struct MorphologyParams {
    /**
     * @brief Type of morphological operation.
     * Options: "noop", "dilation", "erosion", "closing", "opening".
     */
    std::string operation = "dilation";

    // /**
    //  * @brief Method to use for offset computation.
    //  * Options: "ours" (Voronoi Power-based) or "brute_force".
    //  */
    // std::string method = "ours";

    /**
     * @brief Size of each dexel/grid cell in world units (e.g., millimeters).
     */
    double dexelSize = 1.0;

    /**
     * @brief Radius of the morphological operation.
     * Units depend on `radiusInMM`.
     */
    double radius = 8.0;

    /**
     * @brief Whether the radius is interpreted in millimeters.
     * If false, radius is in dexel units.
     */
    bool radiusInMM = false;

    /**
     * @brief Number of dexels (resolution) along the longest axis of the bounding box.
     */
    int numDexels = 256;

    /**
     * @brief Padding distance to be added around the input mesh, in world units.
     * Default is 0.0 (no padding).
     */
    double padding = 0.0;
};

class MorphologyProcessor {
public:
    // For C-style API / WASM usage
    static MorphologyProcessor* Create(const MorphologyParams& params);
    static void Delete(MorphologyProcessor* instance);

    // WASM-compatible run method using raw pointers
    bool Run(const float* inVertices, size_t numInVertices,
             const unsigned int* inIndices, size_t numInIndices,
             float** outVertices, size_t* numOutVertices,
             unsigned int** outIndices, size_t* numOutIndices);
public:
    MorphologyProcessor(const MorphologyParams& params);

private:
    bool CreateOffsetOperator();
    bool ApplyOperation();

    MorphologyParams m_params;
    std::vector<float> m_inVertices;
    std::vector<unsigned int> m_inIndices;
    std::vector<float>* m_outVertices;
    std::vector<unsigned int>* m_outIndices;

    vor3d::CompressedVolume m_inputVolume;
    vor3d::CompressedVolume m_outputVolume;
    std::unique_ptr<vor3d::VoronoiMorpho> m_offsetOp;
};

} // namespace morpho3d
