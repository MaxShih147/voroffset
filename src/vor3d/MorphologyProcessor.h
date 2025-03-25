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

struct MorphologyParams {
    double dexelSize = 1.0;
    double padding = 0.0;
    int numDexels = 64;
    double radius = 8.0;
    bool radiusInMM = false;
    std::string method = "ours";       // "ours" or "brute_force"
    std::string operation = "dilation"; // "noop", "erosion", etc.
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
