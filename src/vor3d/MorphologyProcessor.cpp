#include "MorphologyProcessor.h"
// #include "vor3d/DexelConversion.h"
// #include "vor3d/DexelMeshExport.h"
// #include "vor3d/VoronoiMorphoBruteForce.h"
// #include "vor3d/VoronoiMorphoVorPower.h"

// #include <vor3d/CompressedVolume.h>
#include <vor3d/Dexelize.h>
#include <vor3d/VoronoiVorPower.h>
#include <vor3d/VoronoiBruteForce.h>
#include <vor3d/Timer.h>

#include <cstring>
#include <iostream>

using namespace morpho3d;

MorphologyProcessor::MorphologyProcessor(const MorphologyParams& params)
    : m_params(params), m_outVertices(nullptr), m_outIndices(nullptr) {}

MorphologyProcessor* MorphologyProcessor::Create(const MorphologyParams& params) {
    return new MorphologyProcessor(params);
}

void MorphologyProcessor::Delete(MorphologyProcessor* instance) {
    delete instance;
}

bool MorphologyProcessor::Run(const float* inVertices, size_t numInVertices,
                              const unsigned int* inIndices, size_t numInIndices,
                              float** outVertices, size_t* numOutVertices,
                              unsigned int** outIndices, size_t* numOutIndices) {
    // Copy input into internal vectors
    m_inVertices.assign(inVertices, inVertices + numInVertices);
    m_inIndices.assign(inIndices, inIndices + numInIndices);
    m_outVertices = new std::vector<float>();
    m_outIndices = new std::vector<unsigned int>();

    // Convert to Dexel volume
    m_inputVolume = voroffset3d::CreateDexelsFromMeshBuffers(
        m_inVertices, m_inIndices, m_params.dexelSize, m_params.padding, m_params.numDexels);

    if (m_params.radiusInMM) {
        m_params.radius /= m_inputVolume.spacing();
    }

    if (!CreateOffsetOperator()) return false;
    if (!ApplyOperation()) return false;

    // Dump result to mesh
    voroffset3d::DumpDexelsToVoxelsMC(m_outputVolume, *m_outVertices, *m_outIndices);

    // Allocate and copy output
    *numOutVertices = m_outVertices->size();
    *numOutIndices = m_outIndices->size();
    *outVertices = new float[*numOutVertices];
    *outIndices = new unsigned int[*numOutIndices];
    std::memcpy(*outVertices, m_outVertices->data(), *numOutVertices * sizeof(float));
    std::memcpy(*outIndices, m_outIndices->data(), *numOutIndices * sizeof(unsigned int));

    return true;
}

bool MorphologyProcessor::CreateOffsetOperator() {
    if (m_params.method == "ours") {
        m_offsetOp = std::make_unique<vor3d::VoronoiMorphoVorPower>();
    } else if (m_params.method == "brute_force") {
        m_offsetOp = std::make_unique<vor3d::VoronoiMorphoBruteForce>();
    } else {
        std::cerr << "Error: Unknown offset method: " << m_params.method << std::endl;
        return false;
    }
    return true;
}

bool MorphologyProcessor::ApplyOperation() {
    double timeFirst = 0, timeSecond = 0;
    const auto& op = m_params.operation;
    if (op == "noop") {
        m_outputVolume = m_inputVolume;
    } else if (op == "erosion") {
        m_offsetOp->erosion(m_inputVolume, m_outputVolume, m_params.radius, timeFirst, timeSecond);
    } else if (op == "dilation") {
        m_offsetOp->dilation(m_inputVolume, m_outputVolume, m_params.radius, timeFirst, timeSecond);
    } else if (op == "closing") {
        vor3d::CompressedVolume tmp;
        m_offsetOp->dilation(m_inputVolume, tmp, m_params.radius, timeFirst, timeSecond);
        m_offsetOp->erosion(tmp, m_outputVolume, m_params.radius, timeFirst, timeSecond);
    } else if (op == "opening") {
        vor3d::CompressedVolume tmp;
        m_offsetOp->erosion(m_inputVolume, tmp, m_params.radius, timeFirst, timeSecond);
        m_offsetOp->dilation(tmp, m_outputVolume, m_params.radius, timeFirst, timeSecond);
    } else {
        std::cerr << "Error: Unknown operation: " << op << std::endl;
        return false;
    }
    return true;
}
