#include "MorphologyProcessor.h"
#include "OpenBottomShellProcessor.h"
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

static void buildRotationMatrixXYZ(float rx, float ry, float rz, float R[9]) {
    float cx = cosf(rx), sx = sinf(rx);
    float cy = cosf(ry), sy = sinf(ry);
    float cz = cosf(rz), sz = sinf(rz);

    // R = Rz * Ry * Rx
    R[0] = cz * cy;
    R[1] = cz * sy * sx - sz * cx;
    R[2] = cz * sy * cx + sz * sx;

    R[3] = sz * cy;
    R[4] = sz * sy * sx + cz * cx;
    R[5] = sz * sy * cx - cz * sx;

    R[6] = -sy;
    R[7] = cy * sx;
    R[8] = cy * cx;
}

static void applyRotation(float* v, size_t vertCount, const float R[9]) {
    for (size_t i = 0; i < vertCount; ++i) {
        float x = v[i * 3 + 0];
        float y = v[i * 3 + 1];
        float z = v[i * 3 + 2];
        v[i * 3 + 0] = R[0] * x + R[1] * y + R[2] * z;
        v[i * 3 + 1] = R[3] * x + R[4] * y + R[5] * z;
        v[i * 3 + 2] = R[6] * x + R[7] * y + R[8] * z;
    }
}

static void applyRotationInverse(float* v, size_t vertCount, const float R[9]) {
    // inverse of rotation = transpose
    for (size_t i = 0; i < vertCount; ++i) {
        float x = v[i * 3 + 0];
        float y = v[i * 3 + 1];
        float z = v[i * 3 + 2];
        v[i * 3 + 0] = R[0] * x + R[3] * y + R[6] * z;
        v[i * 3 + 1] = R[1] * x + R[4] * y + R[7] * z;
        v[i * 3 + 2] = R[2] * x + R[5] * y + R[8] * z;
    }
}


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

    m_params.radiusInMM = true;
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



bool MorphologyProcessor::RunOpenBottom(const float* inVertices, size_t numInVertices,
    const unsigned int* inIndices, size_t numInIndices,
    float** outVertices, size_t* numOutVertices,
    unsigned int** outIndices, size_t* numOutIndices)
{

    // Copy input into internal vectors
    m_inVertices.assign(inVertices, inVertices + numInVertices);
    m_inIndices.assign(inIndices, inIndices + numInIndices);
    m_outVertices = new std::vector<float>();
    m_outIndices = new std::vector<unsigned int>();

    // Preserve params (Run() historically mutates m_params.radius).
    const double radiusParam = m_params.radius;
    const bool   radiusInMM  = m_params.radiusInMM;

    // Convert to Dexel volume
    m_inputVolume = voroffset3d::CreateDexelsFromMeshBuffers(
        m_inVertices, m_inIndices, m_params.dexelSize, m_params.padding, m_params.numDexels);


    const double spacing = m_inputVolume.spacing();
    const double radiusVox = radiusInMM ? (radiusParam / spacing) : radiusParam;
    const float  radiusMM  = static_cast<float>(radiusInMM ? radiusParam : (radiusParam * spacing));

    // Use voxel units for morphology operator.
    m_params.radius = radiusVox;

    if (!CreateOffsetOperator()) return false;
    if (!ApplyOperation()) return false;

    // Open-bottom modification: extend the lowest solid span downward in each dexel ray.
    // This is applied only in this branch and never touches the normal hollow path.
    voroffset3d::OpenBottomShellProcessor::ApplyToDexelVolume(
        m_outputVolume,
        inVertices,
        (numInVertices / 3),
        static_cast<float>(m_inputVolume.spacing()),
        radiusMM
    );

    // Dump result to mesh
    voroffset3d::DumpDexelsToVoxelsMC(m_outputVolume, *m_outVertices, *m_outIndices);

    // Allocate and copy output
    *numOutVertices = m_outVertices->size();
    *numOutIndices = m_outIndices->size();
    *outVertices = new float[*numOutVertices];
    *outIndices = new unsigned int[*numOutIndices];
    std::memcpy(*outVertices, m_outVertices->data(), *numOutVertices * sizeof(float));
    std::memcpy(*outIndices, m_outIndices->data(), *numOutIndices * sizeof(unsigned int));

    // Restore params to avoid compounding conversions across calls.
    m_params.radius = radiusParam;
    m_params.radiusInMM = radiusInMM;

    return true;
}

bool MorphologyProcessor::CreateOffsetOperator() {
    m_offsetOp = std::make_unique<vor3d::VoronoiMorphoVorPower>();
    // m_offsetOp = std::make_unique<vor3d::VoronoiMorphoBruteForce>();
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
