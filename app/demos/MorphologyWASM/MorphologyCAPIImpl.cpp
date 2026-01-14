#include "MorphologyCAPI.h"
#include <vor3d/MorphologyProcessor.h>

#include <vector>
#include <cmath>
#include <cstring>

using namespace morpho3d;

// ----------------- 3x3 rotation helpers (ZYX Euler) -----------------

struct Mat3f {
    // column-major: [ 0 3 6 ; 1 4 7 ; 2 5 8 ]
    float m[9];
};

static Mat3f matMul(const Mat3f& A, const Mat3f& B)
{
    Mat3f C;
    // C = A * B, column-major
    for (int col = 0; col < 3; ++col) {
        const float b0 = B.m[3 * col + 0];
        const float b1 = B.m[3 * col + 1];
        const float b2 = B.m[3 * col + 2];

        C.m[3 * col + 0] = A.m[0] * b0 + A.m[3] * b1 + A.m[6] * b2;
        C.m[3 * col + 1] = A.m[1] * b0 + A.m[4] * b1 + A.m[7] * b2;
        C.m[3 * col + 2] = A.m[2] * b0 + A.m[5] * b1 + A.m[8] * b2;
    }
    return C;
}

// three.js 的 Euler 'ZYX':
// rotation.set(rx, ry, rz, 'ZYX') 對 column vector v 的效果等價於
// v' = Rz(rz) * Ry(ry) * Rx(rx) * v
static Mat3f makeEulerZYX(float rx, float ry, float rz)
{
    const float sx = std::sin(rx), cx = std::cos(rx);
    const float sy = std::sin(ry), cy = std::cos(ry);
    const float sz = std::sin(rz), cz = std::cos(rz);

    // row-major ZYX:
    // R =
    // [ cy*cz,            cz*sx*sy - cx*sz,   sx*sz + cx*cz*sy ]
    // [ cy*sz,            cx*cz + sx*sy*sz,   cx*sy*sz - cz*sx ]
    // [   -sy,            cy*sx,              cx*cy           ]

    Mat3f R;

    // 轉成 column-major 儲存：
    // column 0 = [R00, R10, R20]
    R.m[0] = cy * cz;
    R.m[1] = cy * sz;
    R.m[2] = -sy;

    // column 1 = [R01, R11, R21]
    R.m[3] = cz * sx * sy - cx * sz;
    R.m[4] = cx * cz + sx * sy * sz;
    R.m[5] = cy * sx;

    // column 2 = [R02, R12, R22]
    R.m[6] = sx * sz + cx * cz * sy;
    R.m[7] = cx * sy * sz - cz * sx;
    R.m[8] = cx * cy;

    return R;
}


static Mat3f transpose(const Mat3f& R)
{
    Mat3f Rt;
    Rt.m[0] = R.m[0]; Rt.m[1] = R.m[3]; Rt.m[2] = R.m[6];
    Rt.m[3] = R.m[1]; Rt.m[4] = R.m[4]; Rt.m[5] = R.m[7];
    Rt.m[6] = R.m[2]; Rt.m[7] = R.m[5]; Rt.m[8] = R.m[8];
    return Rt;
}

// numVerts = 頂點數量（不是 float 數量）
static void applyRotation(const Mat3f& R,
    const float* inVerts,
    float* outVerts,
    std::size_t numVerts)
{
    for (std::size_t i = 0; i < numVerts; ++i) {
        const float x = inVerts[3 * i + 0];
        const float y = inVerts[3 * i + 1];
        const float z = inVerts[3 * i + 2];

        // column-major: v' = R * v
        outVerts[3 * i + 0] = R.m[0] * x + R.m[3] * y + R.m[6] * z;
        outVerts[3 * i + 1] = R.m[1] * x + R.m[4] * y + R.m[7] * z;
        outVerts[3 * i + 2] = R.m[2] * x + R.m[5] * y + R.m[8] * z;
    }
}

extern "C" {

MorphologyHandle Morphology_Create(MorphologyCParams params) {
    MorphologyParams cppParams;
    cppParams.operation   = params.operation ? params.operation : "dilation";
    cppParams.dexelSize   = params.dexelSize;
    cppParams.radius      = params.radius;
    cppParams.radiusInMM  = params.radiusInMM;
    cppParams.numDexels   = params.numDexels;
    cppParams.padding     = params.padding;

    return static_cast<MorphologyHandle>(MorphologyProcessor::Create(cppParams));
}

void Morphology_Delete(MorphologyHandle handle) {
    auto* processor = static_cast<MorphologyProcessor*>(handle);
    MorphologyProcessor::Delete(processor);
}


bool Morphology_Run_OpenBottom(MorphologyHandle handle,
    const float* inVertices, size_t numInVertices,
    const unsigned int* inIndices, size_t numInIndices,
    float rotX_rad, float rotY_rad, float rotZ_rad,
    float** outVertices, size_t* numOutVertices,
    unsigned int** outIndices, size_t* numOutIndices)
{
    auto* processor = static_cast<MorphologyProcessor*>(handle);
    if (!processor) return false;

    if (!inVertices || numInVertices == 0 ||
        !inIndices || numInIndices == 0) {
        return false;
    }

    // numInVertices 是 float 數量，3 * vertex_count
    const std::size_t inFloatCount = static_cast<std::size_t>(numInVertices);
    const std::size_t inVertexCount = inFloatCount / 3;
    if (inVertexCount == 0) {
        return false;
    }



    // 1) 建立 ZYX Euler 旋轉矩陣及其反矩陣
    const Mat3f R = makeEulerZYX(rotX_rad, rotY_rad, rotZ_rad);
    const Mat3f Rinv = transpose(R);


    // 2) 將輸入頂點轉到「開底座標系」
    std::vector<float> rotatedIn(inFloatCount);
    applyRotation(R, inVertices, rotatedIn.data(), inVertexCount);

    // 3) 在旋轉後空間執行原本的 RunOpenBottom
    float* tmpOutVerts = nullptr;
    unsigned int* tmpOutIdx = nullptr;
    size_t tmpOutVertCount = 0;
    size_t tmpOutIdxCount = 0;

    bool ok = processor->RunOpenBottom(
        rotatedIn.data(), numInVertices,    // 注意：仍然傳入 float 數量
        inIndices, numInIndices,
        &tmpOutVerts, &tmpOutVertCount,
        &tmpOutIdx, &tmpOutIdxCount
    );

    if (!ok) {
        if (tmpOutVerts) delete[] tmpOutVerts;
        if (tmpOutIdx)   delete[] tmpOutIdx;
        return false;
    }

    // 4) 將開底後的 cavity 頂點從「開底座標系」轉回「原始座標系」
    if (tmpOutVerts && tmpOutVertCount > 0) {
        const std::size_t outVertexCount = static_cast<std::size_t>(tmpOutVertCount) / 3;
        if (outVertexCount > 0) {
            applyRotation(Rinv, tmpOutVerts, tmpOutVerts, outVertexCount);
        }
    }


    // 5) 將結果交還呼叫端（所有記憶體仍由 MorphologyProcessor 分配）
    *outVertices = tmpOutVerts;
    *numOutVertices = tmpOutVertCount;
    *outIndices = tmpOutIdx;
    *numOutIndices = tmpOutIdxCount;

    return true;
}
bool Morphology_Run(MorphologyHandle handle,
    const float* inVertices, size_t numInVertices,
    const unsigned int* inIndices, size_t numInIndices,
    float** outVertices, size_t* numOutVertices,
    unsigned int** outIndices, size_t* numOutIndices) {

    auto* processor = static_cast<MorphologyProcessor*>(handle);
    if (!processor) return false;

    return processor->Run(
        inVertices, numInVertices,
        inIndices, numInIndices,
        outVertices, numOutVertices,
        outIndices, numOutIndices
    );
}

} // extern "C"