#pragma once

#include <vor3d/MorphologyProcessor.h>
#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration of opaque processor handle
typedef void* MorphologyHandle;

// C-style parameters struct (mirroring MorphologyParams)
typedef struct {
    const char* operation;  // "dilation", "erosion", etc.
    const char* method;     // "ours", "brute_force"
    double dexelSize;
    double radius;
    bool radiusInMM;
    int numDexels;
    double padding;
} MorphologyCParams;

// Create and destroy
MorphologyHandle Morphology_Create(MorphologyCParams params);
void Morphology_Delete(MorphologyHandle handle);

// Run morphology operation
bool Morphology_Run(MorphologyHandle handle,
                    const float* inVertices, size_t numInVertices,
                    const unsigned int* inIndices, size_t numInIndices,
                    float** outVertices, size_t* numOutVertices,
                    unsigned int** outIndices, size_t* numOutIndices);

bool Morphology_Run_OpenBottom(MorphologyHandle handle,
    const float* inVertices, size_t numInVertices,
    const unsigned int* inIndices, size_t numInIndices,
    float rotX_rad, float rotY_rad, float rotZ_rad,
    float** outVertices, size_t* numOutVertices,
    unsigned int** outIndices, size_t* numOutIndices);


#ifdef __cplusplus
}
#endif