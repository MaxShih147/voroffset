#include "MorphologyCAPI.h"
#include <vor3d/MorphologyProcessor.h>

using namespace morpho3d;

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