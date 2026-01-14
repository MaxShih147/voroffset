#pragma once

#include <cstddef>

namespace voroffset3d {

    class CompressedVolume;

    /**
     * Open-bottom by extending the lowest solid segment downward in each dexel ray,
     * then letting the normal voxel/MC pipeline generate the mesh.
     *
     * Assumption: CompressedVolume stores solid spans along Z for each (x,y) cell:
     *   ray = [begin0, end0, begin1, end1, ...] in "voxel units" (z / spacing).
     * This matches CreateDexelsFromMeshBuffers() which writes z/spacing.
     *
     * This does NOT touch the input shell mesh. It only modifies the inner "cavity" volume (eroded solid)
     * so that later boolean subtraction can produce an open bottom.
     */
    class OpenBottomShellProcessor {
    public:
        /**
         * Apply open-bottom on a dexel volume.
         *
         * @param vol            [in/out] dexel volume (m_outputVolume) representing the inner solid to subtract.
         * @param shellVertices  [in]     original shell vertices (x,y,z interleaved, float)
         * @param shellVertCount [in]     number of vertices (not float count)
         * @param spacingMM      [in]     voxel spacing in mm (m_inputVolume.spacing())
         * @param offsetRadiusMM [in]     erosion radius in mm (same meaning as hollow thickness parameter)
         * @param extraDownMM    [in]     optional extra extension in mm (default 0). Use for safety margin.
         * @param minExtendVox   [in]     minimum extend layers in voxels (default 2). You asked 1~2, default picks 2.
         * @return true if something was modified.
         */
        static bool ApplyToDexelVolume(
            CompressedVolume& vol,
            const float* shellVertices,
            std::size_t shellVertCount,
            float spacingMM,
            float offsetRadiusMM,
            float extraDownMM = 0.0f,
            int minExtendVox = 2
        );
    };

} // namespace voroffset3d
