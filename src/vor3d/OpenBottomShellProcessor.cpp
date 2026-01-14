#include "OpenBottomShellProcessor.h"

#include <vor3d/CompressedVolume.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace voroffset3d
{
    // 同時算出 shell 的 minZ / maxZ（world space）
    static inline void ComputeShellZBounds(
        const float* v, std::size_t vertCount,
        float& outMinZ, float& outMaxZ)
    {
        if (!v || vertCount == 0) {
            outMinZ = outMaxZ = 0.0f;
            return;
        }
        outMinZ = 1e30f;
        outMaxZ = -1e30f;
        for (std::size_t i = 0; i < vertCount; ++i) {
            float z = v[i * 3 + 2];
            if (z < outMinZ) outMinZ = z;
            if (z > outMaxZ) outMaxZ = z;
        }
    }

    static inline void ComputeShellAABB3D(
        const float* v, std::size_t vertCount,
        float& minX, float& maxX,
        float& minY, float& maxY,
        float& minZ, float& maxZ)
    {
        if (!v || vertCount == 0) {
            minX = maxX = minY = maxY = minZ = maxZ = 0.0f;
            return;
        }
        minX = minY = minZ = 1e30f;
        maxX = maxY = maxZ = -1e30f;
        for (std::size_t i = 0; i < vertCount; ++i) {
            const float x = v[i * 3 + 0];
            const float y = v[i * 3 + 1];
            const float z = v[i * 3 + 2];
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
            if (z < minZ) minZ = z;
            if (z > maxZ) maxZ = z;
        }
    }


    bool OpenBottomShellProcessor::ApplyToDexelVolume(
        CompressedVolume& vol,
        const float* shellVertices,
        std::size_t shellVertCount,
        float spacingMM,
        float offsetRadiusMM,
        float extraDownMM,
        int minExtendVox
    )
    {
        if (spacingMM <= 0.0f) return false;
        if (!shellVertices || shellVertCount == 0) return false;

        // 1. 算 shell 的 Z 範圍
        float shellMinZ = 0.0f;
        float shellMaxZ = 0.0f;
        ComputeShellZBounds(shellVertices, shellVertCount, shellMinZ, shellMaxZ);

        // ---- 決定「往下挖」多少（mm） ----
        float baseExtendMM =
            (offsetRadiusMM > 0.0f)
            ? offsetRadiusMM
            : (spacingMM * static_cast<float>(minExtendVox));

        float extendMM = std::max(extraDownMM, baseExtendMM);
        float minExtendByVox = spacingMM * static_cast<float>(minExtendVox);
        if (extendMM < minExtendByVox) extendMM = minExtendByVox;

        // 目標 cavity 底部（world Z）
        const float zBottomWorld = shellMinZ - extendMM;

        // 2. 目前 volume 的 origin / extent
        Eigen::Vector3d origin = vol.origin();
        Eigen::Vector3d extent = vol.extent();
        const float originZ_old = static_cast<float>(origin[2]);
        const float spacing = spacingMM; // vol.spacing() 應該等於 spacingMM

        // ---- Step 2.1: 如有需要，整個 volume origin 下移（這次 log 顯示不一定會動到） ----
        const float safetyVoxBelow = 1.0f; // 在 cavity 底以下再留 1 voxel 空間
        const float desiredOriginZ = zBottomWorld - safetyVoxBelow * spacing;

        if (desiredOriginZ < originZ_old) {
            const float originZ_new = desiredOriginZ;
            const float deltaZ = originZ_old - originZ_new;      // > 0
            const float deltaSpan = deltaZ / spacing;               // voxel 單位

            auto gridSize = vol.gridSize();
            const int W = gridSize[0];
            const int H = gridSize[1];

            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    auto& ray = vol.at(x, y);
                    for (double& v : ray) {
                        v += static_cast<double>(deltaSpan);
                    }
                }
            }

            origin[2] = static_cast<double>(originZ_new);
            extent[2] += static_cast<double>(deltaZ);
            vol.setOrigin(origin);
            vol.setExtent(extent);
        }

        // 重新讀 origin（有可能沒變、有可能剛被改過）
        origin = vol.origin();
        const float originZ = static_cast<float>(origin[2]);

        // ---- Step 2.2: 決定「底部」 voxel index ----
        float targetBegin = (zBottomWorld - originZ) / spacing;
        if (targetBegin < 0.0f) targetBegin = 0.0f;

        // ---- Step 2.3: 決定「頂部」 voxel index，用來避免 cavity 在 +Z 超出 shell ----
        //    這裡我們希望 cavity 頂部略微留一點 margin 在 shell 內部。
        const float topMargin = spacing * 0.5f;      // 可調：0.0f 表示 clamp 到 shellMaxZ
        float zTopWorld = shellMaxZ - topMargin;
        if (zTopWorld < shellMinZ) {
            // 理論上不會發生，但保險：至少不要比底還低
            zTopWorld = shellMinZ;
        }
        float targetEndTop = (zTopWorld - originZ) / spacing;
        if (targetEndTop < 0.0f) {
            targetEndTop = 0.0f;
        }

        auto gridSize2 = vol.gridSize();
        const int W2 = gridSize2[0];
        const int H2 = gridSize2[1];

        bool changed = false;

        for (int y = 0; y < H2; ++y) {
            for (int x = 0; x < W2; ++x) {
                auto& ray = vol.at(x, y);
                if (ray.size() < 2) continue;

                std::size_t n = ray.size();
                if (n & 1) n -= 1;
                if (n < 2) continue;

                // 找出這條 ray 中「最下面」的 begin
                double minBegin = ray[0];
                std::size_t minIdx = 0;
                for (std::size_t i = 0; i + 1 < n; i += 2) {
                    if (ray[i] < minBegin) {
                        minBegin = ray[i];
                        minIdx = i;
                    }
                }

                // 讓最低 segment 的起點不高於 targetBegin（往下拉）
                const double desiredBegin =
                    std::min(minBegin, static_cast<double>(targetBegin));
                if (desiredBegin < minBegin) {
                    ray[minIdx] = desiredBegin;
                    changed = true;
                }

                // 另外，把所有 segment 的 end 都 clamp 不超過 targetEndTop（避免 +Z 超出 shell）
                for (std::size_t i = 0; i + 1 < n; i += 2) {
                    double begin = ray[i];
                    double end = ray[i + 1];
                    double clampedEnd = std::min(end, static_cast<double>(targetEndTop));
                    if (clampedEnd < end) {
                        ray[i + 1] = clampedEnd;
                        changed = true;
                    }
                    // 若 clamp 造成 end < begin，這條 segment 理論上應該被視為空；這邊先不特別處理，
                    // 因為這種情況只會出現在非常邊界的少數 voxel。
                }
            }
        }

        return changed;
    }

} // namespace voroffset3d
