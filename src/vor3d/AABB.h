// AABB.h
#pragma once

#include "Box.h"
#include <vector>
#include <functional>
#include <cstddef>

/**
 * @brief Axis-Aligned Bounding Box hierarchy for triangle meshes.
 * Supports recursive construction and spatial traversal.
 */
class AABB {
public:
    /**
     * @brief Initialize AABB hierarchy using per-triangle bounding boxes.
     * @param itemCount Number of leaf items (typically number of triangles)
     * @param GetBox A function that sets the bounding box for each triangle index
     */
    void Initialize(
        size_t itemCount,
        std::function<void(Box&, size_t)> GetBox
    );

    /**
     * @brief Access the internal bounding boxes.
     */
    const std::vector<Box>& Boxes() const;

    /**
     * @brief Traverse the tree and apply a callback to all triangle indices whose boxes intersect queryBox.
     * @param queryBox The query region
     * @param callback Function called with triangle indices
     */
    void Traverse(const Box& queryBox, std::function<void(size_t)> callback) const;

private:
    size_t count_ = 0;
    std::vector<Box> boxes_;

    size_t GetMaxNodeIndex(size_t nodeIndex, size_t b, size_t e) const;
    void BuildRecursive(size_t nodeIndex, size_t b, size_t e, std::function<void(Box&, size_t)> GetBox);
    bool BoxesIntersect(const Box& a, const Box& b) const;
};
