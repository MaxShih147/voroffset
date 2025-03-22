// AABB.h
#pragma once

#include "vor3d/Box.h"
#include <vector>
#include <functional>
#include <cmath>

class AABB {
public:
    void Initialize(
        size_t itemCount,
        std::function<void(Box&, size_t)> GetBox
    ) {
        itemCount_ = itemCount;
        boxes_.resize(GetMaxNodeIndex(1, 0, itemCount) + 1);
        BuildRecursive(1, 0, itemCount, GetBox);
    }

    const std::vector<Box>& Boxes() const { return boxes_; }

private:
    size_t itemCount_;
    std::vector<Box> boxes_;

    size_t GetMaxNodeIndex(size_t nodeIndex, size_t b, size_t e) const {
        if (b + 1 == e) return nodeIndex;
        size_t m = b + (e - b) / 2;
        size_t left = 2 * nodeIndex;
        size_t right = 2 * nodeIndex + 1;
        return std::max(
            GetMaxNodeIndex(left, b, m),
            GetMaxNodeIndex(right, m, e)
        );
    }

    void BuildRecursive(
        size_t nodeIndex,
        size_t b, size_t e,
        std::function<void(Box&, size_t)> GetBox
    ) {
        if (b + 1 == e) {
            GetBox(boxes_[nodeIndex], b);
            return;
        }
        size_t m = b + (e - b) / 2;
        size_t left = 2 * nodeIndex;
        size_t right = 2 * nodeIndex + 1;
        BuildRecursive(left, b, m, GetBox);
        BuildRecursive(right, m, e, GetBox);
        boxes_[nodeIndex] = Box::Union(boxes_[left], boxes_[right]);
    }
};
