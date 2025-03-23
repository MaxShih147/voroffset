// AABB.cpp
#include "vor3d/AABB.h"

void AABB::Initialize(
    size_t itemCount,
    std::function<void(Box&, size_t)> GetBox
) {
    count_ = itemCount;
    boxes_.resize(GetMaxNodeIndex(1, 0, count_) + 1);
    BuildRecursive(1, 0, count_, GetBox);
}

const std::vector<Box>& AABB::Boxes() const {
    return boxes_;
}

void AABB::Traverse(const Box& queryBox, std::function<void(size_t)> callback) const {
    std::function<void(size_t, size_t, size_t)> recurse = [&](size_t nodeIndex, size_t b, size_t e) {
        if (nodeIndex >= boxes_.size()) return;

        const Box& nodeBox = boxes_[nodeIndex];
        if (!BoxesIntersect(queryBox, nodeBox)) return;

        if (b + 1 == e) {
            callback(b);
            return;
        }

        size_t m = b + (e - b) / 2;
        recurse(2 * nodeIndex, b, m);
        recurse(2 * nodeIndex + 1, m, e);
    };

    recurse(1, 0, count_);
}

bool AABB::BoxesIntersect(const Box& a, const Box& b) const {
    return (a.minCorner.x <= b.maxCorner.x && a.maxCorner.x >= b.minCorner.x) &&
           (a.minCorner.y <= b.maxCorner.y && a.maxCorner.y >= b.minCorner.y) &&
           (a.minCorner.z <= b.maxCorner.z && a.maxCorner.z >= b.minCorner.z);
}

size_t AABB::GetMaxNodeIndex(size_t nodeIndex, size_t b, size_t e) const {
    if (b + 1 == e) return nodeIndex;
    size_t m = b + (e - b) / 2;
    size_t left = 2 * nodeIndex;
    size_t right = 2 * nodeIndex + 1;
    return std::max(GetMaxNodeIndex(left, b, m), GetMaxNodeIndex(right, m, e));
}

void AABB::BuildRecursive(size_t nodeIndex, size_t b, size_t e, std::function<void(Box&, size_t)> GetBox) {
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
