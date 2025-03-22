#pragma once

#include "vor3d/Geometry.h"

class Box {
public:
    Vertex minCorner;
    Vertex maxCorner;

    Box() = default;
    Box(const Vertex& _min, const Vertex& _max) : minCorner(_min), maxCorner(_max) {}

    static Box Union(const Box& a, const Box& b) {
        Vertex minCorner(
            std::min(a.minCorner.x, b.minCorner.x),
            std::min(a.minCorner.y, b.minCorner.y),
            std::min(a.minCorner.z, b.minCorner.z)
        );
        Vertex maxCorner(
            std::max(a.maxCorner.x, b.maxCorner.x),
            std::max(a.maxCorner.y, b.maxCorner.y),
            std::max(a.maxCorner.z, b.maxCorner.z)
        );
        return Box(minCorner, maxCorner);
    }
};