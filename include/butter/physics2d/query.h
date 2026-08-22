#pragma once

#include "butter/physics2d/shapes.h"
#include <limits>
#include <optional>

namespace butter::physics2d {

struct RaycastHit {
    std::size_t body_index{0};
    Vec2 point{};
    Vec2 normal{};
    float distance{0};
};

inline std::optional<RaycastHit> ray_aabb(Vec2 origin, Vec2 direction, float max_distance, const AABB& box, std::size_t index) {
    float tmin = 0, tmax = max_distance; Vec2 normal{};
    for (int axis = 0; axis < 2; ++axis) {
        const float o = origin[axis], d = direction[axis], lo = box.min[axis], hi = box.max[axis];
        if (std::abs(d) < 1.0e-8f) { if (o < lo || o > hi) return std::nullopt; continue; }
        float t1 = (lo - o) / d, t2 = (hi - o) / d; Vec2 n1{}, n2{}; n1[axis] = -1; n2[axis] = 1;
        if (t1 > t2) { std::swap(t1, t2); std::swap(n1, n2); }
        if (t1 > tmin) { tmin = t1; normal = n1; }
        tmax = std::min(tmax, t2); if (tmin > tmax) return std::nullopt;
    }
    if (tmax < 0) return std::nullopt; const float distance = std::max(0.0f, tmin);
    return RaycastHit{index, origin + direction * distance, normal, distance};
}

} // namespace butter::physics2d
