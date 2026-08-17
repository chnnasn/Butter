#pragma once

#include <limits>

#include "vec3.h"

namespace butter::math {

struct AABB {
    Vec3 min{std::numeric_limits<float>::infinity(),
             std::numeric_limits<float>::infinity(),
             std::numeric_limits<float>::infinity()};
    Vec3 max{-std::numeric_limits<float>::infinity(),
             -std::numeric_limits<float>::infinity(),
             -std::numeric_limits<float>::infinity()};

    constexpr AABB() = default;
    constexpr AABB(const Vec3& min, const Vec3& max) : min(min), max(max) {}

    static constexpr AABB empty() { return {}; }

    bool is_empty() const {
        return min.x > max.x || min.y > max.y || min.z > max.z;
    }

    Vec3 center() const { return (min + max) * 0.5f; }
    Vec3 extents() const { return max - min; }
    Vec3 half_extents() const { return extents() * 0.5f; }

    bool contains(const Vec3& point) const {
        return point.x >= min.x && point.x <= max.x &&
               point.y >= min.y && point.y <= max.y &&
               point.z >= min.z && point.z <= max.z;
    }

    bool overlaps(const AABB& other) const {
        return min.x <= other.max.x && max.x >= other.min.x &&
               min.y <= other.max.y && max.y >= other.min.y &&
               min.z <= other.max.z && max.z >= other.min.z;
    }

    AABB expanded(float amount) const {
        const Vec3 delta{amount, amount, amount};
        return {min - delta, max + delta};
    }

    void encapsulate(const Vec3& point) {
        min = min.component_min(point);
        max = max.component_max(point);
    }

    void encapsulate(const AABB& other) {
        if (other.is_empty()) return;
        min = min.component_min(other.min);
        max = max.component_max(other.max);
    }

    static AABB union_of(const AABB& a, const AABB& b) {
        if (a.is_empty()) return b;
        if (b.is_empty()) return a;
        return {a.min.component_min(b.min), a.max.component_max(b.max)};
    }
};

} // namespace butter::math
