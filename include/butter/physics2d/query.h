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

inline std::optional<RaycastHit> ray_aabb(Vec2 origin, Vec2 direction, float max_distance,
                                          const AABB &box, std::size_t index) {
    float tmin = 0, tmax = max_distance;
    Vec2 normal{};
    for (int axis = 0; axis < 2; ++axis) {
        const float o = origin[axis], d = direction[axis], lo = box.min[axis], hi = box.max[axis];
        if (std::abs(d) < 1.0e-8f) {
            if (o < lo || o > hi)
                return std::nullopt;
            continue;
        }
        float t1 = (lo - o) / d, t2 = (hi - o) / d;
        Vec2 n1{}, n2{};
        n1[axis] = -1;
        n2[axis] = 1;
        if (t1 > t2) {
            std::swap(t1, t2);
            std::swap(n1, n2);
        }
        if (t1 > tmin) {
            tmin = t1;
            normal = n1;
        }
        tmax = std::min(tmax, t2);
        if (tmin > tmax)
            return std::nullopt;
    }
    if (tmax < 0)
        return std::nullopt;
    const float distance = std::max(0.0f, tmin);
    return RaycastHit{index, origin + direction * distance, normal, distance};
}

// Exact circle/convex raycast. Rays starting inside a shape are ignored.
inline std::optional<RaycastHit> ray_shape(Vec2 origin, Vec2 direction, float max_distance,
                                           const Shape &shape, const Transform &transform,
                                           std::size_t index = 0) {
    direction = direction.normalized();
    if (direction.length_squared() == 0 || max_distance < 0)
        return std::nullopt;
    if (const auto *circle = std::get_if<Circle>(&shape)) {
        const Vec2 offset = origin - transform.position;
        const float c = offset.length_squared() - circle->radius * circle->radius;
        if (c < 0)
            return std::nullopt;
        const float b = offset.dot(direction), disc = b * b - c;
        if (disc < 0)
            return std::nullopt;
        const float distance = -b - std::sqrt(disc);
        if (distance < 0 || distance > max_distance)
            return std::nullopt;
        const Vec2 point = origin + direction * distance;
        return RaycastHit{index, point, (point - transform.position).normalized(), distance};
    }
    if (const auto *mesh = std::get_if<Mesh>(&shape)) {
        std::optional<RaycastHit> best;
        for (const auto &triangle : mesh->triangles)
            if (auto hit =
                    ray_shape(origin, direction, max_distance, Shape{triangle}, transform, index))
                if (!best || hit->distance < best->distance)
                    best = hit;
        return best;
    }
    const auto vertices = world_vertices(shape, transform);
    if (vertices.size() < 3)
        return std::nullopt;
    float area = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i)
        area += vertices[i].cross(vertices[(i + 1) % vertices.size()]);
    float lower = 0, upper = max_distance;
    Vec2 normal{};
    bool inside = true;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 edge = vertices[(i + 1) % vertices.size()] - vertices[i];
        const Vec2 n = (area >= 0 ? Vec2{edge.y, -edge.x} : Vec2{-edge.y, edge.x}).normalized();
        const float separation = n.dot(origin - vertices[i]);
        if (separation > 0)
            inside = false;
        const float denom = n.dot(direction);
        if (std::abs(denom) < 1.0e-8f) {
            if (separation > 0)
                return std::nullopt;
            continue;
        }
        const float t = -separation / denom;
        if (denom < 0 && t >= lower) {
            lower = t;
            normal = n;
        }
        if (denom > 0)
            upper = std::min(upper, t);
        if (lower > upper)
            return std::nullopt;
    }
    if (inside || lower < 0 || lower > max_distance)
        return std::nullopt;
    return RaycastHit{index, origin + direction * lower, normal, lower};
}

} // namespace butter::physics2d
