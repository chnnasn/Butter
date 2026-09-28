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
    if (const auto *capsule = std::get_if<Capsule>(&shape)) {
        // Exact, and the reason it has to be spelled out: the shape fallback
        // below is a convex polygon raycast, and a capsule's outline as a vertex
        // list is an *inscribed* eighteen-gon. A ray grazing the round end would
        // then miss a capsule it geometrically hits, and the ray query would
        // disagree with the discrete test about the same shape.
        //
        // A capsule is the set of points within `radius` of the segment, so the
        // surface is two caps and a barrel. Solving all three exactly and taking
        // the nearest accepted root is the whole algorithm; the caps only own a
        // hit that lands on their own side of the segment, and the barrel only
        // owns one that lands between the two caps.
        const CapsuleSegment segment = capsule_segment(*capsule, transform);
        const Vec2 edge = segment.b - segment.a;
        const float length2 = edge.length_squared();
        const float radius = capsule->radius;
        // Rays that begin inside a shape are ignored, exactly as for a circle.
        if ((origin - closest_on_segment(origin, segment.a, segment.b)).length_squared() <
            radius * radius)
            return std::nullopt;
        float best = max_distance + 1;
        auto consider = [&](float t) {
            if (t >= 0 && t <= max_distance)
                best = std::min(best, t);
        };
        auto cap = [&](Vec2 center, bool upper) {
            const Vec2 offset = origin - center;
            const float b = offset.dot(direction), c = offset.length_squared() - radius * radius;
            const float disc = b * b - c;
            if (disc < 0)
                return;
            const float t = -b - std::sqrt(disc);
            if (t < 0 || t > max_distance)
                return;
            const float along = (origin + direction * t - segment.a).dot(edge);
            if (upper ? along >= length2 : along <= 0)
                consider(t);
        };
        cap(segment.a, false);
        cap(segment.b, true);
        // Barrel: the ray against the infinite cylinder, then discard any root
        // that lies outside the segment's own span. |direction| is one, so the
        // quadratic is a*t^2 + 2*b*t + c with the coefficients below.
        if (length2 > 1.0e-12f) {
            const Vec2 w = origin - segment.a;
            const float dw = direction.dot(edge), ew = w.dot(edge);
            const float a = length2 - dw * dw;
            if (a > 1.0e-12f) {
                const float b = length2 * w.dot(direction) - dw * ew;
                const float c = length2 * w.length_squared() - ew * ew - length2 * radius * radius;
                const float disc = b * b - a * c;
                if (disc >= 0) {
                    const float t = (-b - std::sqrt(disc)) / a;
                    if (t >= 0 && t <= max_distance) {
                        const float along = ew + dw * t;
                        if (along >= 0 && along <= length2)
                            consider(t);
                    }
                }
            }
        }
        if (best > max_distance)
            return std::nullopt;
        const Vec2 point = origin + direction * best;
        const Vec2 closest = closest_on_segment(point, segment.a, segment.b);
        const Vec2 outward = point - closest;
        const Vec2 normal =
            outward.length_squared() > 1.0e-12f
                ? outward.normalized()
                : (length2 > 1.0e-12f ? (segment.a - segment.b).normalized() : Vec2{0, 1});
        return RaycastHit{index, point, normal, best};
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
    // Stack geometry: a heap-allocated world vertex array would dominate the
    // cost of a box/polygon raycast.
    const shape_detail::WorldVertices storage(shape, transform);
    const auto vertices = storage.view();
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
