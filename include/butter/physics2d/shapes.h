#pragma once

#include "butter/math/vec2.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <variant>
#include <utility>
#include <vector>

namespace butter::physics2d {

using math::Vec2;

struct Circle { float radius{0.5f}; };
struct Box { Vec2 half_extents{0.5f, 0.5f}; };
struct Polygon { std::vector<Vec2> vertices; };
struct Capsule { float radius{0.25f}; float half_length{0.5f}; };
struct Mesh { std::vector<Polygon> triangles; };
using Convex = Polygon;
using Shape = std::variant<Circle, Box, Polygon, Capsule, Mesh>;

struct Transform {
    Vec2 position{};
    float angle{0};
};

inline Vec2 rotate(Vec2 v, float angle) {
    const float c = std::cos(angle), s = std::sin(angle);
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}

struct Contact {
    Vec2 normal{}; // from A to B
    Vec2 point{};
    float penetration{0};
};

struct AABB {
    Vec2 min{}, max{};
    bool overlaps(const AABB& other) const {
        return min.x <= other.max.x && max.x >= other.min.x && min.y <= other.max.y && max.y >= other.min.y;
    }
};

inline std::vector<Vec2> world_vertices(const Shape& shape, const Transform& t) {
    std::vector<Vec2> result;
    if (const auto* box = std::get_if<Box>(&shape)) {
        result = {{-box->half_extents.x, -box->half_extents.y},
                  { box->half_extents.x, -box->half_extents.y},
                  { box->half_extents.x,  box->half_extents.y},
                  {-box->half_extents.x,  box->half_extents.y}};
    } else if (const auto* polygon = std::get_if<Polygon>(&shape)) {
        result = polygon->vertices;
    } else if (const auto* capsule = std::get_if<Capsule>(&shape)) {
        constexpr int segments = 8;
        for (int i = 0; i <= segments; ++i) {
            const float a = 3.1415926535f * float(i) / segments;
            result.push_back({capsule->radius * std::cos(a), capsule->half_length + capsule->radius * std::sin(a)});
            result.push_back({capsule->radius * std::cos(a), -capsule->half_length - capsule->radius * std::sin(a)});
        }
    }
    for (auto& v : result) v = t.position + rotate(v, t.angle);
    return result;
}

inline Vec2 center_of(const Shape& shape, const Transform& t) {
    if (const auto* circle = std::get_if<Circle>(&shape)) (void)circle;
    if (const auto* polygon = std::get_if<Polygon>(&shape)) {
        Vec2 sum{};
        for (const auto& v : polygon->vertices) sum += v;
        if (!polygon->vertices.empty()) return t.position + rotate(sum / float(polygon->vertices.size()), t.angle);
    }
    return t.position;
}

inline AABB compute_aabb(const Shape& shape, const Transform& t) {
    if (const auto* circle = std::get_if<Circle>(&shape)) {
        const Vec2 r{circle->radius, circle->radius}; return {t.position - r, t.position + r};
    }
    if (const auto* box = std::get_if<Box>(&shape)) {
        float c=std::abs(std::cos(t.angle)), s=std::abs(std::sin(t.angle));
        Vec2 r{c*box->half_extents.x+s*box->half_extents.y,s*box->half_extents.x+c*box->half_extents.y};
        return {t.position-r,t.position+r};
    }
    if (const auto* capsule = std::get_if<Capsule>(&shape)) {
        const Vec2 axis = rotate({0, capsule->half_length}, t.angle), r{capsule->radius, capsule->radius};
        const Vec2 lo = t.position - axis - r, hi = t.position + axis + r;
        return {{std::min(lo.x, hi.x), std::min(lo.y, hi.y)}, {std::max(lo.x, hi.x), std::max(lo.y, hi.y)}};
    }
    if (const auto* mesh = std::get_if<Mesh>(&shape)) {
        AABB result{{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}, {-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()}};
        for (const auto& triangle : mesh->triangles) { const AABB part = compute_aabb(Shape{triangle}, t); result.min.x = std::min(result.min.x, part.min.x); result.min.y = std::min(result.min.y, part.min.y); result.max.x = std::max(result.max.x, part.max.x); result.max.y = std::max(result.max.y, part.max.y); }
        return result;
    }
    // Polygon bounds do not need an allocated world-space vertex array.
    const auto& vertices = std::get<Polygon>(shape).vertices;
    if (vertices.empty()) return {t.position, t.position};
    const Vec2 first = t.position + rotate(vertices[0], t.angle);
    AABB result{first, first};
    for (const auto& local : vertices) {
        const Vec2 v = t.position + rotate(local, t.angle);
        result.min.x = std::min(result.min.x, v.x); result.min.y = std::min(result.min.y, v.y);
        result.max.x = std::max(result.max.x, v.x); result.max.y = std::max(result.max.y, v.y);
    }
    return result;
}

inline bool circle_circle(const Circle& a, const Transform& ta, const Circle& b,
                          const Transform& tb, Contact& c) {
    const Vec2 delta = tb.position - ta.position;
    const float r = a.radius + b.radius;
    const float d2 = delta.length_squared();
    if (d2 >= r * r) return false;
    const float d = std::sqrt(std::max(d2, 1.0e-12f));
    c.normal = d > 1.0e-6f ? delta / d : Vec2{1, 0};
    c.penetration = r - d;
    c.point = ta.position + c.normal * (a.radius - c.penetration * 0.5f);
    return true;
}

inline bool circle_polygon(const Circle& circle, const Transform& tc, const Shape& polygon,
                           const Transform& tp, Contact& c) {
    const auto vertices = world_vertices(polygon, tp);
    if (vertices.size() < 3) return false;
    float best_dist2 = std::numeric_limits<float>::max(); Vec2 closest{};
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 a = vertices[i], b = vertices[(i + 1) % vertices.size()];
        const Vec2 edge = b - a; const float denom = edge.length_squared();
        const float t = denom > 1.0e-12f ? std::clamp((tc.position - a).dot(edge) / denom, 0.0f, 1.0f) : 0.0f;
        // Project from the circle center on face interiors. Reconstructing from a
        // distant endpoint loses tangential precision and injects spurious torque.
        Vec2 point = a + edge * t;
        if (t > 0 && t < 1) {
            const Vec2 normal{-edge.y, edge.x};
            point = tc.position - normal * ((tc.position - a).dot(normal) / denom);
        }
        const float d2 = (tc.position - point).length_squared();
        if (d2 < best_dist2) { best_dist2 = d2; closest = point; }
    }
    const Vec2 delta = tc.position - closest; const float distance = std::sqrt(best_dist2);
    bool positive=false, negative=false;
    for (std::size_t i=0;i<vertices.size();++i) {
        const float side=(vertices[(i+1)%vertices.size()]-vertices[i]).cross(tc.position-vertices[i]);
        positive=positive || side>1.0e-6f; negative=negative || side<-1.0e-6f;
    }
    const bool inside=!(positive && negative);
    if (!inside && distance >= circle.radius) return false;
    if (inside && distance>1.0e-6f) {
        c.normal=delta/distance; c.penetration=circle.radius+distance; c.point=closest; return true;
    }
    c.normal = distance > 1.0e-6f ? -delta / distance : (tp.position - tc.position).normalized();
    if (c.normal.length_squared() < 1.0e-6f) c.normal = {0, 1};
    c.penetration = circle.radius - distance; c.point = closest; return true;
}

inline Vec2 closest_on_segment(Vec2 p, Vec2 a, Vec2 b) {
    const Vec2 d = b - a; const float dd = d.length_squared();
    return a + d * (dd > 1.0e-12f ? std::clamp((p - a).dot(d) / dd, 0.0f, 1.0f) : 0.0f);
}

inline bool capsule_circle(const Capsule& capsule, const Transform& tc, const Circle& circle,
                           const Transform& ts, Contact& c) {
    const Vec2 axis = rotate({0, capsule.half_length}, tc.angle);
    const Vec2 closest = closest_on_segment(ts.position, tc.position - axis, tc.position + axis);
    const Vec2 delta = ts.position - closest; const float r = capsule.radius + circle.radius;
    const float d2 = delta.length_squared(); if (d2 >= r * r) return false;
    const float d = std::sqrt(std::max(d2, 1.0e-12f)); c.normal = d > 1.0e-6f ? delta / d : Vec2{1, 0}; c.penetration = r - d; c.point = closest; return true;
}

inline bool polygon_polygon(const Shape& sa, const Transform& ta, const Shape& sb,
                            const Transform& tb, Contact& c) {
    const auto va = world_vertices(sa, ta), vb = world_vertices(sb, tb);
    if (va.size() < 3 || vb.size() < 3) return false;
    float best = std::numeric_limits<float>::max(); Vec2 best_axis{};
    const Vec2 delta = center_of(sb, tb) - center_of(sa, ta);
    auto project = [](const std::vector<Vec2>& v, Vec2 axis, float& lo, float& hi) {
        lo = hi = v[0].dot(axis);
        for (const auto& p : v) { const float d = p.dot(axis); lo = std::min(lo, d); hi = std::max(hi, d); }
    };
    auto axes_from = [&](const std::vector<Vec2>& v) {
        for (std::size_t i = 0; i < v.size(); ++i) {
            const Vec2 edge = v[(i + 1) % v.size()] - v[i];
            Vec2 axis{-edge.y, edge.x}; axis.normalize();
            float alo, ahi, blo, bhi; project(va, axis, alo, ahi); project(vb, axis, blo, bhi);
            const float overlap = std::min(ahi, bhi) - std::max(alo, blo);
            if (overlap <= 0) return false;
            if (overlap < best) { best = overlap; best_axis = delta.dot(axis) < 0 ? -axis : axis; }
        }
        return true;
    };
    if (!axes_from(va) || !axes_from(vb)) return false;
    c.normal = best_axis; c.penetration = best;
    // Center of the overlapping support faces, rather than the center of A.
    const Vec2 tangent{-best_axis.y,best_axis.x};
    float alo,ahi,blo,bhi,atl,ath,btl,bth;
    project(va,best_axis,alo,ahi);project(vb,best_axis,blo,bhi);
    project(va,tangent,atl,ath);project(vb,tangent,btl,bth);
    c.point=best_axis*((ahi+blo)*0.5f)+tangent*((std::max(atl,btl)+std::min(ath,bth))*0.5f);
    return true;
}

inline bool test(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb, Contact& c) {
    if (const auto* ca = std::get_if<Circle>(&a)) {
        if (const auto* cb = std::get_if<Circle>(&b)) return circle_circle(*ca, ta, *cb, tb, c);
        if (std::holds_alternative<Box>(b) || std::holds_alternative<Polygon>(b)) return circle_polygon(*ca, ta, b, tb, c);
    } else if (const auto* cap = std::get_if<Capsule>(&a)) {
        if (const auto* cb = std::get_if<Circle>(&b)) return capsule_circle(*cap, ta, *cb, tb, c);
    } else if (const auto* cb = std::get_if<Circle>(&b)) {
        if (std::holds_alternative<Box>(a) || std::holds_alternative<Polygon>(a)) {
            if (!circle_polygon(*cb, tb, a, ta, c)) return false;
            c.normal = -c.normal;
            return true;
        }
    }
    if (const auto* ma = std::get_if<Mesh>(&a)) { for (const auto& tri : ma->triangles) if (test(Shape{tri}, ta, b, tb, c)) return true; return false; }
    if (const auto* mb = std::get_if<Mesh>(&b)) { if (test(b, tb, a, ta, c)) { c.normal = -c.normal; return true; } return false; }
    return polygon_polygon(a, ta, b, tb, c);
}

// Clip the incident edge against the reference face. Feature IDs are edge/endpoint
// IDs (including generated intersections), independent of world coordinates.
struct ManifoldPoint { Vec2 point{}; float separation{}; unsigned feature{}; };
struct Manifold { Vec2 normal{}; ManifoldPoint points[2]{}; int count{}; };
inline Manifold contact_manifold(const Shape& a, const Transform& ta, const Shape& b,
                                 const Transform& tb, Vec2 normal, Vec2 fallback,
                                 float separation, float margin = 0.002f) {
    Manifold m; m.normal = normal;
    const auto va = world_vertices(a, ta), vb = world_vertices(b, tb);
    if (va.size() < 3 || vb.size() < 3) {
        m.points[0] = {fallback, separation, 0}; m.count = 1; return m;
    }
    auto face = [](const std::vector<Vec2> &v, Vec2 n) {
        int best = 0;
        float alignment = -2;
        Vec2 center{};
        for (auto p : v)
            center += p;
        center /= float(v.size());
        for (int i = 0; i < int(v.size()); ++i) {
            auto edge = v[(i + 1) % v.size()] - v[i];
            Vec2 out = Vec2{edge.y, -edge.x}.normalized();
            if (out.dot(v[i] - center) < 0)
                out = -out;
            if (out.dot(n) > alignment) {
                alignment = out.dot(n);
                best = i;
            }
        }
        return std::make_pair(best, alignment);
    };
    auto af = face(va, normal), bf = face(vb, -normal);
    bool flip = bf.second > af.second + 0.001f;
    const auto &ref = flip ? vb : va;
    const auto &inc = flip ? va : vb;
    const int ri = flip ? bf.first : af.first;
    Vec2 r0 = ref[ri], r1 = ref[(ri + 1) % ref.size()];
    Vec2 tangent = (r1 - r0).normalized();
    Vec2 n = {tangent.y, -tangent.x};
    if (n.dot(flip ? -normal : normal) < 0)
        n = -n;
    int ii = face(inc, -n).first;
    Vec2 p0 = inc[ii], p1 = inc[(ii + 1) % inc.size()];
    float lo = 0, hi = 1;
    const float length = (r1 - r0).length();
    float start = (p0 - r0).dot(tangent), delta = (p1 - p0).dot(tangent);
    if (std::abs(delta) < 1e-8f) {
        if (start < 0 || start > length)
            return m;
    } else {
        float u = -start / delta, v = (length - start) / delta;
        lo = std::max(0.0f, std::min(u, v));
        hi = std::min(1.0f, std::max(u, v));
    }
    if (lo > hi)
        return m;
    for (int k = 0; k < 2; ++k) {
        float t = k ? hi : lo;
        if (k && hi - lo < 1e-6f)
            break;
        Vec2 p = p0 + (p1 - p0) * t;
        float d = (p - r0).dot(n);
        if (d > margin)
            continue;
        unsigned endpoint = t < 1e-6f ? 0 : t > 1 - 1e-6f ? 1 : 2 + k;
        m.points[m.count++] = {p - n * (d * 0.5f), d,
                               unsigned((flip ? 1 : 0) << 24 | ri << 16 | ii << 8) | endpoint};
    }
    return m;
}

} // namespace butter::physics2d
