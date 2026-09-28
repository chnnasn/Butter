#pragma once

#include "butter/math/vec2.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <cstddef>
#include <optional>
#include <span>
#include <unordered_map>
#include <variant>
#include <utility>
#include <vector>

namespace butter::physics2d {

using math::Vec2;

struct AABB {
    Vec2 min{}, max{};
    bool overlaps(const AABB& other) const {
        return min.x <= other.max.x && max.x >= other.min.x && min.y <= other.max.y && max.y >= other.min.y;
    }
};

struct Circle { float radius{0.5f}; };
struct Box { Vec2 half_extents{0.5f, 0.5f}; };
struct Polygon { std::vector<Vec2> vertices; };
struct Capsule { float radius{0.25f}; float half_length{0.5f}; };

// Diagnostic counter: how many individual triangles the narrow phase actually
// looked at. The whole point of the mesh's spatial index is that this number
// stops tracking the size of the mesh, so it is exposed rather than inferred
// from a stopwatch. One instance per program, like any inline function local.
inline std::size_t& mesh_triangle_tests() {
    static std::size_t count = 0;
    return count;
}

// A triangle soup.
//
// Two things live here beside the triangles, and both exist so that a mesh costs
// a query rather than a scan:
//
//   * a bounding-volume hierarchy over the triangles, built lazily on first use,
//     so that a narrow-phase test or a ray only descends into the part of the
//     mesh that could possibly answer;
//   * the mesh's *boundary* -- the edges that belong to exactly one triangle --
//     with a smoothed outward normal at each endpoint. A shape resting on a
//     terrain made of adjacent triangles has to see one continuous surface, not
//     a sequence of independent triangles, or it catches on every seam; the
//     smoothed normal is what the neighbouring edge agrees with, and the contact
//     normal at a shared vertex is the direction from that vertex to the shape,
//     which both edges compute identically.
//
// The index is derived data, so it is `mutable` and rebuilt on demand. It is
// built from `triangles` and therefore cannot be kept in step with a *direct*
// write to that vector from outside; the engine invalidates it from the same
// per-step scan that notices an edited outline (`body_hash`/`fixture_hash`), and
// `invalidate_acceleration()` is the hook for anyone else.
struct Mesh {
    // One triangle in the hierarchy. `count` is non-zero only on a leaf, where
    // `begin`/`count` index `order`.
    struct Node {
        AABB box{};
        std::uint32_t begin{};
        std::uint32_t count{};
        std::uint32_t left{};
        std::uint32_t right{};
    };
    // One edge of the boundary. `normal_a`/`normal_b` are the *smoothed* outward
    // normals at the two endpoints: the average of the outward normals of every
    // boundary edge meeting there. On a straight run that is the owner
    // triangle's face normal, and at a seam it is the direction neither triangle
    // can call its own -- which is exactly what a shape crossing the seam needs.
    struct Edge {
        Vec2 a{}, b{};
        Vec2 normal_a{}, normal_b{};
        // The owning triangle's own outward face normal, which survives the
        // smoothing below. It is what a contact direction is compared against to
        // tell a body outside the surface from one that has gone inside it: the
        // smoothed normals are what two triangles *agree* on, and the face normal
        // is what this edge alone faces.
        Vec2 outward{};
        std::uint32_t triangle{};
        std::uint32_t next_in_triangle{0xffffffffu};
    };
    struct Accel {
        std::vector<Node> nodes;
        std::vector<std::uint32_t> order;
        std::vector<AABB> bounds;
        std::vector<Vec2> centers;
        std::vector<Edge> boundary;
        // First entry of each triangle's own chain of boundary edges, or
        // `kNone` when every edge of the triangle is interior.
        std::vector<std::uint32_t> triangle_edges;
        AABB local{};
        float radius{};
        bool empty{true};
    };

    std::vector<Polygon> triangles;
    mutable Accel accel{};
    mutable bool accel_valid{false};

    // Hand every triangle whose local bounds overlap `box` to `visit_triangle`,
    // in the mesh's own (local) frame. Shape queries transform their question
    // into that frame rather than the mesh into the world, so this is a
    // descent over a hierarchy that never has to be updated.
    template <class F>
    void visit(const AABB& box, F&& visit_triangle) const {
        const Accel& tree = acceleration();
        if (tree.nodes.empty())
            return;
        std::uint32_t stack[64];
        std::size_t depth = 0;
        stack[depth++] = 0;
        while (depth) {
            const Node& node = tree.nodes[stack[--depth]];
            if (!node.box.overlaps(box))
                continue;
            if (node.count) {
                for (std::uint32_t k = 0; k < node.count; ++k)
                    visit_triangle(tree.order[node.begin + k]);
                continue;
            }
            stack[depth++] = node.left;
            stack[depth++] = node.right;
        }
    }
    const Accel& acceleration() const {
        if (!accel_valid) {
            rebuild_acceleration();
            accel_valid = true;
        }
        return accel;
    }
    void invalidate_acceleration() const { accel_valid = false; }
    void rebuild_acceleration() const;
};

using Convex = Polygon;
using Shape = std::variant<Circle, Box, Polygon, Capsule, Mesh>;

struct Transform {
    Vec2 position{};
    float angle{0};
};

// Cached sine/cosine pair. Evaluating std::cos/std::sin is by far the most
// repeated transcendental work in the solver: a single narrow-phase test or a
// single constraint point otherwise costs several calls. Hoisting one Rot per
// shape/body and reusing it for every vertex keeps the trigonometry count
// proportional to the number of transforms, not to the number of vertices.
struct Rot {
    float c{1}, s{0};
    Rot() = default;
    inline explicit Rot(float angle) : c(std::cos(angle)), s(std::sin(angle)) {}
    inline Rot(float cosine, float sine) : c(cosine), s(sine) {}
};

inline Vec2 rotate(Vec2 v, float angle) {
    const float c = std::cos(angle), s = std::sin(angle);
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}

inline Vec2 rotate(Vec2 v, Rot r) { return {r.c * v.x - r.s * v.y, r.s * v.x + r.c * v.y}; }
inline Rot inverse(Rot r) { return {r.c, -r.s}; }

struct Contact {
    Vec2 normal{}; // from A to B
    Vec2 point{};
    float penetration{0};
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
    const Rot rot(t.angle);
    for (auto& v : result) v = t.position + rotate(v, rot);
    return result;
}

inline Vec2 center_of(const Shape& shape, const Transform& t) {
    if (const auto* polygon = std::get_if<Polygon>(&shape)) {
        Vec2 sum{};
        for (const auto& v : polygon->vertices) sum += v;
        if (!polygon->vertices.empty())
            return t.position + rotate(sum / float(polygon->vertices.size()), Rot(t.angle));
    }
    return t.position;
}

// ---------------------------------------------------------------------------
// Mesh spatial index
// ---------------------------------------------------------------------------
// Two structures are derived from the triangle list, and both are about turning
// "walk every triangle" into "descend to the triangles that matter":
//
//   * a median-split bounding-volume hierarchy. Every query (narrow phase, ray,
//     point-in-shape) first transforms its question into the mesh's own frame,
//     so the hierarchy is built once in local space and never touched again.
//   * the boundary edge list. A terrain authored as a strip of triangles has its
//     *interior* edges (the vertical cuts between neighbouring triangles) hidden
//     inside the material; only the outer edges are a surface. Contacts are
//     generated against the boundary, which is why a shape sliding along a mesh
//     floor never sees the cut behind it.
//
// The smoothed normal at a boundary vertex is the average of the outward normals
// of every boundary edge meeting there. Two edges of a seam therefore agree
// about the surface at the vertex they share, and a shape crossing the seam
// slides over it instead of catching on the edge of whichever triangle it has
// just left.
namespace mesh_detail {
inline std::int64_t quantize(float value) {
    return std::int64_t(std::llround(double(value) * 4096.0));
}
struct PointKey {
    std::int64_t x{}, y{};
    friend bool operator==(const PointKey&, const PointKey&) = default;
    friend bool operator!=(const PointKey&, const PointKey&) = default;
};
struct PointKeyHash {
    std::size_t operator()(const PointKey& key) const {
        std::uint64_t h = std::uint64_t(key.x) * 0x9e3779b97f4a7c15ull;
        h ^= std::uint64_t(key.y) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        return std::size_t(h);
    }
};
// An undirected edge key: the two endpoint keys in a fixed order, so a shared
// edge is found however the two triangles happen to wind it.
struct EdgeKey {
    PointKey low{}, high{};
    friend bool operator==(const EdgeKey&, const EdgeKey&) = default;
};
struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& key) const {
        const PointKeyHash hash;
        std::size_t h = hash(key.low);
        h ^= hash(key.high) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        return h;
    }
};
inline PointKey point_key(Vec2 v) { return {quantize(v.x), quantize(v.y)}; }
inline EdgeKey edge_key(Vec2 a, Vec2 b) {
    const PointKey ka = point_key(a), kb = point_key(b);
    return ka.x < kb.x || (ka.x == kb.x && ka.y < kb.y) ? EdgeKey{ka, kb} : EdgeKey{kb, ka};
}
// The unit perpendicular of an edge, pointed away from the triangle's own
// centroid. Which way the outline is wound is an authoring choice, so the
// outside is decided geometrically rather than from the winding.
inline Vec2 outward_normal(Vec2 a, Vec2 b, Vec2 centroid) {
    const Vec2 edge = b - a;
    if (edge.length_squared() <= 0)
        return {};
    Vec2 normal{edge.y, -edge.x};
    if (normal.dot((a + b) * 0.5f - centroid) < 0)
        normal = -normal;
    return normal.normalized();
}
} // namespace mesh_detail

inline void Mesh::rebuild_acceleration() const {
    Accel& out = accel;
    out.nodes.clear();
    out.order.clear();
    out.bounds.clear();
    out.centers.clear();
    out.boundary.clear();
    out.triangle_edges.clear();
    out.radius = 0;
    out.local = {};
    out.empty = triangles.empty();
    if (out.empty)
        return;

    const std::size_t count = triangles.size();
    out.bounds.resize(count);
    out.centers.resize(count);
    out.order.resize(count);
    Vec2 lower{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    Vec2 upper{-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
    for (std::size_t i = 0; i < count; ++i) {
        const auto& vertices = triangles[i].vertices;
        AABB box{};
        Vec2 center{};
        if (!vertices.empty()) {
            box.min = box.max = vertices.front();
            for (const Vec2 v : vertices) {
                box.min.x = std::min(box.min.x, v.x);
                box.min.y = std::min(box.min.y, v.y);
                box.max.x = std::max(box.max.x, v.x);
                box.max.y = std::max(box.max.y, v.y);
                center += v;
                out.radius = std::max(out.radius, v.length());
            }
            center /= float(vertices.size());
        }
        out.bounds[i] = box;
        out.centers[i] = center;
        lower.x = std::min(lower.x, box.min.x);
        lower.y = std::min(lower.y, box.min.y);
        upper.x = std::max(upper.x, box.max.x);
        upper.y = std::max(upper.y, box.max.y);
        out.order[i] = std::uint32_t(i);
    }
    out.local = {lower, upper};

    // Median split on the longest axis of the centroids: the hierarchy is built
    // once and never updated, so a balanced-by-construction tree is worth the
    // nth_element passes, which is what makes the descent depth logarithmic and
    // the fixed-size stack in `visit` sufficient.
    struct Frame {
        std::uint32_t node{}, begin{}, end{};
    };
    out.nodes.reserve(2 * count);
    out.nodes.emplace_back();
    std::vector<Frame> pending;
    pending.push_back({0, 0, std::uint32_t(count)});
    constexpr float kHuge = std::numeric_limits<float>::max();
    while (!pending.empty()) {
        const Frame frame = pending.back();
        pending.pop_back();
        AABB box{};
        bool first = true;
        Vec2 low{kHuge, kHuge}, high{-kHuge, -kHuge};
        for (std::uint32_t k = frame.begin; k < frame.end; ++k) {
            const std::uint32_t triangle = out.order[k];
            const AABB& bounds = out.bounds[triangle];
            // The accumulator starts at the first triangle, not at the origin.
            // Starting from a default AABB makes every node box the union of its
            // real box with (0, 0), so a mesh that does not straddle the origin
            // gets boxes that all reach back to it and the descent stops
            // pruning: a query anywhere visited half the mesh.
            if (first) {
                box = bounds;
                first = false;
            } else {
                box.min.x = std::min(box.min.x, bounds.min.x);
                box.min.y = std::min(box.min.y, bounds.min.y);
                box.max.x = std::max(box.max.x, bounds.max.x);
                box.max.y = std::max(box.max.y, bounds.max.y);
            }
            const Vec2 center = out.centers[triangle];
            low.x = std::min(low.x, center.x);
            low.y = std::min(low.y, center.y);
            high.x = std::max(high.x, center.x);
            high.y = std::max(high.y, center.y);
        }
        out.nodes[frame.node].box = box;
        const std::uint32_t size = frame.end - frame.begin;
        if (size <= 4) {
            out.nodes[frame.node].count = size;
            out.nodes[frame.node].begin = frame.begin;
            continue;
        }
        const bool vertical = high.y - low.y > high.x - low.x;
        const std::uint32_t middle = frame.begin + size / 2;
        std::nth_element(out.order.begin() + frame.begin, out.order.begin() + middle,
                         out.order.begin() + frame.end,
                         [&](std::uint32_t a, std::uint32_t b) {
                             return vertical ? out.centers[a].y < out.centers[b].y
                                             : out.centers[a].x < out.centers[b].x;
                         });
        const std::uint32_t left = std::uint32_t(out.nodes.size());
        out.nodes.emplace_back();
        const std::uint32_t right = std::uint32_t(out.nodes.size());
        out.nodes.emplace_back();
        out.nodes[frame.node].count = 0;
        out.nodes[frame.node].left = left;
        out.nodes[frame.node].right = right;
        pending.push_back({right, middle, frame.end});
        pending.push_back({left, frame.begin, middle});
    }

    // Boundary edges: an edge used by exactly one triangle. `count > 2` is a
    // non-manifold edge; it is treated as boundary too, because the alternative
    // is deciding which of the three faces owns it.
    struct RawEdge {
        Vec2 a{}, b{};
        std::uint32_t triangle{};
    };
    std::vector<RawEdge> raw;
    std::unordered_map<mesh_detail::EdgeKey, int, mesh_detail::EdgeKeyHash> uses;
    for (std::size_t t = 0; t < count; ++t) {
        const auto& vertices = triangles[t].vertices;
        if (vertices.size() < 2)
            continue;
        for (std::size_t e = 0; e < vertices.size(); ++e) {
            const Vec2 a = vertices[e], b = vertices[(e + 1) % vertices.size()];
            if (a == b)
                continue;
            ++uses[mesh_detail::edge_key(a, b)];
            raw.push_back({a, b, std::uint32_t(t)});
        }
    }
    for (const RawEdge& edge : raw) {
        const auto found = uses.find(mesh_detail::edge_key(edge.a, edge.b));
        if (found == uses.end() || found->second == 2)
            continue; // Interior: the material continues past it.
        const Vec2 normal =
            mesh_detail::outward_normal(edge.a, edge.b, out.centers[edge.triangle]);
        if (normal.length_squared() <= 0)
            continue;
        out.boundary.push_back({edge.a, edge.b, normal, normal, normal, edge.triangle, 0xffffffffu});
    }
    // The smoothed normal at each boundary vertex, then the per-triangle chains
    // that let a candidate triangle offer its own boundary edges.
    std::unordered_map<mesh_detail::PointKey, Vec2, mesh_detail::PointKeyHash> smoothed;
    for (const Edge& edge : out.boundary) {
        const auto a = smoothed.find(mesh_detail::point_key(edge.a));
        if (a == smoothed.end())
            smoothed.emplace(mesh_detail::point_key(edge.a), edge.normal_a);
        else
            a->second += edge.normal_a;
        const auto b = smoothed.find(mesh_detail::point_key(edge.b));
        if (b == smoothed.end())
            smoothed.emplace(mesh_detail::point_key(edge.b), edge.normal_b);
        else
            b->second += edge.normal_b;
    }
    out.triangle_edges.assign(count, 0xffffffffu);
    for (std::size_t i = 0; i < out.boundary.size(); ++i) {
        Edge& edge = out.boundary[i];
        const Vec2 sum_a = smoothed[mesh_detail::point_key(edge.a)];
        const Vec2 sum_b = smoothed[mesh_detail::point_key(edge.b)];
        // Two boundary edges that meet at a cusp can cancel out. Falling back to
        // the edge's own normal keeps the surface defined rather than producing a
        // zero normal.
        edge.normal_a = sum_a.length_squared() > 1.0e-12f ? sum_a.normalized() : edge.normal_a;
        edge.normal_b = sum_b.length_squared() > 1.0e-12f ? sum_b.normalized() : edge.normal_b;
        edge.next_in_triangle = out.triangle_edges[edge.triangle];
        out.triangle_edges[edge.triangle] = std::uint32_t(i);
    }
}

// ---------------------------------------------------------------------------
// Shape signatures and mass properties
// ---------------------------------------------------------------------------

// Cheap, allocation-free content signature. Contact caches use it to notice
// public edits to a shape without snapshotting (deep-copying) the geometry, and
// without the O(n) vertex comparisons a structural equality needs.
inline std::uint64_t shape_hash(const Shape& shape) {
    std::uint64_t h = 1469598103934665603ull; // FNV-1a
    auto mix = [&h](std::uint64_t value) {
        h ^= value;
        h *= 1099511628211ull;
    };
    auto mix_float = [&mix](float value) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &value, sizeof(bits));
        mix(bits);
    };
    auto mix_point = [&mix_float](Vec2 v) {
        mix_float(v.x);
        mix_float(v.y);
    };
    mix(shape.index());
    if (const auto* circle = std::get_if<Circle>(&shape))
        mix_float(circle->radius);
    else if (const auto* box = std::get_if<Box>(&shape))
        mix_point(box->half_extents);
    else if (const auto* polygon = std::get_if<Polygon>(&shape)) {
        mix(polygon->vertices.size());
        for (auto v : polygon->vertices)
            mix_point(v);
    } else if (const auto* capsule = std::get_if<Capsule>(&shape)) {
        mix_float(capsule->radius);
        mix_float(capsule->half_length);
    } else {
        const auto& triangles = std::get<Mesh>(shape).triangles;
        mix(triangles.size());
        for (const auto& triangle : triangles) {
            mix(triangle.vertices.size());
            for (auto v : triangle.vertices)
                mix_point(v);
        }
    }
    return h;
}

// Mass, centroid and rotational inertia about the centroid, in local space.
struct MassData {
    float mass{};
    Vec2 center{};
    float inertia{};
};

inline MassData combine(const MassData& a, const MassData& b) {
    const float mass = a.mass + b.mass;
    if (mass <= 0)
        return {};
    MassData result;
    result.mass = mass;
    result.center = (a.center * a.mass + b.center * b.mass) / mass;
    // Parallel axis theorem: shift each part's centroid inertia onto the
    // combined centroid.
    result.inertia = a.inertia + a.mass * (a.center - result.center).length_squared() +
                     b.inertia + b.mass * (b.center - result.center).length_squared();
    return result;
}

inline MassData transformed(const MassData& data, const Transform& local) {
    MassData result = data;
    result.center = local.position + rotate(data.center, Rot(local.angle));
    return result;
}

inline MassData mass_data(const Shape& shape, float density) {
    MassData result;
    if (density <= 0)
        return result;
    if (const auto* circle = std::get_if<Circle>(&shape)) {
        result.mass = density * 3.14159265358979f * circle->radius * circle->radius;
        result.inertia = 0.5f * result.mass * circle->radius * circle->radius;
        return result;
    }
    if (const auto* box = std::get_if<Box>(&shape)) {
        const float hx = box->half_extents.x, hy = box->half_extents.y;
        result.mass = density * 4.0f * hx * hy;
        result.inertia = result.mass * (hx * hx + hy * hy) / 3.0f;
        return result;
    }
    if (const auto* capsule = std::get_if<Capsule>(&shape)) {
        const float r = capsule->radius, hl = capsule->half_length;
        const float rect_mass = density * 4.0f * r * hl;
        const float rect_inertia = rect_mass * (r * r + hl * hl) / 3.0f;
        const float cap_mass = density * 3.14159265358979f * r * r;
        const float cap_inertia = 0.5f * cap_mass * r * r;
        result.mass = rect_mass + cap_mass;
        result.inertia = rect_inertia + cap_inertia;
        return result;
    }
    if (const auto* mesh = std::get_if<Mesh>(&shape)) {
        for (const auto& triangle : mesh->triangles)
            result = combine(result, mass_data(Shape{triangle}, density));
        return result;
    }
    // Convex polygon: exact area/centroid/inertia integrals over the fan of
    // signed triangles around the local origin, then shifted to the centroid.
    const auto& vertices = std::get<Polygon>(shape).vertices;
    if (vertices.size() < 3)
        return result;
    float area = 0, moment = 0;
    Vec2 centroid{};
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 a = vertices[i], b = vertices[(i + 1) % vertices.size()];
        const float cross = a.cross(b);
        area += cross;
        centroid += (a + b) * cross;
        moment += cross * (a.dot(a) + a.dot(b) + b.dot(b));
    }
    area *= 0.5f;
    if (std::abs(area) < 1.0e-12f)
        return result; // Degenerate (collinear) outline carries no mass.
    const float magnitude = std::abs(area);
    result.center = centroid / (6.0f * area);
    result.mass = density * magnitude;
    result.inertia = density * magnitude * (moment / (12.0f * area)) - result.mass * result.center.length_squared();
    if (result.inertia < 0)
        result.inertia = 0;
    return result;
}

// Distance from the local origin to the farthest point of the shape. Used for
// conservative bounds (CCD) and motion thresholds.
inline float shape_radius(const Shape& shape) {
    if (const auto* circle = std::get_if<Circle>(&shape))
        return circle->radius;
    if (const auto* box = std::get_if<Box>(&shape))
        return box->half_extents.length();
    if (const auto* capsule = std::get_if<Capsule>(&shape))
        return capsule->half_length + capsule->radius;
    float result = 0;
    auto consider = [&result](const std::vector<Vec2>& vertices) {
        for (auto v : vertices)
            result = std::max(result, v.length());
    };
    if (const auto* polygon = std::get_if<Polygon>(&shape))
        consider(polygon->vertices);
    else
        // A mesh's reach is part of its spatial index, so asking for it does not
        // walk the triangles.
        return std::get<Mesh>(shape).acceleration().radius;
    return result;
}

namespace shape_detail {
// Most convex fixtures have only a few vertices. Keep their temporary world
// geometry on the stack, with an unrestricted fallback for larger polygons.
class WorldVertices {
    std::array<Vec2, 16> local_;
    std::vector<Vec2> overflow_;
    std::size_t size_{};
public:
    WorldVertices(const Shape& shape, const Transform& t) {
        if (const auto* box = std::get_if<Box>(&shape)) {
            size_ = 4;
            local_[0] = {-box->half_extents.x, -box->half_extents.y};
            local_[1] = { box->half_extents.x, -box->half_extents.y};
            local_[2] = { box->half_extents.x,  box->half_extents.y};
            local_[3] = {-box->half_extents.x,  box->half_extents.y};
        } else if (const auto* polygon = std::get_if<Polygon>(&shape)) {
            size_ = polygon->vertices.size();
            if (size_ <= local_.size())
                std::copy(polygon->vertices.begin(), polygon->vertices.end(), local_.begin());
            else
                overflow_ = polygon->vertices;
        } else if (std::holds_alternative<Capsule>(shape)) {
            overflow_ = world_vertices(shape, t);
            size_ = overflow_.size();
            return;
        }
        // One Rot for the whole outline instead of one trig pair per vertex.
        const Rot rot(t.angle);
        auto* data = size_ <= local_.size() ? local_.data() : overflow_.data();
        for (std::size_t i = 0; i < size_; ++i)
            data[i] = t.position + rotate(data[i], rot);
    }
    std::span<const Vec2> view() const {
        return {size_ <= local_.size() ? local_.data() : overflow_.data(), size_};
    }
};
} // namespace shape_detail

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
        // The triangles' own bounds were computed once, in local space, when the
        // index was built. `local` therefore already encloses the whole mesh, and
        // enclosing a rotated box needs its four corners and nothing more -- this
        // used to re-derive a bounds for every triangle on every call, which on a
        // large mesh dominated the narrow phase.
        const Mesh::Accel& tree = mesh->acceleration();
        if (tree.empty)
            return {t.position, t.position};
        const Rot rot(t.angle);
        const Vec2 corners[4] = {{tree.local.min.x, tree.local.min.y},
                                 {tree.local.max.x, tree.local.min.y},
                                 {tree.local.max.x, tree.local.max.y},
                                 {tree.local.min.x, tree.local.max.y}};
        AABB result{};
        for (int k = 0; k < 4; ++k) {
            const Vec2 v = t.position + rotate(corners[k], rot);
            if (k == 0) {
                result.min = result.max = v;
                continue;
            }
            result.min.x = std::min(result.min.x, v.x);
            result.min.y = std::min(result.min.y, v.y);
            result.max.x = std::max(result.max.x, v.x);
            result.max.y = std::max(result.max.y, v.y);
        }
        return result;
    }
    // Polygon bounds do not need an allocated world-space vertex array.
    const auto& vertices = std::get<Polygon>(shape).vertices;
    if (vertices.empty()) return {t.position, t.position};
    const Rot rot(t.angle);
    const Vec2 first = t.position + rotate(vertices[0], rot);
    AABB result{first, first};
    for (const auto& local : vertices) {
        const Vec2 v = t.position + rotate(local, rot);
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
    const shape_detail::WorldVertices storage(polygon, tp);
    const auto vertices = storage.view();
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

// ---------------------------------------------------------------------------
// Capsule: a segment plus a radius, handled exactly
//
// A capsule is precisely the set of points within `radius` of a segment, and
// `capsule_circle` below always meant exactly that. Nothing else did: the
// generic vertex machinery could not express a capsule at all, so
// `world_vertices` tessellated its outline into an eighteen-gon and the
// remaining narrow-phase pairs ran SAT over that tessellation. Two things
// followed from it, and both are bugs rather than approximations:
//
//   * The result depended on argument order. `test(Capsule, Circle)` reached the
//     exact routine; `test(Circle, Capsule)` fell through to polygon baseline
//     SAT, and a Circle has no vertex list, so the pair reported *no contact at
//     all*. Which of the two a broad-phase pair got came down to which body
//     happened to be created first.
//   * The contact normal came from wherever the tessellation put a vertex, so a
//     capsule resting on a flat floor crept sideways and its resting depth
//     depended on the segment count.
//
// The helpers below keep the flat sides flat and the round ends round, so a
// capsule's discrete test, its ray query and its sweep all describe the same
// shape.
// ---------------------------------------------------------------------------
struct CapsuleSegment {
    Vec2 a{}, b{};
};
inline CapsuleSegment capsule_segment(const Capsule& capsule, const Transform& t) {
    const Vec2 axis = rotate({0, capsule.half_length}, t.angle);
    return {t.position - axis, t.position + axis};
}
// Support along `n` is the far end of the segment plus the radius: the round end
// the direction actually points at, not a tessellated vertex near it.
inline Vec2 capsule_support(const Capsule& capsule, const Transform& t, Vec2 n) {
    const CapsuleSegment segment = capsule_segment(capsule, t);
    return (segment.a.dot(n) > segment.b.dot(n) ? segment.a : segment.b) + n * capsule.radius;
}
// Closest pair between two segments (Ericson, Real-Time Collision Detection
// 5.1.9). For two segments that do not intersect the closest pair always
// contains an endpoint of one of them, which is what makes the distance to a
// polygon boundary below exact rather than sampled.
inline std::pair<Vec2, Vec2> segment_closest_points(CapsuleSegment s, CapsuleSegment o) {
    const Vec2 d1 = s.b - s.a, d2 = o.b - o.a, r = s.a - o.a;
    const float a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r), c = d1.dot(r);
    if (a <= 1.0e-12f && e <= 1.0e-12f)
        return {s.a, o.a};
    float u = 0, v = 0;
    if (a <= 1.0e-12f) {
        v = std::clamp(f / e, 0.0f, 1.0f);
    } else if (e <= 1.0e-12f) {
        u = std::clamp(c / a, 0.0f, 1.0f);
    } else {
        const float b = d1.dot(d2), denom = a * e - b * b;
        u = denom > 1.0e-12f ? std::clamp((b * f - c * e) / denom, 0.0f, 1.0f) : 0.0f;
        v = (b * u + f) / e;
        if (v < 0) {
            v = 0;
            u = std::clamp(-c / a, 0.0f, 1.0f);
        } else if (v > 1) {
            v = 1;
            u = std::clamp((b - c) / a, 0.0f, 1.0f);
        }
    }
    return {s.a + d1 * u, o.a + d2 * v};
}
// `normal` must point from A to B, like every other routine in this file.
inline bool capsule_circle(const Capsule& capsule, const Transform& tc, const Circle& circle,
                           const Transform& ts, Contact& c) {
    const CapsuleSegment segment = capsule_segment(capsule, tc);
    const Vec2 closest = closest_on_segment(ts.position, segment.a, segment.b);
    const Vec2 delta = ts.position - closest; const float r = capsule.radius + circle.radius;
    const float d2 = delta.length_squared(); if (d2 >= r * r) return false;
    const float d = std::sqrt(std::max(d2, 1.0e-12f)); c.normal = d > 1.0e-6f ? delta / d : Vec2{1, 0}; c.penetration = r - d; c.point = closest; return true;
}

inline bool capsule_capsule(const Capsule& a, const Transform& ta, const Capsule& b,
                            const Transform& tb, Contact& c) {
    const auto [pa, pb] = segment_closest_points(capsule_segment(a, ta), capsule_segment(b, tb));
    const Vec2 delta = pb - pa;
    const float r = a.radius + b.radius;
    const float d2 = delta.length_squared();
    if (d2 >= r * r) return false;
    const float d = std::sqrt(std::max(d2, 1.0e-12f));
    // Two parallel capsules overlap along a whole strip, so the closest pair is
    // not unique and its direction carries no information. The line between the
    // segment midpoints is the same answer for every point of the strip, which
    // is what keeps a resting stack of parallel capsules from jittering.
    if (d > 1.0e-6f)
        c.normal = delta / d;
    else {
        c.normal = (tb.position - ta.position).normalized();
        if (c.normal.length_squared() < 1.0e-6f) c.normal = {0, 1};
    }
    c.penetration = r - d;
    c.point = (pa + pb) * 0.5f;
    return true;
}

// Capsule against a convex Box or Polygon.
//
// The exact statement is "the segment is closer to the polygon than the
// radius", so the test is the distance from the segment to the polygon's
// boundary. Taking the minimum over edges -- rather than accepting the first
// edge that reports a hit -- is also what stops a body from snagging on the
// shared edge of a terrain mesh: both neighbouring triangles measure the same
// distance to the same segment, so neither can claim the contact alone.
inline bool capsule_convex(const Capsule& capsule, const Transform& tc, const Shape& convex,
                           const Transform& tp, Contact& c) {
    const shape_detail::WorldVertices storage(convex, tp);
    const auto vertices = storage.view();
    if (vertices.size() < 3) return false;
    const CapsuleSegment segment = capsule_segment(capsule, tc);
    float winding = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i)
        winding += vertices[i].cross(vertices[(i + 1) % vertices.size()]);
    // Outward normal of edge `i`, orients by the polygon's winding so a
    // clockwise outline is handled as well as a counter-clockwise one.
    auto outward = [&](std::size_t i) {
        const Vec2 edge = vertices[(i + 1) % vertices.size()] - vertices[i];
        return (winding >= 0 ? Vec2{edge.y, -edge.x} : Vec2{-edge.y, edge.x}).normalized();
    };
    // Nearest boundary point to the segment, over every edge. For two segments
    // that do not cross, the closest pair contains an endpoint of one of them,
    // so this is the exact distance rather than a sample of it.
    float best = std::numeric_limits<float>::max();
    Vec2 segment_point{}, boundary_point{};
    std::size_t nearest_edge = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto [p, q] =
            segment_closest_points(segment, CapsuleSegment{vertices[i], vertices[(i + 1) % vertices.size()]});
        const float d2 = (q - p).length_squared();
        if (d2 < best) {
            best = d2;
            segment_point = p;
            boundary_point = q;
            nearest_edge = i;
        }
    }
    // Both endpoints on the inner side of every edge means the whole segment is
    // inside, which is the only case where the push-out direction reverses.
    bool inside = true;
    for (std::size_t i = 0; i < vertices.size() && inside; ++i) {
        const Vec2 out = outward(i);
        for (Vec2 end : {segment.a, segment.b})
            if (out.dot(end - vertices[i]) > 0) {
                inside = false;
                break;
            }
    }
    const float distance = std::sqrt(std::max(best, 0.0f));
    if (!inside && distance >= capsule.radius)
        return false;
    c.point = distance > 1.0e-6f ? boundary_point : segment_point;
    if (distance > 1.0e-6f) {
        // Outside: A is pushed away from the boundary. Inside: A is pushed back
        // towards it, so the normal points deeper in. Same convention as
        // circle_polygon's two branches, which the solver already relies on.
        const Vec2 delta = inside ? segment_point - boundary_point : boundary_point - segment_point;
        c.normal = delta / distance;
    } else {
        // The segment crosses the boundary or lies along it, so the closest pair
        // is a single point and carries no direction. The edge normal is what
        // still means something.
        c.normal = -outward(nearest_edge);
        if (c.normal.length_squared() < 0.5f) c.normal = {0, 1};
    }
    c.penetration = inside ? capsule.radius + distance : capsule.radius - distance;
    return true;
}

inline bool polygon_polygon(const Shape& sa, const Transform& ta, const Shape& sb,
                            const Transform& tb, Contact& c) {
    const shape_detail::WorldVertices storage_a(sa, ta), storage_b(sb, tb);
    const auto va = storage_a.view(), vb = storage_b.view();
    if (va.size() < 3 || vb.size() < 3) return false;
    float best = std::numeric_limits<float>::max(); Vec2 best_axis{};
    const Vec2 delta = center_of(sb, tb) - center_of(sa, ta);
    auto project = [](std::span<const Vec2> v, Vec2 axis, float& lo, float& hi) {
        lo = hi = v[0].dot(axis);
        for (const auto& p : v) { const float d = p.dot(axis); lo = std::min(lo, d); hi = std::max(hi, d); }
    };
    auto axes_from = [&](std::span<const Vec2> v) {
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

// ---------------------------------------------------------------------------
// Mesh contacts
// ---------------------------------------------------------------------------
// Contacts against a mesh are generated against its *boundary* rather than
// against its triangles. The distinction is the whole of the seam problem: the
// interior edges of a triangle soup are cuts through the material, and treating
// them as geometry is exactly what makes a body catch on the seam between two
// triangles that are supposed to be one floor.
//
// Only the outer edges are tested, then, and the normal a contact gets is one
// of three things: the face of the outline that is touching (an outline face
// contact), the edge's own face (an edge-interior contact), or the direction
// from an endpoint of the edge to the closest point of the outline (a
// vertex contact). That last case is the seam normal: the neighbouring edge
// computes the same direction for the same vertex, so a shape sliding over the
// seam sees one surface instead of two disagreeing ones.
namespace mesh_detail {
inline Transform into_local(const Transform& other, const Transform& mesh) {
    return {rotate(other.position - mesh.position, inverse(Rot(mesh.angle))),
            other.angle - mesh.angle};
}
inline Vec2 direction_into_local(Vec2 v, const Transform& mesh) {
    return rotate(v, inverse(Rot(mesh.angle)));
}
inline Vec2 point_into_local(Vec2 p, const Transform& mesh) {
    return rotate(p - mesh.position, inverse(Rot(mesh.angle)));
}
inline Vec2 point_into_world(Vec2 p, const Transform& mesh) {
    return mesh.position + rotate(p, Rot(mesh.angle));
}
inline Vec2 closest_on_outline(std::span<const Vec2> vertices, Vec2 p) {
    Vec2 best = vertices.front();
    float best_distance = std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 point =
            closest_on_segment(p, vertices[i], vertices[(i + 1) % vertices.size()]);
        const float distance = (point - p).length_squared();
        if (distance < best_distance) {
            best_distance = distance;
            best = point;
        }
    }
    return best;
}
inline Vec2 closest_on_segment_pair(Vec2 p1, Vec2 q1, Vec2 p2, Vec2 q2, Vec2& on_first) {
    const Vec2 d1 = q1 - p1, d2 = q2 - p2, offset = p1 - p2;
    const float a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(offset);
    float s = 0, t = 0;
    if (a <= 1.0e-12f && e <= 1.0e-12f) {
        on_first = p1;
        return p2;
    }
    if (a <= 1.0e-12f) {
        t = std::clamp(f / e, 0.0f, 1.0f);
    } else {
        const float c = d1.dot(offset);
        if (e <= 1.0e-12f) {
            s = std::clamp(-c / a, 0.0f, 1.0f);
        } else {
            const float b = d1.dot(d2);
            const float denominator = a * e - b * b;
            s = denominator > 1.0e-12f ? std::clamp((b * f - c * e) / denominator, 0.0f, 1.0f)
                                       : 0.0f;
            t = (b * s + f) / e;
            if (t < 0) {
                t = 0;
                s = std::clamp(-c / a, 0.0f, 1.0f);
            } else if (t > 1) {
                t = 1;
                s = std::clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    on_first = p1 + d1 * s;
    return p2 + d2 * t;
}
// A circle against one boundary edge. The closest point of a segment is either
// its interior -- where the normal is the segment's own, exact -- or one of its
// endpoints, and there the normal is the direction from that vertex to the
// circle. Both triangles that share the vertex report the same direction, so a
// rolling ball does not trip over the seam.
inline bool circle_edge(const Circle& circle, const Transform& tc, Vec2 a, Vec2 b, Contact& out) {
    const Vec2 closest = closest_on_segment(tc.position, a, b);
    const Vec2 delta = tc.position - closest;
    const float distance = delta.length();
    if (distance >= circle.radius)
        return false;
    out.normal = distance > 1.0e-6f ? delta / distance : Vec2{0, 1};
    out.penetration = circle.radius - distance;
    out.point = closest;
    return true;
}
inline bool capsule_edge(const Capsule& capsule, const Transform& tc, Vec2 a, Vec2 b,
                         Contact& out) {
    const CapsuleSegment segment = capsule_segment(capsule, tc);
    Vec2 on_capsule{};
    const Vec2 on_edge = closest_on_segment_pair(segment.a, segment.b, a, b, on_capsule);
    const Vec2 delta = on_capsule - on_edge;
    const float distance = delta.length();
    if (distance >= capsule.radius)
        return false;
    out.normal = distance > 1.0e-6f ? delta / distance : Vec2{0, 1};
    out.penetration = capsule.radius - distance;
    out.point = on_edge;
    return true;
}
// An outline (a box or a polygon) against one boundary edge, in one frame.
//
// The candidate axes are the outline's own faces, the edge's face, and the
// direction from each endpoint of the edge to the closest point of the outline.
// The outline's faces are needed to notice separation; the edge's face is the
// surface response; and the two vertex directions are what the two triangles of
// a seam agree on. Axes are oriented from the edge towards the outline, so the
// winning axis *is* the contact normal pointing out of the mesh.
inline bool outline_edge(std::span<const Vec2> vertices, Vec2 a, Vec2 b, Vec2 outward,
                         Contact& out) {
    const Vec2 edge = b - a;
    if (vertices.size() < 2 || edge.length_squared() <= 1.0e-12f)
        return false;
    Vec2 outline_center{};
    for (const Vec2 v : vertices)
        outline_center += v;
    outline_center /= float(vertices.size());
    const Vec2 segment_center = (a + b) * 0.5f;
    // Which way is "out of the material". The offset between the two centres is
    // the honest answer; it degenerates only when the outline is centred exactly
    // on the edge, and then the edge's own outward normal is the reference.
    const Vec2 offset = outline_center - segment_center;
    const Vec2 reference = offset.length_squared() > 1.0e-8f ? offset : outward;
    const Vec2 segment[2] = {a, b};
    const auto project = [](std::span<const Vec2> points, Vec2 axis) {
        float low = points.front().dot(axis), high = low;
        for (const Vec2 p : points) {
            const float value = p.dot(axis);
            low = std::min(low, value);
            high = std::max(high, value);
        }
        return std::pair<float, float>{low, high};
    };
    float best_gap = -std::numeric_limits<float>::max();
    Vec2 best_axis{};
    Vec2 deepest = vertices.front();
    auto consider = [&](Vec2 axis) {
        if (axis.length_squared() <= 1.0e-12f)
            return;
        axis = axis.normalized();
        if (axis.dot(reference) < 0)
            axis = -axis;
        const auto outline_span = project(vertices, axis);
        const auto segment_span = project(segment, axis);
        const float gap = outline_span.first - segment_span.second;
        if (gap <= best_gap)
            return;
        best_gap = gap;
        best_axis = axis;
        deepest = vertices.front();
        for (const Vec2 v : vertices)
            if (v.dot(axis) < deepest.dot(axis))
                deepest = v;
    };
    consider(Vec2{edge.y, -edge.x});
    consider(closest_on_outline(vertices, a) - a);
    consider(closest_on_outline(vertices, b) - b);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec2 face = vertices[(i + 1) % vertices.size()] - vertices[i];
        consider(Vec2{face.y, -face.x});
    }
    if (best_gap > 0)
        return false;
    out.normal = best_axis;
    out.penetration = -best_gap;
    const Vec2 on_segment = closest_on_segment(deepest, a, b);
    out.point = on_segment + best_axis * (best_gap * 0.5f);
    return true;
}
inline bool edge_contact(const Mesh::Edge& edge, const Shape& other, const Transform& transform,
                         Contact& out) {
    bool touching;
    if (const auto* circle = std::get_if<Circle>(&other))
        touching = circle_edge(*circle, transform, edge.a, edge.b, out);
    else if (const auto* capsule = std::get_if<Capsule>(&other))
        touching = capsule_edge(*capsule, transform, edge.a, edge.b, out);
    else if (std::holds_alternative<Box>(other) || std::holds_alternative<Polygon>(other)) {
        const shape_detail::WorldVertices storage(other, transform);
        touching = outline_edge(storage.view(), edge.a, edge.b, edge.normal_a, out);
    } else {
        return false;
    }
    if (!touching)
        return false;
    // The edge's own outward normal is the only direction that is a face of the
    // mesh, and for a body outside the material the primitive tests above
    // already agree with it: the direction from the edge to a body outside is
    // the direction out of the material. For a body that has gone *inside* they
    // disagree -- the direction from the edge to its centre points further in --
    // and following it would drive the body deeper through the surface rather
    // than ejecting it through the nearest face. The contact is with this edge
    // either way; only the direction needs correcting.
    //
    // Only an edge's *interior* has the edge's own face as its surface. At an
    // endpoint the surface direction is the vertex's, which is what the
    // primitive test computed and what the neighbouring edges agree on -- that
    // agreement is what makes a seam smooth -- so a contact at an endpoint is
    // left exactly as it was found.
    if (out.normal.dot(edge.outward) < 0) {
        const Vec2 touch = closest_on_segment(out.point, edge.a, edge.b);
        const float span = (edge.b - edge.a).length();
        if ((touch - edge.a).length() > 1.0e-4f * span &&
            (touch - edge.b).length() > 1.0e-4f * span)
            out.normal = edge.outward;
    }
    return true;
}
// The mesh's surface contacts with `other` (already in the mesh's own frame),
// grouped by the surface they are on.
//
// A mesh floor is a strip of collinear boundary edges, and the seam between two
// of them is not a feature: the same surface continues past it. Reporting only
// the deepest edge gives a body straddling a seam a manifold on one side of the
// seam, so its support sits off-centre and it tips towards the unsupported half.
// The group is therefore the deepest edge *plus* every other boundary edge the
// body touches that lies on the deepest edge's own line -- collinear, so a
// contact on one of them is a contact with the same surface -- and the body is
// supported across all of them.
struct EdgeHit {
    std::uint32_t edge{0xffffffffu};
    Contact contact{};
    bool found{false};
};
constexpr int kSeamEdges = 4;
struct EdgeGroup {
    std::uint32_t edge[kSeamEdges]{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu};
    Contact contact{};
    int count{0};
};
inline EdgeHit deepest_edge(const Mesh& mesh, const Shape& other, const Transform& local,
                            int& touching) {
    EdgeHit best{};
    const Mesh::Accel& tree = mesh.acceleration();
    const AABB query = compute_aabb(other, local);
    mesh.visit(query, [&](std::uint32_t triangle) {
        ++mesh_triangle_tests();
        for (std::uint32_t e = tree.triangle_edges[triangle]; e != 0xffffffffu;
             e = tree.boundary[e].next_in_triangle) {
            Contact candidate;
            if (!edge_contact(tree.boundary[e], other, local, candidate))
                continue;
            ++touching;
            if (!best.found || candidate.penetration > best.contact.penetration) {
                best.found = true;
                best.edge = e;
                best.contact = candidate;
            }
        }
    });
    return best;
}
inline EdgeGroup seam_edges(const Mesh& mesh, const Shape& other, const Transform& local) {
    EdgeGroup group;
    const Mesh::Accel& tree = mesh.acceleration();
    if (tree.boundary.empty())
        return group;
    int touching = 0;
    const EdgeHit best = deepest_edge(mesh, other, local, touching);
    if (!best.found)
        return group;
    group.edge[0] = best.edge;
    group.contact = best.contact;
    group.count = 1;
    // One edge is the answer, and it is the common case -- a body resting inside
    // a single triangle. Only a body straddling a seam pays for the second
    // descent, and then only to collect the collinear run it is standing on.
    if (touching < 2 || kSeamEdges < 2)
        return group;
    const Vec2 origin = tree.boundary[best.edge].a;
    const AABB query = compute_aabb(other, local);
    mesh.visit(query, [&](std::uint32_t triangle) {
        ++mesh_triangle_tests();
        for (std::uint32_t e = tree.triangle_edges[triangle]; e != 0xffffffffu;
             e = tree.boundary[e].next_in_triangle) {
            if (e == best.edge || group.count >= kSeamEdges)
                continue;
            const Mesh::Edge& edge = tree.boundary[e];
            // On the same line: both endpoints are the same distance from the
            // deepest edge's line as its own endpoints are, which is zero.
            const Vec2 normal = best.contact.normal;
            if (std::abs((edge.a - origin).dot(normal)) > 1.0e-3f ||
                std::abs((edge.b - origin).dot(normal)) > 1.0e-3f)
                continue;
            Contact candidate;
            if (!edge_contact(edge, other, local, candidate))
                continue;
            if (candidate.normal.dot(normal) < 0.999f)
                continue;
            group.edge[group.count++] = e;
        }
    });
    return group;
}
} // namespace mesh_detail

inline bool mesh_against(const Mesh& mesh, const Transform& tm, const Shape& other,
                         const Transform& to, Contact& c);

inline bool test(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb, Contact& c) {
    if (const auto* ca = std::get_if<Circle>(&a)) {
        if (const auto* cb = std::get_if<Circle>(&b)) return circle_circle(*ca, ta, *cb, tb, c);
        // Both orders of every mixed pair have to give the same answer. They did
        // not: a circle on the left of a capsule reached the exact routine while
        // a capsule on the left of a circle fell through to polygon SAT, which a
        // Circle cannot satisfy, and reported no contact at all. Which one a
        // broad-phase pair got came down to body creation order.
        if (const auto* cb = std::get_if<Capsule>(&b)) {
            if (!capsule_circle(*cb, tb, *ca, ta, c)) return false;
            c.normal = -c.normal;
            return true;
        }
        if (std::holds_alternative<Box>(b) || std::holds_alternative<Polygon>(b)) return circle_polygon(*ca, ta, b, tb, c);
    } else if (const auto* cap = std::get_if<Capsule>(&a)) {
        if (const auto* cb = std::get_if<Circle>(&b)) return capsule_circle(*cap, ta, *cb, tb, c);
        if (const auto* cb = std::get_if<Capsule>(&b)) return capsule_capsule(*cap, ta, *cb, tb, c);
        if (std::holds_alternative<Box>(b) || std::holds_alternative<Polygon>(b))
            return capsule_convex(*cap, ta, b, tb, c);
    } else if (const auto* cb = std::get_if<Circle>(&b)) {
        if (std::holds_alternative<Box>(a) || std::holds_alternative<Polygon>(a)) {
            if (!circle_polygon(*cb, tb, a, ta, c)) return false;
            c.normal = -c.normal;
            return true;
        }
    } else if (const auto* cb = std::get_if<Capsule>(&b)) {
        if (std::holds_alternative<Box>(a) || std::holds_alternative<Polygon>(a)) {
            if (!capsule_convex(*cb, tb, a, ta, c)) return false;
            c.normal = -c.normal;
            return true;
        }
    }
    // A mesh is answered through its spatial index, never by walking the
    // triangles: the boundary contacts first, and the triangles themselves only
    // where the mesh has no boundary left to offer (a closed shell) or where the
    // shape is inside one.
    if (std::holds_alternative<Mesh>(a)) return mesh_against(std::get<Mesh>(a), ta, b, tb, c);
    if (std::holds_alternative<Mesh>(b)) {
        if (!mesh_against(std::get<Mesh>(b), tb, a, ta, c)) return false;
        c.normal = -c.normal;
        return true;
    }
    return polygon_polygon(a, ta, b, tb, c);
}

// The mesh as the "A" side: `c.normal` comes back pointing from the mesh into
// `other`, and `c.point` is in world space.
inline bool mesh_against(const Mesh& mesh, const Transform& tm, const Shape& other,
                         const Transform& to, Contact& c) {
    const Mesh::Accel& tree = mesh.acceleration();
    if (tree.empty)
        return false;
    const Transform local = mesh_detail::into_local(to, tm);
    const mesh_detail::EdgeGroup group = mesh_detail::seam_edges(mesh, other, local);
    if (group.count == 0) {
        // The boundary offered nothing. Either the mesh is a closed shell, whose
        // material is its triangles, or the shape is inside the surface rather
        // than on it; both are answered by the triangles themselves, through the
        // same hierarchy.
        const AABB query = compute_aabb(other, local);
        float deepest = 0;
        bool found = false;
        mesh.visit(query, [&](std::uint32_t triangle) {
            ++mesh_triangle_tests();
            const Shape& face = mesh.triangles[triangle];
            Contact candidate;
            if (!test(other, local, face, {}, candidate))
                return;
            if (!found || candidate.penetration > deepest) {
                deepest = candidate.penetration;
                found = true;
                c = candidate;
            }
        });
        if (!found)
            return false;
        c.normal = -c.normal; // test() answered "other -> triangle".
    } else {
        c = group.contact;
    }
    c.normal = rotate(c.normal, Rot(tm.angle));
    c.point = mesh_detail::point_into_world(c.point, tm);
    return true;
}

// Clip the incident edge against the reference face. Feature IDs are edge/endpoint
// IDs (including generated intersections), independent of world coordinates.
struct ManifoldPoint { Vec2 point{}; float separation{}; unsigned feature{}; };
struct Manifold { Vec2 normal{}; ManifoldPoint points[2]{}; int count{}; };
inline Manifold mesh_manifold(const Mesh& mesh, const Transform& tm, const Shape& other,
                              const Transform& to, Vec2 normal, Vec2 fallback, float separation,
                              float margin);

inline Manifold contact_manifold(const Shape& a, const Transform& ta, const Shape& b,
                                 const Transform& tb, Vec2 normal, Vec2 fallback,
                                 float separation, float margin = 0.002f) {
    // A mesh has to be clipped against the one edge the contact is actually on,
    // for two reasons: the mesh as a whole is not a shape a face clipper can
    // use, and clipping against the full outline would put the second point of a
    // box resting on a long floor somewhere in the middle of the mesh rather
    // than under the box. `normal` is the A-to-B direction either way, which the
    // mesh path re-expresses in the mesh's own frame.
    if (std::holds_alternative<Mesh>(a))
        return mesh_manifold(std::get<Mesh>(a), ta, b, tb, normal, fallback, separation, margin);
    if (std::holds_alternative<Mesh>(b)) {
        Manifold flipped =
            mesh_manifold(std::get<Mesh>(b), tb, a, ta, -normal, fallback, separation, margin);
        flipped.normal = -flipped.normal;
        return flipped;
    }
    Manifold m; m.normal = normal;
    // A capsule's flat side *is* its segment. Clipping the tessellated outline
    // instead would place the second contact point on whichever facet the
    // eighteen-gon happened to have there, which is exactly the sliding that
    // treating a capsule as a segment is meant to remove. Every other shape
    // still comes from `WorldVertices`; a capsule never builds one, so a capsule
    // contact costs no heap allocation.
    std::array<Vec2, 2> capsule_a{}, capsule_b{};
    std::optional<shape_detail::WorldVertices> storage_a, storage_b;
    auto vertices_of = [](const Shape& shape, const Transform& t, auto& storage,
                          std::array<Vec2, 2>& capsule_storage) -> std::span<const Vec2> {
        if (const auto* capsule = std::get_if<Capsule>(&shape)) {
            const CapsuleSegment segment = capsule_segment(*capsule, t);
            capsule_storage = {segment.a, segment.b};
            return {capsule_storage.data(), capsule_storage.size()};
        }
        storage.emplace(shape, t);
        return storage->view();
    };
    const auto va = vertices_of(a, ta, storage_a, capsule_a);
    const auto vb = vertices_of(b, tb, storage_b, capsule_b);
    if (va.size() < 2 || vb.size() < 2) {
        m.points[0] = {fallback, separation, 0}; m.count = 1; return m;
    }
    auto face = [](std::span<const Vec2> v, Vec2 n) {
        int best = 0;
        float alignment = -2;
        if (v.size() == 2) {
            // A segment has one face and no "outward" side of its own: the
            // centre-relative test below is degenerate for it, so the requested
            // direction is what decides. Without this the flat side of a capsule
            // scored a negative alignment against a box face and lost a contest
            // it should have won, or won one it should have lost, depending on
            // the floating-point sign of a zero.
            Vec2 out = Vec2{(v[1] - v[0]).y, -(v[1] - v[0]).x}.normalized();
            if (out.dot(n) < 0)
                out = -out;
            return std::make_pair(0, out.dot(n));
        }
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
    // A capsule's face list is its *medial axis*, not its surface: the barrel is
    // `radius` away from the material in every direction, and the segment ends
    // are the centres of the round caps. Measuring the gap against the medial
    // axis without accounting for that made a barrel lying flat on a floor
    // report a gap of exactly `radius`, fail the `d > margin` test below, and
    // produce *no constraint at all* -- the pair was silent while the capsule
    // sank through a floor it was geometrically resting on.
    const auto capsule_radius = [](const Shape& shape) {
        const auto* capsule = std::get_if<Capsule>(&shape);
        return capsule ? capsule->radius : 0.0f;
    };
    const float radius_ref = flip ? capsule_radius(b) : capsule_radius(a);
    const float radius_inc = flip ? capsule_radius(a) : capsule_radius(b);
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
        // The gap between the two *surfaces*, and the point halfway across it.
        // With no capsule involved both radii are zero and this is the original
        // face-plane-to-point measure.
        const float d = (p - r0).dot(n) - radius_ref - radius_inc;
        if (d > margin)
            continue;
        unsigned endpoint = t < 1e-6f ? 0 : t > 1 - 1e-6f ? 1 : 2 + k;
        m.points[m.count++] = {p - n * (d * 0.5f + radius_inc), d,
                               unsigned((flip ? 1 : 0) << 24 | ri << 16 | ii << 8) | endpoint};
    }
    return m;
}

// The manifold of `other` against the one boundary edge of the mesh that the
// contact is on. Everything is done in the mesh's own frame, because that is
// where the index and the boundary live, and converted back at the end.
inline Manifold mesh_manifold(const Mesh& mesh, const Transform& tm, const Shape& other,
                              const Transform& to, Vec2 normal, Vec2 fallback, float separation,
                              float margin) {
    Manifold result;
    result.normal = normal;
    const Mesh::Accel& tree = mesh.acceleration();
    if (tree.empty) {
        result.points[0] = {fallback, separation, 0};
        result.count = 1;
        return result;
    }
    const Transform local = mesh_detail::into_local(to, tm);
    const mesh_detail::EdgeGroup group = mesh_detail::seam_edges(mesh, other, local);
    if (group.count == 0) {
        result.points[0] = {fallback, separation, 0};
        result.count = 1;
        return result;
    }
    // The clipper wants a polygon, and a boundary edge *is* one: a two-vertex
    // outline, which the face selection below already treats as a segment with a
    // single face. Clipping against it is what turns a box resting on a mesh
    // floor into two contact points instead of one, and it keeps them under the
    // box instead of at the middle of the terrain.
    //
    // Every edge of the group is clipped and the points are pooled, so a body
    // standing across a seam collects the support from both of its sides rather
    // than from whichever side the descent happened to reach first.
    const Vec2 local_normal = mesh_detail::direction_into_local(normal, tm);
    const Vec2 local_fallback = mesh_detail::point_into_local(fallback, tm);
    ManifoldPoint pooled[2 * mesh_detail::kSeamEdges];
    int total = 0;
    for (int g = 0; g < group.count && total + 2 <= 2 * mesh_detail::kSeamEdges; ++g) {
        const Mesh::Edge& edge = tree.boundary[group.edge[g]];
        const Polygon segment{{edge.a, edge.b}};
        const Manifold clipped = contact_manifold(other, local, Shape{segment}, {},
                                                  -local_normal, // other -> the mesh's surface
                                                  local_fallback, separation, margin);
        for (int k = 0; k < clipped.count && total < 2 * mesh_detail::kSeamEdges; ++k) {
            pooled[total] = clipped.points[k];
            // The tag distinguishes points produced by different edges. Bits
            // 25-31 are free in the clipper's feature id (bit 24 is its flip
            // bit, and the rest is face and endpoint indices), and without it
            // two edges of the same run hand out the same id for different
            // points, which is what the warm start matches on.
            pooled[total].feature = (pooled[total].feature & 0x01ffffffu) |
                                    (std::uint32_t(g) << 25);
            ++total;
        }
    }
    if (total == 0) {
        result.points[0] = {fallback, separation, 0};
        result.count = 1;
        return result;
    }
    // Keep the pair that is furthest apart: on a flat surface that is the two
    // ends of the support, which is what a two-point manifold is for.
    int first = 0, second = -1;
    float widest = -1;
    for (int i = 0; i < total; ++i)
        for (int j = i + 1; j < total; ++j) {
            const float spread = (pooled[i].point - pooled[j].point).length_squared();
            if (spread > widest) {
                widest = spread;
                first = i;
                second = j;
            }
        }
    result.points[0] = pooled[first];
    result.count = 1;
    if (second >= 0) {
        result.points[1] = pooled[second];
        result.count = 2;
    }
    for (int k = 0; k < result.count; ++k)
        result.points[k].point = mesh_detail::point_into_world(result.points[k].point, tm);
    return result;
}

} // namespace butter::physics2d
