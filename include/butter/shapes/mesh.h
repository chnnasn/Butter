#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include "butter/shapes/collider.h"

namespace butter {

// A lightweight indexed triangle mesh collider. It is typically used as a
// static surface (the narrow phase is triangle based), but can also be
// attached to a dynamic body; it is intentionally not a rendering mesh.
class MeshCollider : public Collider {
public:
    struct Triangle {
        std::uint32_t a{0};
        std::uint32_t b{0};
        std::uint32_t c{0};
    };

    MeshCollider() = default;
    MeshCollider(std::vector<math::Vec3> vertices,
                 std::vector<Triangle> triangles)
        : vertices_(std::move(vertices)), triangles_(std::move(triangles)) {}

    Type type() const override { return Type::Mesh; }

    MeshCollider& set_vertices(std::vector<math::Vec3> vertices) {
        vertices_ = std::move(vertices);
        return *this;
    }
    MeshCollider& set_triangles(std::vector<Triangle> triangles) {
        triangles_ = std::move(triangles);
        return *this;
    }
    MeshCollider& add_triangle(std::uint32_t a, std::uint32_t b,
                               std::uint32_t c) {
        triangles_.push_back({a, b, c});
        return *this;
    }

    const std::vector<math::Vec3>& vertices() const { return vertices_; }
    const std::vector<Triangle>& triangles() const { return triangles_; }

    math::AABB compute_aabb(const math::Transform& transform) const override {
        math::AABB aabb;
        const math::Transform world = transform * offset;
        for (const auto& vertex : vertices_) {
            aabb.encapsulate(world.transform_point(vertex));
        }
        return aabb;
    }

    math::Mat3 compute_inertia_tensor(float mass) const override {
        const math::AABB local = compute_aabb(math::Transform::identity());
        if (local.is_empty()) return math::Mat3::scale(mass * 0.1f);
        const math::Vec3 h = local.half_extents();
        return math::Mat3::diagonal(
            mass * (h.y * h.y + h.z * h.z) / 3.0f,
            mass * (h.x * h.x + h.z * h.z) / 3.0f,
            mass * (h.x * h.x + h.y * h.y) / 3.0f);
    }

private:
    std::vector<math::Vec3> vertices_;
    std::vector<Triangle> triangles_;
};

} // namespace butter
