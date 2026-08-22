#pragma once

#include <utility>
#include <vector>

#include "butter/shapes/collider.h"

namespace butter {

// A point set interpreted as its convex hull by the support mapping.  The
// caller does not need to provide hull faces; GJK/EPA only needs extreme
// points in a direction.
class ConvexCollider : public Collider {
public:
    ConvexCollider() = default;

    explicit ConvexCollider(std::vector<math::Vec3> points)
        : points_(std::move(points)) {}

    Type type() const override { return Type::Convex; }

    ConvexCollider& add_point(const math::Vec3& point) {
        points_.push_back(point);
        return *this;
    }

    const std::vector<math::Vec3>& points() const { return points_; }

    math::AABB compute_aabb(const math::Transform& transform) const override {
        math::AABB aabb;
        const math::Transform world = transform * offset;
        for (const auto& point : points_) {
            aabb.encapsulate(world.transform_point(point));
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
    std::vector<math::Vec3> points_;
};

} // namespace butter
