#pragma once

#include <cmath>

#include "butter/shapes/collider.h"

namespace butter {

class BoxCollider : public Collider {
public:
    BoxCollider(float hx = 0.5f, float hy = 0.5f, float hz = 0.5f)
        : half_extents_(hx, hy, hz) {}

    explicit BoxCollider(const math::Vec3& half_extents)
        : half_extents_(half_extents) {}

    Type type() const override { return Type::Box; }

    math::Vec3 half_extents() const { return half_extents_; }
    BoxCollider& set_half_extents(const math::Vec3& extents) {
        half_extents_ = extents;
        return *this;
    }
    BoxCollider& set_half_extents(float hx, float hy, float hz) {
        half_extents_ = {hx, hy, hz};
        return *this;
    }

    math::AABB compute_aabb(const math::Transform& transform) const override {
        const math::Vec3 center = transform.transform_point(offset.position);
        const math::Mat3 rotation = (transform.rotation * offset.rotation).to_mat3();
        const math::Vec3& h = half_extents_;

        // Mat3 is row-major and maps local coordinates to world coordinates.
        // Each world component therefore uses one *row* of the rotation
        // matrix.  Iterating columns underestimates/overestimates the AABB
        // for non-symmetric rotations and can make the broad phase miss a
        // real contact.
        const math::Vec3 extent{
            std::abs(rotation.m[0][0]) * h.x +
                std::abs(rotation.m[0][1]) * h.y +
                std::abs(rotation.m[0][2]) * h.z,
            std::abs(rotation.m[1][0]) * h.x +
                std::abs(rotation.m[1][1]) * h.y +
                std::abs(rotation.m[1][2]) * h.z,
            std::abs(rotation.m[2][0]) * h.x +
                std::abs(rotation.m[2][1]) * h.y +
                std::abs(rotation.m[2][2]) * h.z,
        };
        return {center - extent, center + extent};
    }

    math::Mat3 compute_inertia_tensor(float mass) const override {
        const float hx = half_extents_.x;
        const float hy = half_extents_.y;
        const float hz = half_extents_.z;
        return math::Mat3::diagonal(
            mass * (hy * hy + hz * hz) / 3.0f,
            mass * (hx * hx + hz * hz) / 3.0f,
            mass * (hx * hx + hy * hy) / 3.0f);
    }

private:
    math::Vec3 half_extents_;
};

} // namespace butter
