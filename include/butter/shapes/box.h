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
        const math::Vec3 center = transform.position + offset.position;
        const math::Mat3 rotation = (transform.rotation * offset.rotation).to_mat3();
        const math::Vec3& h = half_extents_;

        math::Vec3 extent;
        for (int col = 0; col < 3; ++col) {
            extent[col] = std::abs(rotation.m[0][col]) * h.x +
                          std::abs(rotation.m[1][col]) * h.y +
                          std::abs(rotation.m[2][col]) * h.z;
        }
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
