#pragma once

#include <cmath>

#include "butter/shapes/collider.h"

namespace butter {

class CapsuleCollider : public Collider {
public:
    CapsuleCollider(float radius = 0.5f, float height = 1.0f)
        : radius_(radius), half_height_(height * 0.5f) {}

    Type type() const override { return Type::Capsule; }

    float radius() const { return radius_; }
    float height() const { return half_height_ * 2.0f; }
    float half_height() const { return half_height_; }

    CapsuleCollider& set_radius(float r) {
        radius_ = r;
        return *this;
    }
    CapsuleCollider& set_height(float h) {
        half_height_ = h * 0.5f;
        return *this;
    }

    math::AABB compute_aabb(const math::Transform& transform) const override {
        const math::Vec3 center = transform.position + offset.position;
        const math::Mat3 rotation = (transform.rotation * offset.rotation).to_mat3();
        // Local capsule axis is Y.
        const math::Vec3 axis{rotation.m[0][1], rotation.m[1][1], rotation.m[2][1]};
        const math::Vec3 extent{
            std::abs(axis.x) * half_height_ + radius_,
            std::abs(axis.y) * half_height_ + radius_,
            std::abs(axis.z) * half_height_ + radius_,
        };
        return {center - extent, center + extent};
    }

    math::Mat3 compute_inertia_tensor(float mass) const override {
        const float h = half_height_ * 2.0f;
        const float r = radius_;
        const float cylinderish = mass * (3.0f * r * r + h * h) / 12.0f;
        const float along_axis = 0.4f * mass * r * r;
        return math::Mat3::diagonal(cylinderish, along_axis, cylinderish);
    }

private:
    float radius_;
    float half_height_;
};

} // namespace butter
