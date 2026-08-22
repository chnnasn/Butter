#pragma once

#include "butter/shapes/collider.h"

namespace butter {

class SphereCollider : public Collider {
public:
    explicit SphereCollider(float radius = 0.5f) : radius_(radius) {}

    Type type() const override { return Type::Sphere; }

    float radius() const { return radius_; }
    SphereCollider& set_radius(float r) {
        radius_ = r;
        return *this;
    }

    math::AABB compute_aabb(const math::Transform& transform) const override {
        const math::Vec3 center = transform.transform_point(offset.position);
        const math::Vec3 delta{radius_, radius_, radius_};
        return {center - delta, center + delta};
    }

    math::Mat3 compute_inertia_tensor(float mass) const override {
        const float i = 0.4f * mass * radius_ * radius_;
        return math::Mat3::diagonal(i, i, i);
    }

private:
    float radius_;
};

} // namespace butter
