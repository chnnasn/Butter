#pragma once

#include "butter/constraints/joint.h"

namespace butter {

// A compact point-to-point hinge. It keeps the body anchor points together and
// damps their relative linear motion. Full angular limits are intentionally
// left as a future extension.
class HingeJoint : public Joint {
public:
    math::Vec3 anchor_a{};
    math::Vec3 anchor_b{};
    float stiffness{30.0f};
    float damping{8.0f};

    HingeJoint(Body& a, Body& b, const math::Vec3& world_anchor)
        : anchor_a(world_anchor - a.position), anchor_b(world_anchor - b.position) {
        body_a = &a;
        body_b = &b;
    }

    HingeJoint(Body& a, Body& b) : HingeJoint(a, b, (a.position + b.position) * 0.5f) {}

    void solve(float dt) override {
        if (!is_valid()) return;
        const float inv_mass = body_a->inverse_mass() + body_b->inverse_mass();
        if (inv_mass <= 0) return;

        const math::Vec3 world_a = body_a->position + body_a->rotation.rotate(anchor_a);
        const math::Vec3 world_b = body_b->position + body_b->rotation.rotate(anchor_b);
        const math::Vec3 delta = world_b - world_a;
        const float distance = delta.length();
        if (distance < 1.0e-6f) return;
        const math::Vec3 normal = delta / distance;

        const float relative_velocity = (body_b->velocity - body_a->velocity).dot(normal);
        const float impulse_magnitude = (-distance * stiffness - relative_velocity * damping) * dt;

        const math::Vec3 impulse = normal * impulse_magnitude;
        body_a->apply_impulse(-impulse);
        body_b->apply_impulse(impulse);
    }
};

} // namespace butter
