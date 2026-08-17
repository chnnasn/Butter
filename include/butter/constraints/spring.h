#pragma once

#include "butter/constraints/joint.h"

namespace butter {

class SpringJoint : public Joint {
public:
    float rest_length{1.0f};
    float stiffness{40.0f};
    float damping{2.0f};

    SpringJoint(Body& a, Body& b, float rest_length)
        : rest_length(rest_length) {
        body_a = &a;
        body_b = &b;
    }

    SpringJoint(Body& a, Body& b)
        : rest_length((a.position - b.position).length()) {
        body_a = &a;
        body_b = &b;
    }

    void solve(float dt) override {
        if (!is_valid()) return;

        math::Vec3 delta = body_b->position - body_a->position;
        const float distance = delta.length();
        if (distance < 1.0e-6f) return;
        const math::Vec3 normal = delta / distance;

        const float relative_velocity = (body_b->velocity - body_a->velocity).dot(normal);
        const float force = -stiffness * (distance - rest_length) - damping * relative_velocity;

        const math::Vec3 impulse = normal * (force * dt);
        body_a->apply_impulse(impulse);
        body_b->apply_impulse(-impulse);
    }
};

} // namespace butter
