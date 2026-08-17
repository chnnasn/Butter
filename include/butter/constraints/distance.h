#pragma once

#include "butter/constraints/joint.h"

namespace butter {

class DistanceJoint : public Joint {
public:
    float target_distance{1.0f};
    float stiffness{20.0f};
    float damping{4.0f};

    DistanceJoint(Body& a, Body& b, float distance)
        : target_distance(distance) {
        body_a = &a;
        body_b = &b;
    }

    DistanceJoint(Body& a, Body& b)
        : target_distance((a.position - b.position).length()) {
        body_a = &a;
        body_b = &b;
    }

    void solve(float dt) override {
        if (!is_valid()) return;
        const float inv_mass = body_a->inverse_mass() + body_b->inverse_mass();
        if (inv_mass <= 0) return;

        math::Vec3 delta = body_b->position - body_a->position;
        const float distance = delta.length();
        if (distance < 1.0e-6f) return;
        const math::Vec3 normal = delta / distance;

        const float relative_velocity = (body_b->velocity - body_a->velocity).dot(normal);
        const float correction = (distance - target_distance) * stiffness;
        const float impulse_magnitude = (-correction - relative_velocity * damping) * dt;

        const math::Vec3 impulse = normal * impulse_magnitude;
        body_a->apply_impulse(-impulse);
        body_b->apply_impulse(impulse);
    }
};

} // namespace butter
