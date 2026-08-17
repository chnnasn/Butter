#pragma once

#include "butter/core/material.h"
#include "butter/math/aabb.h"
#include "butter/math/mat3.h"
#include "butter/math/transform.h"

namespace butter {

struct ContactInfo {
    math::Vec3 point{};
    math::Vec3 normal{0, 1, 0};
    float penetration{0};
    float impulse{0};
    bool hit{false};

    void reset() {
        point = {};
        normal = {0, 1, 0};
        penetration = 0;
        impulse = 0;
        hit = false;
    }
};

class Collider {
public:
    Material material;
    math::Transform offset{};

    enum class Type { Sphere, Box, Capsule, Convex, Mesh };

    virtual ~Collider() = default;

    virtual Type type() const = 0;
    virtual math::AABB compute_aabb(const math::Transform& transform) const = 0;
    virtual math::Mat3 compute_inertia_tensor(float mass) const = 0;

    static bool test(const Collider& a, const math::Transform& ta,
                     const Collider& b, const math::Transform& tb,
                     ContactInfo& contact);

    Collider& set_friction(float f) {
        material.friction = f;
        return *this;
    }

    Collider& set_bounciness(float b) {
        material.restitution = b;
        return *this;
    }

    Collider& set_offset(const math::Transform& t) {
        offset = t;
        return *this;
    }
};

} // namespace butter
