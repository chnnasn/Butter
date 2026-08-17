#pragma once

#include "quat.h"
#include "vec3.h"

namespace butter::math {

struct Transform {
    Vec3 position{0, 0, 0};
    Quat rotation{};

    static constexpr Transform identity() { return {}; }

    constexpr Transform() = default;
    constexpr Transform(const Vec3& position) : position(position) {}
    constexpr Transform(const Vec3& position, const Quat& rotation)
        : position(position), rotation(rotation) {}

    Vec3 transform_point(const Vec3& point) const {
        return position + rotation.rotate(point);
    }

    Vec3 transform_direction(const Vec3& direction) const {
        return rotation.rotate(direction);
    }

    Vec3 world_to_local(const Vec3& world_point) const {
        return rotation.conjugate().rotate(world_point - position);
    }

    Vec3 local_to_world(const Vec3& local_point) const {
        return position + rotation.rotate(local_point);
    }

    Transform inverse() const {
        const Quat inv_rot = rotation.conjugate();
        return {inv_rot.rotate(-position), inv_rot};
    }

    Transform operator*(const Transform& other) const {
        return {position + rotation.rotate(other.position),
                (rotation * other.rotation).normalized()};
    }
};

} // namespace butter::math
