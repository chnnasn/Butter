#pragma once

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "butter/core/body.h"
#include "butter/shapes/box.h"
#include "butter/shapes/capsule.h"
#include "butter/shapes/convex.h"
#include "butter/shapes/sphere.h"

namespace butter {

class World;

class BodyBuilder {
public:
    explicit BodyBuilder(World& world) : world_(world) {}

    BodyBuilder(BodyBuilder&&) = default;
    BodyBuilder& operator=(BodyBuilder&&) = default;
    BodyBuilder(const BodyBuilder&) = delete;
    BodyBuilder& operator=(const BodyBuilder&) = delete;

    BodyBuilder& dynamic() { type_ = BodyType::Dynamic; return *this; }
    BodyBuilder& kinematic() { type_ = BodyType::Kinematic; return *this; }
    BodyBuilder& static_body() { type_ = BodyType::Static; return *this; }

    BodyBuilder& at(float x, float y, float z) {
        position_ = {x, y, z};
        return *this;
    }
    BodyBuilder& at(const math::Vec3& pos) {
        position_ = pos;
        return *this;
    }

    BodyBuilder& rotated(float x, float y, float z) {
        rotation_ = math::Quat::from_euler(x, y, z);
        return *this;
    }
    BodyBuilder& rotated(const math::Quat& rotation) {
        rotation_ = rotation;
        return *this;
    }

    BodyBuilder& mass(float m) {
        mass_ = m;
        return *this;
    }
    BodyBuilder& density(float d) {
        density_ = d;
        material_.density = d;
        return *this;
    }

    BodyBuilder& moving(float vx, float vy, float vz) {
        velocity_ = {vx, vy, vz};
        return *this;
    }
    BodyBuilder& moving(const math::Vec3& velocity) {
        velocity_ = velocity;
        return *this;
    }

    BodyBuilder& box(float hx, float hy, float hz) {
        colliders_.push_back(std::make_unique<BoxCollider>(hx, hy, hz));
        return *this;
    }
    BodyBuilder& box(const math::Vec3& half_extents) {
        colliders_.push_back(std::make_unique<BoxCollider>(half_extents));
        return *this;
    }

    BodyBuilder& sphere(float radius) {
        colliders_.push_back(std::make_unique<SphereCollider>(radius));
        return *this;
    }

    BodyBuilder& capsule(float radius, float height) {
        colliders_.push_back(std::make_unique<CapsuleCollider>(radius, height));
        return *this;
    }

    BodyBuilder& convex(std::vector<math::Vec3> points) {
        colliders_.push_back(std::make_unique<ConvexCollider>(std::move(points)));
        return *this;
    }

    BodyBuilder& friction(float f) {
        material_.friction = f;
        return *this;
    }
    BodyBuilder& bounciness(float b) {
        material_.restitution = b;
        return *this;
    }
    BodyBuilder& trigger(bool is_trigger = true) {
        material_.is_trigger = is_trigger;
        return *this;
    }

    template <typename F>
    BodyBuilder& on_collision(F&& callback) {
        collision_callback_ = std::forward<F>(callback);
        return *this;
    }

    template <typename F>
    BodyBuilder& on_sleep(F&& callback) {
        sleep_callback_ = std::forward<F>(callback);
        return *this;
    }

    // Defined in world.h after World is complete.
    Body& build();
    operator Body&();

private:
    World& world_;
    BodyType type_{BodyType::Dynamic};
    math::Vec3 position_{0, 0, 0};
    math::Quat rotation_{};
    math::Vec3 velocity_{0, 0, 0};
    float mass_{1.0f};
    float density_{1.0f};
    Material material_;
    std::vector<std::unique_ptr<Collider>> colliders_;
    std::function<void(const CollisionEvent&)> collision_callback_;
    std::function<void()> sleep_callback_;

    friend class World;
};

} // namespace butter
