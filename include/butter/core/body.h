#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "butter/core/events.h"
#include "butter/core/material.h"
#include "butter/math/mat3.h"
#include "butter/math/quat.h"
#include "butter/math/vec3.h"
#include "butter/shapes/collider.h"

namespace butter {

enum class BodyType {
    Static,
    Dynamic,
    Kinematic,
};

class Body {
public:
    math::Vec3 position{};
    math::Quat rotation{};
    math::Vec3 velocity{};
    math::Vec3 angular_velocity{};

    std::function<void()> on_sleep;
    std::function<void()> on_wake;
    std::function<void(const CollisionEvent&)> on_collision;

    // Collision filtering: a pair collides only when both bodies accept each
    // other's group. Particles can use this to ignore other particles while
    // still colliding with the central bodies.
    std::uint32_t collision_group{1};
    std::uint32_t collision_mask{0xFFFFFFFF};

    // Damping makes flying debris lose energy to air resistance and settle
    // instead of sliding forever. Zero keeps the classic vacuum behavior.
    float linear_damping{0.0f};
    float angular_damping{0.0f};

    Body() = default;
    explicit Body(BodyType type) : type_(type) {}
    virtual ~Body() = default;

    float mass() const { return mass_; }
    float inverse_mass() const { return inverse_mass_; }
    BodyType body_type() const { return type_; }
    bool is_dynamic() const { return type_ == BodyType::Dynamic; }
    bool is_static() const { return type_ == BodyType::Static; }
    bool is_kinematic() const { return type_ == BodyType::Kinematic; }
    bool is_sleeping() const { return sleeping_; }
    bool is_destroyed() const { return destroyed_; }
    bool gravity_enabled() const { return gravity_enabled_; }

    const std::string& name() const { return name_; }
    Body& set_name(std::string name) {
        name_ = std::move(name);
        return *this;
    }

    Body& set_type(BodyType type) {
        type_ = type;
        update_mass_properties();
        return *this;
    }

    Body& set_mass(float m) {
        mass_ = m > 0 ? m : 0.0f;
        update_mass_properties();
        return *this;
    }

    Body& set_gravity_enabled(bool enabled) {
        gravity_enabled_ = enabled;
        return *this;
    }

    Body& set_sleeping(bool sleeping) {
        if (sleeping_ == sleeping) return *this;
        sleeping_ = sleeping;
        if (sleeping_) {
            if (on_sleep) on_sleep();
        } else {
            if (on_wake) on_wake();
        }
        return *this;
    }

    Body& wake_up() { return set_sleeping(false); }
    Body& sleep() { return set_sleeping(true); }

    void apply_force(const math::Vec3& force) {
        wake_up();
        force_accumulator_ += force;
    }
    void apply_impulse(const math::Vec3& impulse) {
        wake_up();
        if (inverse_mass_ > 0) {
            velocity += impulse * inverse_mass_;
        }
    }

    void apply_impulse_at_point(const math::Vec3& impulse, const math::Vec3& point) {
        wake_up();
        if (inverse_mass_ <= 0) return;

        velocity += impulse * inverse_mass_;
        const math::Vec3 r = point - position;
        angular_velocity += inverse_inertia_tensor_ * r.cross(impulse);
    }

    void apply_torque(const math::Vec3& torque) {
        wake_up();
        torque_accumulator_ += torque;
    }

    void apply_force_at_point(const math::Vec3& force, const math::Vec3& point) {
        apply_force(force);
        apply_torque((point - position).cross(force));
    }

    void add_velocity(const math::Vec3& delta) { velocity += delta; }
    void add_angular_velocity(const math::Vec3& delta) { angular_velocity += delta; }

    template <typename T, typename... Args>
    T& add_collider(Args&&... args) {
        auto collider = std::make_unique<T>(std::forward<Args>(args)...);
        T* ptr = collider.get();
        colliders_.push_back(std::move(collider));
        update_mass_properties();
        return *ptr;
    }

    Collider& add_collider(std::unique_ptr<Collider> collider) {
        Collider* ptr = collider.get();
        colliders_.push_back(std::move(collider));
        update_mass_properties();
        return *ptr;
    }

    std::vector<std::unique_ptr<Collider>>& colliders() { return colliders_; }
    const std::vector<std::unique_ptr<Collider>>& colliders() const { return colliders_; }

    void clear_forces() {
        force_accumulator_ = {};
        torque_accumulator_ = {};
    }

    void mark_for_destruction() { destroyed_ = true; }

    const math::Mat3& inertia_tensor() const { return inertia_tensor_; }
    const math::Mat3& inverse_inertia_tensor() const { return inverse_inertia_tensor_; }

    void update_mass_properties() {
        if (mass_ <= 0.0f || type_ == BodyType::Static) {
            inverse_mass_ = 0.0f;
            inertia_tensor_ = math::Mat3::zero();
            inverse_inertia_tensor_ = math::Mat3::zero();
            return;
        }

        inverse_mass_ = 1.0f / mass_;
        math::Mat3 total = math::Mat3::zero();
        if (colliders_.empty()) {
            total = math::Mat3::scale(mass_ * 0.1f);
        } else {
            const float share = mass_ / static_cast<float>(colliders_.size());
            for (const auto& collider : colliders_) {
                const math::Mat3 c = collider->compute_inertia_tensor(share);
                for (int row = 0; row < 3; ++row) {
                    for (int col = 0; col < 3; ++col) {
                        total.m[row][col] += c.m[row][col];
                    }
                }
            }
        }

        inertia_tensor_ = total;
        inverse_inertia_tensor_ = math::Mat3::zero();
        for (int i = 0; i < 3; ++i) {
            if (total.m[i][i] != 0.0f) {
                inverse_inertia_tensor_.m[i][i] = 1.0f / total.m[i][i];
            }
        }
    }

private:
    BodyType type_{BodyType::Dynamic};
    float mass_{1.0f};
    float inverse_mass_{1.0f};
    math::Mat3 inertia_tensor_{};
    math::Mat3 inverse_inertia_tensor_{};
    math::Vec3 force_accumulator_{};
    math::Vec3 torque_accumulator_{};
    bool sleeping_{false};
    bool gravity_enabled_{true};
    bool destroyed_{false};
    std::string name_;
    std::vector<std::unique_ptr<Collider>> colliders_;

    friend class World;
    friend class BodyBuilder;
};

inline float CollisionEvent::relative_velocity() const {
    return (body_a.velocity - body_b.velocity).dot(normal);
}

} // namespace butter
