#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#include "butter/constraints/distance.h"
#include "butter/constraints/hinge.h"
#include "butter/constraints/joint.h"
#include "butter/constraints/spring.h"
#include "butter/core/body_builder.h"
#include "butter/query/overlap.h"
#include "butter/query/raycast.h"
#include "butter/shapes/collision_test.h"

namespace butter {

namespace detail {

inline bool ray_sphere(const math::Vec3& origin, const math::Vec3& direction,
                       const math::Vec3& center, float radius, float max_distance,
                       float& out_t, math::Vec3& out_normal) {
    const math::Vec3 oc = origin - center;
    const float b = oc.dot(direction);
    const float c = oc.length_squared() - radius * radius;
    const float discriminant = b * b - c;
    if (discriminant < 0.0f) return false;

    const float sqrt_discriminant = std::sqrt(discriminant);
    float t = -b - sqrt_discriminant;
    if (t < 0.0f) t = -b + sqrt_discriminant;
    if (t < 0.0f || t > max_distance) return false;

    out_t = t;
    out_normal = (origin + direction * t - center).normalized();
    return true;
}

inline bool ray_box(const math::Vec3& origin, const math::Vec3& direction,
                    const math::Vec3& center, const math::Quat& rotation,
                    const math::Vec3& half_extents, float max_distance,
                    float& out_t, math::Vec3& out_normal) {
    const math::Quat inv_rotation = rotation.conjugate();
    const math::Vec3 local_origin = inv_rotation.rotate(origin - center);
    const math::Vec3 local_direction = inv_rotation.rotate(direction);

    float t_min = 0.0f;
    float t_max = max_distance;
    int normal_axis = 0;
    float normal_sign = 1.0f;

    for (int axis = 0; axis < 3; ++axis) {
        const float d = local_direction[axis];
        const float o = local_origin[axis];
        const float h = half_extents[axis];

        if (std::abs(d) < 1.0e-8f) {
            if (o < -h || o > h) return false;
            continue;
        }

        float t1 = (-h - o) / d;
        float t2 = (h - o) / d;
        if (t1 > t2) std::swap(t1, t2);

        if (t1 > t_min) {
            t_min = t1;
            normal_axis = axis;
            normal_sign = (d < 0.0f) ? 1.0f : -1.0f;
        }
        if (t2 < t_max) t_max = t2;
        if (t_min > t_max) return false;
    }

    out_t = t_min;
    math::Vec3 local_normal{};
    local_normal[normal_axis] = normal_sign;
    out_normal = rotation.rotate(local_normal);
    return true;
}

inline bool ray_capsule(const math::Vec3& origin, const math::Vec3& direction,
                        const math::Vec3& center, const math::Quat& rotation,
                        float half_height, float radius, float max_distance,
                        float& out_t, math::Vec3& out_normal) {
    const math::Vec3 a = center + rotation.rotate({0, -half_height, 0});
    const math::Vec3 b = center + rotation.rotate({0, half_height, 0});

    float best_t = std::numeric_limits<float>::infinity();
    math::Vec3 best_normal;
    bool hit = false;

    const math::Vec3 endpoints[2] = {a, b};
    for (const math::Vec3& endpoint : endpoints) {
        float t;
        math::Vec3 normal;
        if (ray_sphere(origin, direction, endpoint, radius, max_distance, t, normal) &&
            t < best_t) {
            best_t = t;
            best_normal = normal;
            hit = true;
        }
    }

    // Mid-segment sphere catches the common case; true capsule/cylinder tests
    // are deliberately kept simple in this reference implementation.
    const math::Vec3 mid = (a + b) * 0.5f;
    float t;
    math::Vec3 normal;
    if (ray_sphere(origin, direction, mid, radius + half_height, max_distance, t, normal) &&
        t < best_t) {
        best_t = t;
        best_normal = normal;
        hit = true;
    }

    if (!hit) return false;
    out_t = best_t;
    out_normal = best_normal;
    return true;
}

} // namespace detail

class World {
public:
    struct Config {
        math::Vec3 gravity{0, -9.81f, 0};
        float fixed_timestep{1.0f / 60.0f};
        int solver_iterations{8};
        int velocity_iterations{4};
        int position_iterations{2};
        bool enable_sleeping{true};
        float sleep_threshold{0.05f};
    };

    Config config;
    math::Vec3 gravity;

    std::function<void(const CollisionEvent&)> on_collision;
    std::function<void(const TriggerEvent&)> on_trigger;
    std::function<void(Body&)> on_body_added;
    std::function<void(Body&)> on_body_removed;
    std::function<void(float)> on_pre_step;
    std::function<void(float)> on_post_step;

    explicit World(const Config& configuration = {})
        : config(configuration), gravity(configuration.gravity) {}

    BodyBuilder create_body() { return BodyBuilder(*this); }

    Body& add_body(BodyType type = BodyType::Dynamic) {
        auto body = std::make_unique<Body>(type);
        return add_body(std::move(body));
    }

    std::vector<Body*> create_bodies(int count) {
        std::vector<Body*> result;
        result.reserve(static_cast<std::size_t>(count));
        for (int i = 0; i < count; ++i) {
            result.push_back(&add_body());
        }
        return result;
    }

    void destroy(Body& body) {
        body.mark_for_destruction();
    }

    void step(float dt) {
        accumulator_ += dt;
        while (accumulator_ >= config.fixed_timestep) {
            step_fixed(config.fixed_timestep);
            accumulator_ -= config.fixed_timestep;
        }
    }

    void step() { step(config.fixed_timestep); }

    auto bodies() {
        return bodies_
            | std::views::filter([](const auto& body) { return !body->is_destroyed(); })
            | std::views::transform([](auto& body) -> Body& { return *body; });
    }

    auto dynamic_bodies() {
        return bodies() | std::views::filter([](Body& body) { return body.is_dynamic(); });
    }

    QueryBuilder raycast() { return QueryBuilder(*this); }

    std::optional<RaycastHit> raycast(const math::Vec3& origin,
                                      const math::Vec3& direction,
                                      float max_distance = std::numeric_limits<float>::max()) {
        return raycast_first(origin, direction, max_distance, {}, false);
    }

    std::optional<RaycastHit> raycast_first(const math::Vec3& origin,
                                            const math::Vec3& direction,
                                            float max_distance,
                                            const std::vector<Body*>& ignore_list,
                                            bool dynamic_only) {
        const math::Vec3 dir = direction.normalized();
        RaycastHit best;
        best.distance = max_distance;

        for (auto& body_ptr : bodies_) {
            Body* body = body_ptr.get();
            if (body->is_destroyed()) continue;
            if (dynamic_only && !body->is_dynamic()) continue;
            if (std::find(ignore_list.begin(), ignore_list.end(), body) != ignore_list.end()) continue;

            const math::Transform transform{body->position, body->rotation};
            for (const auto& collider : body->colliders()) {
                float t = 0;
                math::Vec3 normal;
                bool hit = ray_cast_shape(*collider, transform, origin, dir, max_distance, t, normal);
                if (hit && t < best.distance) {
                    best.body = body;
                    best.distance = t;
                    best.point = origin + dir * t;
                    best.normal = normal;
                    best.hit = true;
                }
            }
        }

        if (!best.hit) return std::nullopt;
        return best;
    }

    std::vector<RaycastHit> raycast_all(const math::Vec3& origin,
                                        const math::Vec3& direction,
                                        float max_distance,
                                        const std::vector<Body*>& ignore_list,
                                        bool dynamic_only) {
        std::vector<RaycastHit> hits;
        const math::Vec3 dir = direction.normalized();

        for (auto& body_ptr : bodies_) {
            Body* body = body_ptr.get();
            if (body->is_destroyed()) continue;
            if (dynamic_only && !body->is_dynamic()) continue;
            if (std::find(ignore_list.begin(), ignore_list.end(), body) != ignore_list.end()) continue;

            const math::Transform transform{body->position, body->rotation};
            for (const auto& collider : body->colliders()) {
                float t = 0;
                math::Vec3 normal;
                if (ray_cast_shape(*collider, transform, origin, dir, max_distance, t, normal)) {
                    RaycastHit hit;
                    hit.body = body;
                    hit.distance = t;
                    hit.point = origin + dir * t;
                    hit.normal = normal;
                    hit.hit = true;
                    hits.push_back(hit);
                }
            }
        }

        std::sort(hits.begin(), hits.end(), [](const RaycastHit& a, const RaycastHit& b) {
            return a.distance < b.distance;
        });
        return hits;
    }

    std::vector<Body*> overlap_test(const Collider& collider,
                                    const math::Transform& transform) {
        std::vector<Body*> result;
        for (auto& body_ptr : bodies_) {
            Body* body = body_ptr.get();
            if (body->is_destroyed()) continue;
            const math::Transform body_transform{body->position, body->rotation};
            for (const auto& body_collider : body->colliders()) {
                ContactInfo contact;
                if (Collider::test(collider, transform, *body_collider, body_transform, contact)) {
                    result.push_back(body);
                    break;
                }
            }
        }
        return result;
    }

    std::vector<Body*> query_aabb(const math::AABB& aabb) {
        std::vector<Body*> result;
        for (auto& body_ptr : bodies_) {
            Body* body = body_ptr.get();
            if (body->is_destroyed()) continue;
            const math::Transform transform{body->position, body->rotation};
            for (const auto& collider : body->colliders()) {
                if (collider->compute_aabb(transform).overlaps(aabb)) {
                    result.push_back(body);
                    break;
                }
            }
        }
        return result;
    }

    std::vector<Body*> query_sphere(const math::Vec3& center, float radius) {
        const math::AABB aabb{center - math::Vec3{radius, radius, radius},
                              center + math::Vec3{radius, radius, radius}};
        return query_aabb(aabb);
    }

    int body_count() const { return static_cast<int>(bodies_.size()); }
    int dynamic_body_count() const {
        int count = 0;
        for (const auto& body : bodies_) {
            if (!body->is_destroyed() && body->is_dynamic()) ++count;
        }
        return count;
    }
    float simulation_time() const { return simulation_time_; }

    template <typename T, typename... Args>
    T& add_joint(Args&&... args) {
        auto joint = std::make_unique<T>(std::forward<Args>(args)...);
        T* ptr = joint.get();
        joints_.push_back(std::move(joint));
        return *ptr;
    }

    std::vector<std::unique_ptr<Joint>>& joints() { return joints_; }
    const std::vector<std::unique_ptr<Joint>>& joints() const { return joints_; }

private:
    struct ContactPair {
        Body* a{nullptr};
        Body* b{nullptr};
        ContactInfo info;
        float restitution{0};
        float friction{0.5f};
        bool is_trigger{false};
    };

    std::vector<std::unique_ptr<Body>> bodies_;
    std::vector<std::unique_ptr<Joint>> joints_;
    std::vector<ContactPair> contacts_;
    float accumulator_{0};
    float simulation_time_{0};

    Body& add_body(std::unique_ptr<Body> body) {
        Body* ptr = body.get();
        bodies_.push_back(std::move(body));
        if (on_body_added) on_body_added(*ptr);
        return *ptr;
    }

    static bool ray_cast_shape(const Collider& collider,
                               const math::Transform& transform,
                               const math::Vec3& origin,
                               const math::Vec3& direction,
                               float max_distance,
                               float& out_t,
                               math::Vec3& out_normal) {
        const math::Vec3 center = transform.position + collider.offset.position;
        const math::Quat rotation = (transform.rotation * collider.offset.rotation).normalized();

        switch (collider.type()) {
            case Collider::Type::Sphere: {
                const auto& sphere = static_cast<const SphereCollider&>(collider);
                return detail::ray_sphere(origin, direction, center, sphere.radius(),
                                          max_distance, out_t, out_normal);
            }
            case Collider::Type::Box: {
                const auto& box = static_cast<const BoxCollider&>(collider);
                return detail::ray_box(origin, direction, center, rotation,
                                       box.half_extents(), max_distance, out_t, out_normal);
            }
            case Collider::Type::Capsule: {
                const auto& capsule = static_cast<const CapsuleCollider&>(collider);
                return detail::ray_capsule(origin, direction, center, rotation,
                                           capsule.half_height(), capsule.radius(),
                                           max_distance, out_t, out_normal);
            }
            default: {
                // Convex/mesh fallback: test against the world-space AABB.
                return detail::ray_box(origin, direction, center, math::Quat::identity(),
                                       collider.compute_aabb(transform).half_extents(),
                                       max_distance, out_t, out_normal);
            }
        }
    }

    void step_fixed(float dt) {
        if (on_pre_step) on_pre_step(dt);

        update_bodies(dt);
        detect_collisions();
        solve_joints(dt);
        solve_collisions();

        if (on_post_step) on_post_step(dt);
        simulation_time_ += dt;
        cleanup();
    }

    void update_bodies(float dt) {
        for (auto& body_ptr : bodies_) {
            Body& body = *body_ptr;
            if (body.is_destroyed() || !body.is_dynamic()) continue;

            if (body.is_sleeping()) {
                if (body.force_accumulator_.length_squared() > 0.0f ||
                    body.torque_accumulator_.length_squared() > 0.0f) {
                    body.wake_up();
                } else {
                    continue;
                }
            }

            if (body.gravity_enabled()) {
                body.apply_force(gravity * body.mass());
            }

            body.velocity += body.force_accumulator_ * (body.inverse_mass() * dt);
            body.angular_velocity += body.torque_accumulator_ * dt;

            if (body.linear_damping > 0.0f) {
                body.velocity *= 1.0f / (1.0f + body.linear_damping * dt);
            }
            if (body.angular_damping > 0.0f) {
                body.angular_velocity *= 1.0f / (1.0f + body.angular_damping * dt);
            }

            body.position += body.velocity * dt;

            const math::Quat spin{
                1.0f,
                body.angular_velocity.x * dt * 0.5f,
                body.angular_velocity.y * dt * 0.5f,
                body.angular_velocity.z * dt * 0.5f,
            };
            body.rotation = (spin * body.rotation).normalized();

            body.clear_forces();

            if (config.enable_sleeping &&
                body.velocity.length_squared() < config.sleep_threshold * config.sleep_threshold) {
                body.set_sleeping(true);
            }
        }
    }

    void detect_collisions() {
        contacts_.clear();
        const int count = static_cast<int>(bodies_.size());
        for (int i = 0; i < count; ++i) {
            Body* a = bodies_[i].get();
            if (a->is_destroyed()) continue;
            const math::Transform transform_a{a->position, a->rotation};

            for (int j = i + 1; j < count; ++j) {
                Body* b = bodies_[j].get();
                if (b->is_destroyed()) continue;
                if (!a->is_dynamic() && !b->is_dynamic()) continue;
                if ((a->collision_mask & b->collision_group) == 0 ||
                    (b->collision_mask & a->collision_group) == 0) continue;

                const math::Transform transform_b{b->position, b->rotation};
                for (const auto& collider_a : a->colliders()) {
                    for (const auto& collider_b : b->colliders()) {
                        ContactInfo contact;
                        if (!Collider::test(*collider_a, transform_a, *collider_b, transform_b, contact)) {
                            continue;
                        }

                        ContactPair pair;
                        pair.a = a;
                        pair.b = b;
                        pair.info = contact;
                        pair.restitution = std::max(collider_a->material.restitution,
                                                    collider_b->material.restitution);
                        pair.friction = std::sqrt(collider_a->material.friction *
                                                  collider_b->material.friction);
                        pair.is_trigger = collider_a->material.is_trigger ||
                                          collider_b->material.is_trigger;
                        contacts_.push_back(pair);
                        break;
                    }
                    if (!contacts_.empty() && contacts_.back().a == a && contacts_.back().b == b) break;
                }
            }
        }
    }

    void solve_collisions() {
        for (int iteration = 0; iteration < config.velocity_iterations; ++iteration) {
            for (auto& pair : contacts_) {
                if (pair.is_trigger) continue;
                solve_contact_velocity(pair);
            }
        }

        for (int iteration = 0; iteration < config.position_iterations; ++iteration) {
            for (auto& pair : contacts_) {
                if (pair.is_trigger) continue;
                solve_contact_position(pair);
            }
        }

        for (auto& pair : contacts_) {
            if (pair.is_trigger) {
                if (on_trigger) {
                    TriggerEvent event{*pair.a, *pair.b, true};
                    on_trigger(event);
                }
                continue;
            }

            CollisionEvent event{*pair.a, *pair.b, pair.info.point, pair.info.normal,
                                 pair.info.impulse, pair.info.penetration};
            if (on_collision) on_collision(event);
            if (pair.a->on_collision) pair.a->on_collision(event);
            if (pair.b->on_collision) pair.b->on_collision(event);
        }
    }

    void solve_contact_velocity(ContactPair& pair) {
        Body& a = *pair.a;
        Body& b = *pair.b;
        const float inverse_mass = a.inverse_mass() + b.inverse_mass();
        if (inverse_mass <= 0.0f) return;

        const math::Vec3 normal = pair.info.normal;
        const math::Vec3 contact = pair.info.point;
        const math::Vec3 ra = contact - a.position;
        const math::Vec3 rb = contact - b.position;

        // Relative velocity at the contact point includes angular velocity, so
        // a spinning crate that touches the ground can actually be slowed by
        // friction instead of sliding forever.
        const math::Vec3 relative = (b.velocity + b.angular_velocity.cross(rb)) -
                                    (a.velocity + a.angular_velocity.cross(ra));
        const float normal_speed = relative.dot(normal);

        if (normal_speed > 0.0f && pair.info.penetration <= 0.0f) return;

        const math::Vec3 ra_cross_n = ra.cross(normal);
        const math::Vec3 rb_cross_n = rb.cross(normal);
        const float effective_mass = inverse_mass +
            ra_cross_n.dot(a.inverse_inertia_tensor() * ra_cross_n) +
            rb_cross_n.dot(b.inverse_inertia_tensor() * rb_cross_n);
        if (effective_mass <= 0.0f) return;

        const float restitution = pair.restitution;
        float impulse_magnitude = -normal_speed / effective_mass;
        if (normal_speed < 0.0f) {
            impulse_magnitude *= (1.0f + restitution);
        }

        math::Vec3 tangent = relative - normal * normal_speed;
        const float tangent_length = tangent.length();
        if (tangent_length > 1.0e-6f) {
            tangent = tangent / tangent_length;
            float friction_impulse = -relative.dot(tangent) / effective_mass;
            const float max_friction = pair.friction * std::abs(impulse_magnitude);
            friction_impulse = std::clamp(friction_impulse, -max_friction, max_friction);
            const math::Vec3 impulse = normal * impulse_magnitude + tangent * friction_impulse;
            a.apply_impulse_at_point(-impulse, contact);
            b.apply_impulse_at_point(impulse, contact);
            pair.info.impulse = impulse_magnitude;
        } else {
            const math::Vec3 impulse = normal * impulse_magnitude;
            a.apply_impulse_at_point(-impulse, contact);
            b.apply_impulse_at_point(impulse, contact);
            pair.info.impulse = impulse_magnitude;
        }
    }

    void solve_contact_position(ContactPair& pair) {
        Body& a = *pair.a;
        Body& b = *pair.b;
        const float inverse_mass = a.inverse_mass() + b.inverse_mass();
        if (inverse_mass <= 0.0f) return;

        constexpr float slop = 0.01f;
        constexpr float percent = 0.4f;
        const float correction = std::max(pair.info.penetration - slop, 0.0f) /
                                 inverse_mass * percent;
        const math::Vec3 delta = pair.info.normal * correction;
        a.position -= delta * a.inverse_mass();
        b.position += delta * b.inverse_mass();
    }

    void solve_joints(float dt) {
        for (auto& joint : joints_) {
            if (joint && joint->enabled) {
                joint->solve(dt);
            }
        }
    }

    void cleanup() {
        for (auto it = bodies_.begin(); it != bodies_.end();) {
            if ((*it)->is_destroyed()) {
                Body* body = it->get();
                if (on_body_removed) on_body_removed(*body);
                it = bodies_.erase(it);
            } else {
                ++it;
            }
        }
    }

    friend class BodyBuilder;
    friend class QueryBuilder;
};

inline Body& BodyBuilder::build() {
    auto body = std::make_unique<Body>();
    body->type_ = type_;
    body->position = position_;
    body->rotation = rotation_;
    body->velocity = velocity_;
    body->gravity_enabled_ = true;

    for (auto& collider : colliders_) {
        collider->material = material_;
        body->add_collider(std::move(collider));
    }
    colliders_.clear();

    body->set_mass(mass_);
    body->on_collision = collision_callback_;
    body->on_sleep = sleep_callback_;

    return world_.add_body(std::move(body));
}

inline BodyBuilder::operator Body&() {
    return build();
}

inline std::optional<RaycastHit> QueryBuilder::first() {
    return world_.raycast_first(origin_, direction_, max_distance_, ignore_list_, dynamic_only_);
}

inline std::vector<RaycastHit> QueryBuilder::all() {
    return world_.raycast_all(origin_, direction_, max_distance_, ignore_list_, dynamic_only_);
}

} // namespace butter
