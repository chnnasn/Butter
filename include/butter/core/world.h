#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
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
    // Transform the ray into the capsule's local frame.  A capsule is the
    // Minkowski sum of a finite Y-axis segment and a sphere: solve the exact
    // finite cylinder interval, then test both spherical caps.  The previous
    // midpoint-sphere approximation could report hits in the empty corners
    // around a long capsule and miss its cylindrical side.
    const math::Quat inverse_rotation = rotation.conjugate();
    const math::Vec3 local_origin = inverse_rotation.rotate(origin - center);
    const float direction_length = direction.length();
    if (direction_length <= 1.0e-12f) return false;
    const math::Vec3 local_direction =
        inverse_rotation.rotate(direction / direction_length);
    const float h = std::max(0.0f, half_height);
    const float r = std::max(0.0f, radius);
    const float max_t = max_distance;

    float best_t = std::numeric_limits<float>::infinity();
    math::Vec3 best_normal_local{};
    bool hit = false;
    const auto record = [&](float t, const math::Vec3& normal_local) {
        if (!std::isfinite(t) || t < 0.0f || t > max_t || t >= best_t) return;
        best_t = t;
        best_normal_local = normal_local;
        hit = true;
    };

    // Infinite cylinder around the local Y axis, clipped to the finite
    // segment between the two cap centres.
    const float quadratic_a = local_direction.x * local_direction.x +
                              local_direction.z * local_direction.z;
    const float quadratic_b = 2.0f * (local_origin.x * local_direction.x +
                                      local_origin.z * local_direction.z);
    const float quadratic_c = local_origin.x * local_origin.x +
                              local_origin.z * local_origin.z - r * r;
    if (quadratic_a > 1.0e-12f) {
        const float discriminant = quadratic_b * quadratic_b -
                                   4.0f * quadratic_a * quadratic_c;
        if (discriminant >= 0.0f) {
            const float root = std::sqrt(std::max(0.0f, discriminant));
            float t0 = (-quadratic_b - root) / (2.0f * quadratic_a);
            float t1 = (-quadratic_b + root) / (2.0f * quadratic_a);
            if (t0 > t1) std::swap(t0, t1);
            const float roots[2] = {t0, t1};
            for (const float t : roots) {
                if (t < 0.0f || t > max_t) continue;
                const float y = local_origin.y + local_direction.y * t;
                if (y < -h - 1.0e-6f || y > h + 1.0e-6f) continue;
                const math::Vec3 point = local_origin + local_direction * t;
                record(t, math::Vec3{point.x, 0.0f, point.z}.normalized());
            }
        }
    }

    // Spherical caps use the ray-sphere quadratic, clipped to the outward
    // hemisphere. A full-sphere test without this clip would report the
    // hidden half of a cap inside the cylindrical section as an early hit.
    const math::Vec3 cap_centres[2] = {{0, -h, 0}, {0, h, 0}};
    for (const math::Vec3& cap : cap_centres) {
        const math::Vec3 offset = local_origin - cap;
        const float b = offset.dot(local_direction);
        const float c = offset.length_squared() - r * r;
        const float discriminant = b * b - c;
        if (discriminant < 0.0f) continue;
        const float root = std::sqrt(std::max(0.0f, discriminant));
        const float roots[2] = {-b - root, -b + root};
        const float cap_sign = cap.y >= 0.0f ? 1.0f : -1.0f;
        for (const float t : roots) {
            if (t < 0.0f || t > max_t) continue;
            const math::Vec3 point = local_origin + local_direction * t;
            if (h > 1.0e-6f && cap_sign * (point.y - cap.y) < -1.0e-6f)
                continue;
            const math::Vec3 normal = (point - cap).normalized();
            record(t, normal);
        }
    }

    if (!hit) return false;
    out_t = best_t;
    out_normal = rotation.rotate(best_normal_local).normalized();
    return true;
}

} // namespace detail

class World {
public:
    enum class SolverMode { Impulse, PBD };

    struct Config {
        math::Vec3 gravity{0, -9.81f, 0};
        float fixed_timestep{1.0f / 60.0f};
        int solver_iterations{8};
        int velocity_iterations{4};
        int position_iterations{2};
        bool enable_sleeping{true};
        float sleep_threshold{0.05f};
        SolverMode solver_mode{SolverMode::Impulse};
        float pbd_compliance{0.0f};

        // Dynamic spatial-hash broad-phase settings.  The hash is rebuilt at
        // each fixed step so public position/rotation fields can be edited by
        // callers without an explicit proxy update call.  A body spanning
        // more than broadphase_max_cells_per_body cells is kept in a small
        // "large body" list and tested against all proxies.
        bool enable_broadphase{true};
        float broadphase_cell_size{2.0f};
        float broadphase_fat_margin{0.05f};
        int broadphase_max_cells_per_body{1024};
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
        std::vector<int> candidate_indices;
        if (config.enable_broadphase) {
            candidate_indices = broadphase_query_indices(aabb);
        } else {
            candidate_indices.reserve(bodies_.size());
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                candidate_indices.push_back(static_cast<int>(i));
            }
        }

        for (const int index : candidate_indices) {
            if (index < 0 || index >= static_cast<int>(bodies_.size())) continue;
            Body* body = bodies_[static_cast<std::size_t>(index)].get();
            if (!body || body->is_destroyed()) continue;
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

    // Number of body pairs emitted by the most recent broad-phase pass.  It
    // is useful for diagnostics and lets applications verify that the
    // broad-phase is actually culling narrow-phase work.
    std::size_t broadphase_candidate_count() const {
        return broadphase_candidate_count_;
    }

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
    struct BroadphaseCell {
        int x{0};
        int y{0};
        int z{0};

        bool operator==(const BroadphaseCell& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct BroadphaseCellHash {
        std::size_t operator()(const BroadphaseCell& cell) const noexcept {
            // SplitMix-style integer mixing gives a good distribution for
            // both positive and negative grid coordinates.
            const auto mix = [](std::uint64_t value) {
                value += 0x9e3779b97f4a7c15ull;
                value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
                value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
                return value ^ (value >> 31);
            };
            const std::uint64_t x = static_cast<std::uint32_t>(cell.x);
            const std::uint64_t y = static_cast<std::uint32_t>(cell.y);
            const std::uint64_t z = static_cast<std::uint32_t>(cell.z);
            return static_cast<std::size_t>(mix(x) ^ (mix(y) << 1) ^ (mix(z) >> 1));
        }
    };

    struct ContactPair {
        Body* a{nullptr};
        Body* b{nullptr};
        // A contact's collider order is not necessarily the logical trigger
        // direction (the second collider can be the trigger), so retain the
        // direction explicitly for event tracking.
        Body* trigger_body{nullptr};
        Body* other_body{nullptr};
        ContactInfo info;
        float restitution{0};
        float friction{0.5f};
        bool is_trigger{false};
    };

    struct TriggerPairKey {
        Body* trigger{nullptr};
        Body* other{nullptr};

        friend bool operator==(const TriggerPairKey& lhs,
                               const TriggerPairKey& rhs) noexcept {
            return lhs.trigger == rhs.trigger && lhs.other == rhs.other;
        }
    };

    struct TriggerPairHash {
        std::size_t operator()(const TriggerPairKey& key) const noexcept {
            const auto trigger_hash = std::hash<const void*>{}(key.trigger);
            const auto other_hash = std::hash<const void*>{}(key.other);
            return trigger_hash ^ (other_hash + static_cast<std::size_t>(0x9e3779b9u) +
                                   (trigger_hash << 6) + (trigger_hash >> 2));
        }
    };

    std::vector<std::unique_ptr<Body>> bodies_;
    std::vector<std::unique_ptr<Joint>> joints_;
    std::vector<ContactPair> contacts_;
    // The broad-phase is intentionally stored as body indices instead of
    // Body* values.  Indices remain stable for the duration of a step, while
    // cleanup only erases bodies after the step has finished.
    std::unordered_map<BroadphaseCell, std::vector<int>, BroadphaseCellHash>
        broadphase_cells_;
    std::vector<int> broadphase_large_bodies_;
    std::vector<math::AABB> broadphase_bounds_;
    std::size_t broadphase_candidate_count_{0};
    // Contacts are rebuilt every frame (and repeatedly during PBD
    // projection), so trigger events must be compared with the prior fixed
    // step rather than emitted directly from each detect pass.
    std::unordered_set<TriggerPairKey, TriggerPairHash> active_trigger_pairs_;
    // Velocities immediately after force integration and before PBD contact
    // projection. They are used as the physically predicted state for the
    // restitution/friction pass, so positional correction is not mistaken for
    // a real bounce.
    std::vector<math::Vec3> pbd_predicted_velocities_;
    std::vector<math::Vec3> pbd_predicted_angular_velocities_;
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
        const math::Vec3 center = transform.transform_point(collider.offset.position);
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
                // Use the AABB centre (rather than the body origin), since a
                // point cloud or mesh is not required to be origin-centred.
                const math::AABB bounds = collider.compute_aabb(transform);
                if (bounds.is_empty()) return false;
                return detail::ray_box(origin, direction, bounds.center(),
                                       math::Quat::identity(), bounds.half_extents(),
                                       max_distance, out_t, out_normal);
            }
        }
    }

    void step_fixed(float dt) {
        if (on_pre_step) on_pre_step(dt);

        update_bodies(dt);
        solve_joints(dt);
        if (config.solver_mode == SolverMode::PBD) {
            // The PBD solver rebuilds contacts before every projection pass;
            // avoid a discarded broad-phase/narrow-phase pass here.
            solve_collisions_pbd(dt);
        } else {
            detect_collisions();
            solve_collisions();
        }

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

            if (config.enable_sleeping) {
                const float speed_limit = config.sleep_threshold * config.sleep_threshold;
                const bool quiet = body.velocity.length_squared() < speed_limit &&
                                   body.angular_velocity.length_squared() < speed_limit;
                // A body must remain quiet for a settling period. Sleeping on
                // the first low-speed gravity step freezes objects in mid-air.
                if (quiet) {
                    body.sleep_time_ += dt;
                    if (body.sleep_time_ >= 0.5f) body.set_sleeping(true);
                } else {
                    body.sleep_time_ = 0.0f;
                }
            }
        }
    }

    static int broadphase_cell_coordinate(float value, float cell_size) {
        const double coordinate = std::floor(static_cast<double>(value) /
                                              static_cast<double>(cell_size));
        // Avoid implementation-defined float-to-int overflow for malformed
        // user transforms. Such a body is handled as a large proxy instead.
        constexpr double min_int = static_cast<double>(std::numeric_limits<int>::min());
        constexpr double max_int = static_cast<double>(std::numeric_limits<int>::max());
        if (!std::isfinite(coordinate) || coordinate < min_int || coordinate > max_int) {
            return 0;
        }
        return static_cast<int>(coordinate);
    }

    static std::uint64_t broadphase_pair_key(int a, int b) {
        if (a > b) std::swap(a, b);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
               static_cast<std::uint32_t>(b);
    }

    static bool body_has_trigger(const Body& body) {
        for (const auto& collider : body.colliders()) {
            if (collider && collider->material.is_trigger) return true;
        }
        return false;
    }

    bool world_has_trigger_colliders() const {
        for (const auto& body : bodies_) {
            if (body && !body->is_destroyed() && body_has_trigger(*body)) return true;
        }
        return false;
    }

    void rebuild_broadphase() {
        broadphase_cells_.clear();
        broadphase_large_bodies_.clear();

        const int count = static_cast<int>(bodies_.size());
        broadphase_bounds_.assign(static_cast<std::size_t>(count), math::AABB{});

        const float configured_cell_size = config.broadphase_cell_size;
        const float cell_size = std::isfinite(configured_cell_size) &&
                                        configured_cell_size > 1.0e-4f
                                    ? configured_cell_size
                                    : 2.0f;
        const float margin = std::isfinite(config.broadphase_fat_margin)
                                 ? std::max(config.broadphase_fat_margin, 0.0f)
                                 : 0.0f;
        const int configured_limit = config.broadphase_max_cells_per_body;
        const std::int64_t max_cells = configured_limit > 0
                                           ? configured_limit
                                           : 1024;

        for (int i = 0; i < count; ++i) {
            Body* body = bodies_[static_cast<std::size_t>(i)].get();
            if (!body || body->is_destroyed()) continue;

            const math::Transform transform{body->position, body->rotation};
            math::AABB bounds;
            for (const auto& collider : body->colliders()) {
                if (collider) bounds.encapsulate(collider->compute_aabb(transform));
            }
            broadphase_bounds_[static_cast<std::size_t>(i)] = bounds;
            if (bounds.is_empty()) continue;

            // A non-finite public transform must not be allowed to produce a
            // gigantic/undefined grid range. Treat it as a large proxy; the
            // exact narrow phase can then reject it safely (or report a
            // contact if the shape implementation supports it).
            const bool finite_bounds = std::isfinite(bounds.min.x) &&
                                       std::isfinite(bounds.min.y) &&
                                       std::isfinite(bounds.min.z) &&
                                       std::isfinite(bounds.max.x) &&
                                       std::isfinite(bounds.max.y) &&
                                       std::isfinite(bounds.max.z);
            if (!finite_bounds) {
                broadphase_large_bodies_.push_back(i);
                continue;
            }

            const math::AABB fat_bounds = bounds.expanded(margin);
            const int min_x = broadphase_cell_coordinate(fat_bounds.min.x, cell_size);
            const int min_y = broadphase_cell_coordinate(fat_bounds.min.y, cell_size);
            const int min_z = broadphase_cell_coordinate(fat_bounds.min.z, cell_size);
            const int max_x = broadphase_cell_coordinate(fat_bounds.max.x, cell_size);
            const int max_y = broadphase_cell_coordinate(fat_bounds.max.y, cell_size);
            const int max_z = broadphase_cell_coordinate(fat_bounds.max.z, cell_size);

            const std::int64_t span_x = static_cast<std::int64_t>(max_x) - min_x + 1;
            const std::int64_t span_y = static_cast<std::int64_t>(max_y) - min_y + 1;
            const std::int64_t span_z = static_cast<std::int64_t>(max_z) - min_z + 1;
            const bool invalid_range = span_x <= 0 || span_y <= 0 || span_z <= 0;
            bool too_many_cells = invalid_range;
            if (!too_many_cells) {
                // Compare against the limit before multiplying. This avoids
                // signed 64-bit overflow for extreme user coordinates.
                too_many_cells = span_x > max_cells || span_y > max_cells ||
                                 span_z > max_cells ||
                                 span_x > max_cells / std::max<std::int64_t>(span_y, 1) ||
                                 span_x * span_y > max_cells / std::max<std::int64_t>(span_z, 1);
            }
            if (too_many_cells) {
                // Huge ground planes and malformed/infinite proxies should
                // not explode the hash table. They are few and are tested
                // against the ordinary candidate list below.
                broadphase_large_bodies_.push_back(i);
                continue;
            }

            for (std::int64_t x = min_x; x <= static_cast<std::int64_t>(max_x); ++x) {
                for (std::int64_t y = min_y; y <= static_cast<std::int64_t>(max_y); ++y) {
                    for (std::int64_t z = min_z; z <= static_cast<std::int64_t>(max_z); ++z) {
                        broadphase_cells_[BroadphaseCell{static_cast<int>(x),
                                                         static_cast<int>(y),
                                                         static_cast<int>(z)}].push_back(i);
                    }
                }
            }
        }
    }

    std::vector<std::pair<int, int>> broadphase_candidates() {
        std::vector<std::pair<int, int>> candidates;
        const int count = static_cast<int>(bodies_.size());

        // A cell can contain the same pair in many neighbouring cells. Keep a
        // compact integer key and sort/unique once instead of allocating a
        // hash set for every pair insertion.
        std::vector<std::uint64_t> keys;
        for (const auto& entry : broadphase_cells_) {
            const auto& occupants = entry.second;
            for (std::size_t i = 0; i < occupants.size(); ++i) {
                for (std::size_t j = i + 1; j < occupants.size(); ++j) {
                    const int a = occupants[i];
                    const int b = occupants[j];
                    if (a != b) keys.push_back(broadphase_pair_key(a, b));
                }
            }
        }

        // Large proxies are intentionally paired with every body. The final
        // exact-AABB test still rejects distant pairs and the sorted key pass
        // removes duplicates with regular grid cells.
        for (const int large : broadphase_large_bodies_) {
            for (int other = 0; other < count; ++other) {
                if (large != other) keys.push_back(broadphase_pair_key(large, other));
            }
        }

        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
        candidates.reserve(keys.size());
        for (const std::uint64_t key : keys) {
            const int a = static_cast<int>(static_cast<std::uint32_t>(key >> 32));
            const int b = static_cast<int>(static_cast<std::uint32_t>(key & 0xFFFFFFFFu));
            candidates.emplace_back(a, b);
        }
        broadphase_candidate_count_ = candidates.size();
        return candidates;
    }

    std::vector<int> broadphase_query_indices(const math::AABB& query) {
        std::vector<int> indices;
        if (query.is_empty()) return indices;
        rebuild_broadphase();

        const float configured_cell_size = config.broadphase_cell_size;
        const float cell_size = std::isfinite(configured_cell_size) &&
                                        configured_cell_size > 1.0e-4f
                                    ? configured_cell_size
                                    : 2.0f;
        const int configured_limit = config.broadphase_max_cells_per_body;
        const std::int64_t max_cells = configured_limit > 0
                                           ? configured_limit
                                           : 1024;
        const bool finite_query = std::isfinite(query.min.x) &&
                                  std::isfinite(query.min.y) &&
                                  std::isfinite(query.min.z) &&
                                  std::isfinite(query.max.x) &&
                                  std::isfinite(query.max.y) &&
                                  std::isfinite(query.max.z);

        std::vector<std::uint64_t> unique_keys;
        if (finite_query) {
            const int min_x = broadphase_cell_coordinate(query.min.x, cell_size);
            const int min_y = broadphase_cell_coordinate(query.min.y, cell_size);
            const int min_z = broadphase_cell_coordinate(query.min.z, cell_size);
            const int max_x = broadphase_cell_coordinate(query.max.x, cell_size);
            const int max_y = broadphase_cell_coordinate(query.max.y, cell_size);
            const int max_z = broadphase_cell_coordinate(query.max.z, cell_size);
            const std::int64_t span_x = static_cast<std::int64_t>(max_x) - min_x + 1;
            const std::int64_t span_y = static_cast<std::int64_t>(max_y) - min_y + 1;
            const std::int64_t span_z = static_cast<std::int64_t>(max_z) - min_z + 1;
            bool too_many_cells = span_x <= 0 || span_y <= 0 || span_z <= 0;
            if (!too_many_cells) {
                too_many_cells = span_x > max_cells || span_y > max_cells ||
                                 span_z > max_cells ||
                                 span_x > max_cells / std::max<std::int64_t>(span_y, 1) ||
                                 span_x * span_y > max_cells / std::max<std::int64_t>(span_z, 1);
            }
            if (!too_many_cells) {
                for (std::int64_t x = min_x; x <= static_cast<std::int64_t>(max_x); ++x) {
                    for (std::int64_t y = min_y; y <= static_cast<std::int64_t>(max_y); ++y) {
                        for (std::int64_t z = min_z; z <= static_cast<std::int64_t>(max_z); ++z) {
                            const auto it = broadphase_cells_.find(
                                BroadphaseCell{static_cast<int>(x), static_cast<int>(y),
                                               static_cast<int>(z)});
                            if (it == broadphase_cells_.end()) continue;
                            for (const int body_index : it->second) {
                                unique_keys.push_back(static_cast<std::uint32_t>(body_index));
                            }
                        }
                    }
                }
            } else {
                // A huge query (for example an unbounded editor selection)
                // is cheaper and safer to answer by scanning all proxies.
                unique_keys.reserve(bodies_.size());
                for (std::size_t i = 0; i < bodies_.size(); ++i) {
                    unique_keys.push_back(static_cast<std::uint32_t>(i));
                }
            }
        } else {
            unique_keys.reserve(bodies_.size());
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                unique_keys.push_back(static_cast<std::uint32_t>(i));
            }
        }

        // Large proxies are not in ordinary cells, so they always need to be
        // considered for a query.
        for (const int body_index : broadphase_large_bodies_) {
            unique_keys.push_back(static_cast<std::uint32_t>(body_index));
        }
        std::sort(unique_keys.begin(), unique_keys.end());
        unique_keys.erase(std::unique(unique_keys.begin(), unique_keys.end()), unique_keys.end());
        indices.reserve(unique_keys.size());
        for (const std::uint64_t key : unique_keys) {
            indices.push_back(static_cast<int>(static_cast<std::uint32_t>(key)));
        }
        return indices;
    }

    void detect_collisions() {
        contacts_.clear();
        const int count = static_cast<int>(bodies_.size());
        if (count < 2) {
            broadphase_candidate_count_ = 0;
            return;
        }

        std::vector<std::pair<int, int>> candidates;
        if (config.enable_broadphase) {
            rebuild_broadphase();
            candidates = broadphase_candidates();
        } else {
            broadphase_candidate_count_ =
                static_cast<std::size_t>(count) * static_cast<std::size_t>(count - 1) / 2;
            candidates.reserve(broadphase_candidate_count_);
            for (int i = 0; i < count; ++i) {
                for (int j = i + 1; j < count; ++j) candidates.emplace_back(i, j);
            }
        }

        for (const auto& candidate : candidates) {
            const int i = candidate.first;
            const int j = candidate.second;
            if (i < 0 || j < 0 || i >= count || j >= count || i == j) continue;
            Body* a = bodies_[static_cast<std::size_t>(i)].get();
            Body* b = bodies_[static_cast<std::size_t>(j)].get();
            if (!a || !b || a->is_destroyed() || b->is_destroyed()) continue;
            // Static/kinematic trigger volumes still need overlap events even
            // when the other object is static.  Ordinary static-static solid
            // pairs remain culled as before.
            if (!a->is_dynamic() && !b->is_dynamic() &&
                !body_has_trigger(*a) && !body_has_trigger(*b)) {
                continue;
            }
            if (config.enable_broadphase &&
                !broadphase_bounds_[static_cast<std::size_t>(i)].overlaps(
                    broadphase_bounds_[static_cast<std::size_t>(j)])) {
                continue;
            }
            if ((a->collision_mask & b->collision_group) == 0 ||
                (b->collision_mask & a->collision_group) == 0) continue;

            const math::Transform transform_a{a->position, a->rotation};
            const math::Transform transform_b{b->position, b->rotation};
            const bool body_a_has_trigger = body_has_trigger(*a);
            const bool body_b_has_trigger = body_has_trigger(*b);
            bool solid_contact_added = false;
            bool trigger_a_added = false;
            bool trigger_b_added = false;
            for (const auto& collider_a : a->colliders()) {
                if (!collider_a) continue;
                for (const auto& collider_b : b->colliders()) {
                    if (!collider_b) continue;
                    ContactInfo contact;
                    if (!Collider::test(*collider_a, transform_a, *collider_b,
                                        transform_b, contact)) {
                        continue;
                    }

                    const bool collider_a_is_trigger = collider_a->material.is_trigger;
                    const bool collider_b_is_trigger = collider_b->material.is_trigger;

                    // A body pair may contain both a solid collider and a
                    // trigger collider.  Keep one solid contact for physics
                    // and one directional contact per trigger body instead of
                    // letting whichever collider happens to be iterated first
                    // suppress the other behavior.
                    if (!collider_a_is_trigger && !collider_b_is_trigger) {
                        if (solid_contact_added) continue;
                        ContactPair pair;
                        pair.a = a;
                        pair.b = b;
                        pair.info = contact;
                        pair.restitution = std::max(collider_a->material.restitution,
                                                    collider_b->material.restitution);
                        pair.friction = std::sqrt(std::max(
                            0.0f, collider_a->material.friction * collider_b->material.friction));
                        pair.is_trigger = false;
                        contacts_.push_back(pair);
                        solid_contact_added = true;
                        if (!body_a_has_trigger && !body_b_has_trigger) break;
                        continue;
                    }

                    if (collider_a_is_trigger && !trigger_a_added) {
                        ContactPair pair;
                        pair.a = a;
                        pair.b = b;
                        pair.info = contact;
                        pair.restitution = std::max(collider_a->material.restitution,
                                                    collider_b->material.restitution);
                        pair.friction = std::sqrt(std::max(
                            0.0f, collider_a->material.friction * collider_b->material.friction));
                        pair.is_trigger = true;
                        pair.trigger_body = a;
                        pair.other_body = b;
                        contacts_.push_back(pair);
                        trigger_a_added = true;
                    }

                    if (collider_b_is_trigger && !trigger_b_added) {
                        ContactPair pair;
                        pair.a = a;
                        pair.b = b;
                        pair.info = contact;
                        pair.restitution = std::max(collider_a->material.restitution,
                                                    collider_b->material.restitution);
                        pair.friction = std::sqrt(std::max(
                            0.0f, collider_a->material.friction * collider_b->material.friction));
                        pair.is_trigger = true;
                        pair.trigger_body = b;
                        pair.other_body = a;
                        contacts_.push_back(pair);
                        trigger_b_added = true;
                    }

                    if (solid_contact_added &&
                        (!body_a_has_trigger || trigger_a_added) &&
                        (!body_b_has_trigger || trigger_b_added)) {
                        break;
                    }
                }
                if (solid_contact_added &&
                    (!body_a_has_trigger || trigger_a_added) &&
                    (!body_b_has_trigger || trigger_b_added)) {
                    break;
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

        // Position correction can move a dynamic body into or out of a
        // trigger volume that was detected before the solver ran.  Refresh
        // the trigger contacts for the final geometry, while retaining the
        // solved contact list so collision callbacks still receive the
        // velocity impulse accumulated above.
        auto solved_contacts = std::move(contacts_);
        // Avoid a second full broad/narrow-phase pass for worlds that have no
        // trigger colliders at all.  If a prior trigger pair is still active,
        // keep the refresh so its exit can be observed after destruction or a
        // material change.
        if (world_has_trigger_colliders() || !active_trigger_pairs_.empty()) {
            detect_collisions();
        } else {
            contacts_.clear();
        }

        for (auto& pair : solved_contacts) {
            if (pair.is_trigger) continue;

            CollisionEvent event{*pair.a, *pair.b, pair.info.point, pair.info.normal,
                                 pair.info.impulse, pair.info.penetration};
            if (on_collision) on_collision(event);
            if (pair.a->on_collision) pair.a->on_collision(event);
            if (pair.b->on_collision) pair.b->on_collision(event);
        }

        dispatch_trigger_events();
    }

    // Position Based Dynamics contact projection.  Forces are integrated to a
    // predicted state first, then contacts are treated as positional
    // constraints. The velocity pass starts from the force-integrated
    // prediction (not the raw projection displacement, which can inject an
    // artificial bounce) and applies restitution plus Coulomb friction for
    // approaching support contacts.
    void solve_collisions_pbd(float dt) {
        if (dt <= 0.0f) return;

        pbd_predicted_velocities_.clear();
        pbd_predicted_angular_velocities_.clear();
        pbd_predicted_velocities_.reserve(bodies_.size());
        pbd_predicted_angular_velocities_.reserve(bodies_.size());
        for (const auto& body : bodies_) {
            pbd_predicted_velocities_.push_back(body->velocity);
            pbd_predicted_angular_velocities_.push_back(body->angular_velocity);
        }

        const auto predicted_velocity = [&](const Body& body) -> const math::Vec3& {
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                if (bodies_[i].get() == &body) return pbd_predicted_velocities_[i];
            }
            // The body can only be absent when user code destroys it from an
            // event callback.  Returning the current value keeps this helper
            // safe for that edge case.
            return body.velocity;
        };
        const auto predicted_angular_velocity = [&](const Body& body) -> const math::Vec3& {
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                if (bodies_[i].get() == &body) return pbd_predicted_angular_velocities_[i];
            }
            return body.angular_velocity;
        };

        for (int iteration = 0; iteration < config.position_iterations; ++iteration) {
            // Contact penetration changes after every projection. Rebuild the
            // contact set instead of applying the same stale penetration many
            // times.  This is essential for a stack: moving the bottom crate
            // changes the penetration of every crate above it.
            detect_collisions();
            for (auto& pair : contacts_) {
                if (pair.is_trigger) continue;
                Body& a = *pair.a;
                Body& b = *pair.b;
                if (a.is_destroyed() || b.is_destroyed() ||
                    pair.info.penetration <= 0.0f) continue;

                const math::Vec3 normal = pair.info.normal;
                const math::Vec3 contact = pair.info.point;
                const math::Vec3 ra = contact - a.position;
                const math::Vec3 rb = contact - b.position;
                const math::Vec3 ra_cross_n = ra.cross(normal);
                const math::Vec3 rb_cross_n = rb.cross(normal);
                const float effective_mass = a.inverse_mass() + b.inverse_mass() +
                    ra_cross_n.dot(a.inverse_inertia_tensor() * ra_cross_n) +
                    rb_cross_n.dot(b.inverse_inertia_tensor() * rb_cross_n);
                if (effective_mass <= 0.0f) continue;

                constexpr float slop = 0.005f;
                const float c = std::max(pair.info.penetration - slop, 0.0f);
                if (c <= 0.0f) continue;
                const float alpha = config.pbd_compliance / (dt * dt);
                const float correction = c / (effective_mass + alpha);
                const math::Vec3 delta = pair.info.normal * correction;
                a.position -= delta * a.inverse_mass();
                b.position += delta * b.inverse_mass();

                // Keep contact projection translational in this first stable
                // rigid-body path. Applying an unscaled angular correction
                // here injects energy into a stack; angular motion is still
                // integrated and friction is solved below at the contact.
            }
        }

        // Get the contacts at the final projected state.  The previous loop
        // intentionally rebuilds contacts before each projection, so this
        // extra pass is needed for events and velocity constraints to see the
        // actual final geometry.
        detect_collisions();

        // Start the velocity pass from the force-integrated prediction. Using
        // (projected_position - previous_position) here turns a large stack
        // correction into an artificial upward launch when a pair is no
        // longer overlapping on the final pass. Contact constraints below
        // remove approaching normal velocity and apply restitution/friction.
        for (std::size_t i = 0; i < bodies_.size(); ++i) {
            Body& body = *bodies_[i];
            if (body.is_dynamic() && !body.is_destroyed()) {
                body.velocity = pbd_predicted_velocities_[i];
                body.angular_velocity = pbd_predicted_angular_velocities_[i];
            }
        }

        // Apply velocity-level restitution and friction after projection.  A
        // few iterations are important for a stack because correcting one
        // contact changes the relative velocity at neighbouring contacts.
        const int velocity_iterations = std::max(0, config.velocity_iterations);
        for (int iteration = 0; iteration < velocity_iterations; ++iteration) {
            for (auto& pair : contacts_) {
                if (pair.is_trigger || !pair.a || !pair.b) continue;
                Body& a = *pair.a;
                Body& b = *pair.b;
                if (a.is_destroyed() || b.is_destroyed()) continue;

                const math::Vec3 normal = pair.info.normal.normalized();
                // Face-to-face vertical support is the stable part of the
                // compact rigid-body PBD solver. Lateral face/corner contacts
                // are still projected positionally above, but feeding their
                // approximate SAT normal into a velocity impulse can create
                // a horizontal ratchet in a dense crate stack.
                if (std::abs(normal.y) < 0.7f) continue;
                const math::Vec3 contact = pair.info.point;
                const math::Vec3 ra = contact - a.position;
                const math::Vec3 rb = contact - b.position;
                const math::Vec3 ra_cross_n = ra.cross(normal);
                const math::Vec3 rb_cross_n = rb.cross(normal);
                const float normal_mass = a.inverse_mass() + b.inverse_mass() +
                    ra_cross_n.dot(a.inverse_inertia_tensor() * ra_cross_n) +
                    rb_cross_n.dot(b.inverse_inertia_tensor() * rb_cross_n);
                if (normal_mass <= 0.0f) continue;

                const math::Vec3 relative =
                    (b.velocity + b.angular_velocity.cross(rb)) -
                    (a.velocity + a.angular_velocity.cross(ra));
                const float current_normal_speed = relative.dot(normal);

                const math::Vec3 pre_relative =
                    (predicted_velocity(b) + predicted_angular_velocity(b).cross(rb)) -
                    (predicted_velocity(a) + predicted_angular_velocity(a).cross(ra));
                const float pre_normal_speed = pre_relative.dot(normal);

                // Do not manufacture support from a post-projection contact.
                // Only a genuinely approaching pair may receive a normal
                // impulse; this prevents side contacts in a resting stack
                // from accumulating lateral drift.
                if (pre_normal_speed >= -1.0e-4f) continue;

                const float target_normal_speed =
                    -pair.restitution * pre_normal_speed;

                const float normal_impulse =
                    (target_normal_speed - current_normal_speed) / normal_mass;
                if (std::abs(normal_impulse) > 1.0e-7f) {
                    const math::Vec3 impulse = normal * normal_impulse;
                    a.apply_impulse_at_point(-impulse, contact);
                    b.apply_impulse_at_point(impulse, contact);
                    pair.info.impulse = normal_impulse;
                }

                // Estimate the support impulse from the predicted approach as
                // well as the explicit restitution impulse.  Position
                // projection supplies support without a velocity impulse, but
                // friction still needs that normal load to cap tangential
                // impulse (otherwise a resting crate slides forever).
                const float support_impulse =
                    std::max(0.0f, -pre_normal_speed) / normal_mass;
                const float normal_load = std::max(
                    support_impulse, std::abs(normal_impulse));

                const math::Vec3 after_normal =
                    (b.velocity + b.angular_velocity.cross(rb)) -
                    (a.velocity + a.angular_velocity.cross(ra));
                const float after_normal_speed = after_normal.dot(normal);
                math::Vec3 tangent = after_normal - normal * after_normal_speed;
                const float tangent_length = tangent.length();
                if (tangent_length <= 1.0e-6f || normal_load <= 0.0f) continue;
                tangent /= tangent_length;

                const math::Vec3 ra_cross_t = ra.cross(tangent);
                const math::Vec3 rb_cross_t = rb.cross(tangent);
                const float tangent_mass = a.inverse_mass() + b.inverse_mass() +
                    ra_cross_t.dot(a.inverse_inertia_tensor() * ra_cross_t) +
                    rb_cross_t.dot(b.inverse_inertia_tensor() * rb_cross_t);
                if (tangent_mass <= 0.0f) continue;

                float friction_impulse = -after_normal.dot(tangent) / tangent_mass;
                const float max_friction = pair.friction * normal_load;
                friction_impulse = std::clamp(friction_impulse,
                                              -max_friction, max_friction);
                if (std::abs(friction_impulse) <= 1.0e-7f) continue;
                const math::Vec3 impulse = tangent * friction_impulse;
                a.apply_impulse_at_point(-impulse, contact);
                b.apply_impulse_at_point(impulse, contact);
            }
        }

        // Report contacts consistently with impulse mode.  The normal impulse
        // is not an exact XPBD lambda, but it is useful diagnostic information
        // and reflects the velocity correction applied above.
        for (auto& pair : contacts_) {
            if (pair.is_trigger) continue;
            CollisionEvent event{*pair.a, *pair.b, pair.info.point, pair.info.normal,
                                 pair.info.impulse, pair.info.penetration};
            if (on_collision) on_collision(event);
            if (pair.a->on_collision) pair.a->on_collision(event);
            if (pair.b->on_collision) pair.b->on_collision(event);
        }

        dispatch_trigger_events();
    }

    void dispatch_trigger_events() {
        std::unordered_set<TriggerPairKey, TriggerPairHash> current;
        current.reserve(contacts_.size());
        for (const auto& pair : contacts_) {
            if (!pair.is_trigger || !pair.trigger_body || !pair.other_body) continue;
            if (pair.trigger_body->is_destroyed() || pair.other_body->is_destroyed()) {
                continue;
            }
            current.insert(TriggerPairKey{pair.trigger_body, pair.other_body});
        }

        // Publish the new state before callbacks.  This makes a callback that
        // calls step() or destroy() observe a coherent state and prevents a
        // re-entrant dispatch from turning one overlap into duplicate enters.
        auto previous = std::move(active_trigger_pairs_);
        active_trigger_pairs_ = current;

        // Copy the callable so a callback is free to replace/clear
        // `world.on_trigger` without invalidating the remainder of this
        // dispatch pass.
        const auto callback = on_trigger;
        if (!callback) return;

        // Report exits first, then enters/stays.  Exits are important when a
        // trigger and an object exchange places in the same fixed step.
        for (const auto& key : previous) {
            if (current.find(key) != current.end()) continue;
            if (!key.trigger || !key.other) continue;
            callback(TriggerEvent{*key.trigger, *key.other, false, false, true});
        }

        for (const auto& key : current) {
            if (!key.trigger || !key.other) continue;
            // Only transitions are events.  Keeping persistent overlaps in
            // active_trigger_pairs_ is what prevents a fresh enter every
            // fixed step while preserving the original `is_enter == false`
            // convention for exits.
            if (previous.find(key) != previous.end()) continue;
            callback(TriggerEvent{*key.trigger, *key.other, true, false, false});
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

    void clear_trigger_pairs_for_body(Body& body) {
        std::vector<TriggerPairKey> removed;
        for (auto it = active_trigger_pairs_.begin();
             it != active_trigger_pairs_.end();) {
            if (it->trigger == &body || it->other == &body) {
                removed.push_back(*it);
                it = active_trigger_pairs_.erase(it);
            } else {
                ++it;
            }
        }

        // Bodies remain alive until the caller erases their unique_ptr, so it
        // is safe for the exit callback to inspect either reference.  Remove
        // all keys before invoking user code to make callbacks re-entrancy
        // safe and to avoid duplicate exits if another body is destroyed from
        // an exit handler.
        const auto callback = on_trigger;
        if (callback) {
            for (const auto& key : removed) {
                if (key.trigger && key.other) {
                    callback(TriggerEvent{*key.trigger, *key.other,
                                          false, false, true});
                }
            }
        }
    }

    void cleanup() {
        // Gather pointers first.  Trigger/body-removed callbacks are user
        // code and may add another body, which can reallocate `bodies_`; doing
        // that while erasing through an iterator would invalidate the loop.
        std::vector<Body*> removed;
        removed.reserve(bodies_.size());
        std::unordered_set<Body*> processed;
        std::size_t next = 0;
        for (;;) {
            // Callbacks can mark another body for destruction (or even add a
            // body and destroy it immediately), so extend the queue before
            // each callback and process newly discovered entries as well.
            for (const auto& body : bodies_) {
                if (body && body->is_destroyed() && processed.insert(body.get()).second) {
                    removed.push_back(body.get());
                }
            }
            if (next >= removed.size()) break;

            Body* body = removed[next++];
            clear_trigger_pairs_for_body(*body);
            if (on_body_removed) on_body_removed(*body);
        }

        bodies_.erase(std::remove_if(bodies_.begin(), bodies_.end(),
                                     [](const auto& body) {
                                         return !body || body->is_destroyed();
                                     }),
                       bodies_.end());
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
