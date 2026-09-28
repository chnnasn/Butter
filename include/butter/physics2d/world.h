#pragma once

#include "butter/core/material.h"
#include "butter/physics2d/broadphase.h"
#include "butter/physics2d/ccd.h"
#include "butter/physics2d/query.h"
#include "butter/physics2d/shapes.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace butter::physics2d {

enum class BodyType { Static, Dynamic, Kinematic };
enum class SolverMode { Impulse, PBD };

struct Body;
struct Fixture {
    Body *body{};
    Shape shape{Circle{}};
    Transform local{};
    Material material{};
    float restitution_threshold{1};
    bool trigger{false};
    std::uint32_t collision_group{1}, collision_mask{0xffffffffu};
    std::uint64_t user_data{};
    // Revision counters. Every field the engine has to react to is grouped by
    // what it invalidates: geometry (shape/offset), material, and filtering.
    // The world polls these counters once per step and only does the expensive
    // work for the objects that moved. Direct field writes still work; they
    // bypass the counter and are caught by the exact snapshot comparison, which
    // is the documented cost of the compatibility path.
    std::uint32_t geometry_version{0};
    std::uint32_t material_version{0};
    std::uint32_t filter_version{0};
    // Per-step memo of the content hash. A shape that participates in many
    // contacts is hashed once per step instead of once per contact.
    mutable std::uint64_t shape_hash_cache{0};
    mutable std::uint32_t shape_hash_epoch{0};

    Fixture &set_shape(Shape value);
    Fixture &set_local(Transform value);
    Fixture &set_material(const Material &value);
    Fixture &set_friction(float value);
    Fixture &set_restitution(float value);
    Fixture &set_density(float value);
    Fixture &set_filter(std::uint32_t group, std::uint32_t mask);
    Fixture &set_trigger(bool value = true);
};

struct Body {
    BodyType type{BodyType::Dynamic};
    Transform transform{};
    // Index of this body inside its World. Engine-managed; kept current so the
    // solver can map a body back to a slot in O(1) while building islands.
    std::uint32_t slot{0};
    // Revision counters, grouped by what they invalidate. See Fixture for the
    // contract: the controlled setters below bump them, direct writes do not
    // and are picked up by the per-step snapshot comparison instead.
    std::uint32_t transform_version{0};
    std::uint32_t geometry_version{0};
    std::uint32_t material_version{0};
    std::uint32_t filter_version{0};
    mutable std::uint64_t shape_hash_cache{0};
    mutable std::uint32_t shape_hash_epoch{0};
    Vec2 velocity{};
    float angular_velocity{0};
    float mass{1};
    float inverse_mass{1};
    float inertia{1};
    float inverse_inertia{1};
    // Centroid offset in body space. Constraints use the centroid as their
    // lever-arm origin and integration rotates about it, exactly like Box2D's
    // sweep localCenter. Zero for symmetric single-shape bodies.
    Vec2 local_center{};
    float linear_damping{0.02f};
    float angular_damping{0.02f};
    Material material{};
    Shape shape{Circle{}};
    bool trigger{false};
    std::uint32_t collision_group{1}, collision_mask{0xffffffffu};
    bool bullet{false}; // Also sweep against other dynamic bodies.
    bool sleeping{false};
    int sleep_counter{0};
    std::uint64_t user_data{};
    bool use_default_shape{true};
    bool fixed_rotation{false};
    // Fixture-driven mass data is recomputed on every fixture change while
    // this is set. Clear it to keep hand-assigned mass properties.
    bool auto_mass{true};
    Vec2 force{};
    float torque{};
    std::vector<std::unique_ptr<Fixture>> fixtures;
    bool is_dynamic() const { return type == BodyType::Dynamic && inverse_mass > 0 && !sleeping; }
    Vec2 center_of_mass() const {
        return transform.position + rotate(local_center, Rot(transform.angle));
    }
    Vec2 center_of_mass(Rot rotation) const { return transform.position + rotate(local_center, rotation); }
    // Rotate in place about the centroid: the body origin has to compensate so
    // the centroid stays where the physics put it.
    void rotate_by(float delta_angle) {
        if (local_center.x == 0 && local_center.y == 0) {
            transform.angle += delta_angle;
            return;
        }
        const Vec2 before = rotate(local_center, Rot(transform.angle));
        transform.angle += delta_angle;
        transform.position += before - rotate(local_center, Rot(transform.angle));
    }
    // Assign explicit mass properties, bypassing fixture-derived data.
    Body& set_mass_data(const MassData& data) {
        auto_mass = false;
        if (type != BodyType::Dynamic) {
            mass = inverse_mass = inertia = inverse_inertia = 0;
            local_center = {};
            return *this;
        }
        mass = data.mass;
        inertia = data.inertia;
        local_center = data.center;
        inverse_mass = mass > 0 ? 1.0f / mass : 0.0f;
        inverse_inertia = inertia > 0 ? 1.0f / inertia : 0.0f;
        return *this;
    }
    // Recompute mass, centroid and inertia from every fixture's shape, density
    // and local offset. Static/kinematic bodies carry no mass.
    Body& reset_mass_data() {
        if (type != BodyType::Dynamic) {
            mass = inverse_mass = inertia = inverse_inertia = 0;
            local_center = {};
            return *this;
        }
        MassData total{};
        for (const auto& fixture : fixtures)
            total = combine(total, transformed(mass_data(fixture->shape, fixture->material.density),
                                               fixture->local));
        auto_mass = true;
        if (total.mass <= 0) {
            // No mass-bearing geometry: stay dynamic but immovable. Explicit
            // set_mass_data() is required to revive the body.
            mass = inertia = 0;
            inverse_mass = inverse_inertia = 0;
            local_center = {};
            return *this;
        }
        set_mass_data(total);
        auto_mass = true;
        return *this;
    }
    void wake() {
        sleeping = false;
        sleep_counter = 0;
    }
    // Controlled mutation paths. They keep the fluent style but also record what
    // changed, so the engine re-derives data only for objects that were actually
    // touched. A controlled edit also wakes the body: a body held still by a
    // joint would otherwise stay asleep and ignore the new value.
    Body &set_position(Vec2 value) {
        transform.position = value;
        ++transform_version;
        wake();
        return *this;
    }
    Body &set_angle(float value) {
        transform.angle = value;
        ++transform_version;
        wake();
        return *this;
    }
    Body &set_transform(Transform value) {
        transform = value;
        ++transform_version;
        wake();
        return *this;
    }
    Body &translate(Vec2 delta) {
        transform.position += delta;
        ++transform_version;
        wake();
        return *this;
    }
    Body &set_shape(Shape value) {
        shape = std::move(value);
        ++geometry_version;
        wake();
        return *this;
    }
    Body &set_type(BodyType value) {
        type = value;
        ++geometry_version;
        wake();
        return *this;
    }
    Body &set_material(const Material &value) {
        material = value;
        ++material_version;
        wake();
        return *this;
    }
    Body &set_filter(std::uint32_t group, std::uint32_t mask) {
        collision_group = group;
        collision_mask = mask;
        ++filter_version;
        wake();
        return *this;
    }
};

inline Fixture &Fixture::set_shape(Shape value) {
    shape = std::move(value);
    ++geometry_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_local(Transform value) {
    local = value;
    ++geometry_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_material(const Material &value) {
    material = value;
    ++material_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_friction(float value) {
    material.friction = value;
    ++material_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_restitution(float value) {
    material.restitution = value;
    ++material_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_density(float value) {
    material.density = value;
    ++material_version;
    if (body)
        body->reset_mass_data();
    return *this;
}
inline Fixture &Fixture::set_filter(std::uint32_t group, std::uint32_t mask) {
    collision_group = group;
    collision_mask = mask;
    ++filter_version;
    if (body)
        body->wake();
    return *this;
}
inline Fixture &Fixture::set_trigger(bool value) {
    trigger = value;
    ++filter_version;
    if (body)
        body->wake();
    return *this;
}

// Snapshot of everything a joint's velocity-phase geometry is derived from.
//
// The velocity phase only writes velocities: it never moves a body, and
// integration runs after it. The world anchors, the lever arms, the constraint
// axis and the effective masses a joint builds are therefore identical on
// every solver iteration, and rebuilding them meant paying two sine/cosine
// pairs, a square root and the effective-mass assembly per joint per iteration
// to arrive at numbers that had not changed.
//
// It is deliberately a comparison of the inputs themselves rather than a
// revision counter. The engine's own position integration writes `transform`
// directly and so does not move `transform_version`, and `set_mass_data()`
// moves no counter at all. Comparing the values that were actually used is
// correct for every write path.
//
// The two derived flags belong in the key just as much as the geometry does.
// `movable_*` comes from `is_dynamic()` (which reads the sleeping flag) and the
// inverse inertia is suppressed by `fixed_rotation`; a body that *wakes* or has
// its rotation unlocked changes neither its position nor its angle, so a
// transform-only key would hand the solver a stale entry that declares both
// endpoints immovable and leave the joint inert for the rest of the run.
struct JointGeometryKey {
    Vec2 position_a{}, position_b{};
    Vec2 anchor_a{}, anchor_b{};
    Vec2 local_center_a{}, local_center_b{};
    float angle_a{}, angle_b{};
    float inverse_mass_a{}, inverse_mass_b{};
    float inverse_inertia_a{}, inverse_inertia_b{};
    bool movable_a{}, movable_b{};
    bool fixed_a{}, fixed_b{};

    void capture(const Body *a, const Body *b, Vec2 local_anchor_a, Vec2 local_anchor_b) {
        position_a = a ? a->transform.position : Vec2{};
        position_b = b ? b->transform.position : Vec2{};
        angle_a = a ? a->transform.angle : 0;
        angle_b = b ? b->transform.angle : 0;
        anchor_a = local_anchor_a;
        anchor_b = local_anchor_b;
        local_center_a = a ? a->local_center : Vec2{};
        local_center_b = b ? b->local_center : Vec2{};
        inverse_mass_a = a ? a->inverse_mass : 0;
        inverse_mass_b = b ? b->inverse_mass : 0;
        inverse_inertia_a = a ? a->inverse_inertia : 0;
        inverse_inertia_b = b ? b->inverse_inertia : 0;
        movable_a = a && a->is_dynamic();
        movable_b = b && b->is_dynamic();
        fixed_a = a && a->fixed_rotation;
        fixed_b = b && b->fixed_rotation;
    }
    bool matches(const Body *a, const Body *b, Vec2 local_anchor_a, Vec2 local_anchor_b) const {
        if (!a || !b)
            return false;
        return position_a == a->transform.position && position_b == b->transform.position &&
               angle_a == a->transform.angle && angle_b == b->transform.angle &&
               anchor_a == local_anchor_a && anchor_b == local_anchor_b &&
               local_center_a == a->local_center && local_center_b == b->local_center &&
               inverse_mass_a == a->inverse_mass && inverse_mass_b == b->inverse_mass &&
               inverse_inertia_a == a->inverse_inertia && inverse_inertia_b == b->inverse_inertia &&
               movable_a == a->is_dynamic() && movable_b == b->is_dynamic() &&
               fixed_a == a->fixed_rotation && fixed_b == b->fixed_rotation;
    }
};

struct DistanceJoint {
    Body *a{};
    Body *b{};
    float length{1};
    float stiffness{1}; // Projection fraction for legacy distance joints.
    Vec2 anchor_a{}, anchor_b{};
    float spring_stiffness{0}, damping{0};
    bool collide_connected{true};
    std::uint64_t user_data{};
    virtual ~DistanceJoint() = default;
    // Warm starting hook, run once per step in the velocity phase.
    virtual void warm_start() {}
    // Velocity phase, run before positions are integrated. Solving here is what
    // keeps a joint from injecting energy: the position update then uses a
    // velocity that already satisfies the constraint, exactly like the contact
    // solver does.
    virtual void solve(float dt) {
        // Every geometric term below is a function of the two transforms only,
        // and the velocity phase never touches a transform: integration runs
        // after it. The cache therefore holds for `warm_start()` and for all
        // `solver_iterations` calls, and only the velocity terms are re-read.
        const AnchorGeometry &g = anchor_geometry();
        if (!g.usable)
            return;
        const Vec2 va =
            a->velocity + Vec2{-a->angular_velocity * g.ra.y, a->angular_velocity * g.ra.x};
        const Vec2 vb =
            b->velocity + Vec2{-b->angular_velocity * g.rb.y, b->angular_velocity * g.rb.x};
        float impulse = -(vb - va).dot(g.direction) / g.effective_inverse_mass;
        if (spring_stiffness > 0 && dt > 0) {
            const float gamma = 1.0f / (dt * (damping + dt * spring_stiffness));
            impulse = -((vb - va).dot(g.direction) +
                        (g.distance - length) * dt * spring_stiffness * gamma) /
                      (g.effective_inverse_mass + gamma);
        }
        if (a->is_dynamic()) {
            a->velocity -= g.direction * impulse * a->inverse_mass;
            a->angular_velocity -= g.lever_a * impulse * a->inverse_inertia;
        }
        if (b->is_dynamic()) {
            b->velocity += g.direction * impulse * b->inverse_mass;
            b->angular_velocity += g.lever_b * impulse * b->inverse_inertia;
        }
    }
    // Position phase, run after positions are integrated.
    virtual void solve_position(float dt) {
        if (!a || !b)
            return;
        if (spring_stiffness > 0)
            return;
        const Rot rot_a(a->transform.angle), rot_b(b->transform.angle);
        const Vec2 anchor_world_a = a->transform.position + rotate(anchor_a, rot_a);
        const Vec2 anchor_world_b = b->transform.position + rotate(anchor_b, rot_b);
        const Vec2 ra = anchor_world_a - a->center_of_mass(rot_a);
        const Vec2 rb = anchor_world_b - b->center_of_mass(rot_b);
        const Vec2 delta = anchor_world_b - anchor_world_a;
        const float d = delta.length();
        if (d < 1.0e-6f)
            return;
        const Vec2 n = delta / d;
        const float ca = ra.cross(n), cb = rb.cross(n);
        const float inv = a->inverse_mass + b->inverse_mass + ca * ca * a->inverse_inertia +
                          cb * cb * b->inverse_inertia;
        if (inv <= 0)
            return;
        const float correction = (d - length) * stiffness / inv;
        if (a->is_dynamic()) {
            a->transform.position += n * correction * a->inverse_mass;
            a->rotate_by(ca * correction * a->inverse_inertia);
        }
        if (b->is_dynamic()) {
            b->transform.position -= n * correction * b->inverse_mass;
            b->rotate_by(-cb * correction * b->inverse_inertia);
        }
        // The position pass moves the bodies, so anything built from the old
        // transforms -- including the parent's velocity cache -- is stale now.
        invalidate_geometry_cache();
    }

  protected:
    // ---- Velocity-phase geometry -------------------------------------------
    // `solve()` runs once per solver iteration, but nothing it reads from the
    // bodies changes during that phase: the velocity pass only writes
    // velocities, and integration happens afterwards. The world anchors, the
    // lever arms, the axis and the effective inverse mass are therefore
    // identical on every iteration. Before this cache existed a hinge chain
    // paid two sine/cosine pairs and a square root per joint per iteration to
    // arrive at the same numbers again.
    //
    struct AnchorGeometry {
        Vec2 ra{}, rb{}, direction{};
        float lever_a{}, lever_b{};
        float distance{};
        float effective_inverse_mass{};
        bool usable{false};
    };
    mutable AnchorGeometry anchor_geometry_{};
    mutable bool anchor_cache_valid_{false};
    mutable JointGeometryKey anchor_cache_key_{};

    const AnchorGeometry &anchor_geometry() const {
        if (anchor_cache_valid_ && anchor_cache_key_.matches(a, b, anchor_a, anchor_b))
            return anchor_geometry_;
        AnchorGeometry built;
        if (a && b) {
            const Rot rot_a(a->transform.angle), rot_b(b->transform.angle);
            const Vec2 anchor_world_a = a->transform.position + rotate(anchor_a, rot_a);
            const Vec2 anchor_world_b = b->transform.position + rotate(anchor_b, rot_b);
            // Lever arms are measured from each centroid, not the body origin.
            built.ra = anchor_world_a - a->center_of_mass(rot_a);
            built.rb = anchor_world_b - b->center_of_mass(rot_b);
            const Vec2 delta = anchor_world_b - anchor_world_a;
            built.distance = delta.length();
            if (built.distance >= 1.0e-6f) {
                built.direction = delta / built.distance;
                built.lever_a = built.ra.cross(built.direction);
                built.lever_b = built.rb.cross(built.direction);
                built.effective_inverse_mass =
                    a->inverse_mass + b->inverse_mass +
                    built.lever_a * built.lever_a * a->inverse_inertia +
                    built.lever_b * built.lever_b * b->inverse_inertia;
                // Mirrors the original `if (inv <= 0) return`, including its
                // behaviour on a non-finite value.
                built.usable = !(built.effective_inverse_mass <= 0);
            }
        }
        anchor_geometry_ = built;
        anchor_cache_key_.capture(a, b, anchor_a, anchor_b);
        anchor_cache_valid_ = true;
        return anchor_geometry_;
    }
    void invalidate_geometry_cache() const { anchor_cache_valid_ = false; }
};
struct SpringJoint : DistanceJoint {};

// A revolute (pin) joint. The two local anchors are forced to stay coincident,
// which removes both relative translation degrees of freedom while leaving
// relative rotation free. `anchor_a`/`anchor_b` are offsets from each body
// origin, exactly like DistanceJoint, so the same anchor helpers apply.
//
// The old placeholder (a zero-length distance joint) only removed the radial
// degree of freedom, so a chain was free to swing sideways and stretch. This
// version solves the full 2x2 point-to-point block, so a hinge actually holds,
// and adds the motor and angle limits that a revolute joint is expected to
// have.
struct HingeJoint : DistanceJoint {
    // Relative angle of the reference pose. The motor drives `angle()` towards
    // zero-based values and the limits are measured from here, so a hinge that
    // was just created sits at angle 0.
    float reference_angle{0};
    float motor_speed{0};      // Target relative angular velocity, rad/s.
    float max_motor_torque{0}; // <= 0 leaves the motor off.
    bool enable_limit{false};
    float lower_angle{0}, upper_angle{0};
    // Bias applied to the limit velocity constraint. Full Baumgarte (1) makes
    // a chainjitter; a mild bias plus the positional pass below converges
    // without injecting energy.
    float limit_bias{0.2f};
    // Largest anchor error the position pass corrects in one go.
    float max_position_correction{0.2f};
    // Accumulated impulses, kept across steps for warm starting.
    Vec2 impulse{};
    float motor_impulse{0}, lower_impulse{0}, upper_impulse{0};

    // Fluent configuration, matching the rest of the library's style. Changing
    // a joint is an external edit, so it wakes the bodies it connects: a joint
    // that is holding a body still would otherwise leave it asleep and ignore
    // the new setting entirely.
    HingeJoint &motor(float speed, float max_torque) {
        motor_speed = speed;
        max_motor_torque = max_torque < 0 ? -max_torque : max_torque;
        wake_bodies();
        return *this;
    }
    HingeJoint &limit(float lower, float upper) {
        enable_limit = true;
        lower_angle = lower < upper ? lower : upper;
        upper_angle = lower < upper ? upper : lower;
        wake_bodies();
        return *this;
    }
    HingeJoint &no_limit() {
        enable_limit = false;
        wake_bodies();
        return *this;
    }
    void wake_bodies() {
        if (a && a->type == BodyType::Dynamic)
            a->wake();
        if (b && b->type == BodyType::Dynamic)
            b->wake();
    }
    // Current relative angle, including the reference pose offset.
    float angle() const {
        if (!a || !b)
            return 0;
        return b->transform.angle - a->transform.angle - reference_angle;
    }
    void warm_start() override {
        // The motor and limit branches below touch `a` and `b` directly, so the
        // null check has to stay here rather than hiding behind `usable`.
        if (!a || !b)
            return;
        // A feature that is switched off must not hand its stale impulse to the
        // warm start, or the body takes one last kick after the user disabled it.
        if (max_motor_torque <= 0)
            motor_impulse = 0;
        if (!enable_limit)
            lower_impulse = upper_impulse = 0;
        const Geometry &g = hinge_geometry();
        if (impulse.x != 0 || impulse.y != 0) {
            if (g.movable_a) {
                a->velocity -= impulse * a->inverse_mass;
                a->angular_velocity -= g.ra.cross(impulse) * g.ia;
            }
            if (g.movable_b) {
                b->velocity += impulse * b->inverse_mass;
                b->angular_velocity += g.rb.cross(impulse) * g.ib;
            }
        }
        if (motor_impulse != 0) {
            a->angular_velocity -= g.ia * motor_impulse;
            b->angular_velocity += g.ib * motor_impulse;
        }
        if (lower_impulse != 0 || upper_impulse != 0) {
            a->angular_velocity += g.ia * (upper_impulse - lower_impulse);
            b->angular_velocity -= g.ib * (upper_impulse - lower_impulse);
        }
    }
    void solve(float dt) override {
        const Geometry &g = hinge_geometry();
        if (!g.usable)
            return;
        // Neither endpoint can take a velocity, so this iteration cannot change
        // anything. Zeroing the accumulators here is the point: an impulse that
        // was accumulated while the pair was movable must not survive into the
        // step where one of them wakes up, or the joint spends it as a kick.
        if (!g.movable_a && !g.movable_b) {
            impulse = {};
            motor_impulse = lower_impulse = upper_impulse = 0;
            return;
        }
        const float inv_dt = dt > 0 ? 1.0f / dt : 0.0f;

        // --- Motor: drive the relative angular velocity, capped by torque. ---
        if (max_motor_torque > 0 && g.axial_mass > 0) {
            const float cdot = b->angular_velocity - a->angular_velocity - motor_speed;
            float delta = -g.axial_mass * cdot;
            const float cap = max_motor_torque * dt;
            const float previous = motor_impulse;
            motor_impulse = std::clamp(previous + delta, -cap, cap);
            delta = motor_impulse - previous;
            if (g.movable_a)
                a->angular_velocity -= g.ia * delta;
            if (g.movable_b)
                b->angular_velocity += g.ib * delta;
        } else {
            motor_impulse = 0;
        }

        // --- Limits: two one-sided angular constraints. ---
        if (enable_limit && g.axial_mass > 0) {
            const float relative = angle();
            // Lower limit: stops the angle from falling below lower_angle.
            float c = relative - lower_angle;
            float cdot = b->angular_velocity - a->angular_velocity;
            float delta = -g.axial_mass * (cdot + std::max(c, 0.0f) * inv_dt * limit_bias);
            float previous = lower_impulse;
            lower_impulse = std::max(previous + delta, 0.0f);
            delta = lower_impulse - previous;
            if (g.movable_a)
                a->angular_velocity -= g.ia * delta;
            if (g.movable_b)
                b->angular_velocity += g.ib * delta;
            // Upper limit: stops the angle from rising above upper_angle.
            c = upper_angle - relative;
            cdot = a->angular_velocity - b->angular_velocity;
            delta = -g.axial_mass * (cdot + std::max(c, 0.0f) * inv_dt * limit_bias);
            previous = upper_impulse;
            upper_impulse = std::max(previous + delta, 0.0f);
            delta = upper_impulse - previous;
            if (g.movable_a)
                a->angular_velocity += g.ia * delta;
            if (g.movable_b)
                b->angular_velocity -= g.ib * delta;
        } else {
            lower_impulse = upper_impulse = 0;
        }

        // --- Point-to-point: the 2x2 block over both translation axes. ---
        // Solving the anchor as a block instead of one axis at a time is what
        // makes the hinge hold sideways instead of only radially.
        if (g.det <= 1.0e-12f)
            return;
        const Vec2 cdot = point_velocity(*b, g.rb) - point_velocity(*a, g.ra);
        const float ix = -(g.k22 * cdot.x - g.k12 * cdot.y) / g.det;
        const float iy = -(-g.k12 * cdot.x + g.k11 * cdot.y) / g.det;
        const Vec2 delta{ix, iy};
        impulse += delta;
        if (g.movable_a) {
            a->velocity -= delta * a->inverse_mass;
            a->angular_velocity -= g.ra.cross(delta) * a->inverse_inertia;
        }
        if (g.movable_b) {
            b->velocity += delta * b->inverse_mass;
            b->angular_velocity += g.rb.cross(delta) * b->inverse_inertia;
        }
    }

    void solve_position(float dt) override {
        if (!a || !b)
            return;
        const Geometry g = geometry();
        if (g.det <= 1.0e-12f || (!g.movable_a && !g.movable_b))
            return;
        Vec2 separation = g.point_b - g.point_a;
        const float error = separation.length();
        if (error <= 1.0e-7f)
            return;
        // Box2D clamps the position error before solving; without it a heavily
        // loaded chain overshoots and oscillates instead of converging.
        const float clamped = std::min(error, max_position_correction);
        separation *= clamped / error;
        // The same 2x2 block inverts the position error, so the rotational
        // coupling of an offset anchor is handled exactly.
        const float px = -(g.k22 * separation.x - g.k12 * separation.y) / g.det * stiffness;
        const float py = -(-g.k12 * separation.x + g.k11 * separation.y) / g.det * stiffness;
        const Vec2 correction{px, py};
        if (g.movable_a) {
            a->transform.position -= correction * g.ma;
            a->rotate_by(-g.ra.cross(correction) * g.ia);
        }
        if (g.movable_b) {
            b->transform.position += correction * g.mb;
            b->rotate_by(g.rb.cross(correction) * g.ib);
        }
    }

  private:
    // Everything the velocity and position phases share: the anchor geometry,
    // the effective 2x2 mass block and its inverse determinant.
    struct Geometry {
        Vec2 point_a{}, point_b{}, ra{}, rb{};
        float ma{}, mb{}, ia{}, ib{};
        float k11{}, k12{}, k22{}, det{};
        float axial_mass{};
        bool movable_a{}, movable_b{};
        // False only when the joint has no endpoints to solve at all. The
        // degenerate *geometry* case is deliberately not folded in here: the
        // motor and the limits are angular constraints that do not depend on
        // the 2x2 block, so they still have to run when its determinant
        // vanishes. `solve()` keeps its own separate `det` test.
        bool usable{false};
    };
    // The velocity phase freezes the transforms, so the whole 2x2 block -- not
    // just the anchors -- is built once per step and reused by the warm start
    // and by every solver iteration. The key mirrors the one DistanceJoint
    // uses; consult JointGeometryKey for why it compares inputs instead of a
    // revision counter, and why the movability flags have to be part of it.
    mutable Geometry velocity_geometry_{};
    mutable bool velocity_cache_valid_{false};
    mutable JointGeometryKey velocity_cache_key_{};
    const Geometry &hinge_geometry() const {
        if (velocity_cache_valid_ && velocity_cache_key_.matches(a, b, anchor_a, anchor_b))
            return velocity_geometry_;
        if (a && b) {
            velocity_geometry_ = geometry();
            velocity_cache_key_.capture(a, b, anchor_a, anchor_b);
        } else {
            velocity_geometry_ = Geometry{};
        }
        velocity_cache_valid_ = true;
        return velocity_geometry_;
    }
    Geometry geometry() const {
        Geometry g;
        g.usable = true;
        const Rot rot_a(a->transform.angle), rot_b(b->transform.angle);
        g.point_a = a->transform.position + rotate(anchor_a, rot_a);
        g.point_b = b->transform.position + rotate(anchor_b, rot_b);
        // Solve at the midpoint between the two anchors. They are coincident
        // while the joint holds, so the choice only matters mid-correction.
        const Vec2 point = (g.point_a + g.point_b) * 0.5f;
        // Lever arms from the centroids, so an offset centroid cannot fake a
        // lighter body.
        g.ra = point - a->center_of_mass(rot_a);
        g.rb = point - b->center_of_mass(rot_b);
        g.ma = mass_inv(*a);
        g.mb = mass_inv(*b);
        g.ia = inertia_inv(*a);
        g.ib = inertia_inv(*b);
        g.movable_a = a->is_dynamic();
        g.movable_b = b->is_dynamic();
        g.k11 = g.ma + g.mb + g.ia * g.ra.y * g.ra.y + g.ib * g.rb.y * g.rb.y;
        g.k12 = -g.ia * g.ra.x * g.ra.y - g.ib * g.rb.x * g.rb.y;
        g.k22 = g.ma + g.mb + g.ia * g.ra.x * g.ra.x + g.ib * g.rb.x * g.rb.x;
        g.det = g.k11 * g.k22 - g.k12 * g.k12;
        g.axial_mass = (g.ia + g.ib) > 0 ? 1.0f / (g.ia + g.ib) : 0.0f;
        return g;
    }
    static float mass_inv(const Body &body) { return body.is_dynamic() ? body.inverse_mass : 0; }
    static float inertia_inv(const Body &body) {
        return body.is_dynamic() && !body.fixed_rotation ? body.inverse_inertia : 0;
    }
    // A sleeping body contributes no velocity, exactly like the contact solver.
    static Vec2 point_velocity(const Body &body, Vec2 r) {
        if (body.type == BodyType::Static || body.sleeping)
            return {};
        return body.velocity + Vec2{-body.angular_velocity * r.y, body.angular_velocity * r.x};
    }
};

inline Transform fixture_transform(const Fixture &fixture, Rot body_rot) {
    const auto &t = fixture.body->transform;
    if (fixture.local.angle == 0)
        return {t.position + rotate(fixture.local.position, body_rot), t.angle};
    return {t.position + rotate(fixture.local.position, body_rot), t.angle + fixture.local.angle};
}

inline Transform fixture_transform(const Fixture &fixture) {
    return fixture_transform(fixture, Rot(fixture.body->transform.angle));
}

class World;
class BodyBuilder {
  public:
    explicit BodyBuilder(World &world) : world_(world) {}
    BodyBuilder &bullet(bool value = true) {
        bullet_ = value;
        return *this;
    }
    BodyBuilder &dynamic() {
        type_ = BodyType::Dynamic;
        return *this;
    }
    BodyBuilder &kinematic() {
        type_ = BodyType::Kinematic;
        mass_ = 0;
        return *this;
    }
    BodyBuilder &static_body() {
        type_ = BodyType::Static;
        mass_ = 0;
        return *this;
    }
    BodyBuilder &at(float x, float y) {
        position_ = {x, y};
        return *this;
    }
    BodyBuilder &angle(float radians) {
        angle_ = radians;
        return *this;
    }
    BodyBuilder &mass(float value) {
        mass_ = std::max(value, 0.0001f);
        return *this;
    }
    // Derive mass from the shape's area at this density, instead of taking the
    // (default 1) explicit mass.
    BodyBuilder &density(float value) {
        density_ = std::max(value, 0.0f);
        return *this;
    }
    BodyBuilder &circle(float radius) {
        shape_ = Circle{radius};
        return *this;
    }
    BodyBuilder &box(float hx, float hy) {
        shape_ = Box{{hx, hy}};
        return *this;
    }
    BodyBuilder &polygon(std::vector<Vec2> vertices) {
        shape_ = Polygon{std::move(vertices)};
        return *this;
    }
    BodyBuilder &convex(std::vector<Vec2> vertices) { return polygon(std::move(vertices)); }
    BodyBuilder &capsule(float radius, float half_length) {
        shape_ = Capsule{radius, half_length};
        return *this;
    }
    BodyBuilder &mesh(std::vector<Polygon> triangles) {
        shape_ = Mesh{std::move(triangles)};
        return *this;
    }
    BodyBuilder &velocity(float x, float y) {
        velocity_ = {x, y};
        return *this;
    }
    BodyBuilder &angular_velocity(float value) {
        angular_velocity_ = value;
        return *this;
    }
    BodyBuilder &collision_filter(std::uint32_t group, std::uint32_t mask) {
        group_ = group;
        mask_ = mask;
        return *this;
    }
    BodyBuilder &friction(float value) {
        material_.friction = value;
        return *this;
    }
    BodyBuilder &restitution(float value) {
        material_.restitution = value;
        return *this;
    }
    BodyBuilder &trigger(bool value = true) {
        trigger_ = value;
        return *this;
    }
    Body &build();

  private:
    World &world_;
    BodyType type_{BodyType::Dynamic};
    Vec2 position_{};
    float angle_{};
    float mass_{1};
    Vec2 velocity_{};
    float angular_velocity_{0};
    Shape shape_{Circle{}};
    Material material_{};
    bool trigger_{false};
    bool bullet_{false};
    float density_{-1};
    std::uint32_t group_{1}, mask_{0xffffffffu};
};

class World {
  public:
    // Upper bound `Config::substeps` is clamped to. Past a few windows the
    // per-window detection cost outgrows anything the extra position passes
    // buy, and an unbounded value is a foot-gun rather than a knob.
    static constexpr int kMaxSubsteps = 64;
    struct StepStatistics {
        double ccd_ms{}, detection_ms{}, velocity_ms{}, position_ms{};
        std::size_t constraints{}, warm_started_points{}, position_clamps{};
        double sleeping_ms{}, cache_ms{};
        std::size_t position_corrections{}, projection_candidates{}, projection_sweeps{};
        std::size_t projection_fast_rejections{};
        std::size_t cached_manifolds{}, sleeping_contacts{}, active_constraints{};
        std::size_t awake_bodies{}, sleep_groups{}, moving_groups{}, settling_groups{};
        std::size_t solver_islands{};
        std::size_t stationary_steps{};
        float max_speed{}, max_angular_speed{}, max_penetration{}, max_correction{};
        // Contact-management breakdown. A step walks the candidate pairs up to
        // three times (velocity build, post-integration build, event
        // transitions). These counters say where `detection_ms` goes and how
        // much of the work the reuse paths avoided, so a change to contact
        // management can be judged without guessing.
        //
        //   pairs_ms      candidate enumeration, fixture pairing, filtering and
        //                 the persistent-cache lookup (broadphase_ms is the
        //                 proxy re-filing and partner-list refresh inside it)
        //   narrow_ms     geometry detection (test / separation / manifold)
        //   assemble_ms   constraint assembly and warm-start matching
        //   events_ms     contact/trigger transition bookkeeping
        //   save_ms       flushing the built contacts into the persistent cache
        double pairs_ms{}, narrow_ms{}, assemble_ms{}, events_ms{}, save_ms{};
        double broadphase_ms{};
        std::size_t candidate_pairs{}, fixture_pairs{}, narrow_tests{}, reused_manifolds{};
        // Steps on which the contact set was provably unchanged, so the whole
        // build (and its cache round trip) was skipped.
        std::size_t unchanged_contact_steps{};
        // Steps on which only the contacts of the bodies that moved were
        // re-detected, and every surviving constraint was rewritten in place
        // instead of the whole set being rebuilt.
        std::size_t refreshed_contact_steps{};
    };
    struct ContactDiagnostic {
        Body *a{}, *b{};
        Vec2 velocity_a{}, velocity_b{}, normal{};
        float penetration{};
        Transform before_a{}, before_b{}, after_a{}, after_b{};
        bool after_integration{}; // false denotes a position-solver correction.
    };
    // Optional first-error trace at post-integration detection and position updates.
    std::function<void(const ContactDiagnostic &)> on_contact_diagnostic;
    const StepStatistics &step_statistics() const { return step_statistics_; }
    double simulation_time() const { return simulation_time_; }

    struct Config {
        Vec2 gravity{0, -9.81f};
        int solver_iterations{8};
        float fixed_timestep{1.0f / 60.0f};
        SolverMode solver_mode{SolverMode::Impulse};
        bool enable_broadphase{true};
        // Extra width of the dynamic-tree fat AABBs. Negative selects the
        // automatic margin derived from the CCD tolerance.
        float broadphase_fat_margin{-1};
        float max_position_correction{0.2f};
        CcdSettings ccd{};
        // Integration/solve windows per step. 1 is the classic single-window
        // step and reproduces every existing trajectory bit for bit. Raising it
        // trades solver iterations for sub-stepping, which is the lever that
        // actually buys stacking accuracy per unit of time: see
        // benchmarks/bench_2d_solver_quality.cpp for the measurements.
        //
        // Declared last on purpose. `Config` is an aggregate that callers
        // initialise positionally, so inserting a field anywhere but the end
        // would silently reinterpret their arguments.
        int substeps{1};
        // Re-detect only the contacts of the bodies that moved, rewriting the
        // constraints in place, instead of rebuilding the whole contact set
        // whenever anything did. It is a pure implementation choice and never
        // changes the result -- a scene where most bodies move falls back to
        // the rebuild on its own -- so it exists here to let a test run both
        // ways and compare.
        bool incremental_contacts{true};
        // Contact compliance, in newtons per metre. 0 makes the velocity pass
        // solve the normal constraint against its exact effective mass.
        //
        // A positive value adds the standard soft-constraint relaxation term,
        // `gamma = 1 / (dt^2 * stiffness)`, to the diagonal. What that does here
        // is under-relax the *increment* of the accumulated impulse, because the
        // compliance that actually supports a load is position error, and this
        // engine's position pass owns it: the textbook companion term
        // `-gamma * impulse` makes the constraint a real spring but leaves it
        // unable to carry a static load at all, and the stack simply falls. So
        // the honest measurement is that the *column* and *mass ratio* scenes get
        // worse at equal budget (+72.6% and +406.0% of their own primary metric),
        // and the scene it helps is a joint chain whose links also touch (-66.7%),
        // where a relaxed contact stops fighting the pins. Off by default for that
        // reason, and kept because those numbers are the answer to "should the
        // contact be made compliant", which is worth being able to reproduce
        // rather than assert. See benchmarks/bench_2d_solver_quality.cpp.
        float contact_stiffness{0};
        // Solve each island's joints inside the contact iteration instead of
        // sweeping every joint after every contact. A joint and a contact that
        // share a body then see each other's corrections within the same
        // iteration rather than one sweep later.
        //
        // Measured on a 24-link chain whose links also touch (the scene in
        // tests/test_2d_solver_options.cpp, default budget): the worst joint gap
        // during the swing is 0.00447 m interleaved against 0.00865 m separate,
        // so it halves the transient. The *settled* gap is the same reading
        // either way -- both are at the float resolution of a six-metre world --
        // because the equilibrium the two arrangements converge to is the same
        // and only the path to it differs. A chain is therefore the reason to
        // turn this on, and it is not a stacking lever:
        // benchmarks/bench_2d_solver_quality.cpp shows the tall stack and mass
        // ratio scenes (which have no joints) are unaffected, and the steady
        // reading on its own joint chain is at numerical noise either way.
        //
        // The cost is proportional to the island's joint count, not to the
        // scene: joints are solved `solver_iterations` times instead of once.
        // Off by default because it is a different solve order, and the default
        // order is what every existing trajectory was recorded against.
        bool interleave_joints{false};
        // Fraction of a contact's measured overlap that one position sweep
        // removes. The position pass runs `solver_iterations` times per step, so
        // a resting stack reaches the overlap at which this fraction and the
        // sink gravity adds each step balance: the steady penetration of a stack
        // is proportional to 1 / relaxation. 0.2 is the historical default and
        // reproduces every existing trajectory bit for bit; raising it trades
        // jitter in a badly conditioned scene for proportionally less sink, and
        // costs one multiply per contact per sweep. See
        // benchmarks/bench_2d_solver_quality.cpp for the measurement.
        float position_relaxation{0.2f};
    };
    World() : World(Config{}) {}
    explicit World(const Config &config) : config_(config) {}
    BodyBuilder create_body() { return BodyBuilder(*this); }
    DistanceJoint &add_distance_joint(Body &a, Body &b, float length) {
        require_unlocked();
        auto joint = std::make_unique<DistanceJoint>();
        joint->a = &a;
        joint->b = &b;
        joint->length = length;
        joints_.push_back(std::move(joint));
        return *joints_.back();
    }
    SpringJoint &add_spring_joint(Body &a, Body &b, float length, float stiffness = 0.5f) {
        require_unlocked();
        auto joint = std::make_unique<SpringJoint>();
        joint->a = &a;
        joint->b = &b;
        joint->length = length;
        joint->spring_stiffness = stiffness;
        SpringJoint &ref = *joint;
        joints_.push_back(std::move(joint));
        return ref;
    }
    // Hinge with per-body local anchors (offsets from each body origin).
    HingeJoint &add_hinge_joint(Body &a, Body &b, Vec2 anchor_a = {}, Vec2 anchor_b = {}) {
        require_unlocked();
        auto joint = std::make_unique<HingeJoint>();
        joint->a = &a;
        joint->b = &b;
        joint->anchor_a = anchor_a;
        joint->anchor_b = anchor_b;
        // The limits and the motor are measured from the pose the joint was
        // created in, so a fresh hinge reads angle 0.
        joint->reference_angle = b.transform.angle - a.transform.angle;
        HingeJoint &ref = *joint;
        joints_.push_back(std::move(joint));
        return ref;
    }
    // Hinge pinned at one shared world point: the local anchors are derived so
    // both bodies rotate about that exact point.
    HingeJoint &add_hinge_joint_at(Body &a, Body &b, Vec2 world_point) {
        return add_hinge_joint(a, b, to_local(a, world_point), to_local(b, world_point));
    }
    // Engine integration uses stable heap-owned bodies and fixtures. Mutations
    // are rejected while callbacks/solving are active.
    bool locked() const { return locked_; }
    Body &create_empty_body() {
        require_unlocked();
        auto body = std::make_unique<Body>();
        body->use_default_shape = false;
        return add_body(std::move(body));
    }
    Fixture &add_fixture(Body &body, const Fixture &definition) {
        require_unlocked();
        if (body.use_default_shape)
            throw std::invalid_argument("Explicit fixtures require create_empty_body()");
        auto fixture = std::make_unique<Fixture>(definition);
        fixture->body = &body;
        auto &result = *fixture;
        body.fixtures.push_back(std::move(fixture));
        broadphase_dirty_ = true; // The body's proxy set changed.
        // Box2D-style ResetMassData: density-weighted mass, centroid and
        // inertia of the whole compound outline.
        if (body.auto_mass)
            body.reset_mass_data();
        return result;
    }
    void destroy_fixture(Fixture &fixture) {
        require_unlocked();
        // The fixture owns itself through the body's `unique_ptr` list, so the
        // erase below frees it. Everything after that point has to work off the
        // owner, not off the fixture: reaching back through `fixture.body` would
        // read a freed object.
        Body *const body = fixture.body;
        wake_neighbors(*body);
        // Only this fixture's contacts go. Everything the body is not touching
        // keeps its manifolds and its accumulated impulses, so destroying one
        // object no longer makes every other stack in the world re-solve from
        // zero -- which is what a blanket `contact_cache_.clear()` did. The
        // constraint's key carries the two fixture addresses, so a fixture is
        // identified there exactly.
        const auto owner = reinterpret_cast<std::uintptr_t>(&fixture);
        drop_constraints_for([&](const Constraint &c) {
            return c.key[2] == owner || c.key[3] == owner;
        });
        contact_cache_.erase_if([&](const ContactKey &key, const CachedContact &) {
            return key[2] == owner || key[3] == owner;
        });
        forget_contacts(&fixture);
        auto &fixtures = body->fixtures;
        std::erase_if(fixtures, [&](const auto &p) { return p.get() == &fixture; });
        if (body->auto_mass)
            body->reset_mass_data();
        // The body's proxy set changed, so only this body is re-filed; the
        // obstacle index notices a changed fixture list on its own.
        reindex_body(*body);
    }
    // A handle to a body that survives the body's own destruction. The slot is
    // reused, so the index alone is not an identity: the generation is what
    // makes a stale handle detectable instead of silently naming whichever body
    // took the slot over.
    struct BodyId {
        std::uint32_t index{0xffffffffu};
        std::uint64_t generation{0};
        bool valid() const { return index != 0xffffffffu; }
        friend bool operator==(const BodyId &, const BodyId &) = default;
    };
    BodyId id(const Body &body) const {
        return {body.slot, body.slot < generations_.size() ? generations_[body.slot] : 0};
    }
    // The body a handle names, or nullptr if the handle predates a destroy.
    Body *body(BodyId handle) {
        if (!alive(handle))
            return nullptr;
        return bodies_[handle.index].get();
    }
    const Body *body(BodyId handle) const {
        if (!alive(handle))
            return nullptr;
        return bodies_[handle.index].get();
    }
    bool alive(BodyId handle) const {
        return handle.index < bodies_.size() && handle.index < generations_.size() &&
               generations_[handle.index] == handle.generation;
    }
    void destroy_joint(DistanceJoint &joint) {
        require_unlocked();
        wake_neighbors(*joint.a);
        wake_neighbors(*joint.b);
        std::erase_if(joints_, [&](const auto &p) { return p.get() == &joint; });
    }
    void destroy_body(Body &body) {
        require_unlocked();
        wake_neighbors(body);
        const auto it = std::find_if(bodies_.begin(), bodies_.end(),
                                     [&](const auto &p) { return p.get() == &body; });
        if (it == bodies_.end())
            return;
        const std::size_t index = std::size_t(it - bodies_.begin());
        const std::size_t last = bodies_.size() - 1;
        const Body *dead = &body;
        // 1. Contacts, cached manifolds, accumulated impulses and joints that
        //    touch this body -- and nothing else. The work, and the disturbance,
        //    is confined to this body's neighbourhood.
        drop_constraints_for(
            [&](const Constraint &c) { return c.a == dead || c.b == dead; });
        contact_cache_.erase_if([&](const ContactKey &key, const CachedContact &) {
            return key[0] == reinterpret_cast<std::uintptr_t>(dead) ||
                   key[1] == reinterpret_cast<std::uintptr_t>(dead);
        });
        std::erase_if(joints_, [&](const auto &p) { return p->a == dead || p->b == dead; });
        for (auto &fixture : body.fixtures)
            forget_contacts(fixture.get());
        // 2. Triggers, whose keys are slots: the pairs this body was in are over,
        //    and every pair the body that takes its slot was in moves with it.
        std::unordered_set<std::uint64_t> remapped;
        for (auto key : active_triggers_) {
            const auto [i, j] = unpack_key(key);
            if (i == index || j == index) {
                if (on_trigger)
                    emit_trigger(*bodies_[i], *bodies_[j], false);
                continue;
            }
            const std::size_t a = i == last ? index : i, b = j == last ? index : j;
            remapped.insert(pair_key(std::min(a, b), std::max(a, b)));
        }
        active_triggers_ = std::move(remapped);
        // 3. Broad phase. The destroyed body leaves the tree and every partner
        //    list; the body that moves into the freed slot has its proxies
        //    re-pointed rather than rebuilt, and its own partner list carried
        //    over. Both halves are proportional to the two bodies involved, where
        //    the previous `broadphase_dirty_ = true` re-filed and re-paired the
        //    whole world.
        const std::vector<std::size_t> dead_partners = list_or_empty(neighbors_, index);
        std::vector<std::size_t> moved_partners;
        drop_body_proxies(index);
        drop_obstacle_proxies(index);
        if (index != last) {
            // Re-point the moved body's proxies at the freed slot. The proxies
            // themselves stay where they are in both trees: only the body they
            // name changes, which is what keeps this a local edit.
            if (const std::uint32_t *first = at_or_null(body_first_proxy_, last))
                for (std::uint32_t p = *first; p != kNoFixture; p = proxies_[p].next)
                    proxies_[p].body = std::uint32_t(index);
            if (const std::uint32_t *first = at_or_null(obstacle_first_proxy_, last))
                for (std::uint32_t p = *first; p != kNoFixture; p = obstacle_proxies_[p].next)
                    obstacle_proxies_[p].body = std::uint32_t(index);
            moved_partners = list_or_empty(neighbors_, last);
            std::erase(moved_partners, index);
        }
        for (std::size_t other : dead_partners)
            erase_sorted(neighbors_[other], index);
        if (index != last)
            for (std::size_t other : moved_partners) {
                erase_sorted(neighbors_[other], last);
                insert_sorted(neighbors_[other], index);
            }
        broadphase_candidate_count_ -= dead_partners.size();
        // 4. The body and every per-slot view of it. Anything that is not moved
        //    here would keep the destroyed body's snapshot and report the body
        //    that took the slot over as edited on the next scan.
        bodies_[index] = std::move(bodies_[last]);
        bodies_.pop_back();
        if (index != last)
            bodies_[index]->slot = std::uint32_t(index);
        move_slot(body_first_proxy_, index, last);
        move_slot(body_bounds_, index, last);
        move_slot(sync_transform_, index, last);
        move_slot(sync_type_, index, last);
        move_slot(sync_hash_, index, last);
        move_slot(sync_default_, index, last);
        move_slot(sync_sleeping_, index, last);
        move_slot(sync_fixtures_, index, last);
        move_slot(edit_transform_, index, last);
        move_slot(edit_type_, index, last);
        move_slot(edit_transform_version_, index, last);
        move_slot(edit_geometry_version_, index, last);
        move_slot(edit_fixtures_, index, last);
        move_slot(dirty_stamp_, index, last);
        move_slot(neighbors_, index, last);
        // The views only exist once the broad phase has been synced. Writing the
        // moved body's partner list into a view that was never built is what a
        // body destroyed out of a freshly built scene would do.
        if (index != last && index < neighbors_.size())
            neighbors_[index] = std::move(moved_partners);
        move_slot(pair_stamp_, index, last);
        move_slot(member_stamp_, index, last);
        move_slot(query_stamp_, index, last);
        move_slot(generations_, index, last);
        move_slot(obstacle_first_proxy_, index, last);
        move_slot(obstacle_transform_, index, last);
        move_slot(obstacle_types_, index, last);
        move_slot(obstacle_fixtures_, index, last);
        // Any handle to this slot -- the destroyed body's or the one that just
        // moved in -- is stale from here on. Only a slot that is still there
        // needs the bump: a slot that was the last one is gone, and the create
        // that takes its place bumps it again.
        if (index != last)
            generations_[index] = ++generation_counter_;
        compact_proxies_if_needed();
        // Slot-indexed views that are cheaper to invalidate than to move: the CCD
        // collider list holds Fixture pointers, and the next contact build must
        // not take the "nothing changed" shortcut past a body that is gone.
        ccd_view_invalid_ = true;
        refresh_stamp_.clear();
        moved_bodies_.clear();
        affected_pairs_.clear();
        velocity_geometry_.clear();
        dirty_bodies_.clear();
    }
    void step(float dt = -1.0f) {
        if (dt < 0)
            dt = config_.fixed_timestep;
        if (!std::isfinite(dt) || dt <= 0)
            return;
        require_unlocked();
        struct Lock {
            bool &flag;
            Lock(bool &f) : flag(f) { flag = true; }
            ~Lock() { flag = false; }
        } lock(locked_);
        step_statistics_ = {};
        contact_build_ran_ = false;
        // Everything below may move a body. Nothing has moved yet, so this is
        // also the point at which a public write made between two steps becomes
        // visible to the broad-phase snapshot scan.
        ++motion_epoch_;
        // New step, new hash epoch: every shape is hashed at most once from here
        // until the next step, however many contacts reference it.
        ++hash_epoch_;
        auto cache_start = Clock::now();
        invalidate_edited_contacts();
        step_statistics_.cache_ms += milliseconds(cache_start);
        // Sub-stepping splits the frame into `substeps` integration/solve
        // windows of dt/substeps each. The persistent contact cache and the
        // accumulated impulses carry from one window to the next exactly as
        // they carry between frames, so every window warm starts from the
        // equilibrium its predecessor reached. `Config::substeps` defaults to
        // 1, which is a single window and reproduces the step bit for bit.
        const int substeps = std::clamp(config_.substeps, 1, kMaxSubsteps);
        const float window = dt / float(substeps);
        for (int s = 0; s < substeps; ++s)
            advance(window, s == 0, s + 1 == substeps);
        auto sleep_start = Clock::now();
        update_sleep();
        step_statistics_.sleeping_ms += milliseconds(sleep_start);
        // Record the state this step ended in. Next step anything that differs
        // from it was written from outside the engine.
        capture_edit_state();
        // Only a step that rebuilt contacts can have made the persistent cache
        // stale. A step that reused them leaves the cache exactly as it was,
        // which is what makes a long resting stretch cost nothing.
        if (contact_build_ran_)
            save_constraints(window);
        simulation_time_ += dt;
    }
    // One integration/solve window of `step()`.
    //
    // `first` marks the window that applies this frame's external forces:
    // gravity plus the accumulated force/torque are frame-level inputs, so
    // applying them in every window would multiply them by the substep count.
    // `last` marks the window that publishes contact and trigger transitions,
    // which are still reported once per frame.
    //
    // Everything else -- detection, warm starting, the velocity solve, the
    // position integration, re-detection and the position solve -- runs once
    // per window. That is the whole point of sub-stepping: a stack corrected by
    // one large position pass against a pose the solver has already
    // over-corrected is instead corrected by `substeps` small passes, each one
    // taken from a pose the previous window left in equilibrium.
    void advance(float dt, bool first, bool last) {
        for (auto &p : bodies_) {
            if (first) {
                if (p->sleeping && p->type == BodyType::Dynamic &&
                    (p->force.length_squared() > 0 || p->torque != 0 ||
                     p->velocity.length_squared() > 0 || p->angular_velocity != 0))
                    p->wake();
                if (p->is_dynamic()) {
                    p->velocity += (config_.gravity + p->force * p->inverse_mass) * dt;
                    p->angular_velocity += p->torque * p->inverse_inertia * dt;
                }
                p->force = {};
                p->torque = 0;
            }
            if (p->is_dynamic()) {
                p->velocity *= std::max(0.0f, 1.0f - p->linear_damping * dt);
                p->angular_velocity *= std::max(0.0f, 1.0f - p->angular_damping * dt);
            }
            if (p->fixed_rotation)
                p->angular_velocity = 0;
        }
        // Resting contacts absorb gravity before CCD sees them as new impacts.
        // Joints warm start alongside them: applying the previous window's
        // accumulated impulse up front is what turns a hinge chain from soft to
        // rigid.
        for (auto &joint : joints_)
            joint->warm_start();
        build_constraints(dt, true);
        build_solver_islands();
        solve_velocities(dt);
        // Joint velocities are solved *before* positions are integrated, so the
        // integration never uses a velocity the joint already rejects. Doing it
        // afterwards (the old order) made the position pass do the velocity
        // pass's job and pumped energy into every swing.
        //
        // With `interleave_joints` they were already solved, inside the island
        // sweep above.
        if (!config_.interleave_joints)
            solve_joint_velocities(dt);
        // The first detection still validates public edits and filtering every
        // step. With no awake dynamics or moving kinematics, no geometry can
        // change during integration: retain its constraints and event candidates.
        bool stationary = std::all_of(bodies_.begin(), bodies_.end(), [](const auto &b) {
            return (b->type != BodyType::Dynamic || b->sleeping) &&
                   (b->type != BodyType::Kinematic ||
                    (b->velocity.length_squared() == 0 && b->angular_velocity == 0));
        });
        // `step_continuous` publishes its own counters, so a frame with several
        // windows needs to fold them together rather than let the last window
        // overwrite the rest. With one window this is an identity.
        if (first)
            ccd_frame_totals_ = {};
        if (stationary) {
            CcdStatistics idle;
            idle.advanced_time = dt;
            ccd_frame_totals_.accumulate(idle);
            ccd_statistics_ = ccd_frame_totals_;
            ++step_statistics_.stationary_steps;
        } else {
            auto ccd_start = Clock::now();
            step_continuous(dt);
            step_statistics_.ccd_ms += milliseconds(ccd_start);
            ccd_frame_totals_.accumulate(ccd_statistics_);
            ccd_statistics_ = ccd_frame_totals_;
            // Integration moved bodies, so the broad-phase snapshot scan has to
            // run again before anyone may claim the contacts are unchanged.
            ++motion_epoch_;
            build_constraints(dt, false);
            build_solver_islands();
            solve_velocities(dt);
        }
        // Transitions are a frame-level report: they are published once, from
        // the pose the frame ended in. With a single window this is the same
        // place as before.
        if (last) {
        auto event_start = Clock::now();
        if (contact_build_ran_) {
        current_triggers_.clear();
        current_contacts_.clear();
        for_each_candidate_pair([&](std::size_t i, std::size_t j) {
            Body &a = *bodies_[i];
            Body &b = *bodies_[j];
            visit_pairs(a, b, [&](Fixture &fa, Fixture &fb) {
                if (a.use_default_shape && b.use_default_shape && !fa.trigger && !fb.trigger)
                    return; // Legacy solid bodies have no fixture contact events.
                Contact c;
                if (!allowed(fa, fb))
                    return;
                // Velocity solving has not changed geometry since detection.
                // Solid event membership uses that same result; triggers still
                // need their own overlap test because they create no constraints.
                bool touching =
                    !fa.trigger && !fb.trigger
                        ? std::binary_search(solid_contacts_.begin(), solid_contacts_.end(),
                                             contact_key(fa, fb))
                        : test(fa.shape, fixture_transform(fa), fb.shape, fixture_transform(fb), c);
                if (!touching)
                    return;
                if (!contact_cache_.contains(contact_key(fa, fb))) {
                    if (fa.body->type == BodyType::Dynamic && fa.body->sleeping)
                        fa.body->wake();
                    if (fb.body->type == BodyType::Dynamic && fb.body->sleeping)
                        fb.body->wake();
                }
                if (a.use_default_shape && b.use_default_shape) {
                    if (fa.trigger || fb.trigger)
                        current_triggers_.push_back(pair_key(i, j));
                } else if (!a.use_default_shape && !b.use_default_shape) {
                    auto key = std::make_pair(&fa, &fb);
                    current_contacts_.push_back(key);
                    if (active_contacts_.insert(key).second && on_contact)
                        on_contact(fa, fb, true);
                }
            });
        });
        std::sort(current_contacts_.begin(), current_contacts_.end());
        for (auto it = active_contacts_.begin(); it != active_contacts_.end();) {
            if (!std::binary_search(current_contacts_.begin(), current_contacts_.end(), *it)) {
                auto key = *it;
                it = active_contacts_.erase(it);
                if (on_contact)
                    on_contact(*key.first, *key.second, false);
            } else
                ++it;
        }
        std::sort(current_triggers_.begin(), current_triggers_.end());
        for (auto key : current_triggers_)
            if (active_triggers_.insert(key).second && on_trigger) {
                auto [i, j] = unpack_key(key);
                emit_trigger(*bodies_[i], *bodies_[j], true);
            }
        for (auto it = active_triggers_.begin(); it != active_triggers_.end();) {
            if (!std::binary_search(current_triggers_.begin(), current_triggers_.end(), *it)) {
                auto [i, j] = unpack_key(*it);
                it = active_triggers_.erase(it);
                if (on_trigger)
                    emit_trigger(*bodies_[i], *bodies_[j], false);
            } else
                ++it;
        }
        step_statistics_.detection_ms += milliseconds(event_start);
        step_statistics_.events_ms += milliseconds(event_start);
        } // contacts unchanged: event state, triggers and callbacks are unchanged too
        } // not the last window: publish no transitions yet
        // Kinematics have been integrated by now, so this is the pose the
        // position solver will see. Static obstacles are re-filed once and then
        // never touched again.
        auto obstacle_start = Clock::now();
        sync_obstacles();
        step_statistics_.cache_ms += milliseconds(obstacle_start);
        auto position_start = Clock::now();
        for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
            for (const SolverIsland &island : solver_islands_) {
                // Springs fold their whole response into solve(), so they only
                // need the velocity pass; everything else projects its position
                // error here, after the integration has already happened.
                for (std::size_t k = island.joint_begin; k < island.joint_end; ++k) {
                    auto &joint = *joints_[island_joints_[k]];
                    if (joint.spring_stiffness > 0)
                        continue;
                    auto ta = joint.a->transform, tb = joint.b->transform;
                    joint.solve_position(dt);
                    guard_projection(*joint.a, ta);
                    guard_projection(*joint.b, tb);
                }
                for (std::size_t k = island.constraint_begin; k < island.constraint_end; ++k)
                    solve_position(constraints_[island_constraints_[k]]);
            }
        }
        step_statistics_.position_ms += milliseconds(position_start);
    }
    void set_solver_iterations(int iterations) {
        config_.solver_iterations = std::max(1, iterations);
    }
    std::function<bool(const Fixture &, const Fixture &)> contact_filter;
    std::function<void(Fixture &, Fixture &, bool)> on_contact;
    const CcdStatistics &ccd_statistics() const { return ccd_statistics_; }
    std::size_t body_count() const { return bodies_.size(); }
    std::size_t joint_count() const { return joints_.size(); }
    std::size_t broadphase_candidate_count() const { return broadphase_candidate_count_; }    // Broad phase accelerated queries. Both refresh the tree first, so they see
    // public transform edits without requiring a step.
    std::vector<Body *> query_aabb(const AABB &area) const {
        std::vector<Body *> result;
        if (!config_.enable_broadphase) {
            for (const auto &body : bodies_)
                if (body_aabb(*body).overlaps(area))
                    result.push_back(body.get());
            return result;
        }
        // Forced: these calls can arrive between steps, after a public write.
        sync_broadphase(true);
        if (query_stamp_.size() != bodies_.size()) {
            query_stamp_.assign(bodies_.size(), 0);
            query_stamp_value_ = 0;
        }
        ++query_stamp_value_;
        const std::uint32_t stamp = std::uint32_t(query_stamp_value_);
        tree_.query(area, 0.0f, [&](std::uint32_t payload) {
            // One body can own several proxies, so a hit has to be reported once.
            const std::size_t index = proxies_[payload].body;
            if (index >= bodies_.size() || query_stamp_[index] == stamp)
                return;
            query_stamp_[index] = stamp;
            if (body_aabb(*bodies_[index]).overlaps(area))
                result.push_back(bodies_[index].get());
        });
        return result;
    }
    std::optional<RaycastHit>
    raycast(Vec2 origin, Vec2 direction,
            float max_distance = std::numeric_limits<float>::max()) const {
        direction = direction.normalized();
        if (direction.length_squared() == 0 || max_distance < 0)
            return std::nullopt;
        std::optional<RaycastHit> best;
        auto consider = [&](std::size_t index, const Shape &shape, const Transform &t) {
            if (auto hit = ray_shape(origin, direction, max_distance, shape, t, index))
                if (!best || hit->distance < best->distance)
                    best = hit;
        };
        if (!config_.enable_broadphase) {
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                const auto &body = *bodies_[i];
                if (body.use_default_shape)
                    consider(i, body.shape, body.transform);
                else
                    for (const auto &fixture : body.fixtures)
                        consider(i, fixture->shape, fixture_transform(*fixture));
            }
            return best;
        }
        // Forced: a ray cast can arrive between steps, after a public write.
        sync_broadphase(true);
        tree_.raycast(origin, direction, max_distance, [&](std::uint32_t payload, float entry) {
            // Everything behind the closest exact hit can be pruned. The tree's
            // entry distance is measured on the fat box, so it never exceeds the
            // exact distance of any fixture filed under that proxy.
            if (best && entry > best->distance)
                return true;
            const Proxy &proxy = proxies_[payload];
            const std::size_t index = proxy.body;
            if (index >= bodies_.size())
                return true;
            if (proxy.fixture == kNoFixture)
                consider(index, bodies_[index]->shape, bodies_[index]->transform);
            else {
                const Fixture &fixture = *bodies_[index]->fixtures[proxy.fixture];
                consider(index, fixture.shape, fixture_transform(fixture));
            }
            return true;
        });
        return best;
    }
    std::function<void(Body &, Body &, bool)> on_trigger;

  private:
    CcdStatistics ccd_statistics_{};
    // Per-frame total across the integration windows; see `accumulate`.
    CcdStatistics ccd_frame_totals_{};
    CcdWorkspace ccd_workspace_;
    std::vector<CcdMotion> ccd_motions_;
    std::vector<Fixture> ccd_legacy_;
    // Per-body staleness for the CCD view. Rebuilding every body's collider list
    // on every step is pure bookkeeping for a world that is mostly parked, and
    // it is O(bodies) in a path that otherwise only touches what moves: eight
    // thousand static obstacles turned 1.6 ms of CCD into 789 ms over a
    // 300-frame run. A body's view is rebuilt when the broad phase or the edit
    // registry reports a change it caches, so the cost follows what changed.
    std::vector<unsigned char> ccd_view_stale_;
    // Slot of a default-shape body's stand-in Fixture in `ccd_legacy_`. Slots
    // are stable so the pointers stored in the collider list stay valid across
    // steps; the slot map is only reset when the body set itself changes.
    std::vector<std::uint32_t> ccd_legacy_slot_;
    std::size_t ccd_legacy_slots_{0};
    // `sync_broadphase` is const and is what discovers a changed body set, so
    // this flag is written from a const context.
    mutable bool ccd_view_invalid_{true};
    void step_continuous(float dt) {
        auto &motions = ccd_motions_;
        auto &legacy = ccd_legacy_;
        const std::size_t count = bodies_.size();
        // Make sure the edit registry is fresh: `build_constraints` normally
        // synced the broad phase already this epoch, in which case this is a
        // memoised no-op, but a world with a custom contact filter returns from
        // `contacts_unchanged` before it gets there.
        sync_broadphase();
        const bool rebuild_all =
            ccd_view_invalid_ || motions.size() != count || ccd_view_stale_.size() != count;
        if (rebuild_all) {
            ccd_view_invalid_ = false;
            ccd_legacy_slots_ = 0;
            // `legacy` must not be cleared and repopulated element by element:
            // destroying a Fixture frees the heap payload of its Shape, and for
            // a Polygon or Mesh the next assignment has to allocate it again.
            // This branch only runs when the body set itself changed, so that
            // cost is paid once per insertion, not once per step.
            legacy.clear();
            ccd_legacy_slot_.assign(count, std::uint32_t(-1));
            ccd_view_stale_.assign(count, 1);
            motions.clear();
            motions.resize(count);
        } else {
            // Two sources, both small. The broad phase reports what moved or
            // changed shape; the edit registry reports every public write the
            // engine did not make itself, which is what a restitution or a
            // collision mask written straight onto a fixture is. The view caches
            // both, so a view kept alive by the broad phase alone would keep
            // bouncing off a restitution the user has already overwritten.
            for (std::size_t index : edited_bodies_)
                if (index < ccd_view_stale_.size())
                    ccd_view_stale_[index] = 1;
            for (std::size_t index : dirty_bodies_)
                if (index < ccd_view_stale_.size())
                    ccd_view_stale_[index] = 1;
        }
        // Every default-shape body owns at most one slot, so a reservation of
        // `count` keeps the addresses handed to the collider lists stable.
        legacy.reserve(count);
        for (std::size_t index = 0; index < count; ++index) {
            Body &body = *bodies_[index];
            auto &motion = motions[index];
            // Pointer and scalar fields are a handful of stores and can change
            // without the geometry moving -- `set_mass_data` moves the centre of
            // mass and the inverse mass without bumping any version -- so they
            // are refreshed unconditionally.
            motion.transform = &body.transform;
            motion.velocity = &body.velocity;
            motion.angular_velocity = &body.angular_velocity;
            motion.sleeping = &body.sleeping;
            motion.sleep_counter = &body.sleep_counter;
            motion.local_center = body.local_center;
            motion.inverse_mass = body.type == BodyType::Dynamic ? body.inverse_mass : 0;
            motion.inverse_inertia =
                body.type == BodyType::Dynamic && !body.fixed_rotation ? body.inverse_inertia : 0;
            motion.dynamic = body.type == BodyType::Dynamic;
            motion.kinematic = body.type == BodyType::Kinematic;
            motion.bullet = body.bullet;
            if (!ccd_view_stale_[index])
                continue;
            ccd_view_stale_[index] = 0;
            motion.colliders.clear();
            auto append = [&](Fixture &f) {
                motion.colliders.push_back({&f.shape, f.local, f.material.friction,
                                            f.material.restitution, f.restitution_threshold,
                                            f.collision_group, f.collision_mask, f.trigger, &f});
            };
            if (body.use_default_shape) {
                if (ccd_legacy_slot_[index] == std::uint32_t(-1))
                    ccd_legacy_slot_[index] = std::uint32_t(ccd_legacy_slots_++);
                const std::size_t slot = ccd_legacy_slot_[index];
                if (legacy.size() <= slot)
                    legacy.resize(slot + 1);
                Fixture &f = legacy[slot];
                f.body = &body;
                f.shape = body.shape;
                f.local = {};
                f.material = body.material;
                f.trigger = body.trigger;
                f.collision_group = body.collision_group;
                f.collision_mask = body.collision_mask;
                f.restitution_threshold = 0;
                append(f);
            } else
                for (auto &fixture : body.fixtures)
                    append(*fixture);
        }
        legacy.resize(ccd_legacy_slots_);
        // The continuous pass keeps its own index of the geometry that cannot
        // move, so it does not re-derive the bounds of every wall in the world
        // on every iteration. That index is only rebuilt when one of those
        // bodies actually changed, which the broad-phase scan above has already
        // decided: walking the edited list costs what changed, not what exists.
        // Kinematics are deliberately not counted -- they are on the mover list
        // and their bounds are recomputed per iteration regardless -- and a body
        // that leaves the immovable set changes the count the index validates
        // itself against.
        bool stationary_changed = rebuild_all;
        if (!stationary_changed)
            for (std::size_t index : edited_bodies_)
                if (index < count && bodies_[index]->type == BodyType::Static) {
                    stationary_changed = true;
                    break;
                }
        ccd_statistics_ = advance_continuous(
            motions, dt, config_.ccd,
            [&](const CcdCollider &a, const CcdCollider &b) {
                return allowed(*static_cast<Fixture *>(a.tag), *static_cast<Fixture *>(b.tag));
            },
            [&](const CcdCollider &a, const CcdCollider &b, const SweepHit &) {
                auto *fa = static_cast<Fixture *>(a.tag);
                auto *fb = static_cast<Fixture *>(b.tag);
                if (fa->body->use_default_shape || fb->body->use_default_shape)
                    return;
                if (active_contacts_.insert({fa, fb}).second && on_contact)
                    on_contact(*fa, *fb, true);
            },
            &ccd_workspace_, stationary_changed);
    }
    struct CallbackLock {
        bool &flag;
        bool previous;
        explicit CallbackLock(bool &value) : flag(value), previous(value) { flag = true; }
        ~CallbackLock() { flag = previous; }
    };
    bool locked_{false};
    std::set<std::pair<Fixture *, Fixture *>> active_contacts_;
    std::vector<std::pair<Fixture *, Fixture *>> current_contacts_;
    std::vector<std::uint64_t> current_triggers_;
    void require_unlocked() const {
        if (locked_)
            throw std::logic_error("Cannot mutate a stepping physics world");
    }
    // Remove every constraint `drop` accepts, rebuilding the flat point array
    // around the survivors. A surviving constraint keeps its identity, its
    // manifold, its friction and -- because the accumulated impulses live in the
    // points -- its warm-start state, so a scene that is not touching the
    // destroyed body keeps solving from where it was instead of starting from
    // zero. The solve order is preserved, which is what keeps the trajectories
    // of the surviving contacts comparable.
    template <class Drop> void drop_constraints_for(Drop &&drop) {
        if (constraints_.empty())
            return;
        ContactPoints kept;
        kept.reserve(points_.size());
        std::size_t write = 0;
        for (Constraint &c : constraints_) {
            if (drop(c))
                continue;
            const int begin = int(kept.size());
            for (int k = 0; k < c.count; ++k) {
                const auto p = points_.at(std::size_t(c.point_begin + k));
                kept.push(ConstraintPoint{p.local_a, p.local_b, p.normal_impulse, p.tangent_impulse,
                                          p.target, p.feature});
            }
            c.point_begin = begin;
            constraints_[write++] = c;
        }
        if (write == constraints_.size())
            return;
        constraints_.resize(write);
        points_ = std::move(kept);
        points_live_ = points_.size();
        active_constraints_.clear();
        for (std::size_t i = 0; i < constraints_.size(); ++i)
            if (constraints_[i].a->is_dynamic() || constraints_[i].b->is_dynamic())
                active_constraints_.push_back(i);
        solid_contacts_.clear();
        for (const Constraint &c : constraints_)
            if (c.event_contact)
                solid_contacts_.push_back(c.key);
        std::sort(solid_contacts_.begin(), solid_contacts_.end());
        velocity_geometry_.clear();
    }
    // Move the last slot's entry into `index` and drop it. Every per-slot view
    // of a body has to go through this when a slot is reused, or that slot would
    // keep the destroyed body's snapshot and report the survivor as edited.
    template <class T>
    static void move_slot(std::vector<T> &values, std::size_t index, std::size_t last) {
        if (values.size() <= last)
            return;
        if (index != last)
            values[index] = std::move(values[last]);
        values.pop_back();
    }
    // Give the per-slot views an entry for a newly appended slot. Only valid
    // while the broad-phase index is up to date, which the caller checks.
    void grow_slots(std::size_t count) {
        body_first_proxy_.resize(count, kNoFixture);
        body_bounds_.resize(count, AABB{});
        sync_transform_.resize(count);
        sync_type_.resize(count);
        sync_hash_.resize(count, 0);
        sync_default_.resize(count, 0);
        sync_sleeping_.resize(count, 0);
        sync_fixtures_.resize(count);
        edit_transform_.resize(count);
        edit_type_.resize(count);
        edit_transform_version_.resize(count, 0);
        edit_geometry_version_.resize(count, 0);
        edit_fixtures_.resize(count);
        dirty_stamp_.resize(count, 0);
        neighbors_.resize(count);
        pair_stamp_.resize(count, 0);
        member_stamp_.resize(count, 0);
        query_stamp_.resize(count, 0);
    }
    // `values[i]`, or nullptr when that view has not been built yet. A body can
    // be destroyed before the broad phase has ever been synced -- a world that
    // builds a scene and immediately drops something out of it -- and then there
    // are no per-slot views to move. The next sync builds them from the bodies
    // that are still there.
    template <class T> static const T *at_or_null(const std::vector<T> &values, std::size_t i) {
        return i < values.size() ? &values[i] : nullptr;
    }
    static const std::vector<std::size_t> &
    list_or_empty(const std::vector<std::vector<std::size_t>> &values, std::size_t i) {
        static const std::vector<std::size_t> empty;
        return i < values.size() ? values[i] : empty;
    }
    // Take one body's proxies out of the tree. The tree's node ids come from a
    // free list, so destroying one body's proxies leaves every other proxy's id
    // valid: nothing else has to be touched.
    void drop_body_proxies(std::size_t index) const {
        if (index >= body_first_proxy_.size())
            return;
        for (std::uint32_t p = body_first_proxy_[index]; p != kNoFixture; p = proxies_[p].next)
            tree_.destroy_proxy(proxies_[p].node);
        body_first_proxy_[index] = kNoFixture;
        if (index < body_bounds_.size())
            body_bounds_[index] = AABB{};
    }
    // Take one body's obstacle proxies out of the obstacle tree. Same shape as
    // `drop_body_proxies`: the tree's node ids come from a free list, so only
    // this body's proxies have to be removed and nothing else moves.
    void drop_obstacle_proxies(std::size_t index) const {
        if (index >= obstacle_first_proxy_.size())
            return;
        for (std::uint32_t p = obstacle_first_proxy_[index]; p != kNoFixture;
             p = obstacle_proxies_[p].next)
            obstacle_tree_.destroy_proxy(obstacle_proxies_[p].node);
        obstacle_first_proxy_[index] = kNoFixture;
    }
    // Both proxy stores keep dead entries for destroyed bodies, because the
    // entries of every other body are addressed by index and renumbering them
    // would cost exactly the full rebuild this is here to avoid. Compact once the
    // dead ones outnumber the live ones: the rebuild is paid for by the churn
    // that created them, so the amortised cost of a destroy stays constant.
    void compact_proxies_if_needed() {
        if (proxies_.size() > 2 * bodies_.size() + 64)
            broadphase_dirty_ = true;
        if (obstacle_proxies_.size() > 2 * bodies_.size() + 64)
            obstacle_types_.clear(); // Forces sync_obstacles to rebuild the index.
    }
    // Re-file one body against its current fixture set, locally: only this
    // body's proxies are destroyed and recreated, and only its own candidate
    // list is re-derived. Used when a fixture is added to or removed from a body
    // that is already in the index.
    void reindex_body(Body &body) {
        if (broadphase_dirty_ || body_first_proxy_.size() != bodies_.size())
            return; // Nothing filed yet: the next sync builds it from the fixtures.
        const std::size_t index = body.slot;
        const float margin = fat_margin();
        drop_body_proxies(index);
        create_body_proxies(index, margin);
        refresh_pairs(index, margin);
    }
    void forget_contacts(Fixture *fixture) {
        CallbackLock lock(locked_);
        for (auto it = active_contacts_.begin(); it != active_contacts_.end();) {
            if (it->first == fixture || it->second == fixture) {
                auto pair = *it;
                it = active_contacts_.erase(it);
                if (on_contact)
                    on_contact(*pair.first, *pair.second, false);
            } else
                ++it;
        }
    }
    static Fixture legacy_fixture(Body &b) {
        Fixture f;
        f.body = &b;
        f.shape = b.shape;
        f.material = b.material;
        f.trigger = b.trigger;
        f.collision_group = b.collision_group;
        f.collision_mask = b.collision_mask;
        f.restitution_threshold = 0;
        return f;
    }
    // The legacy single-shape body has no Fixture of its own, so a temporary one
    // stands in for it. Both sides used to build that temporary unconditionally,
    // which meant two throwaway Fixture values (each carrying a Shape variant and
    // a Material) for every candidate pair even though the explicit-fixture case
    // never looks at them. They are now built only for the side that needs one.
    template <class F> static void visit_pairs(Body &a, Body &b, F &&fn) {
        if (a.use_default_shape && b.use_default_shape) {
            Fixture fa = legacy_fixture(a), fb = legacy_fixture(b);
            fn(fa, fb);
        } else if (a.use_default_shape) {
            Fixture fa = legacy_fixture(a);
            for (auto &y : b.fixtures)
                fn(fa, *y);
        } else if (b.use_default_shape) {
            Fixture fb = legacy_fixture(b);
            for (auto &x : a.fixtures)
                fn(*x, fb);
        } else
            for (auto &x : a.fixtures)
                for (auto &y : b.fixtures)
                    fn(*x, *y);
    }
    bool allowed(const Fixture &a, const Fixture &b) const {
        if (!(a.collision_group & b.collision_mask) || !(b.collision_group & a.collision_mask))
            return false;
        for (const auto &j : joints_)
            if (!j->collide_connected &&
                ((j->a == a.body && j->b == b.body) || (j->a == b.body && j->b == a.body)))
                return false;
        return !contact_filter || contact_filter(a, b);
    }
    static AABB body_aabb(const Body &b) {
        if (b.use_default_shape)
            return compute_aabb(b.shape, b.transform);
        if (b.fixtures.empty())
            return {b.transform.position, b.transform.position};
        AABB result =
            compute_aabb(b.fixtures.front()->shape, fixture_transform(*b.fixtures.front()));
        for (const auto &f : b.fixtures) {
            auto bounds = compute_aabb(f->shape, fixture_transform(*f));
            result.min.x = std::min(result.min.x, bounds.min.x);
            result.min.y = std::min(result.min.y, bounds.min.y);
            result.max.x = std::max(result.max.x, bounds.max.x);
            result.max.y = std::max(result.max.y, bounds.max.y);
        }
        return result;
    }
    using Clock = std::chrono::steady_clock;
    static double milliseconds(Clock::time_point start) {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    }
    struct ConstraintPoint {
        Vec2 local_a{}, local_b{};
        float normal_impulse{}, tangent_impulse{}, target{};
        unsigned feature{};
    };
    using ContactKey = std::array<std::uintptr_t, 4>;
    // Structure-of-arrays storage for every live contact point. The solver
    // sweeps the whole point set several times per step; keeping the anchors
    // and accumulators in separate contiguous arrays means that sweep does not
    // pull per-contact metadata (transforms, manifolds, shape pointers) through
    // the cache alongside each point.
    struct ContactPoints {
        struct Point {
            Vec2 &local_a;
            Vec2 &local_b;
            float &normal_impulse;
            float &tangent_impulse;
            float &target;
            unsigned &feature;
        };
        struct ConstPoint {
            const Vec2 &local_a;
            const Vec2 &local_b;
            const float &normal_impulse;
            const float &tangent_impulse;
            const float &target;
            const unsigned &feature;
        };
        std::vector<Vec2> local_a, local_b;
        std::vector<float> normal_impulse, tangent_impulse, target;
        std::vector<unsigned> feature;
        void clear() {
            local_a.clear();
            local_b.clear();
            normal_impulse.clear();
            tangent_impulse.clear();
            target.clear();
            feature.clear();
        }
        std::size_t size() const { return local_a.size(); }
        void reserve(std::size_t count) {
            local_a.reserve(count);
            local_b.reserve(count);
            normal_impulse.reserve(count);
            tangent_impulse.reserve(count);
            target.reserve(count);
            feature.reserve(count);
        }
        void push(const ConstraintPoint &point) {
            local_a.push_back(point.local_a);
            local_b.push_back(point.local_b);
            normal_impulse.push_back(point.normal_impulse);
            tangent_impulse.push_back(point.tangent_impulse);
            target.push_back(point.target);
            feature.push_back(point.feature);
        }
        ConstraintPoint value(std::size_t index) const {
            return {local_a[index],      local_b[index],   normal_impulse[index],
                    tangent_impulse[index], target[index], feature[index]};
        }
        Point at(std::size_t index) {
            return {local_a[index],      local_b[index],          normal_impulse[index],
                    tangent_impulse[index], target[index],        feature[index]};
        }
        ConstPoint at(std::size_t index) const {
            return {local_a[index],      local_b[index],          normal_impulse[index],
                    tangent_impulse[index], target[index],        feature[index]};
        }
    };
    struct Constraint {
        ContactKey key{};
        Body *a{}, *b{};
        Vec2 normal{};
        float friction{};
        float projection_radius_a{}, projection_radius_b{};
        int count{};
        int point_begin{};
        const Shape *shape_a{}, *shape_b{};
        // Content hashes captured when the constraint was built, so the cache
        // never has to re-hash an outline while saving or validating it.
        std::uint64_t shape_hash_a{}, shape_hash_b{};
        // Whether each endpoint used its legacy single shape when the key was
        // built. A change here makes the stored fixture tokens meaningless.
        bool default_a{}, default_b{};
        Transform fixture_a{}, fixture_b{};
        Transform detected_a{}, detected_b{};
        Contact contact{};
        Manifold manifold{};
        bool event_contact{};
    };
    struct CachedContact {
        Constraint c;
        // Warm-start snapshot of the contact points. The live points live in
        // the SoA store, which is rebuilt every detection pass, so the cache
        // keeps its own compact copy (never more than two points).
        ConstraintPoint points[2]{};
        float dt{};
        std::size_t generation{};
        // Content signatures instead of shape snapshots. Storing `Shape` by
        // value deep-copied every outline (a Mesh copied all of its triangles)
        // twice per contact per step, and comparing them walked all vertices.
        std::uint64_t shape_hash_a{}, shape_hash_b{};
        Transform transform_a{}, transform_b{};
    };
    // Open-addressed contact cache.
    //
    // This replaced std::map<ContactKey, CachedContact>. The red-black tree
    // allocated a node for every contact pair, compared four words per level,
    // and chased pointers on the hottest path of the step. The keys are four
    // pointers, so hashing them into a flat table turns every lookup into one
    // probe over a contiguous array.
    class FlatContactCache {
      public:
        enum : std::uint8_t { Empty = 0, Tombstone = 1, Occupied = 2 };
        struct Slot {
            ContactKey key{};
            CachedContact value{};
            std::uint8_t state{Empty};
        };

        std::size_t size() const { return size_; }
        std::size_t capacity() const { return slots_.size(); }
        void clear() {
            slots_.clear();
            size_ = 0;
            tombstones_ = 0;
        }
        bool contains(const ContactKey &key) const { return const_cast<FlatContactCache *>(this)->find(key) != nullptr; }
        CachedContact *find(const ContactKey &key) {
            if (slots_.empty())
                return nullptr;
            auto index = slot_of(key);
            for (std::size_t probe = 0; probe < slots_.size(); ++probe) {
                Slot &slot = slots_[index];
                if (slot.state == Empty)
                    return nullptr;
                if (slot.state == Occupied && slot.key == key)
                    return &slot.value;
                index = (index + 1) & (slots_.size() - 1);
            }
            return nullptr;
        }
        const CachedContact *find(const ContactKey &key) const {
            return const_cast<FlatContactCache *>(this)->find(key);
        }
        // Insert-or-get. The value is default constructed on first use.
        CachedContact &operator[](const ContactKey &key) {
            if (slots_.empty())
                grow(64);
            else if ((size_ + tombstones_ + 1) * 10 >= slots_.size() * 7)
                grow(slots_.size() * 2);
            return insert(key).value;
        }
        template <class Predicate> void erase_if(Predicate predicate) {
            for (Slot &slot : slots_)
                if (slot.state == Occupied && predicate(slot.key, slot.value)) {
                    slot.state = Tombstone;
                    slot.value = {};
                    --size_;
                    ++tombstones_;
                }
        }
        template <class Visitor> void each(Visitor visitor) {
            for (Slot &slot : slots_)
                if (slot.state == Occupied)
                    visitor(slot.key, slot.value);
        }

      private:
        static std::uint64_t hash_key(const ContactKey &key) {
            std::uint64_t h = 1469598103934665603ull;
            for (auto word : key) {
                h ^= std::uint64_t(word);
                h *= 1099511628211ull;
            }
            h ^= h >> 29;
            h *= 0xbf58476d1ce4e5b9ull;
            h ^= h >> 32;
            return h;
        }
        std::size_t slot_of(const ContactKey &key) const {
            return std::size_t(hash_key(key) & (slots_.size() - 1));
        }
        Slot &insert(const ContactKey &key) {
            std::size_t index = slot_of(key);
            Slot *reusable = nullptr;
            for (;;) {
                Slot &slot = slots_[index];
                if (slot.state == Occupied) {
                    if (slot.key == key)
                        return slot;
                } else if (slot.state == Tombstone) {
                    if (!reusable)
                        reusable = &slot;
                } else
                    break;
                index = (index + 1) & (slots_.size() - 1);
            }
            // Prefer a tombstone already scanned: it keeps probe chains short.
            if (reusable) {
                reusable->key = key;
                reusable->value = {};
                reusable->state = Occupied;
                --tombstones_;
                ++size_;
                return *reusable;
            }
            Slot &slot = slots_[index];
            slot.key = key;
            slot.value = {};
            slot.state = Occupied;
            ++size_;
            return slot;
        }
        void grow(std::size_t capacity) {
            std::vector<Slot> old = std::move(slots_);
            slots_.assign(capacity, Slot{});
            size_ = 0;
            tombstones_ = 0;
            for (Slot &slot : old)
                if (slot.state == Occupied) {
                    Slot &target = insert(slot.key);
                    target.value = std::move(slot.value);
                }
        }
        std::vector<Slot> slots_;
        std::size_t size_{};
        std::size_t tombstones_{};
    };
    std::size_t cache_generation_{};
    FlatContactCache contact_cache_;
    std::vector<Constraint> constraints_;
    ContactPoints points_;
    std::vector<std::size_t> active_constraints_;
    // Set when this step actually rebuilt the contact set. When it did not, the
    // contacts, the event state and the persistent cache are all still the ones
    // the previous step left behind.
    bool contact_build_ran_{false};
    // Set by `assemble_pair` when an existing constraint cannot express the new
    // state without the constraint array growing or shrinking -- the pair
    // separated, the manifold vanished, or its event membership flipped. The
    // in-place refresh abandons itself and falls back to a full rebuild, which
    // is the only path that can move the array.
    bool refresh_conflict_{false};
    // Scratch for the in-place contact refresh: which bodies moved this build,
    // which candidate pairs touch them, and how many contact points are live in
    // `points_` (the rest of that store is the waste a narrow manifold leaves
    // behind, and the refresh defers to a rebuild once it is worth reclaiming).
    std::vector<std::uint32_t> refresh_stamp_;
    std::uint32_t refresh_stamp_value_{0};
    std::vector<std::size_t> moved_bodies_;
    std::vector<std::pair<std::size_t, std::size_t>> affected_pairs_;
    std::size_t points_live_{0};
    struct VelocityGeometry {
        Vec2 ra[2], rb[2], tangent;
        float normal_mass[2]{}, tangent_mass[2]{};
        float k01{}, determinant{};
        bool movable{};
    };
    std::vector<VelocityGeometry> velocity_geometry_;
    std::vector<ContactKey> solid_contacts_;
    mutable StepStatistics step_statistics_{};
    double simulation_time_{};
    static float inv_mass(const Body &b) { return b.is_dynamic() ? b.inverse_mass : 0; }
    static float inv_inertia(const Body &b) {
        return b.is_dynamic() && !b.fixed_rotation ? b.inverse_inertia : 0;
    }
    // World point to a body-local offset from its origin (what anchors store).
    static Vec2 to_local(const Body &body, Vec2 world_point) {
        return rotate(world_point - body.transform.position, inverse(Rot(body.transform.angle)));
    }
    static Vec2 point_velocity(const Body &b, Vec2 r) {
        if (b.type == BodyType::Static || b.sleeping)
            return {};
        return b.velocity + Vec2{-b.angular_velocity * r.y, b.angular_velocity * r.x};
    }
    static void apply(Constraint &c, Vec2 ra, Vec2 rb, Vec2 impulse) {
        c.a->velocity -= impulse * inv_mass(*c.a);
        c.b->velocity += impulse * inv_mass(*c.b);
        c.a->angular_velocity -= ra.cross(impulse) * inv_inertia(*c.a);
        c.b->angular_velocity += rb.cross(impulse) * inv_inertia(*c.b);
    }
    void save_constraints(float dt) {
        auto start = Clock::now();
        ++cache_generation_;
        for (auto &c : constraints_) {
            auto &cached = contact_cache_[c.key];
            cached.c = c;
            cached.dt = dt;
            cached.generation = cache_generation_;
            // Copied, not recomputed: the hashes were captured once when the
            // constraint was built.
            cached.shape_hash_a = c.shape_hash_a;
            cached.shape_hash_b = c.shape_hash_b;
            cached.transform_a = c.a->transform;
            cached.transform_b = c.b->transform;
            for (int k = 0; k < c.count; ++k)
                cached.points[k] = points_.value(std::size_t(c.point_begin + k));
        }
        contact_cache_.erase_if([&](const ContactKey &, const CachedContact &item) {
            return item.generation != cache_generation_;
        });
        step_statistics_.cache_ms += milliseconds(start);
    }
    static ContactKey contact_key(const Fixture &a, const Fixture &b) {
        return {reinterpret_cast<std::uintptr_t>(a.body), reinterpret_cast<std::uintptr_t>(b.body),
                a.body->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&a),
                b.body->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&b)};
    }
    static bool same_transform(const Transform &a, const Transform &b) {
        return a.position == b.position && a.angle == b.angle;
    }
    // -----------------------------------------------------------------------
    // Per-step shape hash memo
    //
    // A shape that takes part in many contacts used to be re-hashed for every
    // one of them, several times per step. The hash is now computed at most once
    // per shape per step: the memo lives on the shape's owner and is keyed by a
    // step epoch, so a direct write to the outline is still noticed (the hash is
    // recomputed on the next step) while a big shared floor costs one hash.
    // -----------------------------------------------------------------------
    mutable std::uint32_t hash_epoch_{1};
    static std::uint64_t memo_hash(const Shape &shape, std::uint64_t &cache,
                                   std::uint32_t &epoch, std::uint32_t current) {
        if (epoch != current) {
            cache = shape_hash(shape);
            epoch = current;
        }
        return cache;
    }
    std::uint64_t body_hash(const Body &body) const {
        return memo_hash(body.shape, body.shape_hash_cache, body.shape_hash_epoch, hash_epoch_);
    }
    std::uint64_t fixture_hash(const Fixture &fixture) const {
        return memo_hash(fixture.shape, fixture.shape_hash_cache, fixture.shape_hash_epoch,
                         hash_epoch_);
    }
    // -----------------------------------------------------------------------
    // Edit registry
    //
    // Public fields stay writable, so the engine records what every body looked
    // like at the end of the previous step. At the start of the next step
    // anything that differs was written from outside: the engine has not run yet
    // this step. The controlled setters bump a version counter, which lets the
    // registry skip the field comparison; a direct write is still caught by it.
    // That comparison is the documented cost of the direct-write path.
    // -----------------------------------------------------------------------
    struct FixtureState {
        std::uint64_t hash{0};
        Transform local{};
        float friction{0};
        float restitution{0};
        float threshold{0};
        std::uint32_t group{1};
        std::uint32_t mask{0xffffffffu};
        bool trigger{false};
    };
    mutable std::vector<Transform> edit_transform_;
    mutable std::vector<BodyType> edit_type_;
    mutable std::vector<std::uint32_t> edit_transform_version_, edit_geometry_version_;
    mutable std::vector<std::vector<FixtureState>> edit_fixtures_;
    mutable std::vector<std::uint32_t> dirty_stamp_;
    mutable std::size_t dirty_stamp_value_{0};
    mutable std::vector<std::size_t> dirty_bodies_;
    static FixtureState state_of(const Body &body, const Fixture &fixture, std::uint64_t hash) {
        return {hash,
                fixture.local,
                fixture.material.friction,
                fixture.material.restitution,
                fixture.restitution_threshold,
                fixture.collision_group,
                fixture.collision_mask,
                fixture.trigger};
    }
    static FixtureState state_of(const Body &body, std::uint64_t hash) {
        // The legacy body shape is exposed as a synthetic fixture whose local
        // offset is zero and whose threshold matches legacy_fixture().
        return {hash, {}, body.material.friction, body.material.restitution, 0.0f,
                body.collision_group, body.collision_mask, body.trigger};
    }
    template <class Snapshot>
    static bool state_matches(const Body &body, const Fixture &fixture, const Snapshot &state,
                              std::uint64_t hash) {
        return state.hash == hash && state.local.position == fixture.local.position &&
               state.local.angle == fixture.local.angle &&
               state.friction == fixture.material.friction &&
               state.restitution == fixture.material.restitution &&
               state.threshold == fixture.restitution_threshold &&
               state.group == fixture.collision_group && state.mask == fixture.collision_mask &&
               state.trigger == fixture.trigger;
    }
    static bool state_matches(const Body &body, const FixtureState &state, std::uint64_t hash) {
        return state.hash == hash && state.friction == body.material.friction &&
               state.restitution == body.material.restitution &&
               state.group == body.collision_group && state.mask == body.collision_mask &&
               state.trigger == body.trigger;
    }
    bool body_externally_edited(std::size_t index) const {
        const Body &body = *bodies_[index];
        // A controlled edit is authoritative, so the field comparison is skipped.
        if (body.transform_version != edit_transform_version_[index] ||
            body.geometry_version != edit_geometry_version_[index])
            return true;
        if (body.type != edit_type_[index] ||
            !same_transform(body.transform, edit_transform_[index]))
            return true;
        const auto &snapshot = edit_fixtures_[index];
        if (body.use_default_shape)
            return snapshot.empty() ||
                   !state_matches(body, snapshot.front(), body_hash(body));
        if (snapshot.size() != body.fixtures.size())
            return true;
        for (std::size_t k = 0; k < snapshot.size(); ++k)
            if (!state_matches(body, *body.fixtures[k], snapshot[k], fixture_hash(*body.fixtures[k])))
                return true;
        return false;
    }
    void capture_edit_state() const {
        const std::size_t count = bodies_.size();
        edit_transform_.resize(count);
        edit_type_.resize(count);
        edit_transform_version_.resize(count);
        edit_geometry_version_.resize(count);
        edit_fixtures_.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            const Body &body = *bodies_[i];
            edit_transform_[i] = body.transform;
            edit_type_[i] = body.type;
            edit_transform_version_[i] = body.transform_version;
            edit_geometry_version_[i] = body.geometry_version;
            auto &snapshot = edit_fixtures_[i];
            if (body.use_default_shape) {
                if (snapshot.empty())
                    snapshot.emplace_back();
                snapshot.resize(1);
                snapshot[0] = state_of(body, body_hash(body));
                continue;
            }
            snapshot.resize(body.fixtures.size());
            for (std::size_t k = 0; k < body.fixtures.size(); ++k)
                snapshot[k] = state_of(body, *body.fixtures[k], fixture_hash(*body.fixtures[k]));
        }
    }
    // Rebuild the dirty list for this step. Bodies the engine edited itself
    // (through integration or the position solve) happened after this snapshot,
    // so they are not reported here.
    void refresh_edit_registry() const {
        ++dirty_stamp_value_;
        const std::uint32_t stamp = std::uint32_t(dirty_stamp_value_);
        dirty_bodies_.clear();
        const std::size_t count = bodies_.size();
        if (edit_transform_.size() != count) {
            dirty_stamp_.assign(count, stamp);
            for (std::size_t i = 0; i < count; ++i)
                dirty_bodies_.push_back(i);
            return;
        }
        if (dirty_stamp_.size() != count)
            dirty_stamp_.assign(count, 0);
        for (std::size_t i = 0; i < count; ++i)
            if (body_externally_edited(i)) {
                dirty_stamp_[i] = stamp;
                dirty_bodies_.push_back(i);
            }
    }
    bool dirty(std::size_t index) const {
        return index < dirty_stamp_.size() && dirty_stamp_[index] == dirty_stamp_value_;
    }
    // O(1) re-evaluation of the filter rules of a cached contact. The fixture
    // tokens in the key identify the exact fixtures, so there is no need to walk
    // every fixture pair of the two bodies.
    bool contact_enabled(const Constraint &c) const {
        const auto *fa = reinterpret_cast<const Fixture *>(c.key[2]);
        const auto *fb = reinterpret_cast<const Fixture *>(c.key[3]);
        const bool trigger_a = fa ? fa->trigger : c.a->trigger;
        const bool trigger_b = fb ? fb->trigger : c.b->trigger;
        if (trigger_a || trigger_b)
            return false;
        const std::uint32_t group_a = fa ? fa->collision_group : c.a->collision_group;
        const std::uint32_t mask_a = fa ? fa->collision_mask : c.a->collision_mask;
        const std::uint32_t group_b = fb ? fb->collision_group : c.b->collision_group;
        const std::uint32_t mask_b = fb ? fb->collision_mask : c.b->collision_mask;
        if (!(group_a & mask_b) || !(group_b & mask_a))
            return false;
        for (const auto &j : joints_)
            if (!j->collide_connected &&
                ((j->a == c.a && j->b == c.b) || (j->a == c.b && j->b == c.a)))
                return false;
        if (!contact_filter)
            return true;
        Fixture legacy_a{}, legacy_b{};
        if (!fa) {
            legacy_a = legacy_fixture(*c.a);
            fa = &legacy_a;
        }
        if (!fb) {
            legacy_b = legacy_fixture(*c.b);
            fb = &legacy_b;
        }
        return contact_filter(*fa, *fb);
    }
    // Current content hash of the fixture a cached contact was built from,
    // without re-hashing a shape that has already been hashed this step.
    std::uint64_t cached_shape_hash(const Constraint &c, bool first) const {
        const bool default_shape = first ? c.default_a : c.default_b;
        if (default_shape)
            return body_hash(first ? *c.a : *c.b);
        const auto *fixture =
            reinterpret_cast<const Fixture *>(first ? c.key[2] : c.key[3]);
        return fixture ? fixture_hash(*fixture) : 0;
    }
    void invalidate_edited_contacts() {
        refresh_edit_registry();
        if (!contact_filter) {
            // No user callback: a step in which nothing was written from outside
            // does no cache work at all, and an external edit invalidates its
            // contacts outright. Warm-start impulses belong to the old state, so
            // keeping them across a teleport or an outline edit would be wrong.
            if (dirty_bodies_.empty())
                return;
            contact_cache_.erase_if([&](const ContactKey &, const CachedContact &old) {
                auto &c = old.c;
                if (!dirty(c.a->slot) && !dirty(c.b->slot))
                    return false;
                if (c.a->type == BodyType::Dynamic)
                    c.a->wake();
                if (c.b->type == BodyType::Dynamic)
                    c.b->wake();
                return true;
            });
            return;
        }
        // An arbitrary user callback cannot be introspected, so while one is
        // installed every cached contact is revalidated each step. The per-
        // contact snapshots decide whether the geometry changed and the filter
        // rules are re-run in O(1) through the fixture tokens in the key.
        contact_cache_.erase_if([&](const ContactKey &, const CachedContact &old) {
            auto &c = old.c;
            auto edited = [&](Body *body, bool first, const Transform &transform,
                              std::uint64_t saved, Transform local) {
                const bool default_shape = first ? c.default_a : c.default_b;
                return default_shape != body->use_default_shape ||
                       !same_transform(body->transform, transform) ||
                       cached_shape_hash(c, first) != saved ||
                       (!default_shape &&
                        !same_transform(reinterpret_cast<const Fixture *>(first ? c.key[2]
                                                                               : c.key[3])
                                            ->local,
                                        local));
            };
            const bool changed =
                edited(c.a, true, old.transform_a, old.shape_hash_a, c.fixture_a) ||
                edited(c.b, false, old.transform_b, old.shape_hash_b, c.fixture_b);
            if (!changed && contact_enabled(c))
                return false;
            if (c.a->type == BodyType::Dynamic)
                c.a->wake();
            if (c.b->type == BodyType::Dynamic)
                c.b->wake();
            return true;
        });
    }
    // Content hash of the shape a `visit_pairs` fixture stands for, routed
    // through the per-step memo instead of re-hashing the outline.
    std::uint64_t current_hash(const Fixture &fixture) const {
        return fixture.body->use_default_shape ? body_hash(*fixture.body) : fixture_hash(fixture);
    }
    const Constraint *sleeping_cache(const Fixture &a, const Fixture &b,
                                    const CachedContact *cached) const {
        if (a.body->is_dynamic() || b.body->is_dynamic() || a.body->type == BodyType::Kinematic ||
            b.body->type == BodyType::Kinematic)
            return nullptr;
        if (!cached)
            return nullptr;
        const auto &old = *cached;
        // Bodies and shapes are publicly writable: exact snapshots invalidate
        // the sleep shortcut after teleports, geometry edits or fixture offsets.
        if (old.c.default_a != a.body->use_default_shape ||
            old.c.default_b != b.body->use_default_shape ||
            !same_transform(a.body->transform, old.transform_a) ||
            !same_transform(b.body->transform, old.transform_b) ||
            !same_transform(a.local, old.c.fixture_a) ||
            !same_transform(b.local, old.c.fixture_b) ||
            current_hash(a) != old.shape_hash_a || current_hash(b) != old.shape_hash_b)
            return nullptr;
        return &old.c;
    }
    const Constraint *geometry_cache(const Fixture &a, const Fixture &b,
                                    const CachedContact *cached) const {
        if (!cached)
            return nullptr;
        const auto &old = *cached;
        // Exact detection snapshots only: position correction, integration,
        // public geometry edits and fixture offset changes invalidate reuse.
        if (old.c.default_a != a.body->use_default_shape ||
            old.c.default_b != b.body->use_default_shape ||
            !same_transform(a.body->transform, old.c.detected_a) ||
            !same_transform(b.body->transform, old.c.detected_b) ||
            !same_transform(a.local, old.c.fixture_a) ||
            !same_transform(b.local, old.c.fixture_b) ||
            current_hash(a) != old.shape_hash_a || current_hash(b) != old.shape_hash_b)
            return nullptr;
        return &old.c;
    }
    // -----------------------------------------------------------------------
    // Resting contacts
    //
    // A step in which no body that takes part in a contact moved, no candidate
    // pair was gained or lost, no public field was written and the world's
    // structure did not change cannot have altered a single contact. The build
    // from the previous step is then still exactly right, and the whole of
    // contact management -- pair enumeration, fixture pairing, filtering, the
    // manifold-reuse checks, the warm-start transfer and the cache round trip --
    // is skipped.
    //
    // Bodies that move without touching anything do not force a rebuild. A
    // kinematic platform patrolling a floor far from the stack leaves every
    // resting contact intact; only the cheap part (re-deriving the pair set) is
    // redone, so a body moving *into* range is still noticed.
    // -----------------------------------------------------------------------
    // Record what the contact set was derived from, so the next step can tell
    // whether it is still valid without re-deriving anything.
    void note_build(float dt) const {
        built_body_count_ = bodies_.size();
        built_joint_count_ = joints_.size();
        built_with_filter_ = bool(contact_filter);
        built_dt_ = dt;
    }
    bool contacts_unchanged(float dt) const {
        if (contact_filter || broadphase_dirty_ || !dirty_bodies_.empty())
            return false;
        // Installing or removing the user filter flips every filter answer at
        // once, so the step that does it has to rebuild.
        if (built_with_filter_ != bool(contact_filter))
            return false;
        if (built_body_count_ != bodies_.size() || built_joint_count_ != joints_.size())
            return false;
        // A different timestep rescales the warm-start impulses and the
        // restitution targets, so it is a rebuild.
        if (dt != built_dt_)
            return false;
        // Refresh the pair set first: a pair that appeared or disappeared can
        // add or remove a contact, and a body moving *into* range shows up here.
        if (sync_broadphase())
            return false;
        // A body that moved while sharing no candidate pair with anyone cannot
        // have changed a contact. A neighbouring pair whose gap closed is caught
        // by the check above only if the pair itself is new; when it already
        // existed both endpoints are listed as partners and are rejected here.
        for (std::size_t index : edited_bodies_)
            if (index < neighbors_.size() && !neighbors_[index].empty())
                return false;
        return true;
    }
    // Narrow phase, manifold assembly and impulse transfer for one fixture
    // pair. Returns whether the pair contributes a constraint.
    //
    // `in_place` is the constraint this pair already owns, or nullptr when the
    // pair is being added. Refreshing an existing constraint writes it where it
    // already sits, reusing its point span whenever the new manifold is no
    // wider, so a pair that is merely re-detected relocates nothing. The
    // constraint array is the solve order, and keeping it stable is what lets a
    // refresh produce the same trajectory a full rebuild would.
    //
    // `refresh_conflict_` is set when an existing constraint cannot express the
    // new state in place: the pair separated, the manifold vanished, or its
    // event membership flipped. Those are the cases that need the array to grow
    // or shrink, and the caller answers them with a full rebuild.
    bool assemble_pair(std::size_t i, std::size_t j, Fixture &fa, Fixture &fb, float dt, bool warm,
                       Constraint *in_place) {
        ++step_statistics_.fixture_pairs;
        if (!allowed(fa, fb) || fa.trigger || fb.trigger)
            return false;
        auto write_point = [&](std::size_t index, const ConstraintPoint &source) {
            auto p = points_.at(index);
            p.local_a = source.local_a;
            p.local_b = source.local_b;
            p.normal_impulse = source.normal_impulse;
            p.tangent_impulse = source.tangent_impulse;
            p.target = source.target;
            p.feature = source.feature;
        };
        auto remember = [&](Constraint c, const ConstraintPoint *source, int count) {
            auto radius = [&](std::size_t index) {
                const auto &bounds = body_bounds_[index];
                auto center = bodies_[index]->transform.position;
                Vec2 reach{std::max(std::abs(bounds.min.x - center.x),
                                    std::abs(bounds.max.x - center.x)),
                           std::max(std::abs(bounds.min.y - center.y),
                                    std::abs(bounds.max.y - center.y))};
                // The full body's expanded AABB includes every fixture.
                // Its corner radius encloses any subsequent rotation.
                return reach.length() * 1.001f + 0.001f;
            };
            // Sleeping contacts do not need bounds work. If island
            // propagation wakes them below, -1 uses the full guard.
            c.projection_radius_a = config_.ccd.enabled && c.a->is_dynamic() ? radius(i) : -1;
            c.projection_radius_b = config_.ccd.enabled && c.b->is_dynamic() ? radius(j) : -1;
            c.count = count;
            if (in_place && count <= in_place->count) {
                c.point_begin = in_place->point_begin;
                for (int k = 0; k < count; ++k)
                    write_point(std::size_t(c.point_begin + k), source[k]);
                *in_place = c;
                return;
            }
            c.point_begin = int(points_.size());
            for (int k = 0; k < count; ++k)
                points_.push(source[k]);
            if (in_place)
                *in_place = c;
            else
                constraints_.push_back(c);
        };
        auto *old = contact_cache_.find(contact_key(fa, fb));
        const auto *previous_contact = old;
        auto t_narrow = Clock::now();
        if (auto cached = sleeping_cache(fa, fb, previous_contact)) {
            step_statistics_.narrow_ms += milliseconds(t_narrow);
            ++step_statistics_.reused_manifolds;
            auto t_build2 = Clock::now();
            remember(*cached, previous_contact->points, cached->count);
            (in_place ? *in_place : constraints_.back()).friction =
                std::sqrt(std::max(0.0f, fa.material.friction * fb.material.friction));
            if (in_place && cached->event_contact != in_place->event_contact)
                refresh_conflict_ = true;
            if (cached->event_contact && !in_place)
                solid_contacts_.push_back(contact_key(fa, fb));
            step_statistics_.assemble_ms += milliseconds(t_build2);
            ++step_statistics_.sleeping_contacts;
            return true;
        }
        Contact contact;
        auto ta = fixture_transform(fa), tb = fixture_transform(fb);
        const auto *geometry = geometry_cache(fa, fb, previous_contact);
        bool touching;
        if (geometry) {
            contact = geometry->contact;
            touching = true;
            ++step_statistics_.cached_manifolds;
            ++step_statistics_.reused_manifolds;
        } else {
            ++step_statistics_.narrow_tests;
            touching = test(fa.shape, ta, fb.shape, tb, contact);
        }
        step_statistics_.narrow_ms += milliseconds(t_narrow);
        auto t_build2 = Clock::now();
        bool event_contact = geometry ? geometry->event_contact : touching;
        if (!touching && ccd_detail::supported(fa.shape) && ccd_detail::supported(fb.shape)) {
            auto sep = ccd_detail::separation(fa.shape, ta, fb.shape, tb);
            if (sep.distance > (config_.ccd.enabled ? 2 * config_.ccd.tolerance : 0.002f)) {
                step_statistics_.assemble_ms += milliseconds(t_build2);
                refresh_conflict_ |= in_place != nullptr;
                return false;
            }
            contact = {sep.normal, sep.point, -sep.distance};
            touching = true;
            event_contact = config_.ccd.enabled;
        }
        if (!touching) {
            step_statistics_.assemble_ms += milliseconds(t_build2);
            refresh_conflict_ |= in_place != nullptr;
            return false;
        }
        if (in_place && event_contact != in_place->event_contact)
            refresh_conflict_ = true;
        if (event_contact && !in_place)
            solid_contacts_.push_back(contact_key(fa, fb));
        if (!warm && !in_place && on_contact_diagnostic && contact.penetration > 0)
            on_contact_diagnostic({fa.body, fb.body, fa.body->velocity, fb.body->velocity,
                                   contact.normal, contact.penetration, fa.body->transform,
                                   fb.body->transform, fa.body->transform, fb.body->transform,
                                   true});
        Constraint c;
        c.a = fa.body;
        c.b = fb.body;
        c.normal = contact.normal;
        c.shape_a = fa.body->use_default_shape ? &fa.body->shape : &fa.shape;
        c.shape_b = fb.body->use_default_shape ? &fb.body->shape : &fb.shape;
        c.fixture_a = fa.local;
        c.fixture_b = fb.local;
        c.detected_a = fa.body->transform;
        c.detected_b = fb.body->transform;
        c.contact = contact;
        c.event_contact = event_contact;
        c.key = {reinterpret_cast<std::uintptr_t>(c.a),
                 reinterpret_cast<std::uintptr_t>(c.b),
                 c.a->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&fa),
                 c.b->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&fb)};
        c.default_a = c.a->use_default_shape;
        c.default_b = c.b->use_default_shape;
        c.shape_hash_a = c.default_a ? body_hash(*c.a) : fixture_hash(fa);
        c.shape_hash_b = c.default_b ? body_hash(*c.b) : fixture_hash(fb);
        c.friction = std::sqrt(std::max(0.0f, fa.material.friction * fb.material.friction));
        auto m = geometry ? geometry->manifold
                          : contact_manifold(fa.shape, ta, fb.shape, tb, contact.normal,
                                             contact.point, -contact.penetration);
        c.manifold = m;
        // One rotation/centroid pair per constraint: every contact
        // point shares it, and the solver converts back for free.
        const Rot rot_a(c.a->transform.angle), rot_b(c.b->transform.angle);
        const Vec2 center_a = c.a->center_of_mass(rot_a);
        const Vec2 center_b = c.b->center_of_mass(rot_b);
        const Rot inv_a = inverse(rot_a), inv_b = inverse(rot_b);
        ConstraintPoint built[2]{};
        for (int k = 0; k < m.count; ++k) {
            auto &p = built[k];
            auto &mp = m.points[k];
            p.feature = mp.feature;
            Vec2 pa = mp.point - c.normal * (mp.separation * 0.5f),
                 pb = mp.point + c.normal * (mp.separation * 0.5f);
            const Vec2 ra = pa - center_a, rb = pb - center_b;
            p.local_a = rotate(ra, inv_a);
            p.local_b = rotate(rb, inv_b);
            float speed = (point_velocity(*c.b, rb) - point_velocity(*c.a, ra)).dot(c.normal);
            p.target = -speed > std::max(1.0f, std::min(fa.restitution_threshold,
                                                        fb.restitution_threshold))
                           ? -std::max(fa.material.restitution, fb.material.restitution) * speed
                           : 0;
            if (old && old->c.normal.dot(c.normal) > 0.95f) {
                for (int n = 0; n < old->c.count; ++n) {
                    auto &prev = old->points[n];
                    if (prev.feature == p.feature &&
                        (prev.local_a - p.local_a).length_squared() < 0.04f &&
                        (prev.local_b - p.local_b).length_squared() < 0.04f) {
                        float ratio =
                            old->dt > 0 ? std::clamp(dt / old->dt, 0.0f, 2.0f) : 0;
                        p.normal_impulse = prev.normal_impulse * ratio;
                        p.tangent_impulse = prev.tangent_impulse * ratio;
                        break;
                    }
                }
            }
        }
        if (!m.count) {
            step_statistics_.assemble_ms += milliseconds(t_build2);
            refresh_conflict_ |= in_place != nullptr;
            return false;
        }
        remember(c, built, m.count);
        step_statistics_.assemble_ms += milliseconds(t_build2);
        return true;
    }
    // Bring the contacts of the bodies that moved up to date, in place.
    //
    // A step detects contacts twice: once before integration, and once after
    // it. The second pass exists to bring the manifolds, the anchors and the
    // restitution targets of the pairs that *moved* up to date; everything else
    // about them is already right. Rebuilding the whole contact set for that
    // re-ran the narrow phase, the filter and the fixture pairing for every
    // pair in the world, and moved every constraint in the array -- which is
    // the solve order.
    //
    // This pass visits only the candidate pairs of a body that moved, or that
    // was written from outside, and rewrites their constraints where they sit.
    // It refuses and defers to the full rebuild whenever the contact *set* has
    // to change: a pair that appeared or separated, a manifold that vanished,
    // or an event membership that flipped. Those are the cases where the array
    // has to grow or shrink, and they are rare next to a stack settling.
    bool refresh_moved_contacts(float dt, bool warm) {
        if (!config_.incremental_contacts)
            return false;
        if (contact_filter || broadphase_dirty_ || constraints_.empty())
            return false;
        // Installing or removing the filter flips every filter answer at once.
        if (built_with_filter_ != bool(contact_filter))
            return false;
        if (built_body_count_ != bodies_.size() || built_joint_count_ != joints_.size())
            return false;
        if (dt != built_dt_)
            return false;
        // The candidate pair set is what decides whether a contact can appear
        // or disappear at all; when it changed, the full rebuild is the answer.
        if (sync_broadphase())
            return false;
        // A pair whose manifold narrows keeps its span, so the point store
        // cannot grow much -- but it cannot shrink either, and the way to
        // reclaim it is the rebuild this pass is avoiding. Defer once the waste
        // is larger than the live set.
        if (points_.size() > 4 * points_live_ + 1024)
            return false;
        // Every body whose pose or credentials the engine did not write itself:
        // the edit registry reports the external writes, the broad phase reports
        // what integration moved.
        ++refresh_stamp_value_;
        const std::uint32_t stamp = refresh_stamp_value_;
        if (refresh_stamp_.size() != bodies_.size())
            refresh_stamp_.assign(bodies_.size(), 0);
        auto &moved_list = moved_bodies_;
        moved_list.clear();
        for (std::size_t list = 0; list < 2; ++list) {
            const auto &indices = list ? dirty_bodies_ : edited_bodies_;
            for (std::size_t index : indices)
                if (index < refresh_stamp_.size() && refresh_stamp_[index] != stamp) {
                    refresh_stamp_[index] = stamp;
                    moved_list.push_back(index);
                }
        }
        if (moved_list.empty())
            return false;
        auto &affected = affected_pairs_;
        affected.clear();
        for (std::size_t index : moved_list) {
            if (index >= neighbors_.size())
                continue;
            for (std::size_t partner : neighbors_[index])
                affected.push_back({std::min(index, partner), std::max(index, partner)});
        }
        std::sort(affected.begin(), affected.end());
        affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
        // When nearly every pair moved there is nothing to step over: the
        // rebuild is the cheaper of the two, and it is also what the counters
        // of a fully active scene have always described.
        if (affected.size() * 2 >= constraints_.size())
            return false;
        refresh_conflict_ = false;
        // The constraint array is ordered by the candidate pair each constraint
        // came from, so a pair's constraints are contiguous and both lists can
        // be walked together. The pairs that did not move are stepped over
        // rather than rebuilt: saving and restoring their accumulators through
        // the cache would hand them exactly what they already hold.
        std::size_t live_points = 0;
        std::size_t cursor = 0;
        for (const auto &pair : affected) {
            const std::size_t i = pair.first, j = pair.second;
            while (cursor < constraints_.size() && pair_before(constraints_[cursor], i, j))
                live_points += std::size_t(constraints_[cursor++].count);
            const std::size_t begin = cursor;
            while (cursor < constraints_.size() && constraints_[cursor].a->slot == i &&
                   constraints_[cursor].b->slot == j)
                live_points += std::size_t(constraints_[cursor++].count);
            ++step_statistics_.candidate_pairs;
            const std::size_t end = cursor;
            std::size_t next = begin;
            visit_pairs(*bodies_[i], *bodies_[j], [&](Fixture &fa, Fixture &fb) {
                if (next < end && constraints_[next].key == contact_key(fa, fb)) {
                    assemble_pair(i, j, fa, fb, dt, warm, &constraints_[next]);
                    ++next;
                    return;
                }
                // This fixture pair owns no constraint. If it touches now -- or
                // its filter answer flipped in a way the filter itself did not
                // report -- the array has to grow.
                if (assemble_pair(i, j, fa, fb, dt, warm, nullptr))
                    refresh_conflict_ = true;
                else
                    ++step_statistics_.fixture_pairs;
            });
            // A constraint whose fixture pair was not visited has lost its
            // reason to exist: the array has to shrink.
            if (next != end)
                refresh_conflict_ = true;
            if (refresh_conflict_)
                return false;
        }
        for (; cursor < constraints_.size(); ++cursor)
            live_points += std::size_t(constraints_[cursor].count);
        points_live_ = live_points;
        return true;
    }
    // Whether `c` belongs to a candidate pair ordered before `(i, j)`. The
    // constraint array follows the candidate visit order, which is ascending in
    // the lower body slot and then the higher one.
    static bool pair_before(const Constraint &c, std::size_t i, std::size_t j) {
        const std::size_t a = c.a->slot, b = c.b->slot;
        return a < i || (a == i && b < j);
    }
    // Everything a contact build does once the constraint array is final,
    // whichever pass produced it.
    void finish_contact_build(Clock::time_point start, double cache_before, float dt, bool warm) {
        auto wake_start = Clock::now();
        wake_connected();
        double wake_ms = milliseconds(wake_start);
        step_statistics_.sleeping_ms += wake_ms;
        std::sort(solid_contacts_.begin(), solid_contacts_.end());
        active_constraints_.clear();
        for (std::size_t i = 0; i < constraints_.size(); ++i)
            if (constraints_[i].a->is_dynamic() || constraints_[i].b->is_dynamic())
                active_constraints_.push_back(i);
        step_statistics_.active_constraints = active_constraints_.size();
        if (warm)
            apply_warm_start();
        step_statistics_.constraints = constraints_.size();
        step_statistics_.detection_ms +=
            milliseconds(start) - wake_ms - (step_statistics_.cache_ms - cache_before);
        note_build(dt);
    }
    void build_constraints(float dt, bool warm) {
        auto start = Clock::now();
        double cache_before = step_statistics_.cache_ms;
        if (contacts_unchanged(dt)) {
            ++step_statistics_.unchanged_contact_steps;
            step_statistics_.constraints = constraints_.size();
            step_statistics_.active_constraints = active_constraints_.size();
            // The manifolds, filter answers and accumulated impulses are all
            // still right, but the warm start is a *per-step* application: the
            // solver expects gravity to have been absorbed before it runs, so
            // the retained impulses have to be pushed into the velocities here
            // exactly as a rebuild would have.
            if (warm)
                apply_warm_start();
            return;
        }
        contact_build_ran_ = true;
        auto t_ids = Clock::now();
        // During the second solve, impulses are already in velocities. Transfer
        // accumulators but do not apply them again.
        if (!warm)
            save_constraints(dt);
        step_statistics_.save_ms += milliseconds(t_ids);
        // Only the pairs that moved need re-detecting, and they can be rewritten
        // where they are. When the contact *set* has to change instead, this
        // declines and the rebuild below does the job as before.
        if (refresh_moved_contacts(dt, warm)) {
            ++step_statistics_.refreshed_contact_steps;
            finish_contact_build(start, cache_before, dt, warm);
            return;
        }
        constraints_.clear();
        points_.clear();
        solid_contacts_.clear();
        auto t_pairs = Clock::now();
        refresh_conflict_ = false;
        for_each_candidate_pair([&](std::size_t i, std::size_t j) {
            ++step_statistics_.candidate_pairs;
            visit_pairs(*bodies_[i], *bodies_[j], [&](Fixture &fa, Fixture &fb) {
                assemble_pair(i, j, fa, fb, dt, warm, nullptr);
            });
        });
        step_statistics_.pairs_ms += milliseconds(t_pairs);
        points_live_ = points_.size();
        finish_contact_build(start, cache_before, dt, warm);
    }
    // Re-apply the impulses accumulated by the previous solve to the velocities.
    // This is what lets a resting contact absorb gravity before CCD sees it, so
    // the resting fast path has to run it as well: skipping it would leave the
    // solver to rediscover the whole support impulse from scratch.
    void apply_warm_start() {
        for (auto i : active_constraints_) {
            auto &c = constraints_[i];
            if (!c.a->is_dynamic() && !c.b->is_dynamic())
                continue;
            const Rot rot_a(c.a->transform.angle), rot_b(c.b->transform.angle);
            for (int k = 0; k < c.count; ++k) {
                auto p = points_.at(std::size_t(c.point_begin + k));
                Vec2 ra = rotate(p.local_a, rot_a), rb = rotate(p.local_b, rot_b);
                if (p.normal_impulse > 0)
                    ++step_statistics_.warm_started_points;
                apply(c, ra, rb,
                      c.normal * p.normal_impulse +
                          Vec2{-c.normal.y, c.normal.x} * p.tangent_impulse);
            }
        }
    }
    // Sequential-impulse pass over every joint, run once per step before the
    // position integration. Springs settle in a single pass because their
    // softness is baked into the impulse they compute.
    void solve_joint_velocities(float dt) {
        if (joints_.empty())
            return;
        auto start = Clock::now();
        for (const SolverIsland &island : solver_islands_)
            for (int iteration = 0; iteration < config_.solver_iterations; ++iteration)
                for (std::size_t k = island.joint_begin; k < island.joint_end; ++k) {
                    auto &joint = *joints_[island_joints_[k]];
                    if (joint.spring_stiffness <= 0 || iteration == 0)
                        joint.solve(dt);
                }
        step_statistics_.velocity_ms += milliseconds(start);
    }
    // Islands of mutually touching awake bodies. Two islands share no body, so
    // solving them one after another in any order is equivalent to one global
    // Gauss-Seidel sweep: the result is unchanged, while the work becomes
    // contiguous and an island with nothing awake in it is skipped outright.
    struct SolverIsland {
        std::size_t constraint_begin{}, constraint_end{};
        std::size_t joint_begin{}, joint_end{};
    };
    void build_solver_islands() {
        solver_islands_.clear();
        island_constraints_.clear();
        island_joints_.clear();
        step_statistics_.solver_islands = 0;
        const std::size_t slots = bodies_.size();
        const std::size_t constraint_count = active_constraints_.size();
        const std::size_t joint_count = joints_.size();
        if (slots == 0 || (constraint_count == 0 && joint_count == 0))
            return;
        island_union_.resize(slots);
        for (std::size_t i = 0; i < slots; ++i)
            island_union_[i] = i;
        auto root = [&](std::size_t i) {
            while (island_union_[i] != i) {
                island_union_[i] = island_union_[island_union_[i]];
                i = island_union_[i];
            }
            return i;
        };
        // Union first, then label: the roots handed out during the union pass
        // are not final until every edge has been processed. Only awake dynamic
        // bodies join islands: a static or kinematic body taken no velocity from
        // the solver, so merging through one (a shared floor, say) would collapse
        // the whole world into a single island for no reason.
        for (std::size_t k = 0; k < constraint_count; ++k) {
            const Constraint &c = constraints_[active_constraints_[k]];
            if (c.a->is_dynamic() && c.b->is_dynamic())
                island_union_[root(c.a->slot)] = root(c.b->slot);
        }
        for (const auto &joint : joints_) {
            if (!joint->a || !joint->b)
                continue;
            if (joint->a->is_dynamic() && joint->b->is_dynamic())
                island_union_[root(joint->a->slot)] = root(joint->b->slot);
        }
        // Label every item by an awake dynamic endpoint. `active_constraints_`
        // only holds pairs with one, so this is always defined.
        auto island_slot = [](const Body &a, const Body &b) {
            return a.is_dynamic() ? a.slot : b.slot;
        };
        island_root_.resize(constraint_count + joint_count);
        for (std::size_t k = 0; k < constraint_count; ++k) {
            const Constraint &c = constraints_[active_constraints_[k]];
            island_root_[k] = root(island_slot(*c.a, *c.b));
        }
        for (std::size_t k = 0; k < joint_count; ++k) {
            const auto &joint = *joints_[k];
            island_root_[constraint_count + k] =
                joint.a && joint.b ? root(island_slot(*joint.a, *joint.b)) : 0;
        }
        // Counting sort: count, prefix sum, scatter. Two linear passes, no
        // per-island allocation and no pointer chasing.
        auto prefix = [&](std::size_t begin, std::size_t end, std::size_t offset,
                          std::vector<std::size_t> &counts, std::vector<std::size_t> &starts) {
            counts.assign(slots, 0);
            for (std::size_t k = begin; k < end; ++k)
                ++counts[island_root_[offset + k]];
            starts.assign(slots + 1, 0);
            for (std::size_t i = 0; i < slots; ++i)
                starts[i + 1] = starts[i] + counts[i];
        };
        prefix(0, constraint_count, 0, island_constraint_counts_, island_constraint_starts_);
        prefix(0, joint_count, constraint_count, island_joint_counts_, island_joint_starts_);
        island_constraints_.resize(constraint_count);
        island_cursor_.assign(island_constraint_starts_.begin(),
                             island_constraint_starts_.end());
        // Store constraint indices, not positions in the active list: the
        // solver then reaches a constraint and its cached geometry in one hop.
        for (std::size_t k = 0; k < constraint_count; ++k)
            island_constraints_[island_cursor_[island_root_[k]]++] = active_constraints_[k];
        island_joints_.resize(joint_count);
        island_cursor_.assign(island_joint_starts_.begin(), island_joint_starts_.end());
        for (std::size_t k = 0; k < joint_count; ++k) {
            if (!joints_[k]->a || !joints_[k]->b)
                continue;
            const std::size_t r = island_root_[constraint_count + k];
            island_joints_[island_cursor_[r]++] = k;
        }
        for (std::size_t i = 0; i < slots; ++i) {
            SolverIsland island;
            island.constraint_begin = island_constraint_starts_[i];
            island.constraint_end = island_constraint_starts_[i] + island_constraint_counts_[i];
            island.joint_begin = island_joint_starts_[i];
            island.joint_end = island_joint_starts_[i] + island_joint_counts_[i];
            if (island.constraint_begin == island.constraint_end &&
                island.joint_begin == island.joint_end)
                continue;
            // A constraint only reaches this list while one of its bodies is an
            // awake dynamic, so an island with a constraint is always awake.
            bool awake = island.constraint_end > island.constraint_begin;
            for (std::size_t k = island.joint_begin; !awake && k < island.joint_end; ++k) {
                const auto &joint = *joints_[island_joints_[k]];
                awake = joint.a->is_dynamic() || joint.b->is_dynamic();
            }
            if (!awake)
                continue;
            solver_islands_.push_back(island);
        }
        step_statistics_.solver_islands = solver_islands_.size();
    }
    void solve_velocities(float dt) {
        auto start = Clock::now();
        build_velocity_geometry(dt);
        for (const SolverIsland &island : solver_islands_)
            for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
                for (std::size_t k = island.constraint_begin; k < island.constraint_end; ++k)
                    solve_velocity(island_constraints_[k]);
                if (!config_.interleave_joints)
                    continue;
                // The island's own joints, in the same sweep. Springs fold their
                // whole response into solve(), so they run once per step.
                for (std::size_t k = island.joint_begin; k < island.joint_end; ++k) {
                    auto &joint = *joints_[island_joints_[k]];
                    if (joint.spring_stiffness <= 0 || iteration == 0)
                        joint.solve(dt);
                }
            }
        step_statistics_.velocity_ms += milliseconds(start);
    }
    void build_velocity_geometry(float dt) {
        // Transforms and mass properties are constant throughout this velocity pass.
        // Rebuild after integration; never reuse these values for position corrections.
        //
        // The cache is keyed by constraint index so the solver reaches a
        // constraint and its geometry in one hop. It is only grown while there
        // is something to solve: `constraints_` keeps every sleeping pair too,
        // and value-initialising that much storage on a resting step would cost
        // more than the whole solve.
        //
        // A positive `contact_stiffness` relaxes the normal constraint's
        // increment by `K / (K + gamma)` with `gamma = 1 / (dt^2 * stiffness)`.
        // It is not a spring: see the field's comment in `Config` for why this
        // engine cannot use one. The friction constraint is relaxed with the
        // same term so the ratio between them stays what the friction
        // coefficient says.
        const float gamma =
            config_.contact_stiffness > 0 && dt > 0
                ? 1.0f / (dt * dt * config_.contact_stiffness)
                : 0.0f;
        if (active_constraints_.empty()) {
            velocity_geometry_.clear();
            return;
        }
        velocity_geometry_.resize(constraints_.size());
        for (auto index : active_constraints_) {
            auto &c = constraints_[index];
            auto &g = velocity_geometry_[index];
            float ma = inv_mass(*c.a), mb = inv_mass(*c.b), ia = inv_inertia(*c.a),
                  ib = inv_inertia(*c.b);
            g.movable = ma + mb > 0;
            g.tangent = {-c.normal.y, c.normal.x};
            const Rot rot_a(c.a->transform.angle), rot_b(c.b->transform.angle);
            for (int k = 0; k < c.count; ++k) {
                const auto point = points_.at(std::size_t(c.point_begin + k));
                g.ra[k] = rotate(point.local_a, rot_a);
                g.rb[k] = rotate(point.local_b, rot_b);
                float na = g.ra[k].cross(c.normal), nb = g.rb[k].cross(c.normal),
                      sa = g.ra[k].cross(g.tangent), sb = g.rb[k].cross(g.tangent);
                g.normal_mass[k] = ma + mb + ia * na * na + ib * nb * nb + gamma;
                g.tangent_mass[k] = ma + mb + ia * sa * sa + ib * sb * sb + gamma;
            }
            if (c.count == 2) {
                float a0 = g.ra[0].cross(c.normal), a1 = g.ra[1].cross(c.normal),
                      b0 = g.rb[0].cross(c.normal), b1 = g.rb[1].cross(c.normal);
                g.k01 = ma + mb + ia * a0 * a1 + ib * b0 * b1;
                g.determinant = g.normal_mass[0] * g.normal_mass[1] - g.k01 * g.k01;
            }
        }
    }
    // One constraint, one solver iteration. `index` addresses `constraints_`
    // and the cached velocity geometry directly, so the island lists stay
    // plain integer arrays.
    void solve_velocity(std::size_t index) {
        auto &c = constraints_[index];
        const auto &g = velocity_geometry_[index];
        if (!g.movable)
            return;
        if (c.count == 2) {
            Vec2 ra0 = g.ra[0], rb0 = g.rb[0], ra1 = g.ra[1], rb1 = g.rb[1];
            float k00 = g.normal_mass[0], k11 = g.normal_mass[1], k01 = g.k01;
            float det = g.determinant;
            auto p0 = points_.at(std::size_t(c.point_begin));
            auto p1 = points_.at(std::size_t(c.point_begin + 1));
            if (det > 1e-8f) {
                float v0 =
                    p0.target - (point_velocity(*c.b, rb0) - point_velocity(*c.a, ra0)).dot(c.normal);
                float v1 =
                    p1.target - (point_velocity(*c.b, rb1) - point_velocity(*c.a, ra1)).dot(c.normal);
                float d0 = (k11 * v0 - k01 * v1) / det, d1 = (k00 * v1 - k01 * v0) / det;
                if (p0.normal_impulse + d0 >= 0 && p1.normal_impulse + d1 >= 0) {
                    p0.normal_impulse += d0;
                    p1.normal_impulse += d1;
                    apply(c, ra0, rb0, c.normal * d0);
                    apply(c, ra1, rb1, c.normal * d1);
                }
            }
        }
        for (int k = 0; k < c.count; ++k) {
            auto p = points_.at(std::size_t(c.point_begin + k));
            Vec2 ra = g.ra[k], rb = g.rb[k];
            Vec2 rv = point_velocity(*c.b, rb) - point_velocity(*c.a, ra);
            float denom = g.normal_mass[k];
            float previous = p.normal_impulse;
            p.normal_impulse =
                std::max(0.0f, previous + (p.target - rv.dot(c.normal)) / denom);
            apply(c, ra, rb, c.normal * (p.normal_impulse - previous));
            Vec2 tangent = g.tangent;
            rv = point_velocity(*c.b, rb) - point_velocity(*c.a, ra);
            previous = p.tangent_impulse;
            float limit = c.friction * p.normal_impulse;
            p.tangent_impulse =
                std::clamp(previous - rv.dot(tangent) / g.tangent_mass[k], -limit, limit);
            apply(c, ra, rb, tangent * (p.tangent_impulse - previous));
        }
    }
    void guard_projection(Body &body, Transform before, float radius = -1) {
        if (!config_.ccd.enabled || !body.is_dynamic())
            return;
        const auto after = body.transform;
        if ((after.position - before.position).length_squared() < 1e-16f &&
            std::abs(after.angle - before.angle) < 1e-8f)
            return;
        if (radius >= 0) {
            // Cheap full-rotation envelope first. Most stack corrections are
            // far from static geometry; avoid rebuilding their fixture AABBs.
            // A hit here still goes through the tighter bounds and exact sweep.
            AABB sweep{{std::min(before.position.x, after.position.x) - radius,
                        std::min(before.position.y, after.position.y) - radius},
                       {std::max(before.position.x, after.position.x) + radius,
                        std::max(before.position.y, after.position.y) + radius}};
            if (!obstacle_near(sweep)) {
                ++step_statistics_.projection_fast_rejections;
                return;
            }
        }
        float fraction = 1;
        // Bound the full correction path, not just its endpoints. Every point
        // rotates by at most radius * angle; the endpoint AABB bounds radius
        // about the body pivot even for offset and compound fixtures.
        auto bounds = body_aabb(body);
        Vec2 reach{std::max(std::abs(bounds.min.x - after.position.x),
                            std::abs(bounds.max.x - after.position.x)),
                   std::max(std::abs(bounds.min.y - after.position.y),
                            std::abs(bounds.max.y - after.position.y))};
        // Rotation pivots on the centroid, so the envelope about the body
        // origin has to include the centroid offset.
        reach += Vec2{std::abs(body.local_center.x), std::abs(body.local_center.y)};
        float padding = reach.length() * std::min(2.0f, std::abs(after.angle - before.angle));
        auto delta = before.position - after.position;
        bounds.min -= Vec2{padding, padding};
        bounds.max += Vec2{padding, padding};
        bounds.min.x += std::min(0.0f, delta.x);
        bounds.min.y += std::min(0.0f, delta.y);
        bounds.max.x += std::max(0.0f, delta.x);
        bounds.max.y += std::max(0.0f, delta.y);
        for_each_obstacle(bounds, [&](Body &obstacle) {
            ++step_statistics_.projection_candidates;
            visit_pairs(body, obstacle, [&](Fixture &a, Fixture &b) {
                if (a.trigger || b.trigger || !allowed(a, b))
                    return;
                ++step_statistics_.projection_sweeps;
                if (ccd_detail::supported(a.shape) && ccd_detail::supported(b.shape)) {
                    ShapeSweep projection{before, after, a.local};
                    auto bt = fixture_transform(b);
                    auto initial = ccd_detail::separation(a.shape, projection.at(0), b.shape, bt);
                    // An already overlapping pair is outside sweep_shapes' contract.
                    // Permit depenetration, but never deepen overlap during projection.
                    if (initial.distance <= config_.ccd.tolerance) {
                        ShapeSweep fixed{obstacle.transform, obstacle.transform, b.local};
                        float depth = std::max(config_.ccd.tolerance, -initial.distance);
                        if (!ccd_detail::separated_during_sweep(a.shape, projection, b.shape, fixed,
                                                                initial.normal, depth))
                            fraction = 0;
                        return;
                    }
                }
                auto hit =
                    sweep_shapes(a.shape, {before, after, a.local}, b.shape,
                                 {obstacle.transform, obstacle.transform, b.local}, config_.ccd);
                if (hit)
                    fraction = std::min(fraction, hit->fraction);
            });
        });
        if (fraction < 1) {
            body.transform = {before.position + (after.position - before.position) * fraction,
                              before.angle + (after.angle - before.angle) * fraction};
            ++step_statistics_.position_clamps;
        }
    }
    void solve_position(Constraint &c) {
        float ma = inv_mass(*c.a), mb = inv_mass(*c.b), ia = inv_inertia(*c.a),
              ib = inv_inertia(*c.b);
        if (ma + mb <= 0)
            return;
        const Rot rot_a(c.a->transform.angle), rot_b(c.b->transform.angle);
        const Vec2 center_a = c.a->center_of_mass(rot_a), center_b = c.b->center_of_mass(rot_b);
        for (int k = 0; k < c.count; ++k) {
            const auto p = points_.at(std::size_t(c.point_begin + k));
            // Lever arms are measured from the centroids, so the world contact
            // points have to be reconstructed from them.
            Vec2 ra = rotate(p.local_a, rot_a), rb = rotate(p.local_b, rot_b);
            float separation = (center_b + rb - center_a - ra).dot(c.normal);
            float correction =
                std::clamp(-separation - 0.001f, 0.0f, config_.max_position_correction);
            step_statistics_.max_penetration =
                std::max(step_statistics_.max_penetration, -separation);
            step_statistics_.max_correction = std::max(step_statistics_.max_correction, correction);
            if (correction == 0)
                continue;
            ++step_statistics_.position_corrections;
            auto before_a = c.a->transform, before_b = c.b->transform;
            float na = ra.cross(c.normal), nb = rb.cross(c.normal);
            const float relaxation = config_.solver_mode == SolverMode::PBD
                                         ? 1.0f
                                         : std::clamp(config_.position_relaxation, 0.0f, 1.0f);
            float impulse = correction * relaxation / (ma + mb + ia * na * na + ib * nb * nb);
            c.a->transform.position -= c.normal * (impulse * ma);
            c.b->transform.position += c.normal * (impulse * mb);
            c.a->rotate_by(-na * impulse * ia);
            c.b->rotate_by(nb * impulse * ib);
            guard_projection(*c.a, before_a, c.projection_radius_a);
            guard_projection(*c.b, before_b, c.projection_radius_b);
            if (on_contact_diagnostic)
                on_contact_diagnostic({c.a, c.b, c.a->velocity, c.b->velocity, c.normal,
                                       -separation, before_a, before_b, c.a->transform,
                                       c.b->transform});
        }
    }
    // Static floors do not merge unrelated islands. Contacts and joints connect
    // dynamic bodies; an island sleeps only when every member is quiet.
    std::vector<Body *> island_bodies_, previous_island_bodies_;
    std::vector<std::pair<Body *, Body *>> island_edges_, previous_island_edges_;
    std::vector<std::pair<Body *, std::size_t>> island_index_;
    std::vector<std::size_t> island_parent_;
    std::vector<std::vector<Body *>> island_groups_;
    const std::vector<std::vector<Body *>> &islands() {
        island_bodies_.clear();
        island_edges_.clear();
        for (auto &b : bodies_)
            if (b->type == BodyType::Dynamic)
                island_bodies_.push_back(b.get());
        auto edge = [&](Body *a, Body *b) {
            if (a->type == BodyType::Dynamic && b->type == BodyType::Dynamic)
                island_edges_.emplace_back(a, b);
        };
        for (auto &c : constraints_)
            edge(c.a, c.b);
        for (auto &j : joints_)
            edge(j->a, j->b);
        if (island_bodies_ == previous_island_bodies_ && island_edges_ == previous_island_edges_)
            return island_groups_;
        previous_island_bodies_ = island_bodies_;
        previous_island_edges_ = island_edges_;
        island_index_.clear();
        island_parent_.resize(island_bodies_.size());
        for (std::size_t i = 0; i < island_bodies_.size(); ++i) {
            island_index_.emplace_back(island_bodies_[i], i);
            island_parent_[i] = i;
        }
        auto less = [](const auto &a, const auto &b) {
            return std::less<Body *>{}(a.first, b.first);
        };
        std::sort(island_index_.begin(), island_index_.end(), less);
        auto index = [&](Body *b) {
            return std::lower_bound(island_index_.begin(), island_index_.end(),
                                    std::make_pair(b, std::size_t{}), less)
                ->second;
        };
        auto root = [&](std::size_t i) {
            while (island_parent_[i] != i) {
                island_parent_[i] = island_parent_[island_parent_[i]];
                i = island_parent_[i];
            }
            return i;
        };
        for (auto [a, b] : island_edges_)
            island_parent_[root(index(a))] = root(index(b));
        island_groups_.resize(island_bodies_.size());
        for (auto &g : island_groups_)
            g.clear();
        for (auto [b, i] : island_index_)
            island_groups_[root(i)].push_back(b);
        return island_groups_;
    }
    void wake_neighbors(Body &body) {
        if (body.type == BodyType::Dynamic)
            body.wake();
        auto wake = [&](Body *a, Body *b) {
            if (a != &body && b != &body)
                return;
            if (a->type == BodyType::Dynamic)
                a->wake();
            if (b->type == BodyType::Dynamic)
                b->wake();
        };
        for (auto &c : constraints_)
            wake(c.a, c.b);
        for (auto &j : joints_)
            wake(j->a, j->b);
        wake_connected();
    }
    void wake_connected() {
        for (auto &c : constraints_) {
            if (c.a->type == BodyType::Kinematic &&
                (c.a->velocity.length_squared() > 0 || c.a->angular_velocity != 0) && c.b->sleeping)
                c.b->wake();
            if (c.b->type == BodyType::Kinematic &&
                (c.b->velocity.length_squared() > 0 || c.b->angular_velocity != 0) && c.a->sleeping)
                c.a->wake();
        }
        for (auto &group : islands()) {
            bool awake = false;
            for (auto *b : group)
                awake |= !b->sleeping;
            if (awake)
                for (auto *b : group)
                    if (b->sleeping)
                        b->wake();
        }
    }
    void update_sleep() {
        for (auto &group : islands()) {
            if (group.empty())
                continue;
            ++step_statistics_.sleep_groups;
            bool quiet = true;
            int counter = 31;
            for (auto *b : group) {
                step_statistics_.max_speed =
                    std::max(step_statistics_.max_speed, b->velocity.length());
                step_statistics_.max_angular_speed =
                    std::max(step_statistics_.max_angular_speed, std::abs(b->angular_velocity));
                quiet &=
                    b->velocity.length_squared() < 0.0025f && std::abs(b->angular_velocity) < 0.05f;
                counter = std::min(counter, b->sleep_counter);
            }
            step_statistics_.moving_groups += !quiet;
            step_statistics_.settling_groups += quiet && counter < 31;
            for (auto *b : group) {
                b->sleep_counter = quiet ? counter + 1 : 0;
                if (b->sleep_counter > 30) {
                    b->sleeping = true;
                    b->velocity = {};
                    b->angular_velocity = 0;
                }
                step_statistics_.awake_bodies += !b->sleeping;
            }
        }
    }
    friend class BodyBuilder;
    Config config_;
    std::vector<std::unique_ptr<Body>> bodies_;
    // One generation per slot, bumped whenever the slot changes occupant. Slot
    // indices are reused, so this is what separates "body 7" from "body 7 after
    // body 7 was destroyed" -- a `BodyId` that outlives its body resolves to
    // nullptr instead of to whoever moved in.
    std::vector<std::uint64_t> generations_;
    std::uint64_t generation_counter_{0};
    std::vector<std::unique_ptr<DistanceJoint>> joints_;
    // Solver islands, rebuilt whenever the active constraint set changes. The
    // scratch arrays live here so a steady step allocates nothing.
    std::vector<SolverIsland> solver_islands_;
    std::vector<std::size_t> island_constraints_; // positions in active_constraints_
    std::vector<std::size_t> island_joints_;      // indices into joints_
    std::vector<std::size_t> island_root_, island_union_, island_cursor_;
    std::vector<std::size_t> island_constraint_counts_, island_constraint_starts_;
    std::vector<std::size_t> island_joint_counts_, island_joint_starts_;
    std::unordered_set<std::uint64_t> active_triggers_{};
    mutable std::size_t broadphase_candidate_count_{0};
    static std::uint64_t pair_key(std::size_t i, std::size_t j) {
        return (std::uint64_t(i) << 32) | std::uint64_t(j);
    }
    static std::pair<std::size_t, std::size_t> unpack_key(std::uint64_t key) {
        return {std::size_t(key >> 32), std::size_t(key & 0xffffffffu)};
    }
    void emit_trigger(Body &a, Body &b, bool enter) {
        CallbackLock lock(locked_);
        if (a.trigger)
            on_trigger(a, b, enter);
        else
            on_trigger(b, a, enter);
    }
    // -----------------------------------------------------------------------
    // Broad phase
    //
    // One proxy per fixture over a single dynamic AABB tree, plus a candidate
    // pair set that is maintained incrementally.
    //
    // The old version recomputed every body's AABB on every pass and, once any
    // single body moved, re-queried the tree for every dynamic body and then
    // sorted and deduplicated the whole candidate list. The work was therefore
    // global even when one object moved on a large map. Here a body's proxies
    // are only re-filed when its public state actually changed, and only the
    // pairs that touch a changed body are re-derived, so per-step cost follows
    // the number of moving bodies instead of the world size.
    //
    // Proxies are created per fixture rather than per body. A compound body
    // whose fixtures sit far apart used to file one enormous union box, which
    // produced candidates for everything in the gap between them.
    // -----------------------------------------------------------------------
    static constexpr std::uint32_t kNoFixture = 0xffffffffu;
    struct Proxy {
        std::uint32_t node{DynamicTree::null_node};
        std::uint32_t body{0};
        std::uint32_t fixture{kNoFixture}; // kNoFixture = the legacy body shape
        std::uint32_t next{kNoFixture};    // intrusive list of one body's proxies
        AABB tight{};
    };
    mutable std::vector<Proxy> proxies_;
    mutable DynamicTree tree_;
    mutable std::vector<std::uint32_t> body_first_proxy_;
    // Last state each body was filed under. These snapshots are what make a
    // direct public field write detectable without recomputing every AABB: the
    // world compares them per body and only pays the geometry work for the ones
    // that differ. A controlled setter (Body::set_position and friends) still
    // has to be seen here, so the compare is the documented cost of the
    // direct-write compatibility path; the version counters on Body/Fixture let
    // callers of the setter skip nothing but let the engine skip *derived* work.
    mutable std::vector<Transform> sync_transform_;
    mutable std::vector<BodyType> sync_type_;
    mutable std::vector<std::uint64_t> sync_hash_;
    mutable std::vector<char> sync_default_;
    mutable std::vector<char> sync_sleeping_;
    // One entry per fixture, or a single synthetic entry for the legacy
    // single-shape body. The same table serves the broad phase (content hash and
    // local offset), the position-correction obstacle index and the per-step
    // edit registry, so a body's public state is captured once.
    mutable std::vector<std::vector<FixtureState>> sync_fixtures_;
    mutable std::vector<AABB> body_bounds_; // fattened union bounds, by body slot
    // Persistent candidate pairs. `neighbors_[i]` is the sorted partner list of
    // body i; the candidate set is exactly { (i, j) : j in neighbors_[i], j > i }.
    // Keeping it sorted means the stream handed to the solver is in the same
    // order the old sorted flat list used, so trajectories do not change.
    mutable std::vector<std::vector<std::size_t>> neighbors_;
    mutable std::vector<std::size_t> edited_bodies_, fresh_neighbors_, pair_removed_,
        pair_added_;
    // Parallel to `edited_bodies_`: whether that body's proxies have to be
    // re-filed and whether its partner list has to be re-derived.
    mutable std::vector<char> edited_relist_, edited_refile_;
    mutable std::vector<std::uint32_t> pair_stamp_;
    mutable std::size_t pair_stamp_value_{0};
    // Same idea for "is this body in the list I am currently rebuilding".
    mutable std::vector<std::uint32_t> member_stamp_;
    mutable std::size_t member_stamp_value_{0};
    mutable std::vector<std::uint32_t> query_stamp_;
    mutable std::size_t query_stamp_value_{0};
    mutable bool broadphase_dirty_{true};
    mutable float broadphase_margin_{-1};
    // Bumped whenever the engine itself moves a body (integration, position
    // solve) and at the start of every step, i.e. whenever a public field write
    // may have happened. `sync_broadphase`'s O(bodies) snapshot scan is memoised
    // on it, so several callers inside one pose share a single scan.
    mutable std::uint64_t motion_epoch_{0};
    mutable std::uint64_t synced_motion_{~std::uint64_t(0)};
    // Whether the last scan saw a pair appear or disappear.
    mutable bool broadphase_changed_{true};
    // Shape of the world the current `constraints_` were derived from.
    mutable std::size_t built_body_count_{std::size_t(-1)};
    mutable std::size_t built_joint_count_{std::size_t(-1)};
    mutable bool built_with_filter_{false};
    mutable float built_dt_{-1};
    float fat_margin() const {
        if (config_.broadphase_fat_margin >= 0)
            return config_.broadphase_fat_margin;
        return std::max(0.002f, config_.ccd.enabled ? config_.ccd.tolerance * 2 : 0.0f);
    }
    static bool pair_allowed(const Body &a, const Body &b) {
        return a.type == BodyType::Dynamic || b.type == BodyType::Dynamic;
    }
    static void insert_sorted(std::vector<std::size_t> &list, std::size_t value) {
        list.insert(std::lower_bound(list.begin(), list.end(), value), value);
    }
    static void erase_sorted(std::vector<std::size_t> &list, std::size_t value) {
        const auto it = std::lower_bound(list.begin(), list.end(), value);
        if (it != list.end() && *it == value)
            list.erase(it);
    }
    // A body is re-filed only when one of the snapshots it was filed under no
    // longer matches. For the common single-shape body this is three float
    // compares and one shape hash, i.e. far cheaper than rebuilding its AABB and
    // re-deriving its candidate pairs.
    //
    // The result is split so the caller can tell mere motion from a real
    // geometry or type change. It is tempting to skip the partner-list refresh
    // for a body whose fat box still covers its tight box (move_proxy then
    // returns false), because the tree's coarse traversal cannot reach a
    // different set of leaves. That is *not* sound here: the final leaf test is
    // `fatten(tight, margin)`, and that predicate is exactly the speculative
    // contact condition used by build_constraints. Two neighbours whose gap is
    // near the margin therefore gain and lose candidates while both fat boxes
    // stand still. Only the tight boxes may drive membership, so any body that
    // moved has to re-derive its list.
    struct EditFlags {
        bool moved{};       // transform no longer matches the filed one
        bool geometry{};    // outline, local offset or fixture set changed
        bool type_changed{};// pair_allowed() may now answer differently
        // The sleeping flag is derived state that contact management acts on:
        // is_dynamic() reads it, so `active_constraints_`, island propagation
        // and the velocity solve all change when it flips even though no
        // geometry moved.
        bool activity{};
        bool refile() const { return moved || geometry; }
        bool relist() const { return moved || geometry || type_changed; }
        bool any() const { return refile() || type_changed || activity; }
    };
    EditFlags classify_edit(std::size_t index) const {
        EditFlags flags;
        const Body &body = *bodies_[index];
        flags.type_changed = body.type != sync_type_[index];
        flags.activity = body.sleeping != (sync_sleeping_[index] != 0);
        flags.moved = body.transform.position != sync_transform_[index].position ||
                      body.transform.angle != sync_transform_[index].angle;
        const auto &snapshot = sync_fixtures_[index];
        if (sync_default_[index] != char(body.use_default_shape)) {
            flags.geometry = true;
            return flags;
        }
        if (body.use_default_shape) {
            flags.geometry = snapshot.size() != 1 || snapshot.front().hash != body_hash(body);
            if (flags.geometry)
                invalidate_shape_index(body.shape);
            return flags;
        }
        if (snapshot.size() != body.fixtures.size()) {
            flags.geometry = true;
            return flags;
        }
        for (std::size_t k = 0; k < snapshot.size(); ++k) {
            const Fixture &fixture = *body.fixtures[k];
            if (snapshot[k].hash != fixture_hash(fixture) ||
                snapshot[k].local.position != fixture.local.position ||
                snapshot[k].local.angle != fixture.local.angle) {
                invalidate_shape_index(fixture.shape);
                flags.geometry = true;
                return flags;
            }
        }
        return flags;
    }
    // A mesh's spatial index is derived from its triangles and lives in the
    // shape, so a direct write to those triangles cannot be seen from inside the
    // shape. It is seen from here instead: the same comparison that notices an
    // edited outline -- a hash that no longer matches the snapshot taken at the
    // end of the previous step -- is what drops the index, so the rebuild
    // happens before the edit can be queried. This is the same "noticed on the
    // next step" cost as every other direct write.
    static void invalidate_shape_index(const Shape &shape) {
        if (const auto *mesh = std::get_if<Mesh>(&shape))
            mesh->invalidate_acceleration();
    }
    void snapshot_body(std::size_t index) const {
        const Body &body = *bodies_[index];
        sync_transform_[index] = body.transform;
        sync_type_[index] = body.type;
        sync_sleeping_[index] = char(body.sleeping);
        sync_default_[index] = char(body.use_default_shape);
        sync_hash_[index] = body.use_default_shape ? body_hash(body) : 0;
        auto &snapshot = sync_fixtures_[index];
        if (body.use_default_shape) {
            if (snapshot.empty())
                snapshot.emplace_back();
            snapshot.resize(1);
            snapshot[0] = state_of(body, body_hash(body));
            return;
        }
        snapshot.resize(body.fixtures.size());
        for (std::size_t k = 0; k < body.fixtures.size(); ++k)
            snapshot[k] = state_of(body, *body.fixtures[k], fixture_hash(*body.fixtures[k]));
    }
    void proxy_geometry(std::size_t index, const Proxy &proxy, const Shape *&shape,
                        Transform &transform) const {
        const Body &body = *bodies_[index];
        if (proxy.fixture == kNoFixture) {
            shape = &body.shape;
            transform = body.transform;
            return;
        }
        const Fixture &fixture = *body.fixtures[proxy.fixture];
        shape = &fixture.shape;
        transform = fixture_transform(fixture);
    }
    AABB union_bounds(std::size_t index, float margin) const {
        AABB result{};
        bool first = true;
        for (std::uint32_t p = body_first_proxy_[index]; p != kNoFixture; p = proxies_[p].next) {
            const AABB &tight = proxies_[p].tight;
            result = first ? tight : merged(result, tight);
            first = false;
        }
        if (first) {
            // No geometry at all: a degenerate box at the body origin, which is
            // what the previous per-body union returned as well.
            const Vec2 position = bodies_[index]->transform.position;
            return {position, position};
        }
        return fatten(result, margin);
    }
    void create_body_proxies(std::size_t index, float margin) const {
        const Body &body = *bodies_[index];
        body_first_proxy_[index] = kNoFixture;
        auto make = [&](const Shape &shape, const Transform &t, std::uint32_t fixture) {
            Proxy proxy;
            proxy.tight = compute_aabb(shape, t);
            proxy.body = std::uint32_t(index);
            proxy.fixture = fixture;
            proxy.next = body_first_proxy_[index];
            proxy.node = tree_.create_proxy(proxy.tight, std::uint32_t(proxies_.size()), margin);
            body_first_proxy_[index] = std::uint32_t(proxies_.size());
            proxies_.push_back(proxy);
        };
        if (body.use_default_shape)
            make(body.shape, body.transform, kNoFixture);
        else
            for (std::size_t k = 0; k < body.fixtures.size(); ++k)
                make(body.fixtures[k]->shape, fixture_transform(*body.fixtures[k]),
                     std::uint32_t(k));
        snapshot_body(index);
        body_bounds_[index] = union_bounds(index, margin);
    }
    // Re-file the proxies of one edited body. `move_proxy` keeps the fat box
    // when the tight box still fits inside it, so a slowly moving body is not
    // re-inserted every step. Returns true when at least one proxy was.
    bool refresh_body_proxies(std::size_t index, float margin) const {
        bool refiled = false;
        for (std::uint32_t p = body_first_proxy_[index]; p != kNoFixture; p = proxies_[p].next) {
            Proxy &proxy = proxies_[p];
            const Shape *shape = nullptr;
            Transform transform{};
            proxy_geometry(index, proxy, shape, transform);
            const AABB tight = compute_aabb(*shape, transform);
            const Vec2 displacement =
                (tight.min + tight.max) * 0.5f - (proxy.tight.min + proxy.tight.max) * 0.5f;
            proxy.tight = tight;
            refiled |= tree_.move_proxy(proxy.node, tight, displacement, margin);
        }
        snapshot_body(index);
        body_bounds_[index] = union_bounds(index, margin);
        return refiled;
    }
    bool proxies_overlap(std::size_t a, std::size_t b, float margin) const {
        for (std::uint32_t p = body_first_proxy_[a]; p != kNoFixture; p = proxies_[p].next) {
            const AABB box = fatten(proxies_[p].tight, margin);
            for (std::uint32_t q = body_first_proxy_[b]; q != kNoFixture; q = proxies_[q].next)
                if (box.overlaps(fatten(proxies_[q].tight, margin)))
                    return true;
        }
        return false;
    }
    // Re-derive one body's candidate partners and update the persistent pair set
    // by the difference. Only partner changes touch the pair lists, so a body
    // that keeps the same neighbours costs one tree query and one pass over the
    // two small lists.
    //
    // The difference used to be `std::sort(fresh)` followed by two
    // `std::set_difference` passes. The tree hands partners back in arbitrary
    // order and the stamp above has already made `fresh` duplicate-free, so the
    // sort existed only to feed set_difference. Membership is now decided by two
    // O(1) stamp probes and the body's own list keeps its existing sorted order,
    // patched by the handful of insertions and erasures that actually happened.
    bool refresh_pairs(std::size_t index, float margin) const {
        const std::size_t body_count = bodies_.size();
        if (pair_stamp_.size() != body_count)
            pair_stamp_.assign(body_count, 0);
        if (member_stamp_.size() != body_count)
            member_stamp_.assign(body_count, 0);
        auto &fresh = fresh_neighbors_;
        fresh.clear();
        ++pair_stamp_value_;
        const std::uint32_t fresh_stamp = std::uint32_t(pair_stamp_value_);
        auto q_start = Clock::now();
        if (config_.enable_broadphase) {
            for (std::uint32_t p = body_first_proxy_[index]; p != kNoFixture;
                 p = proxies_[p].next) {
                tree_.query(fatten(proxies_[p].tight, margin), margin,
                            [&](std::uint32_t payload) {
                                const std::size_t other = proxies_[payload].body;
                                if (other == index || pair_stamp_[other] == fresh_stamp ||
                                    !pair_allowed(*bodies_[index], *bodies_[other]))
                                    return;
                                pair_stamp_[other] = fresh_stamp;
                                fresh.push_back(other);
                            });
            }
        } else {
            for (std::size_t other = 0; other < body_count; ++other) {
                if (other != index && pair_allowed(*bodies_[index], *bodies_[other]) &&
                    proxies_overlap(index, other, margin)) {
                    pair_stamp_[other] = fresh_stamp;
                    fresh.push_back(other);
                }
            }
        }
        (void)q_start;
        auto &current = neighbors_[index];
        ++member_stamp_value_;
        const std::uint32_t member_stamp = std::uint32_t(member_stamp_value_);
        for (std::size_t other : current)
            member_stamp_[other] = member_stamp;
        auto &removed = pair_removed_;
        auto &added = pair_added_;
        removed.clear();
        added.clear();
        for (std::size_t other : current)
            if (pair_stamp_[other] != fresh_stamp)
                removed.push_back(other);
        for (std::size_t other : fresh)
            if (member_stamp_[other] != member_stamp)
                added.push_back(other);
        for (std::size_t other : removed) {
            erase_sorted(neighbors_[other], index);
            --broadphase_candidate_count_;
        }
        for (std::size_t other : added) {
            insert_sorted(neighbors_[other], index);
            ++broadphase_candidate_count_;
        }
        for (std::size_t other : removed)
            erase_sorted(current, other);
        for (std::size_t other : added)
            insert_sorted(current, other);
        return !removed.empty() || !added.empty();
    }
    void rebuild_broadphase(float margin) const {
        tree_.reset();
        proxies_.clear();
        const std::size_t count = bodies_.size();
        body_first_proxy_.assign(count, kNoFixture);
        sync_transform_.resize(count);
        sync_type_.resize(count);
        sync_hash_.resize(count);
        sync_default_.assign(count, 0);
        sync_sleeping_.assign(count, 0);
        sync_fixtures_.resize(count);
        body_bounds_.assign(count, AABB{});
        pair_stamp_.assign(count, 0);
        pair_stamp_value_ = 0;
        member_stamp_.assign(count, 0);
        member_stamp_value_ = 0;
        query_stamp_.assign(count, 0);
        query_stamp_value_ = 0;
        neighbors_.resize(count);
        for (auto &list : neighbors_)
            list.clear();
        broadphase_candidate_count_ = 0;
        proxies_.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            create_body_proxies(i, margin);
        for (std::size_t i = 0; i < count; ++i)
            refresh_pairs(i, margin);
        broadphase_dirty_ = false;
    }
    // Refresh the proxy set to match the current body state. Returns true when
    // the candidate pair set changed, which is the only case in which the
    // constraint set can differ from the previous step.
    bool sync_broadphase(bool force = false) const {
        const float margin = fat_margin();
        // A changed margin invalidates the padding stored in every fat box.
        if (margin != broadphase_margin_) {
            broadphase_margin_ = margin;
            broadphase_dirty_ = true;
        }
        const std::size_t count = bodies_.size();
        if (broadphase_dirty_ || body_first_proxy_.size() != count) {
            rebuild_broadphase(margin);
            // Everything was re-filed, so no body can be reported as unchanged.
            edited_bodies_.clear();
            edited_relist_.clear();
            edited_refile_.clear();
            // Body slots may have moved as well, so every cached per-body view
            // that is indexed by slot is stale -- the CCD collider lists
            // included, since they hold pointers to Fixtures.
            ccd_view_invalid_ = true;
            synced_motion_ = motion_epoch_;
            broadphase_changed_ = true;
            return true;
        }
        // The scan below is O(bodies). Nothing moves a body without bumping the
        // motion epoch, so within one pose the answer cannot have changed and a
        // second caller can reuse it instead of walking every body again.
        if (!force && synced_motion_ == motion_epoch_)
            return broadphase_changed_;
        synced_motion_ = motion_epoch_;
        edited_bodies_.clear();
        edited_relist_.clear();
        edited_refile_.clear();
        for (std::size_t i = 0; i < count; ++i) {
            const EditFlags flags = classify_edit(i);
            if (!flags.any())
                continue;
            edited_bodies_.push_back(i);
            // Membership follows the tight boxes, so every body whose geometry
            // moved re-derives its partner list; a pure type flip does too, and
            // a pure sleep flip does not.
            edited_relist_.push_back(char(flags.relist()));
            edited_refile_.push_back(char(flags.refile()));
            sync_sleeping_[i] = char(bodies_[i]->sleeping);
        }
        if (edited_bodies_.empty()) {
            broadphase_changed_ = false;
            return false;
        }
        // Two passes: re-file every edited body first, then re-derive its pairs.
        // A single pass would derive one body's partners from another edited
        // body's stale box.
        for (std::size_t k = 0; k < edited_bodies_.size(); ++k)
            if (edited_refile_[k] && refresh_body_proxies(edited_bodies_[k], margin))
                edited_relist_[k] = 1;
        bool changed = false;
        auto relist_start = Clock::now();
        for (std::size_t k = 0; k < edited_bodies_.size(); ++k)
            if (edited_relist_[k])
                changed |= refresh_pairs(edited_bodies_[k], margin);
        step_statistics_.broadphase_ms += milliseconds(relist_start);
        broadphase_changed_ = changed;
        return changed;
    }
    // Visit the persistent candidate pairs in the same order the old sorted
    // flat list used: ascending i, then ascending j.
    template <class Visitor> void for_each_candidate_pair(Visitor &&visit) const {
        sync_broadphase();
        for (std::size_t i = 0; i < neighbors_.size(); ++i)
            for (std::size_t j : neighbors_[i])
                if (j > i)
                    visit(i, j);
    }
    // Position correction needs the static and kinematic obstacles near a body
    // that is being pushed out of an overlap. That set used to be rebuilt and
    // sorted every step, then scanned linearly with an early break, so a scene
    // with thousands of walls, platforms and terrain blocks paid for a box pile
    // that only touches a handful of them.
    //
    // The replacement is a dedicated index holding only the obstacles. It has
    // to be separate from the broad-phase tree: that one also contains every
    // dynamic body, so a query around a stack would walk all of them before
    // reaching the floor. Proxies are created once per obstacle fixture; a
    // static body is therefore free after the first step. Kinematics are
    // re-filed by sync_obstacles() right before the position solve, i.e. after
    // integration has moved them, so an expired bound is never reused.
    struct ObstacleProxy {
        std::uint32_t node{DynamicTree::null_node};
        std::uint32_t next{kNoFixture};
        std::uint32_t body{0};
        AABB tight{};
    };
    mutable DynamicTree obstacle_tree_;
    mutable std::vector<ObstacleProxy> obstacle_proxies_;
    mutable std::vector<std::uint32_t> obstacle_first_proxy_;
    mutable std::vector<Transform> obstacle_transform_;
    mutable std::vector<BodyType> obstacle_types_;
    mutable std::vector<std::vector<FixtureState>> obstacle_fixtures_;
    mutable std::vector<std::size_t> obstacle_stamp_;
    mutable std::size_t obstacle_stamp_value_{0};
    // Obstacle bounds depend only on the outline, its local offset and the body
    // transform, so the broad phase's already-memoised geometry snapshots are
    // enough to tell whether an obstacle has to be re-filed.
    bool obstacle_state_changed(std::size_t index) const {
        const auto &current = sync_fixtures_[index];
        const auto &saved = obstacle_fixtures_[index];
        if (current.size() != saved.size())
            return true;
        for (std::size_t k = 0; k < current.size(); ++k)
            if (current[k].hash != saved[k].hash ||
                current[k].local.position != saved[k].local.position ||
                current[k].local.angle != saved[k].local.angle)
                return true;
        return false;
    }
    void rebuild_obstacles() const {
        obstacle_tree_.reset();
        obstacle_proxies_.clear();
        const std::size_t count = bodies_.size();
        obstacle_first_proxy_.assign(count, kNoFixture);
        obstacle_transform_.resize(count);
        obstacle_types_.resize(count);
        obstacle_fixtures_.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            const Body &body = *bodies_[i];
            obstacle_types_[i] = body.type;
            obstacle_transform_[i] = body.transform;
            obstacle_fixtures_[i] = sync_fixtures_[i];
            obstacle_first_proxy_[i] = kNoFixture;
            if (body.type == BodyType::Dynamic)
                continue;
            auto make = [&](const Shape &shape, const Transform &t) {
                ObstacleProxy proxy;
                proxy.tight = compute_aabb(shape, t);
                proxy.body = std::uint32_t(i);
                proxy.next = obstacle_first_proxy_[i];
                proxy.node = obstacle_tree_.create_proxy(
                    proxy.tight, std::uint32_t(obstacle_proxies_.size()), 0.0f);
                obstacle_first_proxy_[i] = std::uint32_t(obstacle_proxies_.size());
                obstacle_proxies_.push_back(proxy);
            };
            if (body.use_default_shape)
                make(body.shape, body.transform);
            else
                for (const auto &fixture : body.fixtures)
                    make(fixture->shape, fixture_transform(*fixture));
        }
    }
    void refresh_obstacle_body(std::size_t index) const {
        const Body &body = *bodies_[index];
        const bool was_obstacle = obstacle_types_[index] != BodyType::Dynamic;
        obstacle_types_[index] = body.type;
        obstacle_transform_[index] = body.transform;
        obstacle_fixtures_[index] = sync_fixtures_[index];
        if (!was_obstacle)
            return;
        auto refresh = [&](std::uint32_t p, const Shape &shape, const Transform &t) {
            ObstacleProxy &proxy = obstacle_proxies_[p];
            const AABB tight = compute_aabb(shape, t);
            const Vec2 displacement =
                (tight.min + tight.max) * 0.5f - (proxy.tight.min + proxy.tight.max) * 0.5f;
            proxy.tight = tight;
            obstacle_tree_.move_proxy(proxy.node, tight, displacement, 0.0f);
        };
        std::uint32_t p = obstacle_first_proxy_[index];
        if (body.use_default_shape) {
            if (p != kNoFixture)
                refresh(p, body.shape, body.transform);
            return;
        }
        for (const auto &fixture : body.fixtures) {
            refresh(p, fixture->shape, fixture_transform(*fixture));
            p = obstacle_proxies_[p].next;
        }
    }
    // Bring the obstacle index up to the current transforms. Runs once per step,
    // after integration and before position solving, so a kinematic platform is
    // filed at the pose the position solver actually sees.
    void sync_obstacles() const {
        const std::size_t count = bodies_.size();
        bool structural = obstacle_types_.size() != count;
        if (!structural)
            for (std::size_t i = 0; i < count; ++i)
                if (obstacle_types_[i] != bodies_[i]->type) {
                    structural = true;
                    break;
                }
        if (structural) {
            rebuild_obstacles();
            return;
        }
        for (std::size_t i = 0; i < count; ++i) {
            if (obstacle_types_[i] == BodyType::Dynamic)
                continue;
            if (same_transform(bodies_[i]->transform, obstacle_transform_[i]) &&
                !obstacle_state_changed(i))
                continue;
            refresh_obstacle_body(i);
        }
    }
    // Visit every non-dynamic body whose fixture bounds can overlap `area`.
    // Candidates are confirmed exactly by the narrow phase below; the index only
    // replaces the linear scan over the whole world.
    template <class Visitor> void for_each_obstacle(const AABB &area, Visitor &&visit) {
        ++obstacle_stamp_value_;
        const std::size_t stamp = obstacle_stamp_value_;
        if (obstacle_stamp_.size() != bodies_.size())
            obstacle_stamp_.assign(bodies_.size(), 0);
        obstacle_tree_.query(area, 0.0f, [&](std::uint32_t payload) {
            const std::size_t index = obstacle_proxies_[payload].body;
            if (index >= bodies_.size() || obstacle_stamp_[index] == stamp)
                return; // one body can own several fixture proxies
            obstacle_stamp_[index] = stamp;
            visit(*bodies_[index]);
        });
    }
    // Whether any static or kinematic obstacle overlaps `area`. Stops at the
    // first hit instead of materialising the candidate list.
    bool obstacle_near(const AABB &area) const {
        bool found = false;
        obstacle_tree_.query_any(area, 0.0f, [&](std::uint32_t) {
            found = true;
            return false;
        });
        return found;
    }
    Body &add_body(std::unique_ptr<Body> body) {
        require_unlocked();
        const std::size_t index = bodies_.size();
        Body &ref = *body;
        ref.slot = std::uint32_t(index);
        bodies_.push_back(std::move(body));
        // A slot is reused as soon as a body is destroyed, so the generation is
        // what makes a handle an identity: it is bumped every time a slot
        // changes occupant, and `alive()` compares it.
        generations_.push_back(++generation_counter_);
        // Appending needs an entry in every per-slot view, but only this body's
        // proxies and candidate list have to be derived: a world that creates
        // and destroys bodies pays for what it added, not for what exists.
        if (broadphase_dirty_ || body_first_proxy_.size() != index) {
            broadphase_dirty_ = true; // Nothing filed yet; the first scan builds all of it.
            return ref;
        }
        grow_slots(index + 1);
        const float margin = fat_margin();
        create_body_proxies(index, margin);
        refresh_pairs(index, margin);
        return ref;
    }
};

inline Body &BodyBuilder::build() {
    auto body = std::make_unique<Body>();
    body->type = type_;
    body->transform = {position_, angle_};
    body->velocity = velocity_;
    body->mass = type_ != BodyType::Dynamic ? 0.0f : mass_;
    body->inverse_mass = body->mass > 0 ? 1.0f / body->mass : 0.0f;
    body->shape = std::move(shape_);
    body->material = material_;
    body->trigger = trigger_;
    body->bullet = bullet_;
    body->angular_velocity = angular_velocity_;
    body->collision_group = group_;
    body->collision_mask = mask_;
    if (body->type != BodyType::Dynamic) {
        body->inertia = body->inverse_inertia = 0;
        body->local_center = {};
        body->auto_mass = false;
        return world_.add_body(std::move(body));
    }
    // Inertia and the centroid come from the actual outline. A shared
    // "point mass at the origin" inertia made rotation depend only on mass:
    // a small circle and a large square span up identically.
    MassData data = mass_data(body->shape, density_ > 0 ? density_ : 1.0f);
    if (density_ > 0 && data.mass > 0) {
        body->mass = data.mass;
        body->inverse_mass = 1.0f / body->mass;
    } else if (data.mass > 0) {
        // Keep the requested mass and scale the outline's inertia with it.
        data.inertia *= body->mass / data.mass;
    }
    body->inertia = data.inertia;
    body->inverse_inertia = body->inertia > 0 ? 1.0f / body->inertia : 0.0f;
    body->local_center = data.center;
    body->auto_mass = false; // Builder mass properties are explicit.
    return world_.add_body(std::move(body));
}

} // namespace butter::physics2d
