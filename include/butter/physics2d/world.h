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
};

struct Body {
    BodyType type{BodyType::Dynamic};
    Transform transform{};
    // Index of this body inside its World. Engine-managed; kept current so the
    // solver can map a body back to a slot in O(1) while building islands.
    std::uint32_t slot{0};
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
        if (!a || !b)
            return;
        // One sine/cosine pair per body per call: this runs once per solver
        // iteration for every joint.
        const Rot rot_a(a->transform.angle), rot_b(b->transform.angle);
        const Vec2 anchor_world_a = a->transform.position + rotate(anchor_a, rot_a);
        const Vec2 anchor_world_b = b->transform.position + rotate(anchor_b, rot_b);
        // Lever arms are measured from each centroid, not from the body origin.
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
        const Vec2 va = a->velocity + Vec2{-a->angular_velocity * ra.y, a->angular_velocity * ra.x};
        const Vec2 vb = b->velocity + Vec2{-b->angular_velocity * rb.y, b->angular_velocity * rb.x};
        float impulse = -(vb - va).dot(n) / inv;
        if (spring_stiffness > 0 && dt > 0) {
            const float gamma = 1.0f / (dt * (damping + dt * spring_stiffness));
            impulse =
                -((vb - va).dot(n) + (d - length) * dt * spring_stiffness * gamma) / (inv + gamma);
        }
        if (a->is_dynamic()) {
            a->velocity -= n * impulse * a->inverse_mass;
            a->angular_velocity -= ca * impulse * a->inverse_inertia;
        }
        if (b->is_dynamic()) {
            b->velocity += n * impulse * b->inverse_mass;
            b->angular_velocity += cb * impulse * b->inverse_inertia;
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
    }
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
        if (!a || !b)
            return;
        // A feature that is switched off must not hand its stale impulse to the
        // warm start, or the body takes one last kick after the user disabled it.
        if (max_motor_torque <= 0)
            motor_impulse = 0;
        if (!enable_limit)
            lower_impulse = upper_impulse = 0;
        const Geometry g = geometry();
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
        if (!a || !b)
            return;
        const Geometry g = geometry();
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
    };
    Geometry geometry() const {
        Geometry g;
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
        // Box2D-style ResetMassData: density-weighted mass, centroid and
        // inertia of the whole compound outline.
        if (body.auto_mass)
            body.reset_mass_data();
        return result;
    }
    void destroy_fixture(Fixture &fixture) {
        require_unlocked();
        wake_neighbors(*fixture.body);
        contact_cache_.clear();
        constraints_.clear();
        forget_contacts(&fixture);
        auto &fixtures = fixture.body->fixtures;
        std::erase_if(fixtures, [&](const auto &p) { return p.get() == &fixture; });
        if (fixture.body->auto_mass)
            fixture.body->reset_mass_data();
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
        contact_cache_.clear();
        constraints_.clear();
        for (auto &fixture : body.fixtures)
            forget_contacts(fixture.get());
        std::erase_if(joints_, [&](const auto &p) { return p->a == &body || p->b == &body; });
        // Legacy trigger keys are indices. Remap them before erasing a body.
        auto it = std::find_if(bodies_.begin(), bodies_.end(),
                               [&](const auto &p) { return p.get() == &body; });
        if (it == bodies_.end())
            return;
        const auto index = std::size_t(it - bodies_.begin());
        std::unordered_set<std::uint64_t> remapped;
        for (auto key : active_triggers_) {
            auto [i, j] = unpack_key(key);
            if (i == index || j == index) {
                if (on_trigger)
                    emit_trigger(*bodies_[i], *bodies_[j], false);
            } else
                remapped.insert(pair_key(i - (i > index), j - (j > index)));
        }
        active_triggers_ = std::move(remapped);
        bodies_.erase(it);
        // Erasing shifts every index after it, so refresh the slots the solver
        // reads. Destroying a body is rare enough that a linear fix-up here is
        // still cheaper than hashing on the hot path.
        for (std::size_t i = index; i < bodies_.size(); ++i)
            bodies_[i]->slot = std::uint32_t(i);
        broadphase_dirty_ = true; // Erasing shifts every index after it.
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
        auto cache_start = Clock::now();
        invalidate_edited_contacts();
        step_statistics_.cache_ms += milliseconds(cache_start);
        for (auto &p : bodies_) {
            if (p->sleeping && p->type == BodyType::Dynamic &&
                (p->force.length_squared() > 0 || p->torque != 0 ||
                 p->velocity.length_squared() > 0 || p->angular_velocity != 0))
                p->wake();
            if (p->is_dynamic()) {
                p->velocity += (config_.gravity + p->force * p->inverse_mass) * dt;
                p->angular_velocity += p->torque * p->inverse_inertia * dt;
                p->velocity *= std::max(0.0f, 1.0f - p->linear_damping * dt);
                p->angular_velocity *= std::max(0.0f, 1.0f - p->angular_damping * dt);
            }
            if (p->fixed_rotation)
                p->angular_velocity = 0;
            p->force = {};
            p->torque = 0;
        }
        // Resting contacts absorb gravity before CCD sees them as new impacts.
        // Joints warm start alongside them: applying last step's accumulated
        // impulse up front is what turns a hinge chain from soft to rigid.
        for (auto &joint : joints_)
            joint->warm_start();
        build_constraints(dt, true);
        build_solver_islands();
        solve_velocities();
        // Joint velocities are solved *before* positions are integrated, so the
        // integration never uses a velocity the joint already rejects. Doing it
        // afterwards (the old order) made the position pass do the velocity
        // pass's job and pumped energy into every swing.
        solve_joint_velocities(dt);
        // The first detection still validates public edits and filtering every
        // step. With no awake dynamics or moving kinematics, no geometry can
        // change during integration: retain its constraints and event candidates.
        bool stationary = std::all_of(bodies_.begin(), bodies_.end(), [](const auto &b) {
            return (b->type != BodyType::Dynamic || b->sleeping) &&
                   (b->type != BodyType::Kinematic ||
                    (b->velocity.length_squared() == 0 && b->angular_velocity == 0));
        });
        if (stationary) {
            ccd_statistics_ = {};
            ccd_statistics_.advanced_time = dt;
            ++step_statistics_.stationary_steps;
        } else {
            auto ccd_start = Clock::now();
            step_continuous(dt);
            step_statistics_.ccd_ms += milliseconds(ccd_start);
            build_constraints(dt, false);
            build_solver_islands();
            solve_velocities();
        }
        auto event_start = Clock::now();
        const auto &pairs = candidate_pairs(false);
        broadphase_candidate_count_ = pairs.size();
        current_triggers_.clear();
        current_contacts_.clear();
        for (const auto [i, j] : pairs) {
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
        }
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
        auto sleep_start = Clock::now();
        update_sleep();
        step_statistics_.sleeping_ms += milliseconds(sleep_start);
        save_constraints(dt);
        simulation_time_ += dt;
    }
    void set_solver_iterations(int iterations) {
        config_.solver_iterations = std::max(1, iterations);
    }
    std::function<bool(const Fixture &, const Fixture &)> contact_filter;
    std::function<void(Fixture &, Fixture &, bool)> on_contact;
    const CcdStatistics &ccd_statistics() const { return ccd_statistics_; }
    std::size_t body_count() const { return bodies_.size(); }
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
        sync_broadphase();
        tree_.query(area, 0.0f, [&](std::uint32_t payload) {
            const auto index = std::size_t(payload);
            if (index < bodies_.size() && body_aabb(*bodies_[index]).overlaps(area))
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
        sync_broadphase();
        tree_.raycast(origin, direction, max_distance, [&](std::uint32_t payload, float entry) {
            // Everything behind the closest exact hit can be pruned.
            if (best && entry > best->distance)
                return true;
            const auto index = std::size_t(payload);
            if (index >= bodies_.size())
                return true;
            const auto &body = *bodies_[index];
            if (body.use_default_shape)
                consider(index, body.shape, body.transform);
            else
                for (const auto &fixture : body.fixtures)
                    consider(index, fixture->shape, fixture_transform(*fixture));
            return true;
        });
        return best;
    }
    std::function<void(Body &, Body &, bool)> on_trigger;

  private:
    CcdStatistics ccd_statistics_{};
    CcdWorkspace ccd_workspace_;
    std::vector<CcdMotion> ccd_motions_;
    std::vector<Fixture> ccd_legacy_;
    void step_continuous(float dt) {
        auto &motions = ccd_motions_;
        auto &legacy = ccd_legacy_;
        // Both arrays are owned by the world and only ever grow, so a steady
        // world refreshes them in place instead of rebuilding them. `legacy` in
        // particular must not be cleared and repopulated: destroying a Fixture
        // frees the heap payload of its Shape, and for a Polygon or Mesh the
        // next assignment has to allocate it again. Assigning into a live
        // element instead reuses that storage.
        motions.resize(bodies_.size());
        legacy.reserve(bodies_.size());
        std::size_t legacy_count = 0;
        for (std::size_t index = 0; index < bodies_.size(); ++index) {
            Body &body = *bodies_[index];
            auto &motion = motions[index];
            motion.colliders.clear();
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
            auto append = [&](Fixture &f) {
                motion.colliders.push_back({&f.shape, f.local, f.material.friction,
                                            f.material.restitution, f.restitution_threshold,
                                            f.collision_group, f.collision_mask, f.trigger, &f});
            };
            if (body.use_default_shape) {
                if (legacy.size() <= legacy_count)
                    legacy.emplace_back(); // Reserve above makes this allocation free.
                Fixture &f = legacy[legacy_count++];
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
        legacy.resize(legacy_count);
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
            &ccd_workspace_);
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
    template <class F> static void visit_pairs(Body &a, Body &b, F &&fn) {
        Fixture fa = legacy_fixture(a), fb = legacy_fixture(b);
        if (a.use_default_shape && b.use_default_shape)
            fn(fa, fb);
        else if (a.use_default_shape) {
            for (auto &y : b.fixtures)
                fn(fa, *y);
        } else if (b.use_default_shape) {
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
    struct VelocityGeometry {
        Vec2 ra[2], rb[2], tangent;
        float normal_mass[2]{}, tangent_mass[2]{};
        float k01{}, determinant{};
        bool movable{};
    };
    std::vector<VelocityGeometry> velocity_geometry_;
    std::vector<ContactKey> solid_contacts_;
    struct ProjectionObstacle {
        Body *body;
        AABB bounds;
    };
    std::vector<ProjectionObstacle> projection_obstacles_;
    StepStatistics step_statistics_{};
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
            cached.shape_hash_a = shape_hash(*c.shape_a);
            cached.shape_hash_b = shape_hash(*c.shape_b);
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
    // A stored content signature is enough to notice a public shape edit; it
    // avoids both the snapshot copy and the O(n) vertex comparison.
    static bool same_shape(const Shape &shape, std::uint64_t hash) {
        return shape_hash(shape) == hash;
    }
    void invalidate_edited_contacts() {
        // Public transforms/shapes have no revision setter. Exact end-of-step
        // snapshots distinguish user edits from normal integration. Destruction
        // APIs clear the cache before releasing any pointed-to object.
        contact_cache_.erase_if([&](const ContactKey &, const CachedContact &old) {
            auto &c = old.c;
            auto edited = [&](Body *b, const Shape *shape, std::uint64_t saved,
                              const Transform &transform, Transform local, std::uintptr_t fixture) {
                return !same_transform(b->transform, transform) || !same_shape(*shape, saved) ||
                       (fixture &&
                        !same_transform(reinterpret_cast<const Fixture *>(fixture)->local, local));
            };
            bool changed =
                edited(c.a, c.shape_a, old.shape_hash_a, old.transform_a, c.fixture_a, c.key[2]) ||
                edited(c.b, c.shape_b, old.shape_hash_b, old.transform_b, c.fixture_b, c.key[3]);
            bool enabled = true;
            visit_pairs(*c.a, *c.b, [&](Fixture &a, Fixture &b) {
                if (contact_key(a, b) == c.key)
                    enabled = !a.trigger && !b.trigger && allowed(a, b);
            });
            if (changed || !enabled) {
                if (c.a->type == BodyType::Dynamic)
                    c.a->wake();
                if (c.b->type == BodyType::Dynamic)
                    c.b->wake();
            }
            return changed || !enabled;
        });
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
        if (!same_transform(a.body->transform, old.transform_a) ||
            !same_transform(b.body->transform, old.transform_b) ||
            !same_transform(a.local, old.c.fixture_a) ||
            !same_transform(b.local, old.c.fixture_b) || !same_shape(a.shape, old.shape_hash_a) ||
            !same_shape(b.shape, old.shape_hash_b))
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
        if (!same_transform(a.body->transform, old.c.detected_a) ||
            !same_transform(b.body->transform, old.c.detected_b) ||
            !same_transform(a.local, old.c.fixture_a) ||
            !same_transform(b.local, old.c.fixture_b) || !same_shape(a.shape, old.shape_hash_a) ||
            !same_shape(b.shape, old.shape_hash_b))
            return nullptr;
        return &old.c;
    }
    void build_constraints(float dt, bool warm) {
        auto start = Clock::now();
        double cache_before = step_statistics_.cache_ms;
        // During the second solve, impulses are already in velocities. Transfer
        // accumulators but do not apply them again.
        if (!warm)
            save_constraints(dt);
        constraints_.clear();
        points_.clear();
        solid_contacts_.clear();
        projection_obstacles_.clear();
        for (auto &b : bodies_)
            if (b->type != BodyType::Dynamic)
                projection_obstacles_.push_back({b.get(), body_aabb(*b)});
        std::sort(projection_obstacles_.begin(), projection_obstacles_.end(),
                  [](const auto &a, const auto &b) { return a.bounds.min.x < b.bounds.min.x; });
        const auto &pairs = candidate_pairs();
        for (auto [i, j] : pairs)
            visit_pairs(*bodies_[i], *bodies_[j], [&](Fixture &fa, Fixture &fb) {
                if (!allowed(fa, fb) || fa.trigger || fb.trigger)
                    return;
                auto remember = [&](Constraint c, const ConstraintPoint *source, int count) {
                    auto radius = [&](std::size_t index) {
                        const auto &bounds = pair_bounds_[index];
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
                    c.point_begin = int(points_.size());
                    for (int k = 0; k < count; ++k)
                        points_.push(source[k]);
                    constraints_.push_back(c);
                };
                auto *old = contact_cache_.find(contact_key(fa, fb));
                const auto *previous_contact = old;
                if (auto cached = sleeping_cache(fa, fb, previous_contact)) {
                    remember(*cached, previous_contact->points, cached->count);
                    constraints_.back().friction =
                        std::sqrt(std::max(0.0f, fa.material.friction * fb.material.friction));
                    if (cached->event_contact)
                        solid_contacts_.push_back(contact_key(fa, fb));
                    ++step_statistics_.sleeping_contacts;
                    return;
                }
                Contact contact;
                auto ta = fixture_transform(fa), tb = fixture_transform(fb);
                const auto *geometry = geometry_cache(fa, fb, previous_contact);
                bool touching;
                if (geometry) {
                    contact = geometry->contact;
                    touching = true;
                    ++step_statistics_.cached_manifolds;
                } else
                    touching = test(fa.shape, ta, fb.shape, tb, contact);
                bool event_contact = geometry ? geometry->event_contact : touching;
                if (!touching && ccd_detail::supported(fa.shape) &&
                    ccd_detail::supported(fb.shape)) {
                    auto sep = ccd_detail::separation(fa.shape, ta, fb.shape, tb);
                    if (sep.distance > (config_.ccd.enabled ? 2 * config_.ccd.tolerance : 0.002f))
                        return;
                    contact = {sep.normal, sep.point, -sep.distance};
                    touching = true;
                    event_contact = config_.ccd.enabled;
                }
                if (!touching)
                    return;
                if (event_contact)
                    solid_contacts_.push_back(contact_key(fa, fb));
                if (!warm && on_contact_diagnostic && contact.penetration > 0)
                    on_contact_diagnostic({fa.body, fb.body, fa.body->velocity, fb.body->velocity,
                                           contact.normal, contact.penetration, fa.body->transform,
                                           fb.body->transform, fa.body->transform,
                                           fb.body->transform, true});
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
                    float speed =
                        (point_velocity(*c.b, rb) - point_velocity(*c.a, ra)).dot(c.normal);
                    p.target =
                        -speed > std::max(1.0f, std::min(fa.restitution_threshold,
                                                         fb.restitution_threshold))
                            ? -std::max(fa.material.restitution, fb.material.restitution) * speed
                            : 0;
                    if (old && old->c.normal.dot(c.normal) > 0.95f) {
                        for (int n = 0; n < old->c.count; ++n) {
                            auto &prev = old->points[n];
                            if (prev.feature == p.feature &&
                                (prev.local_a - p.local_a).length_squared() < 0.04f &&
                                (prev.local_b - p.local_b).length_squared() < 0.04f) {
                                float ratio = old->dt > 0
                                                  ? std::clamp(dt / old->dt, 0.0f, 2.0f)
                                                  : 0;
                                p.normal_impulse = prev.normal_impulse * ratio;
                                p.tangent_impulse = prev.tangent_impulse * ratio;
                                break;
                            }
                        }
                    }
                }
                if (m.count)
                    remember(c, built, m.count);
            });
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
        step_statistics_.constraints = constraints_.size();
        step_statistics_.detection_ms +=
            milliseconds(start) - wake_ms - (step_statistics_.cache_ms - cache_before);
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
    void solve_velocities() {
        auto start = Clock::now();
        build_velocity_geometry();
        for (const SolverIsland &island : solver_islands_)
            for (int iteration = 0; iteration < config_.solver_iterations; ++iteration)
                for (std::size_t k = island.constraint_begin; k < island.constraint_end; ++k)
                    solve_velocity(island_constraints_[k]);
        step_statistics_.velocity_ms += milliseconds(start);
    }
    void build_velocity_geometry() {
        // Transforms and mass properties are constant throughout this velocity pass.
        // Rebuild after integration; never reuse these values for position corrections.
        //
        // The cache is keyed by constraint index so the solver reaches a
        // constraint and its geometry in one hop. It is only grown while there
        // is something to solve: `constraints_` keeps every sleeping pair too,
        // and value-initialising that much storage on a resting step would cost
        // more than the whole solve.
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
                g.normal_mass[k] = ma + mb + ia * na * na + ib * nb * nb;
                g.tangent_mass[k] = ma + mb + ia * sa * sa + ib * sb * sb;
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
            p.normal_impulse = std::max(0.0f, previous + (p.target - rv.dot(c.normal)) / denom);
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
            bool candidate = false;
            for (const auto &entry : projection_obstacles_) {
                if (entry.bounds.min.x > sweep.max.x)
                    break;
                if (sweep.overlaps(entry.bounds)) {
                    candidate = true;
                    break;
                }
            }
            if (!candidate) {
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
        for (const auto &entry : projection_obstacles_) {
            if (entry.bounds.min.x > bounds.max.x)
                break;
            if (!bounds.overlaps(entry.bounds))
                continue;
            ++step_statistics_.projection_candidates;
            auto *obstacle = entry.body;
            visit_pairs(body, *obstacle, [&](Fixture &a, Fixture &b) {
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
                        ShapeSweep fixed{obstacle->transform, obstacle->transform, b.local};
                        float depth = std::max(config_.ccd.tolerance, -initial.distance);
                        if (!ccd_detail::separated_during_sweep(a.shape, projection, b.shape, fixed,
                                                                initial.normal, depth))
                            fraction = 0;
                        return;
                    }
                }
                auto hit =
                    sweep_shapes(a.shape, {before, after, a.local}, b.shape,
                                 {obstacle->transform, obstacle->transform, b.local}, config_.ccd);
                if (hit)
                    fraction = std::min(fraction, hit->fraction);
            });
        }
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
            float impulse = correction * (config_.solver_mode == SolverMode::PBD ? 1.0f : 0.2f) /
                            (ma + mb + ia * na * na + ib * nb * nb);
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
    // Broad phase cache. It is derived entirely from `bodies_`, so it is
    // mutable and refreshed lazily, including from the const query helpers.
    mutable std::vector<AABB> pair_bounds_;
    mutable std::vector<BodyType> pair_types_;
    mutable std::vector<std::pair<std::size_t, std::size_t>> pair_buffer_;
    mutable DynamicTree tree_;
    mutable std::vector<std::uint32_t> proxy_ids_;
    mutable std::vector<AABB> proxy_tight_;
    mutable bool broadphase_dirty_{true};
    mutable float broadphase_margin_{-1};
    float fat_margin() const {
        if (config_.broadphase_fat_margin >= 0)
            return config_.broadphase_fat_margin;
        return std::max(0.002f, config_.ccd.enabled ? config_.ccd.tolerance * 2 : 0.0f);
    }
    // Refresh the proxy set to match the current body transforms. Returns true
    // when any body's exact bounds or type changed, which is the only case in
    // which the candidate pair set can differ from the previous step.
    bool sync_broadphase() const {
        const float margin = fat_margin();
        // A changed margin invalidates the padding stored in every fat box.
        if (margin != broadphase_margin_) {
            broadphase_margin_ = margin;
            broadphase_dirty_ = true;
        }
        const bool structural = broadphase_dirty_ || pair_bounds_.size() != bodies_.size() ||
                                pair_types_.size() != bodies_.size();
        if (structural) {
            tree_.reset();
            proxy_ids_.assign(bodies_.size(), DynamicTree::null_node);
            proxy_tight_.assign(bodies_.size(), AABB{});
            pair_bounds_.assign(bodies_.size(), AABB{});
            pair_types_.assign(bodies_.size(), BodyType::Static);
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                const AABB tight = body_aabb(*bodies_[i]);
                proxy_tight_[i] = tight;
                proxy_ids_[i] = tree_.create_proxy(tight, std::uint32_t(i), margin);
            }
            broadphase_dirty_ = false;
        }
        // Bodies are publicly writable, so the exact bounds are recomputed
        // every pass. A proxy is only re-filed when its tight box leaves the
        // fat box it was stored under, which is the whole point of the tree:
        // resting and static geometry costs nothing.
        bool changed = structural;
        for (std::size_t i = 0; i < bodies_.size(); ++i) {
            const Body &body = *bodies_[i];
            const AABB tight = body_aabb(body);
            const AABB bounds = fatten(tight, margin);
            changed |= bounds.min != pair_bounds_[i].min || bounds.max != pair_bounds_[i].max ||
                       pair_types_[i] != body.type;
            if (proxy_tight_[i].min != tight.min || proxy_tight_[i].max != tight.max) {
                const Vec2 displacement = (tight.min + tight.max) * 0.5f -
                                          (proxy_tight_[i].min + proxy_tight_[i].max) * 0.5f;
                proxy_tight_[i] = tight;
                tree_.move_proxy(proxy_ids_[i], tight, displacement, margin);
            }
            pair_bounds_[i] = bounds;
            pair_types_[i] = body.type;
        }
        return changed;
    }
    const std::vector<std::pair<std::size_t, std::size_t>> &candidate_pairs(bool refresh = true) const {
        if (!refresh)
            return pair_buffer_; // No transforms changed since constraint detection.
        // Even the brute-force oracle path keeps the tree current so both
        // paths agree on when the candidate set has to be recomputed.
        const bool changed = sync_broadphase();
        if (!changed && !pair_buffer_.empty())
            return pair_buffer_;
        pair_buffer_.clear();
        auto append = [&](std::size_t a, std::size_t b) {
            auto i = std::min(a, b), j = std::max(a, b);
            if (i != j &&
                (pair_types_[i] == BodyType::Dynamic || pair_types_[j] == BodyType::Dynamic) &&
                pair_bounds_[i].overlaps(pair_bounds_[j]))
                pair_buffer_.emplace_back(i, j);
        };
        if (!config_.enable_broadphase) {
            for (std::size_t i = 0; i < bodies_.size(); ++i)
                for (std::size_t j = i + 1; j < bodies_.size(); ++j)
                    append(i, j);
        } else {
            // Only dynamic bodies can start a pair, so static geometry is never
            // scanned; it is found through the tree instead.
            const float margin = fat_margin();
            const std::size_t count = bodies_.size();
            for (std::size_t i = 0; i < count; ++i) {
                if (pair_types_[i] != BodyType::Dynamic)
                    continue;
                tree_.query(pair_bounds_[i], margin, [&](std::uint32_t payload) {
                    const auto other = std::size_t(payload);
                    if (other < count)
                        append(i, other);
                });
            }
        }
        std::sort(pair_buffer_.begin(), pair_buffer_.end());
        pair_buffer_.erase(std::unique(pair_buffer_.begin(), pair_buffer_.end()),
                           pair_buffer_.end());
        broadphase_candidate_count_ = pair_buffer_.size();
        return pair_buffer_;
    }
    Body &add_body(std::unique_ptr<Body> body) {
        require_unlocked();
        Body &ref = *body;
        ref.slot = std::uint32_t(bodies_.size());
        bodies_.push_back(std::move(body));
        broadphase_dirty_ = true; // Proxies are indexed by body slot.
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
