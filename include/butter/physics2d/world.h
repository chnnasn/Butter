#pragma once

#include "butter/core/material.h"
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
    Vec2 velocity{};
    float angular_velocity{0};
    float mass{1};
    float inverse_mass{1};
    float inertia{1};
    float inverse_inertia{1};
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
    Vec2 force{};
    float torque{};
    std::vector<std::unique_ptr<Fixture>> fixtures;
    bool is_dynamic() const { return type == BodyType::Dynamic && inverse_mass > 0 && !sleeping; }
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
    void solve(float dt) const {
        if (!a || !b)
            return;
        const Vec2 ra = rotate(anchor_a, a->transform.angle);
        const Vec2 rb = rotate(anchor_b, b->transform.angle);
        const Vec2 delta = b->transform.position + rb - a->transform.position - ra;
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
        if (spring_stiffness > 0)
            return;
        const float correction = (d - length) * stiffness / inv;
        if (a->is_dynamic()) {
            a->transform.position += n * correction * a->inverse_mass;
            a->transform.angle += ca * correction * a->inverse_inertia;
        }
        if (b->is_dynamic()) {
            b->transform.position -= n * correction * b->inverse_mass;
            b->transform.angle -= cb * correction * b->inverse_inertia;
        }
    }
};
struct SpringJoint : DistanceJoint {};
struct HingeJoint : DistanceJoint {};

inline Transform fixture_transform(const Fixture &fixture) {
    const auto &t = fixture.body->transform;
    return {t.position + rotate(fixture.local.position, t.angle), t.angle + fixture.local.angle};
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
    std::uint32_t group_{1}, mask_{0xffffffffu};
};

class World {
  public:
    struct StepStatistics {
        double ccd_ms{}, detection_ms{}, velocity_ms{}, position_ms{};
        std::size_t constraints{}, warm_started_points{}, position_clamps{};
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
        float broadphase_cell_size{2.0f};
        std::size_t broadphase_max_cells_per_body{256};
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
    HingeJoint &add_hinge_joint(Body &a, Body &b, Vec2 anchor_a = {}, Vec2 anchor_b = {}) {
        require_unlocked();
        auto joint = std::make_unique<HingeJoint>();
        joint->a = &a;
        joint->b = &b;
        joint->length =
            (b.transform.position + anchor_b - a.transform.position - anchor_a).length();
        joint->anchor_a = anchor_a;
        joint->anchor_b = anchor_b;
        HingeJoint &ref = *joint;
        joints_.push_back(std::move(joint));
        return ref;
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
        return result;
    }
    void destroy_fixture(Fixture &fixture) {
        require_unlocked();
        contact_cache_.clear();
        constraints_.clear();
        forget_contacts(&fixture);
        auto &fixtures = fixture.body->fixtures;
        std::erase_if(fixtures, [&](const auto &p) { return p.get() == &fixture; });
    }
    void destroy_joint(DistanceJoint &joint) {
        require_unlocked();
        std::erase_if(joints_, [&](const auto &p) { return p.get() == &joint; });
    }
    void destroy_body(Body &body) {
        require_unlocked();
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
        for (auto &p : bodies_) {
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
        build_constraints(dt, true);
        solve_velocities();
        auto ccd_start = Clock::now();
        step_continuous(dt);
        step_statistics_.ccd_ms += milliseconds(ccd_start);
        build_constraints(dt, false);
        solve_velocities();
        auto event_start = Clock::now();
        const auto pairs = candidate_pairs();
        broadphase_candidate_count_ = pairs.size();
        std::unordered_set<std::uint64_t> current_triggers;
        std::set<std::pair<Fixture *, Fixture *>> current_contacts;
        for (const auto [i, j] : pairs) {
            Body &a = *bodies_[i];
            Body &b = *bodies_[j];
            visit_pairs(a, b, [&](Fixture &fa, Fixture &fb) {
                if (a.use_default_shape && b.use_default_shape && !fa.trigger && !fb.trigger)
                    return; // Legacy solid bodies have no fixture contact events.
                Contact c;
                if (!allowed(fa, fb))
                    return;
                bool touching =
                    (!fa.trigger && !fb.trigger && active_contacts_.contains({&fa, &fb}) &&
                     sleeping_cache(fa, fb)) ||
                    test(fa.shape, fixture_transform(fa), fb.shape, fixture_transform(fb), c);
                if (!touching && config_.ccd.enabled && !fa.trigger && !fb.trigger &&
                    ccd_detail::supported(fa.shape) && ccd_detail::supported(fb.shape))
                    touching = ccd_detail::separation(fa.shape, fixture_transform(fa), fb.shape,
                                                      fixture_transform(fb))
                                   .distance <= 2 * config_.ccd.tolerance;
                if (!touching)
                    return;
                if (a.use_default_shape && b.use_default_shape) {
                    if (fa.trigger || fb.trigger)
                        current_triggers.insert(pair_key(i, j));
                } else if (!a.use_default_shape && !b.use_default_shape) {
                    auto key = std::make_pair(&fa, &fb);
                    current_contacts.insert(key);
                    if (!active_contacts_.contains(key) && on_contact)
                        on_contact(fa, fb, true);
                }
            });
        }
        for (auto key : active_contacts_)
            if (!current_contacts.contains(key) && on_contact)
                on_contact(*key.first, *key.second, false);
        active_contacts_ = std::move(current_contacts);
        for (auto key : current_triggers)
            if (!active_triggers_.contains(key) && on_trigger) {
                auto [i, j] = unpack_key(key);
                emit_trigger(*bodies_[i], *bodies_[j], true);
            }
        for (auto key : active_triggers_)
            if (!current_triggers.contains(key) && on_trigger) {
                auto [i, j] = unpack_key(key);
                emit_trigger(*bodies_[i], *bodies_[j], false);
            }
        active_triggers_ = std::move(current_triggers);
        step_statistics_.detection_ms += milliseconds(event_start);
        auto position_start = Clock::now();
        for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
            for (auto &joint : joints_)
                if (joint->spring_stiffness <= 0 || iteration == 0) {
                    auto ta = joint->a->transform, tb = joint->b->transform;
                    joint->solve(dt);
                    guard_projection(*joint->a, ta);
                    guard_projection(*joint->b, tb);
                }
            for (auto &c : constraints_)
                solve_position(c);
        }
        step_statistics_.position_ms += milliseconds(position_start);
        update_sleep();
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
    std::size_t broadphase_candidate_count() const { return broadphase_candidate_count_; }
    std::vector<Body *> query_aabb(const AABB &area) const {
        std::vector<Body *> result;
        for (const auto &body : bodies_)
            if (body_aabb(*body).overlaps(area))
                result.push_back(body.get());
        return result;
    }
    std::optional<RaycastHit>
    raycast(Vec2 origin, Vec2 direction,
            float max_distance = std::numeric_limits<float>::max()) const {
        direction = direction.normalized();
        std::optional<RaycastHit> best;
        for (std::size_t i = 0; i < bodies_.size(); ++i) {
            const auto &body = *bodies_[i];
            auto consider = [&](const Shape &shape, const Transform &t) {
                if (auto hit = ray_shape(origin, direction, max_distance, shape, t, i))
                    if (!best || hit->distance < best->distance)
                        best = hit;
            };
            if (body.use_default_shape)
                consider(body.shape, body.transform);
            else
                for (const auto &fixture : body.fixtures)
                    consider(fixture->shape, fixture_transform(*fixture));
        }
        return best;
    }
    std::function<void(Body &, Body &, bool)> on_trigger;

  private:
    CcdStatistics ccd_statistics_{};
    std::vector<CcdMotion> ccd_motions_;
    std::vector<Fixture> ccd_legacy_;
    void step_continuous(float dt) {
        auto &motions = ccd_motions_;
        auto &legacy = ccd_legacy_;
        motions.resize(bodies_.size());
        legacy.clear();
        legacy.reserve(bodies_.size());
        for (std::size_t index = 0; index < bodies_.size(); ++index) {
            Body &body = *bodies_[index];
            auto &motion = motions[index];
            motion.colliders.clear();
            motion.transform = &body.transform;
            motion.velocity = &body.velocity;
            motion.angular_velocity = &body.angular_velocity;
            motion.sleeping = &body.sleeping;
            motion.sleep_counter = &body.sleep_counter;
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
                legacy.push_back(legacy_fixture(body));
                append(legacy.back());
            } else
                for (auto &fixture : body.fixtures)
                    append(*fixture);
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
            });
    }
    AABB broadphase_bounds(const Body &body) const {
        AABB bounds = body_aabb(body);
        const float margin =
            std::max(0.002f, config_.ccd.enabled ? config_.ccd.tolerance * 2 : 0.0f);
        bounds.min -= Vec2{margin, margin};
        bounds.max += Vec2{margin, margin};
        return bounds;
    }
    struct CallbackLock {
        bool &flag;
        bool previous;
        explicit CallbackLock(bool &value) : flag(value), previous(value) { flag = true; }
        ~CallbackLock() { flag = previous; }
    };
    bool locked_{false};
    std::set<std::pair<Fixture *, Fixture *>> active_contacts_;
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
    struct Constraint {
        ContactKey key{};
        Body *a{}, *b{};
        Vec2 normal{};
        float friction{};
        int count{};
        ConstraintPoint points[2]{};
        const Shape *shape_a{}, *shape_b{};
        Transform fixture_a{}, fixture_b{};
    };
    struct CachedContact {
        Constraint c;
        float dt{};
        std::size_t generation{};
        Shape shape_a{}, shape_b{};
        Transform transform_a{}, transform_b{};
    };
    std::size_t cache_generation_{};
    std::map<ContactKey, CachedContact> contact_cache_;
    std::vector<Constraint> constraints_;
    std::vector<Body *> projection_obstacles_;
    StepStatistics step_statistics_{};
    double simulation_time_{};
    static float inv_mass(const Body &b) { return b.is_dynamic() ? b.inverse_mass : 0; }
    static float inv_inertia(const Body &b) {
        return b.is_dynamic() && !b.fixed_rotation ? b.inverse_inertia : 0;
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
        ++cache_generation_;
        for (auto &c : constraints_) {
            auto &cached = contact_cache_[c.key];
            cached.c = c;
            cached.dt = dt;
            cached.generation = cache_generation_;
            cached.shape_a = *c.shape_a;
            cached.shape_b = *c.shape_b;
            cached.transform_a = c.a->transform;
            cached.transform_b = c.b->transform;
        }
        std::erase_if(contact_cache_, [&](const auto &item) {
            return item.second.generation != cache_generation_;
        });
    }
    static ContactKey contact_key(const Fixture &a, const Fixture &b) {
        return {reinterpret_cast<std::uintptr_t>(a.body), reinterpret_cast<std::uintptr_t>(b.body),
                a.body->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&a),
                b.body->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&b)};
    }
    static bool same_transform(const Transform &a, const Transform &b) {
        return a.position == b.position && a.angle == b.angle;
    }
    static bool same_shape(const Shape &a, const Shape &b) {
        if (a.index() != b.index())
            return false;
        if (auto p = std::get_if<Circle>(&a))
            return p->radius == std::get<Circle>(b).radius;
        if (auto p = std::get_if<Box>(&a))
            return p->half_extents == std::get<Box>(b).half_extents;
        if (auto p = std::get_if<Polygon>(&a))
            return p->vertices == std::get<Polygon>(b).vertices;
        if (auto p = std::get_if<Capsule>(&a)) {
            auto q = std::get<Capsule>(b);
            return p->radius == q.radius && p->half_length == q.half_length;
        }
        auto &x = std::get<Mesh>(a).triangles;
        auto &y = std::get<Mesh>(b).triangles;
        if (x.size() != y.size())
            return false;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i].vertices != y[i].vertices)
                return false;
        return true;
    }
    const Constraint *sleeping_cache(const Fixture &a, const Fixture &b) const {
        if (a.body->is_dynamic() || b.body->is_dynamic() || a.body->type == BodyType::Kinematic ||
            b.body->type == BodyType::Kinematic)
            return nullptr;
        auto it = contact_cache_.find(contact_key(a, b));
        if (it == contact_cache_.end())
            return nullptr;
        auto &old = it->second;
        // Bodies and shapes are publicly writable: exact snapshots invalidate
        // the sleep shortcut after teleports, geometry edits or fixture offsets.
        if (!same_transform(a.body->transform, old.transform_a) ||
            !same_transform(b.body->transform, old.transform_b) ||
            !same_transform(a.local, old.c.fixture_a) ||
            !same_transform(b.local, old.c.fixture_b) || !same_shape(a.shape, old.shape_a) ||
            !same_shape(b.shape, old.shape_b))
            return nullptr;
        return &old.c;
    }
    void build_constraints(float dt, bool warm) {
        auto start = Clock::now();
        // During the second solve, impulses are already in velocities. Transfer
        // accumulators but do not apply them again.
        if (!warm)
            save_constraints(dt);
        constraints_.clear();
        projection_obstacles_.clear();
        for (auto &b : bodies_)
            if (b->type != BodyType::Dynamic)
                projection_obstacles_.push_back(b.get());
        auto pairs = candidate_pairs();
        for (auto [i, j] : pairs)
            visit_pairs(*bodies_[i], *bodies_[j], [&](Fixture &fa, Fixture &fb) {
                if (!allowed(fa, fb) || fa.trigger || fb.trigger)
                    return;
                if (auto cached = sleeping_cache(fa, fb)) {
                    constraints_.push_back(*cached);
                    return;
                }
                Contact contact;
                auto ta = fixture_transform(fa), tb = fixture_transform(fb);
                bool touching = test(fa.shape, ta, fb.shape, tb, contact);
                if (!touching && ccd_detail::supported(fa.shape) &&
                    ccd_detail::supported(fb.shape)) {
                    auto sep = ccd_detail::separation(fa.shape, ta, fb.shape, tb);
                    if (sep.distance > (config_.ccd.enabled ? 2 * config_.ccd.tolerance : 0.002f))
                        return;
                    contact = {sep.normal, sep.point, -sep.distance};
                    touching = true;
                }
                if (!touching)
                    return;
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
                c.key = {reinterpret_cast<std::uintptr_t>(c.a),
                         reinterpret_cast<std::uintptr_t>(c.b),
                         c.a->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&fa),
                         c.b->use_default_shape ? 0 : reinterpret_cast<std::uintptr_t>(&fb)};
                c.friction = std::sqrt(std::max(0.0f, fa.material.friction * fb.material.friction));
                auto m = contact_manifold(fa.shape, ta, fb.shape, tb, contact.normal, contact.point,
                                          -contact.penetration);
                c.count = m.count;
                for (int k = 0; k < c.count; ++k) {
                    auto &p = c.points[k];
                    auto &mp = m.points[k];
                    p.feature = mp.feature;
                    Vec2 pa = mp.point - c.normal * (mp.separation * 0.5f),
                         pb = mp.point + c.normal * (mp.separation * 0.5f);
                    p.local_a = rotate(pa - c.a->transform.position, -c.a->transform.angle);
                    p.local_b = rotate(pb - c.b->transform.position, -c.b->transform.angle);
                    float speed = (point_velocity(*c.b, pb - c.b->transform.position) -
                                   point_velocity(*c.a, pa - c.a->transform.position))
                                      .dot(c.normal);
                    p.target =
                        -speed > std::max(1.0f, std::min(fa.restitution_threshold,
                                                         fb.restitution_threshold))
                            ? -std::max(fa.material.restitution, fb.material.restitution) * speed
                            : 0;
                    auto old = contact_cache_.find(c.key);
                    if (old != contact_cache_.end() && old->second.c.normal.dot(c.normal) > 0.95f) {
                        for (int n = 0; n < old->second.c.count; ++n) {
                            auto &prev = old->second.c.points[n];
                            if (prev.feature == p.feature &&
                                (prev.local_a - p.local_a).length_squared() < 0.04f &&
                                (prev.local_b - p.local_b).length_squared() < 0.04f) {
                                float ratio = old->second.dt > 0
                                                  ? std::clamp(dt / old->second.dt, 0.0f, 2.0f)
                                                  : 0;
                                p.normal_impulse = prev.normal_impulse * ratio;
                                p.tangent_impulse = prev.tangent_impulse * ratio;
                                break;
                            }
                        }
                    }
                }
                if (c.count)
                    constraints_.push_back(c);
            });
        wake_connected();
        if (warm)
            for (auto &c : constraints_) {
                if (!c.a->is_dynamic() && !c.b->is_dynamic())
                    continue;
                for (int k = 0; k < c.count; ++k) {
                    auto &p = c.points[k];
                    Vec2 ra = rotate(p.local_a, c.a->transform.angle),
                         rb = rotate(p.local_b, c.b->transform.angle);
                    if (p.normal_impulse > 0)
                        ++step_statistics_.warm_started_points;
                    apply(c, ra, rb,
                          c.normal * p.normal_impulse +
                              Vec2{-c.normal.y, c.normal.x} * p.tangent_impulse);
                }
            }
        step_statistics_.constraints = constraints_.size();
        step_statistics_.detection_ms += milliseconds(start);
    }
    void solve_velocities() {
        auto start = Clock::now();
        for (int iteration = 0; iteration < config_.solver_iterations; ++iteration)
            for (auto &c : constraints_) {
                float ma = inv_mass(*c.a), mb = inv_mass(*c.b), ia = inv_inertia(*c.a),
                      ib = inv_inertia(*c.b);
                if (ma + mb <= 0)
                    continue;
                if (c.count == 2) {
                    Vec2 ra0 = rotate(c.points[0].local_a, c.a->transform.angle),
                         rb0 = rotate(c.points[0].local_b, c.b->transform.angle);
                    Vec2 ra1 = rotate(c.points[1].local_a, c.a->transform.angle),
                         rb1 = rotate(c.points[1].local_b, c.b->transform.angle);
                    float a0 = ra0.cross(c.normal), a1 = ra1.cross(c.normal),
                          b0 = rb0.cross(c.normal), b1 = rb1.cross(c.normal);
                    float k00 = ma + mb + ia * a0 * a0 + ib * b0 * b0,
                          k11 = ma + mb + ia * a1 * a1 + ib * b1 * b1,
                          k01 = ma + mb + ia * a0 * a1 + ib * b0 * b1;
                    float det = k00 * k11 - k01 * k01;
                    if (det > 1e-8f) {
                        float v0 =
                            c.points[0].target -
                            (point_velocity(*c.b, rb0) - point_velocity(*c.a, ra0)).dot(c.normal);
                        float v1 =
                            c.points[1].target -
                            (point_velocity(*c.b, rb1) - point_velocity(*c.a, ra1)).dot(c.normal);
                        float d0 = (k11 * v0 - k01 * v1) / det, d1 = (k00 * v1 - k01 * v0) / det;
                        if (c.points[0].normal_impulse + d0 >= 0 &&
                            c.points[1].normal_impulse + d1 >= 0) {
                            c.points[0].normal_impulse += d0;
                            c.points[1].normal_impulse += d1;
                            apply(c, ra0, rb0, c.normal * d0);
                            apply(c, ra1, rb1, c.normal * d1);
                        }
                    }
                }
                for (int k = 0; k < c.count; ++k) {
                    auto &p = c.points[k];
                    Vec2 ra = rotate(p.local_a, c.a->transform.angle),
                         rb = rotate(p.local_b, c.b->transform.angle);
                    Vec2 rv = point_velocity(*c.b, rb) - point_velocity(*c.a, ra);
                    float na = ra.cross(c.normal), nb = rb.cross(c.normal),
                          denom = ma + mb + ia * na * na + ib * nb * nb;
                    float previous = p.normal_impulse;
                    p.normal_impulse =
                        std::max(0.0f, previous + (p.target - rv.dot(c.normal)) / denom);
                    apply(c, ra, rb, c.normal * (p.normal_impulse - previous));
                    Vec2 tangent{-c.normal.y, c.normal.x};
                    float sa = ra.cross(tangent), sb = rb.cross(tangent);
                    rv = point_velocity(*c.b, rb) - point_velocity(*c.a, ra);
                    previous = p.tangent_impulse;
                    float limit = c.friction * p.normal_impulse;
                    p.tangent_impulse = std::clamp(
                        previous - rv.dot(tangent) / (ma + mb + ia * sa * sa + ib * sb * sb),
                        -limit, limit);
                    apply(c, ra, rb, tangent * (p.tangent_impulse - previous));
                }
            }
        step_statistics_.velocity_ms += milliseconds(start);
    }
    void guard_projection(Body &body, Transform before) {
        if (!config_.ccd.enabled || !body.is_dynamic())
            return;
        const auto after = body.transform;
        if ((after.position - before.position).length_squared() < 1e-16f &&
            std::abs(after.angle - before.angle) < 1e-8f)
            return;
        float fraction = 1;
        for (auto *obstacle : projection_obstacles_) {
            if (!body_aabb(body).overlaps(body_aabb(*obstacle))) {
                auto bounds = body_aabb(body);
                auto delta = before.position - after.position;
                bounds.min.x += std::min(0.0f, delta.x);
                bounds.min.y += std::min(0.0f, delta.y);
                bounds.max.x += std::max(0.0f, delta.x);
                bounds.max.y += std::max(0.0f, delta.y);
                if (std::abs(after.angle - before.angle) < 1e-8f &&
                    !bounds.overlaps(body_aabb(*obstacle)))
                    continue;
            }
            visit_pairs(body, *obstacle, [&](Fixture &a, Fixture &b) {
                if (a.trigger || b.trigger || !allowed(a, b))
                    return;
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
        for (int k = 0; k < c.count; ++k) {
            auto &p = c.points[k];
            Vec2 ra = rotate(p.local_a, c.a->transform.angle),
                 rb = rotate(p.local_b, c.b->transform.angle);
            float separation =
                (c.b->transform.position + rb - c.a->transform.position - ra).dot(c.normal);
            float correction =
                std::clamp(-separation - 0.001f, 0.0f, config_.max_position_correction);
            if (correction == 0)
                continue;
            auto before_a = c.a->transform, before_b = c.b->transform;
            float na = ra.cross(c.normal), nb = rb.cross(c.normal);
            float impulse = correction * (config_.solver_mode == SolverMode::PBD ? 1.0f : 0.2f) /
                            (ma + mb + ia * na * na + ib * nb * nb);
            c.a->transform.position -= c.normal * (impulse * ma);
            c.a->transform.angle -= na * impulse * ia;
            c.b->transform.position += c.normal * (impulse * mb);
            c.b->transform.angle += nb * impulse * ib;
            guard_projection(*c.a, before_a);
            guard_projection(*c.b, before_b);
            if (on_contact_diagnostic)
                on_contact_diagnostic({c.a, c.b, c.a->velocity, c.b->velocity, c.normal,
                                       -separation, before_a, before_b, c.a->transform,
                                       c.b->transform});
        }
    }
    // Static floors do not merge unrelated islands. Contacts and joints connect
    // dynamic bodies; an island sleeps only when every member is quiet.
    std::vector<std::vector<Body *>> islands() {
        std::map<Body *, std::size_t> index;
        for (auto &b : bodies_)
            if (b->type == BodyType::Dynamic)
                index.emplace(b.get(), index.size());
        std::vector<std::size_t> parent(index.size());
        for (std::size_t i = 0; i < parent.size(); ++i)
            parent[i] = i;
        auto root = [&](std::size_t i) {
            while (parent[i] != i) {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            return i;
        };
        auto join = [&](Body *a, Body *b) {
            if (index.contains(a) && index.contains(b))
                parent[root(index[a])] = root(index[b]);
        };
        for (auto &c : constraints_)
            join(c.a, c.b);
        for (auto &j : joints_)
            join(j->a, j->b);
        std::vector<std::vector<Body *>> groups(parent.size());
        for (auto [b, i] : index)
            groups[root(i)].push_back(b);
        return groups;
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
            bool quiet = true;
            int counter = 31;
            for (auto *b : group) {
                quiet &=
                    b->velocity.length_squared() < 0.0025f && std::abs(b->angular_velocity) < 0.05f;
                counter = std::min(counter, b->sleep_counter);
            }
            for (auto *b : group) {
                b->sleep_counter = quiet ? counter + 1 : 0;
                if (b->sleep_counter > 30) {
                    b->sleeping = true;
                    b->velocity = {};
                    b->angular_velocity = 0;
                }
            }
        }
    }
    friend class BodyBuilder;
    Config config_;
    std::vector<std::unique_ptr<Body>> bodies_;
    std::vector<std::unique_ptr<DistanceJoint>> joints_;
    std::unordered_set<std::uint64_t> active_triggers_{};
    std::size_t broadphase_candidate_count_{0};
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
    std::vector<std::pair<std::size_t, std::size_t>> candidate_pairs() const {
        std::vector<AABB> bounds_cache;
        bounds_cache.reserve(bodies_.size());
        for (auto &b : bodies_)
            bounds_cache.push_back(broadphase_bounds(*b));
        std::unordered_set<std::uint64_t> unique;
        if (!config_.enable_broadphase) {
            for (std::size_t i = 0; i < bodies_.size(); ++i)
                for (std::size_t j = i + 1; j < bodies_.size(); ++j)
                    if (bodies_[i]->type == BodyType::Dynamic ||
                        bodies_[j]->type == BodyType::Dynamic)
                        unique.insert(pair_key(i, j));
        } else {
            std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells;
            std::vector<std::size_t> large;
            const float cell = std::max(config_.broadphase_cell_size, 0.01f);
            auto cell_key = [](int x, int y) {
                return (std::uint64_t(std::uint32_t(x)) << 32) | std::uint32_t(y);
            };
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                const AABB bounds = bounds_cache[i];
                const int min_x = int(std::floor(bounds.min.x / cell)),
                          max_x = int(std::floor(bounds.max.x / cell));
                const int min_y = int(std::floor(bounds.min.y / cell)),
                          max_y = int(std::floor(bounds.max.y / cell));
                const std::size_t count =
                    std::size_t(max_x - min_x + 1) * std::size_t(max_y - min_y + 1);
                if (count > config_.broadphase_max_cells_per_body) {
                    large.push_back(i);
                    continue;
                }
                for (int x = min_x; x <= max_x; ++x)
                    for (int y = min_y; y <= max_y; ++y)
                        cells[cell_key(x, y)].push_back(i);
            }
            for (const auto &[key, list] : cells) {
                (void)key;
                for (std::size_t a = 0; a < list.size(); ++a)
                    for (std::size_t b = a + 1; b < list.size(); ++b) {
                        const auto i = std::min(list[a], list[b]), j = std::max(list[a], list[b]);
                        if (bodies_[i]->type == BodyType::Dynamic ||
                            bodies_[j]->type == BodyType::Dynamic)
                            unique.insert(pair_key(i, j));
                    }
            }
            for (const auto i : large)
                for (std::size_t j = 0; j < bodies_.size(); ++j)
                    if (i != j) {
                        const auto a = std::min(i, j), b = std::max(i, j);
                        if (bodies_[a]->type == BodyType::Dynamic ||
                            bodies_[b]->type == BodyType::Dynamic)
                            unique.insert(pair_key(a, b));
                    }
        }
        std::vector<std::pair<std::size_t, std::size_t>> result;
        result.reserve(unique.size());
        for (const auto key : unique) {
            const auto [i, j] = unpack_key(key);
            if (bounds_cache[i].overlaps(bounds_cache[j]))
                result.emplace_back(i, j);
        }
        std::sort(result.begin(), result.end());
        return result;
    }
    Body &add_body(std::unique_ptr<Body> body) {
        require_unlocked();
        Body &ref = *body;
        bodies_.push_back(std::move(body));
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
    body->inertia = body->mass > 0 ? body->mass : 0;
    body->inverse_inertia = body->inertia > 0 ? 1.0f / body->inertia : 0;
    return world_.add_body(std::move(body));
}

} // namespace butter::physics2d
