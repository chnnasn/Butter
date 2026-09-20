#pragma once

#include "butter/core/material.h"
#include "butter/physics2d/ccd.h"
#include "butter/physics2d/query.h"
#include "butter/physics2d/shapes.h"
#include <algorithm>
#include <cstdint>
#include <functional>
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
        step_continuous(dt);
        const auto pairs = candidate_pairs();
        broadphase_candidate_count_ = pairs.size();
        std::unordered_set<std::uint64_t> current_triggers;
        std::set<std::pair<Fixture *, Fixture *>> current_contacts;
        for (const auto [i, j] : pairs) {
            Body &a = *bodies_[i];
            Body &b = *bodies_[j];
            visit_pairs(a, b, [&](Fixture &fa, Fixture &fb) {
                Contact c;
                if (!allowed(fa, fb))
                    return;
                bool touching =
                    test(fa.shape, fixture_transform(fa), fb.shape, fixture_transform(fb), c);
                if (!touching && config_.ccd.enabled && active_contacts_.contains({&fa, &fb}) &&
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
        for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
            for (auto &joint : joints_)
                if (joint->spring_stiffness <= 0 || iteration == 0)
                    joint->solve(dt);
            for (auto [i, j] : pairs)
                visit_pairs(*bodies_[i], *bodies_[j], [&](Fixture &fa, Fixture &fb) {
                    if (!allowed(fa, fb) || fa.trigger || fb.trigger)
                        return;
                    Contact c;
                    if (test(fa.shape, fixture_transform(fa), fb.shape, fixture_transform(fb), c))
                        solve_contact(fa, fb, c);
                });
        }
        for (auto &body : bodies_)
            if (body->type == BodyType::Dynamic && !body->sleeping) {
                if (body->velocity.length_squared() < 0.0025f &&
                    std::abs(body->angular_velocity) < 0.05f) {
                    if (++body->sleep_counter > 30) {
                        body->sleeping = true;
                        body->velocity = {};
                        body->angular_velocity = 0;
                    }
                } else
                    body->sleep_counter = 0;
            }
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
    void step_continuous(float dt) {
        std::vector<CcdMotion> motions;
        std::vector<Fixture> legacy;
        motions.reserve(bodies_.size());
        legacy.reserve(bodies_.size());
        for (auto &ptr : bodies_) {
            Body &body = *ptr;
            CcdMotion motion;
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
            motions.push_back(std::move(motion));
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
        const float margin = config_.ccd.enabled ? std::max(config_.ccd.tolerance, 1.0e-6f) * 2 : 0;
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
    void solve_contact(Fixture &fa, Fixture &fb, const Contact &c) {
        Body &a = *fa.body;
        Body &b = *fb.body;
        if (a.sleeping && ((b.is_dynamic() && b.velocity.length_squared() > 0.0025f) ||
                           b.type == BodyType::Kinematic))
            a.wake();
        if (b.sleeping && ((a.is_dynamic() && a.velocity.length_squared() > 0.0025f) ||
                           a.type == BodyType::Kinematic))
            b.wake();
        const float inv = a.inverse_mass + b.inverse_mass;
        if (inv <= 0)
            return;
        const Vec2 ra = c.point - a.transform.position, rb = c.point - b.transform.position;
        const Vec2 va = a.velocity + Vec2{-a.angular_velocity * ra.y, a.angular_velocity * ra.x};
        const Vec2 vb = b.velocity + Vec2{-b.angular_velocity * rb.y, b.angular_velocity * rb.x};
        const float rel = (vb - va).dot(c.normal);
        const float projection = config_.solver_mode == SolverMode::PBD ? 1.0f : 0.8f;
        const Vec2 correction =
            c.normal * (std::clamp(c.penetration - 0.001f, 0.0f, config_.max_position_correction) /
                        inv * projection);
        if (a.is_dynamic())
            a.transform.position -= correction * a.inverse_mass;
        if (b.is_dynamic())
            b.transform.position += correction * b.inverse_mass;
        if (rel >= 0)
            return;
        const float e = -rel > std::min(fa.restitution_threshold, fb.restitution_threshold)
                            ? std::max(fa.material.restitution, fb.material.restitution)
                            : 0;
        const float ca = ra.cross(c.normal), cb = rb.cross(c.normal);
        const float denom = inv + ca * ca * a.inverse_inertia + cb * cb * b.inverse_inertia;
        const float impulse = -(1 + e) * rel / denom;
        const Vec2 j = c.normal * impulse;
        if (a.is_dynamic()) {
            a.velocity -= j * a.inverse_mass;
            a.angular_velocity -= ra.cross(j) * a.inverse_inertia;
        }
        if (b.is_dynamic()) {
            b.velocity += j * b.inverse_mass;
            b.angular_velocity += rb.cross(j) * b.inverse_inertia;
        }
        const Vec2 tangent = (vb - va - c.normal * rel).normalized();
        const float ta = ra.cross(tangent), tb = rb.cross(tangent);
        const float jt = -(vb - va).dot(tangent) /
                         (inv + ta * ta * a.inverse_inertia + tb * tb * b.inverse_inertia);
        const float limit = impulse * std::sqrt(fa.material.friction * fb.material.friction);
        const Vec2 friction = tangent * std::clamp(jt, -limit, limit);
        if (a.is_dynamic()) {
            a.velocity -= friction * a.inverse_mass;
            a.angular_velocity -= ra.cross(friction) * a.inverse_inertia;
        }
        if (b.is_dynamic()) {
            b.velocity += friction * b.inverse_mass;
            b.angular_velocity += rb.cross(friction) * b.inverse_inertia;
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
                const AABB bounds = broadphase_bounds(*bodies_[i]);
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
            if (broadphase_bounds(*bodies_[i]).overlaps(broadphase_bounds(*bodies_[j])))
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
