#pragma once

#include "butter/core/material.h"
#include "butter/physics2d/shapes.h"
#include "butter/physics2d/query.h"
#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace butter::physics2d {

enum class BodyType { Static, Dynamic };
enum class SolverMode { Impulse, PBD };

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
    bool sleeping{false};
    int sleep_counter{0};
    bool is_dynamic() const { return type == BodyType::Dynamic && inverse_mass > 0 && !sleeping; }
    void wake() { sleeping = false; sleep_counter = 0; }
};

struct DistanceJoint {
    Body* a{}; Body* b{}; float length{1}; float stiffness{1};
    void solve(float dt) const {
        if (!a || !b) return; const Vec2 delta = b->transform.position - a->transform.position; const float d = delta.length(); if (d < 1.0e-6f) return;
        const float inv = a->inverse_mass + b->inverse_mass; if (inv <= 0) return; const Vec2 correction = delta / d * ((d - length) * stiffness / inv);
        if (a->is_dynamic()) a->transform.position += correction * a->inverse_mass * -1.0f;
        if (b->is_dynamic()) b->transform.position += correction * b->inverse_mass;
        (void)dt;
    }
};

struct SpringJoint : DistanceJoint { float damping{0.1f}; };
struct HingeJoint : DistanceJoint { Vec2 anchor_a{}, anchor_b{}; };

class World;
class BodyBuilder {
public:
    explicit BodyBuilder(World& world) : world_(world) {}
    BodyBuilder& dynamic() { type_ = BodyType::Dynamic; return *this; }
    BodyBuilder& static_body() { type_ = BodyType::Static; mass_ = 0; return *this; }
    BodyBuilder& at(float x, float y) { position_ = {x, y}; return *this; }
    BodyBuilder& angle(float radians) { angle_ = radians; return *this; }
    BodyBuilder& mass(float value) { mass_ = std::max(value, 0.0001f); return *this; }
    BodyBuilder& circle(float radius) { shape_ = Circle{radius}; return *this; }
    BodyBuilder& box(float hx, float hy) { shape_ = Box{{hx, hy}}; return *this; }
    BodyBuilder& polygon(std::vector<Vec2> vertices) { shape_ = Polygon{std::move(vertices)}; return *this; }
    BodyBuilder& convex(std::vector<Vec2> vertices) { return polygon(std::move(vertices)); }
    BodyBuilder& capsule(float radius, float half_length) { shape_ = Capsule{radius, half_length}; return *this; }
    BodyBuilder& mesh(std::vector<Polygon> triangles) { shape_ = Mesh{std::move(triangles)}; return *this; }
    BodyBuilder& velocity(float x, float y) { velocity_ = {x, y}; return *this; }
    BodyBuilder& angular_velocity(float value) { angular_velocity_ = value; return *this; }
    BodyBuilder& collision_filter(std::uint32_t group, std::uint32_t mask) { group_ = group; mask_ = mask; return *this; }
    BodyBuilder& friction(float value) { material_.friction = value; return *this; }
    BodyBuilder& restitution(float value) { material_.restitution = value; return *this; }
    BodyBuilder& trigger(bool value = true) { trigger_ = value; return *this; }
    Body& build();
private:
    World& world_; BodyType type_{BodyType::Dynamic}; Vec2 position_{}; float angle_{};
    float mass_{1}; Vec2 velocity_{}; float angular_velocity_{0}; Shape shape_{Circle{}}; Material material_{}; bool trigger_{false}; std::uint32_t group_{1}, mask_{0xffffffffu};
};

class World {
public:
    struct Config {
        Vec2 gravity{0, -9.81f}; int solver_iterations{8}; float fixed_timestep{1.0f / 60.0f};
        SolverMode solver_mode{SolverMode::Impulse}; bool enable_broadphase{true};
        float broadphase_cell_size{2.0f}; std::size_t broadphase_max_cells_per_body{256};
    };
    explicit World(const Config& config = {}) : config_(config) {}
    BodyBuilder create_body() { return BodyBuilder(*this); }
    DistanceJoint& add_distance_joint(Body& a, Body& b, float length) { joints_.push_back(std::make_unique<DistanceJoint>(DistanceJoint{&a, &b, length})); return *joints_.back(); }
    SpringJoint& add_spring_joint(Body& a, Body& b, float length, float stiffness = 0.5f) { auto joint = std::make_unique<SpringJoint>(); joint->a = &a; joint->b = &b; joint->length = length; joint->stiffness = stiffness; SpringJoint& ref = *joint; joints_.push_back(std::move(joint)); return ref; }
    HingeJoint& add_hinge_joint(Body& a, Body& b, Vec2 anchor_a = {}, Vec2 anchor_b = {}) { auto joint = std::make_unique<HingeJoint>(); joint->a = &a; joint->b = &b; joint->length = (b.transform.position + anchor_b - a.transform.position - anchor_a).length(); joint->anchor_a = anchor_a; joint->anchor_b = anchor_b; HingeJoint& ref = *joint; joints_.push_back(std::move(joint)); return ref; }
    void step(float dt = -1.0f) {
        if (dt < 0) dt = config_.fixed_timestep;
        for (auto& p : bodies_) if (p->is_dynamic()) {
            p->velocity += config_.gravity * dt;
            p->velocity *= std::max(0.0f, 1.0f - p->linear_damping * dt);
            p->angular_velocity *= std::max(0.0f, 1.0f - p->angular_damping * dt);
            p->transform.position += p->velocity * dt;
            p->transform.angle += p->angular_velocity * dt;
        }
        std::unordered_set<std::uint64_t> current_triggers;
        const auto pairs = candidate_pairs(); broadphase_candidate_count_ = pairs.size();
        for (const auto [i, j] : pairs) {
            Body& a = *bodies_[i]; Body& b = *bodies_[j]; if (!a.is_dynamic() && !b.is_dynamic()) continue;
            Contact trigger_contact; if ((a.trigger || b.trigger) && test(a.shape, a.transform, b.shape, b.transform, trigger_contact))
                current_triggers.insert(pair_key(i, j));
        }
        for (const auto key : current_triggers) if (!active_triggers_.contains(key) && on_trigger) {
            const auto [i, j] = unpack_key(key); emit_trigger(*bodies_[i], *bodies_[j], true);
        }
        for (const auto key : active_triggers_) if (!current_triggers.contains(key) && on_trigger) {
            const auto [i, j] = unpack_key(key); if (i < bodies_.size() && j < bodies_.size()) emit_trigger(*bodies_[i], *bodies_[j], false);
        }
        active_triggers_ = std::move(current_triggers);
        for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
            for (const auto& joint : joints_) joint->solve(dt);
            for (const auto [i, j] : pairs) {
                Body& a = *bodies_[i]; Body& b = *bodies_[j]; if (!a.is_dynamic() && !b.is_dynamic()) continue;
                if ((a.collision_mask & b.collision_group) == 0 || (b.collision_mask & a.collision_group) == 0) continue;
                Contact c; if (!test(a.shape, a.transform, b.shape, b.transform, c) || a.trigger || b.trigger) continue;
                const float inv = a.inverse_mass + b.inverse_mass; if (inv <= 0) continue;
                const float projection = config_.solver_mode == SolverMode::PBD ? 1.0f : 0.8f;
                const Vec2 correction = c.normal * (c.penetration / inv * projection);
                if (a.is_dynamic()) a.transform.position -= correction * a.inverse_mass;
                if (b.is_dynamic()) b.transform.position += correction * b.inverse_mass;
                const Vec2 ra = c.point - a.transform.position, rb = c.point - b.transform.position;
                const auto cross = [](Vec2 x, Vec2 y) { return x.cross(y); };
                const Vec2 va = a.velocity + Vec2{-a.angular_velocity * ra.y, a.angular_velocity * ra.x};
                const Vec2 vb = b.velocity + Vec2{-b.angular_velocity * rb.y, b.angular_velocity * rb.x};
                const float rel = (vb - va).dot(c.normal); if (rel >= 0) continue;
                const float e = std::min(a.material.restitution, b.material.restitution);
                const float rn_a = cross(ra, c.normal), rn_b = cross(rb, c.normal);
                const float denom = inv + rn_a * rn_a * a.inverse_inertia + rn_b * rn_b * b.inverse_inertia;
                const float impulse = -(1 + e) * rel / denom; const Vec2 jv = c.normal * impulse;
                if (a.is_dynamic()) a.velocity -= jv * a.inverse_mass;
                if (b.is_dynamic()) b.velocity += jv * b.inverse_mass;
                if (a.is_dynamic()) a.angular_velocity -= cross(ra, jv) * a.inverse_inertia;
                if (b.is_dynamic()) b.angular_velocity += cross(rb, jv) * b.inverse_inertia;
                const Vec2 tangent = (vb - va - c.normal * rel).normalized();
                const float rt_a = cross(ra, tangent), rt_b = cross(rb, tangent);
                const float friction_denom = inv + rt_a * rt_a * a.inverse_inertia + rt_b * rt_b * b.inverse_inertia;
                const float jt = -(vb - va).dot(tangent) / friction_denom;
                const float max_friction = impulse * std::sqrt(a.material.friction * b.material.friction);
                const Vec2 fv = tangent * std::clamp(jt, -max_friction, max_friction);
                if (a.is_dynamic()) { a.velocity -= fv * a.inverse_mass; a.angular_velocity -= cross(ra, fv) * a.inverse_inertia; }
                if (b.is_dynamic()) { b.velocity += fv * b.inverse_mass; b.angular_velocity += cross(rb, fv) * b.inverse_inertia; }
            }
        }
        for (auto& body : bodies_) if (body->type == BodyType::Dynamic && !body->sleeping) {
            if (body->velocity.length_squared() < 0.0025f && std::abs(body->angular_velocity) < 0.05f) {
                if (++body->sleep_counter > 30) { body->sleeping = true; body->velocity = {}; body->angular_velocity = 0; }
            } else body->sleep_counter = 0;
        }
    }
    std::size_t body_count() const { return bodies_.size(); }
    std::size_t broadphase_candidate_count() const { return broadphase_candidate_count_; }
    std::vector<Body*> query_aabb(const AABB& area) const {
        std::vector<Body*> result; for (const auto& body : bodies_) if (compute_aabb(body->shape, body->transform).overlaps(area)) result.push_back(body.get()); return result;
    }
    std::optional<RaycastHit> raycast(Vec2 origin, Vec2 direction, float max_distance = std::numeric_limits<float>::max()) const {
        direction = direction.normalized(); std::optional<RaycastHit> best;
        for (std::size_t i = 0; i < bodies_.size(); ++i) if (auto hit = ray_aabb(origin, direction, max_distance, compute_aabb(bodies_[i]->shape, bodies_[i]->transform), i))
            if (!best || hit->distance < best->distance) best = hit;
        return best;
    }
    std::function<void(Body&, Body&, bool)> on_trigger;
private:
    friend class BodyBuilder; Config config_; std::vector<std::unique_ptr<Body>> bodies_; std::vector<std::unique_ptr<DistanceJoint>> joints_;
    std::unordered_set<std::uint64_t> active_triggers_{}; std::size_t broadphase_candidate_count_{0};
    static std::uint64_t pair_key(std::size_t i, std::size_t j) { return (std::uint64_t(i) << 32) | std::uint64_t(j); }
    static std::pair<std::size_t, std::size_t> unpack_key(std::uint64_t key) { return {std::size_t(key >> 32), std::size_t(key & 0xffffffffu)}; }
    void emit_trigger(Body& a, Body& b, bool enter) { if (a.trigger) on_trigger(a, b, enter); else on_trigger(b, a, enter); }
    std::vector<std::pair<std::size_t, std::size_t>> candidate_pairs() const {
        std::unordered_set<std::uint64_t> unique;
        if (!config_.enable_broadphase) {
            for (std::size_t i = 0; i < bodies_.size(); ++i) for (std::size_t j = i + 1; j < bodies_.size(); ++j)
                if (bodies_[i]->is_dynamic() || bodies_[j]->is_dynamic()) unique.insert(pair_key(i, j));
        } else {
            std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells; std::vector<std::size_t> large;
            const float cell = std::max(config_.broadphase_cell_size, 0.01f);
            auto cell_key = [](int x, int y) { return (std::uint64_t(std::uint32_t(x)) << 32) | std::uint32_t(y); };
            for (std::size_t i = 0; i < bodies_.size(); ++i) {
                const AABB bounds = compute_aabb(bodies_[i]->shape, bodies_[i]->transform);
                const int min_x = int(std::floor(bounds.min.x / cell)), max_x = int(std::floor(bounds.max.x / cell));
                const int min_y = int(std::floor(bounds.min.y / cell)), max_y = int(std::floor(bounds.max.y / cell));
                const std::size_t count = std::size_t(max_x - min_x + 1) * std::size_t(max_y - min_y + 1);
                if (count > config_.broadphase_max_cells_per_body) { large.push_back(i); continue; }
                for (int x = min_x; x <= max_x; ++x) for (int y = min_y; y <= max_y; ++y) cells[cell_key(x, y)].push_back(i);
            }
            for (const auto& [key, list] : cells) { (void)key; for (std::size_t a = 0; a < list.size(); ++a) for (std::size_t b = a + 1; b < list.size(); ++b) {
                const auto i = std::min(list[a], list[b]), j = std::max(list[a], list[b]);
                if (bodies_[i]->is_dynamic() || bodies_[j]->is_dynamic()) unique.insert(pair_key(i, j));
            }}
            for (const auto i : large) for (std::size_t j = 0; j < bodies_.size(); ++j) if (i != j) {
                const auto a = std::min(i, j), b = std::max(i, j);
                if (bodies_[a]->is_dynamic() || bodies_[b]->is_dynamic()) unique.insert(pair_key(a, b));
            }
        }
        std::vector<std::pair<std::size_t, std::size_t>> result; result.reserve(unique.size());
        for (const auto key : unique) { const auto [i, j] = unpack_key(key); if (compute_aabb(bodies_[i]->shape, bodies_[i]->transform).overlaps(compute_aabb(bodies_[j]->shape, bodies_[j]->transform))) result.emplace_back(i, j); }
        return result;
    }
    Body& add_body(std::unique_ptr<Body> body) { Body& ref = *body; bodies_.push_back(std::move(body)); return ref; }
};

inline Body& BodyBuilder::build() {
    auto body = std::make_unique<Body>(); body->type = type_; body->transform = {position_, angle_}; body->velocity = velocity_;
    body->mass = type_ == BodyType::Static ? 0.0f : mass_; body->inverse_mass = body->mass > 0 ? 1.0f / body->mass : 0.0f;
    body->shape = std::move(shape_); body->material = material_; body->trigger = trigger_; body->angular_velocity = angular_velocity_; body->collision_group = group_; body->collision_mask = mask_;
    body->inertia = body->mass > 0 ? body->mass : 0; body->inverse_inertia = body->inertia > 0 ? 1.0f / body->inertia : 0; return world_.add_body(std::move(body));
}

} // namespace butter::physics2d
