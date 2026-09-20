#pragma once

#include "butter/physics2d/shapes.h"
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>

namespace butter::physics2d {

// Convex 2D CCD. Angles are unwrapped: a full revolution must remain 2*pi.
struct ShapeSweep {
    Transform start{}, end{};
    Transform local{};
    Transform at(float fraction) const {
        const float angle = start.angle + (end.angle - start.angle) * fraction;
        return {start.position + (end.position - start.position) * fraction +
                    rotate(local.position, angle),
                angle + local.angle};
    }
};
struct SweepHit {
    float fraction{};
    Vec2 normal{}, point{}; // normal from A to B
    bool converged{true};
};
struct CcdSettings {
    bool enabled{true};
    int max_impacts{32}; // Per moving body, not per world.
    int max_iterations{64};
    float tolerance{0.0001f};
};
enum class CcdFailure { ImpactBudget, NonConvergence, ZeroTimeRepeat };
struct CcdDiagnostic {
    CcdFailure reason{};
    std::size_t body_a{}, body_b{}, collider_a{}, collider_b{};
    float advanced_time{}, remaining_time{};
};
struct CcdStatistics {
    int impacts{};
    std::size_t sweeps{};
    bool limited{};
    float remaining_time{}; // Largest locally clamped interval, not discarded
                            // world time.
    float advanced_time{};
    std::size_t budget_exhaustions{}, non_convergences{}, zero_time_repeats{}, candidates{};
    std::vector<CcdDiagnostic> diagnostics;
};

namespace ccd_detail {
inline bool supported(const Shape &s) {
    return std::holds_alternative<Circle>(s) || std::holds_alternative<Box>(s) ||
           (std::holds_alternative<Polygon>(s) && std::get<Polygon>(s).vertices.size() >= 3);
}
inline float radius(const Shape &s) {
    if (auto *c = std::get_if<Circle>(&s))
        return c->radius;
    if (auto *box = std::get_if<Box>(&s))
        return box->half_extents.length();
    float result = 0;
    for (auto v : world_vertices(s, {}))
        result = std::max(result, v.length());
    return result;
}
inline Vec2 support(const Shape &s, const Transform &t, Vec2 axis) {
    if (auto *c = std::get_if<Circle>(&s))
        return t.position + axis * c->radius;
    const auto vertices = world_vertices(s, t);
    float best = -std::numeric_limits<float>::infinity();
    Vec2 point = t.position;
    for (auto v : vertices)
        if (v.dot(axis) > best) {
            best = v.dot(axis);
            point = v;
        }
    return point;
}
struct Separation {
    float distance{-std::numeric_limits<float>::infinity()};
    Vec2 normal{1, 0}, point{};
};
inline Separation separation(const Shape &a, const Transform &ta, const Shape &b,
                             const Transform &tb) {
    Separation result;
    auto axis_test = [&](Vec2 axis) {
        if (axis.length_squared() < 1.0e-16f)
            return;
        axis = axis.normalized();
        for (Vec2 n : {axis, -axis}) {
            Vec2 pa = support(a, ta, n), pb = support(b, tb, -n);
            float gap = (pb - pa).dot(n);
            if (gap > result.distance)
                result = {gap, n, (pa + pb) * 0.5f};
        }
    };
    const auto va = world_vertices(a, ta), vb = world_vertices(b, tb);
    auto faces = [&](const std::vector<Vec2> &vertices) {
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            Vec2 edge = vertices[(i + 1) % vertices.size()] - vertices[i];
            axis_test({-edge.y, edge.x});
        }
    };
    faces(va);
    faces(vb);
    if (std::holds_alternative<Circle>(a)) {
        if (std::holds_alternative<Circle>(b))
            axis_test(tb.position - ta.position);
        for (auto v : vb)
            axis_test(v - ta.position);
    }
    if (std::holds_alternative<Circle>(b))
        for (auto v : va)
            axis_test(v - tb.position);
    if (!std::isfinite(result.distance))
        axis_test({1, 0});
    const Vec2 tangent{-result.normal.y, result.normal.x};
    auto feature = [&](const Shape &shape, const Transform &t, Vec2 n) {
        Vec2 p = support(shape, t, n);
        float lo = std::numeric_limits<float>::infinity(), hi = -lo;
        if (std::holds_alternative<Circle>(shape))
            return std::make_pair(p.dot(tangent), p.dot(tangent));
        for (auto v : world_vertices(shape, t))
            if (p.dot(n) - v.dot(n) < 1.0e-5f) {
                lo = std::min(lo, v.dot(tangent));
                hi = std::max(hi, v.dot(tangent));
            }
        return std::make_pair(lo, hi);
    };
    const auto [alo, ahi] = feature(a, ta, result.normal);
    const auto [blo, bhi] = feature(b, tb, -result.normal);
    const float along = (std::max(alo, blo) + std::min(ahi, bhi)) * 0.5f;
    result.point = result.normal * result.point.dot(result.normal) + tangent * along;
    return result;
}
inline AABB swept_bounds(const Shape &s, const ShapeSweep &sweep) {
    const float r = radius(s) + sweep.local.position.length();
    return {{std::min(sweep.start.position.x, sweep.end.position.x) - r,
             std::min(sweep.start.position.y, sweep.end.position.y) - r},
            {std::max(sweep.start.position.x, sweep.end.position.x) + r,
             std::max(sweep.start.position.y, sweep.end.position.y) + r}};
}
} // namespace ccd_detail

// Conservative advancement of a separating plane. Rotation is bounded by the
// full angular travel times the radius about the body's pivot, including
// offsets. Initially overlapping shapes belong to the discrete penetration
// solver.
inline std::optional<SweepHit> sweep_shapes(const Shape &a, const ShapeSweep &sa, const Shape &b,
                                            const ShapeSweep &sb,
                                            const CcdSettings &settings = {}) {
    using namespace ccd_detail;
    if (!std::isfinite(settings.tolerance) || settings.tolerance <= 0)
        throw std::invalid_argument("CCD tolerance must be finite and positive");
    if (!supported(a) || !supported(b) || !swept_bounds(a, sa).overlaps(swept_bounds(b, sb)))
        return {};
    const float tolerance = std::max(settings.tolerance, 1.0e-6f);
    const Vec2 relative =
        (sb.end.position - sb.start.position) - (sa.end.position - sa.start.position);
    const float angular =
        std::abs(sa.end.angle - sa.start.angle) *
            ((std::holds_alternative<Circle>(a) ? 0 : radius(a)) + sa.local.position.length()) +
        std::abs(sb.end.angle - sb.start.angle) *
            ((std::holds_alternative<Circle>(b) ? 0 : radius(b)) + sb.local.position.length());
    float fraction = 0;
    Separation distance;
    for (int iteration = 0; iteration < std::max(1, settings.max_iterations); ++iteration) {
        distance = separation(a, sa.at(fraction), b, sb.at(fraction));
        if (fraction == 0 && distance.distance < -tolerance)
            return {};
        const float closing = -relative.dot(distance.normal) + angular;
        if (closing <= 1.0e-8f)
            return {};
        if (distance.distance <= tolerance)
            return SweepHit{fraction, distance.normal, distance.point, true};
        const float advance = (distance.distance - tolerance * 0.5f) / closing;
        if (fraction + advance > 1)
            return {};
        if (advance <= 1.0e-7f || fraction + advance == fraction)
            return SweepHit{fraction, distance.normal, distance.point, false};
        fraction += advance;
    }
    // Do not turn a convergence/budget failure into a false negative.
    distance = separation(a, sa.at(fraction), b, sb.at(fraction));
    return SweepHit{fraction, distance.normal, distance.point, false};
}

// Non-owning views let builder bodies and engine-owned compound fixtures share
// the same TOI solver. The views and their pointed-to state must outlive the
// call.
struct CcdCollider {
    const Shape *shape{};
    Transform local{};
    float friction{0.5f}, restitution{}, restitution_threshold{};
    std::uint32_t group{1}, mask{0xffffffffu};
    bool trigger{};
    void *tag{};
};
struct CcdMotion {
    Transform *transform{};
    Vec2 *velocity{};
    float *angular_velocity{};
    bool *sleeping{};
    int *sleep_counter{};
    float inverse_mass{}, inverse_inertia{};
    bool dynamic{}, kinematic{}, bullet{};
    std::vector<CcdCollider> colliders;
    bool blocked{}; // Local conservative clamp, reset for each advance call.
    bool moves() const { return !blocked && (kinematic || (dynamic && (!sleeping || !*sleeping))); }
    Vec2 linear() const { return moves() ? *velocity : Vec2{}; }
    float angular() const { return moves() ? *angular_velocity : 0; }
    ShapeSweep sweep(float dt, const CcdCollider &f) const {
        return {*transform,
                {transform->position + linear() * dt, transform->angle + angular() * dt},
                f.local};
    }
    void wake() {
        if (dynamic && sleeping && *sleeping) {
            if (sleeping)
                *sleeping = false;
            if (sleep_counter)
                *sleep_counter = 0;
        }
    }
};
using CcdFilter = std::function<bool(const CcdCollider &, const CcdCollider &)>;
using CcdImpact = std::function<void(const CcdCollider &, const CcdCollider &, const SweepHit &)>;

inline CcdStatistics advance_continuous(std::vector<CcdMotion> &bodies, float dt,
                                        const CcdSettings &settings = {},
                                        const CcdFilter &filter = {},
                                        const CcdImpact &impact = {}) {
    CcdStatistics stats;
    if (!std::isfinite(dt) || dt <= 0)
        return stats;
    auto advance = [&](float duration) {
        for (auto &b : bodies)
            if (b.moves()) {
                b.transform->position += b.linear() * duration;
                b.transform->angle += b.angular() * duration;
            }
    };
    for (auto &b : bodies)
        b.blocked = false;
    if (!settings.enabled) {
        advance(dt);
        stats.advanced_time = dt;
        return stats;
    }
    std::vector<int> counts(bodies.size());
    std::map<std::array<std::size_t, 4>, int> zero_hits;
    struct Hit {
        std::size_t i, j, fi, fj;
        SweepHit hit;
    };
    struct Bounds {
        std::size_t i;
        AABB box;
    };
    std::vector<Bounds> bounds;
    bounds.reserve(bodies.size());
    std::vector<Hit> hits;
    hits.reserve(bodies.size());
    float remaining = dt;
    while (remaining > 0) {
        hits.clear();
        bounds.clear();
        float earliest = 1;
        // Sweep-and-prune spatial index; stationary geometry uses tight bounds.
        // Rebuild after impulses because they change swept trajectories.
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            auto &x = bodies[i];
            bool initialized = false;
            AABB total{};
            for (auto &f : x.colliders) {
                auto sweep = x.sweep(remaining, f);
                auto box = std::abs(sweep.end.angle - sweep.start.angle) < 1e-8f
                               ? compute_aabb(*f.shape, sweep.at(0))
                               : ccd_detail::swept_bounds(*f.shape, sweep);
                if (std::abs(sweep.end.angle - sweep.start.angle) < 1e-8f) {
                    auto end = compute_aabb(*f.shape, sweep.at(1));
                    box.min.x = std::min(box.min.x, end.min.x);
                    box.min.y = std::min(box.min.y, end.min.y);
                    box.max.x = std::max(box.max.x, end.max.x);
                    box.max.y = std::max(box.max.y, end.max.y);
                }
                if (!initialized) {
                    total = box;
                    initialized = true;
                } else {
                    total.min.x = std::min(total.min.x, box.min.x);
                    total.min.y = std::min(total.min.y, box.min.y);
                    total.max.x = std::max(total.max.x, box.max.x);
                    total.max.y = std::max(total.max.y, box.max.y);
                }
            }
            if (initialized)
                bounds.push_back({i, total});
        }
        std::sort(bounds.begin(), bounds.end(), [](auto &a, auto &b) {
            return a.box.min.x == b.box.min.x ? a.i < b.i : a.box.min.x < b.box.min.x;
        });
        for (std::size_t bi = 0; bi < bounds.size(); ++bi)
            for (std::size_t bj = bi + 1;
                 bj < bounds.size() && bounds[bj].box.min.x <= bounds[bi].box.max.x; ++bj) {
                if (!bounds[bi].box.overlaps(bounds[bj].box))
                    continue;
                auto i = std::min(bounds[bi].i, bounds[bj].i),
                     j = std::max(bounds[bi].i, bounds[bj].i);
                auto &x = bodies[i];
                auto &y = bodies[j];
                if ((!x.dynamic && !y.dynamic) || (!x.moves() && !y.moves()) ||
                    (x.dynamic && y.dynamic && !x.bullet && !y.bullet))
                    continue;
                ++stats.candidates;
                for (std::size_t fi = 0; fi < x.colliders.size(); ++fi)
                    for (std::size_t fj = 0; fj < y.colliders.size(); ++fj) {
                        auto &fx = x.colliders[fi];
                        auto &fy = y.colliders[fj];
                        if (fx.trigger || fy.trigger || !(fx.group & fy.mask) ||
                            !(fy.group & fx.mask) || (filter && !filter(fx, fy)))
                            continue;
                        ++stats.sweeps;
                        auto hit = sweep_shapes(*fx.shape, x.sweep(remaining, fx), *fy.shape,
                                                y.sweep(remaining, fy), settings);
                        if (!hit)
                            continue;
                        if (hit->fraction == 0) {
                            Vec2 ra = hit->point - x.transform->position,
                                 rb = hit->point - y.transform->position;
                            Vec2 rv = y.linear() + Vec2{-y.angular() * rb.y, y.angular() * rb.x} -
                                      x.linear() - Vec2{-x.angular() * ra.y, x.angular() * ra.x};
                            if (rv.dot(hit->normal) >= -1e-6f &&
                                std::abs(x.angular()) + std::abs(y.angular()) < 1e-6f)
                                continue;
                        }
                        earliest = std::min(earliest, hit->fraction);
                        hits.push_back({i, j, fi, fj, *hit});
                    }
            }
        if (hits.empty()) {
            advance(remaining);
            stats.advanced_time = dt;
            return stats;
        }
        if (earliest > 0)
            zero_hits.clear();
        advance(remaining * earliest);
        remaining *= 1 - earliest;
        // All contacts at this TOI are processed together. Independent bodies
        // never spend each other's budget, even when they share a static floor.
        for (auto &entry : hits) {
            if (entry.hit.fraction > earliest + 1e-6f)
                continue;
            auto *a = &bodies[entry.i];
            auto *b = &bodies[entry.j];
            const auto *fa = &a->colliders[entry.fi];
            const auto *fb = &b->colliders[entry.fj];
            auto *first = &entry.hit;
            bool exhausted = (a->dynamic && counts[entry.i] >= std::max(0, settings.max_impacts)) ||
                             (b->dynamic && counts[entry.j] >= std::max(0, settings.max_impacts));
            bool repeat =
                first->fraction == 0 && ++zero_hits[{entry.i, entry.j, entry.fi, entry.fj}] > 4;
            if (exhausted || !first->converged || repeat) {
                auto reason = !first->converged ? CcdFailure::NonConvergence
                              : exhausted       ? CcdFailure::ImpactBudget
                                                : CcdFailure::ZeroTimeRepeat;
                stats.budget_exhaustions += reason == CcdFailure::ImpactBudget;
                stats.non_convergences += reason == CcdFailure::NonConvergence;
                stats.zero_time_repeats += reason == CcdFailure::ZeroTimeRepeat;
                stats.diagnostics.push_back(
                    {reason, entry.i, entry.j, entry.fi, entry.fj, dt - remaining, remaining});
                stats.limited = true;
                stats.remaining_time = std::max(stats.remaining_time, remaining);
                // Freeze only participants for this interval. Other motion is still
                // swept.
                a->blocked = a->dynamic || a->kinematic;
                b->blocked = b->dynamic || b->kinematic;
                continue;
            }
            a->wake();
            b->wake();
            auto manifold = contact_manifold(*fa->shape, a->sweep(0, *fa).at(0), *fb->shape,
                                             b->sweep(0, *fb).at(0), first->normal, first->point, 0,
                                             settings.tolerance * 2);
            if (!manifold.count) {
                manifold.count = 1;
                manifold.points[0].point = first->point;
            }
            float accumulated[2]{}, friction_sum[2]{}, target[2]{};
            auto apply_impulse = [&](Vec2 ra, Vec2 rb, Vec2 impulse) {
                if (a->dynamic && !a->blocked) {
                    *a->velocity -= impulse * a->inverse_mass;
                    *a->angular_velocity -= ra.cross(impulse) * a->inverse_inertia;
                }
                if (b->dynamic && !b->blocked) {
                    *b->velocity += impulse * b->inverse_mass;
                    *b->angular_velocity += rb.cross(impulse) * b->inverse_inertia;
                }
            };
            auto relative = [&](Vec2 ra, Vec2 rb) {
                return b->linear() + Vec2{-b->angular() * rb.y, b->angular() * rb.x} - a->linear() -
                       Vec2{-a->angular() * ra.y, a->angular() * ra.x};
            };
            const float ma = a->blocked ? 0 : a->inverse_mass,
                        mb = b->blocked ? 0 : b->inverse_mass;
            const float ia = a->blocked ? 0 : a->inverse_inertia,
                        ib = b->blocked ? 0 : b->inverse_inertia;
            for (int k = 0; k < manifold.count; ++k) {
                Vec2 ra = manifold.points[k].point - a->transform->position,
                     rb = manifold.points[k].point - b->transform->position;
                float speed = relative(ra, rb).dot(first->normal);
                target[k] = -speed > std::min(fa->restitution_threshold, fb->restitution_threshold)
                                ? -std::max(fa->restitution, fb->restitution) * speed
                                : 0;
            }
            for (int iteration = 0; iteration < 12; ++iteration) {
                // Coupled normal solve avoids rocking on two-point support faces.
                if (manifold.count == 2) {
                    Vec2 ra0 = manifold.points[0].point - a->transform->position,
                         rb0 = manifold.points[0].point - b->transform->position;
                    Vec2 ra1 = manifold.points[1].point - a->transform->position,
                         rb1 = manifold.points[1].point - b->transform->position;
                    float a0 = ra0.cross(first->normal), a1 = ra1.cross(first->normal),
                          b0 = rb0.cross(first->normal), b1 = rb1.cross(first->normal);
                    float k00 = ma + mb + ia * a0 * a0 + ib * b0 * b0,
                          k11 = ma + mb + ia * a1 * a1 + ib * b1 * b1,
                          k01 = ma + mb + ia * a0 * a1 + ib * b0 * b1;
                    float det = k00 * k11 - k01 * k01;
                    if (det > 1e-8f) {
                        float v0 = target[0] - relative(ra0, rb0).dot(first->normal),
                              v1 = target[1] - relative(ra1, rb1).dot(first->normal);
                        float d0 = (k11 * v0 - k01 * v1) / det, d1 = (k00 * v1 - k01 * v0) / det;
                        if (accumulated[0] + d0 >= 0 && accumulated[1] + d1 >= 0) {
                            accumulated[0] += d0;
                            accumulated[1] += d1;
                            apply_impulse(ra0, rb0, first->normal * d0);
                            apply_impulse(ra1, rb1, first->normal * d1);
                        }
                    }
                }
                for (int k = 0; k < manifold.count; ++k) {
                    Vec2 ra = manifold.points[k].point - a->transform->position,
                         rb = manifold.points[k].point - b->transform->position;
                    float ca = ra.cross(first->normal), cb = rb.cross(first->normal);
                    float denom = ma + mb + ca * ca * ia + cb * cb * ib;
                    if (denom <= 0)
                        continue;
                    float previous = accumulated[k];
                    accumulated[k] = std::max(
                        0.0f, previous + (target[k] - relative(ra, rb).dot(first->normal)) / denom);
                    apply_impulse(ra, rb, first->normal * (accumulated[k] - previous));
                    Vec2 tangent{-first->normal.y, first->normal.x};
                    ca = ra.cross(tangent);
                    cb = rb.cross(tangent);
                    float limit =
                        accumulated[k] * std::sqrt(std::max(0.0f, fa->friction * fb->friction));
                    previous = friction_sum[k];
                    friction_sum[k] =
                        std::clamp(previous - relative(ra, rb).dot(tangent) /
                                                  (ma + mb + ca * ca * ia + cb * cb * ib),
                                   -limit, limit);
                    apply_impulse(ra, rb, tangent * (friction_sum[k] - previous));
                }
            }
            if (a->dynamic)
                ++counts[entry.i];
            if (b->dynamic)
                ++counts[entry.j];
            ++stats.impacts;
            if (impact)
                impact(*fa, *fb, *first);
        }
    }
    stats.advanced_time = dt;
    return stats;
}
} // namespace butter::physics2d
