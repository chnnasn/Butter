#pragma once

#include "butter/physics2d/shapes.h"
#include <cstdint>
#include <functional>
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
    int max_impacts{32};
    int max_iterations{64};
    float tolerance{0.0001f};
};
struct CcdStatistics {
    int impacts{};
    std::size_t sweeps{};
    bool limited{};
    float remaining_time{};
};

namespace ccd_detail {
inline bool supported(const Shape &s) {
    return std::holds_alternative<Circle>(s) || std::holds_alternative<Box>(s) ||
           (std::holds_alternative<Polygon>(s) && std::get<Polygon>(s).vertices.size() >= 3);
}
inline float radius(const Shape &s) {
    if (auto *c = std::get_if<Circle>(&s))
        return c->radius;
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
// full angular travel times the radius about the body's pivot, including offsets.
// Initially overlapping shapes belong to the discrete penetration solver.
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
// the same TOI solver. The views and their pointed-to state must outlive the call.
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
    bool moves() const { return kinematic || (dynamic && (!sleeping || !*sleeping)); }
    Vec2 linear() const { return moves() ? *velocity : Vec2{}; }
    float angular() const { return moves() ? *angular_velocity : 0; }
    ShapeSweep sweep(float dt, const CcdCollider &f) const {
        return {*transform,
                {transform->position + linear() * dt, transform->angle + angular() * dt},
                f.local};
    }
    void wake() {
        if (dynamic) {
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
    if (!settings.enabled) {
        advance(dt);
        return stats;
    }
    float remaining = dt;
    while (remaining > 0) {
        std::optional<SweepHit> first;
        CcdMotion *a = nullptr, *b = nullptr;
        const CcdCollider *fa = nullptr, *fb = nullptr;
        // Enumerate from dynamic bodies: static-only scenes do no pair work.
        for (std::size_t dynamic_index = 0; dynamic_index < bodies.size(); ++dynamic_index)
            if (bodies[dynamic_index].dynamic)
                for (std::size_t other = 0; other < bodies.size(); ++other) {
                    if (other == dynamic_index || (bodies[other].dynamic && other < dynamic_index))
                        continue;
                    const auto i = std::min(dynamic_index, other),
                               j = std::max(dynamic_index, other);
                    auto &x = bodies[i];
                    auto &y = bodies[j];
                    if ((!x.dynamic && !y.dynamic) || (!x.moves() && !y.moves()) ||
                        (x.dynamic && y.dynamic && !x.bullet && !y.bullet))
                        continue;
                    for (const auto &fx : x.colliders)
                        for (const auto &fy : y.colliders) {
                            if (fx.trigger || fy.trigger || !(fx.group & fy.mask) ||
                                !(fy.group & fx.mask))
                                continue;
                            const auto sx = x.sweep(remaining, fx), sy = y.sweep(remaining, fy);
                            if (!ccd_detail::supported(*fx.shape) ||
                                !ccd_detail::supported(*fy.shape) ||
                                !ccd_detail::swept_bounds(*fx.shape, sx)
                                     .overlaps(ccd_detail::swept_bounds(*fy.shape, sy)) ||
                                (filter && !filter(fx, fy)))
                                continue;
                            ++stats.sweeps;
                            auto hit = sweep_shapes(*fx.shape, x.sweep(remaining, fx), *fy.shape,
                                                    y.sweep(remaining, fy), settings);
                            if (!hit)
                                continue;
                            // A touching pair moving apart must not consume the impact budget.
                            if (hit->fraction == 0) {
                                Vec2 ra = hit->point - x.transform->position,
                                     rb = hit->point - y.transform->position;
                                Vec2 rv =
                                    y.linear() + Vec2{-y.angular() * rb.y, y.angular() * rb.x} -
                                    x.linear() - Vec2{-x.angular() * ra.y, x.angular() * ra.x};
                                if (rv.dot(hit->normal) >= -1.0e-6f &&
                                    std::abs(x.angular()) + std::abs(y.angular()) < 1.0e-6f)
                                    continue;
                            }
                            if (!first || hit->fraction < first->fraction) {
                                first = hit;
                                a = &x;
                                b = &y;
                                fa = &fx;
                                fb = &fy;
                            }
                        }
                }
        if (!first) {
            advance(remaining);
            return stats;
        }
        if (stats.impacts >= std::max(0, settings.max_impacts)) {
            stats.limited = true;
            stats.remaining_time = remaining;
            return stats;
        }
        advance(remaining * first->fraction);
        remaining *= 1 - first->fraction;
        if (!first->converged) {
            stats.limited = true;
            stats.remaining_time = remaining;
            return stats;
        }
        a->wake();
        b->wake();
        const Vec2 ra = first->point - a->transform->position,
                   rb = first->point - b->transform->position;
        const Vec2 rv = b->linear() + Vec2{-b->angular() * rb.y, b->angular() * rb.x} -
                        a->linear() - Vec2{-a->angular() * ra.y, a->angular() * ra.x};
        const float ca = ra.cross(first->normal), cb = rb.cross(first->normal);
        const float inv = a->inverse_mass + b->inverse_mass;
        const float denom = inv + ca * ca * a->inverse_inertia + cb * cb * b->inverse_inertia;
        const float speed = rv.dot(first->normal);
        if (denom > 0 && speed < 0) {
            const float restitution =
                -speed > std::min(fa->restitution_threshold, fb->restitution_threshold)
                    ? std::max(fa->restitution, fb->restitution)
                    : 0;
            const float magnitude = -(1 + restitution) * speed / denom;
            const Vec2 impulse = first->normal * magnitude;
            if (a->dynamic) {
                *a->velocity -= impulse * a->inverse_mass;
                *a->angular_velocity -= ra.cross(impulse) * a->inverse_inertia;
            }
            if (b->dynamic) {
                *b->velocity += impulse * b->inverse_mass;
                *b->angular_velocity += rb.cross(impulse) * b->inverse_inertia;
            }
            const Vec2 tangent = (rv - first->normal * speed).normalized();
            const float ta = ra.cross(tangent), tb = rb.cross(tangent);
            const float divisor = inv + ta * ta * a->inverse_inertia + tb * tb * b->inverse_inertia;
            if (divisor > 0) {
                const float limit =
                    magnitude * std::sqrt(std::max(0.0f, fa->friction * fb->friction));
                const Vec2 friction =
                    tangent * std::clamp(-rv.dot(tangent) / divisor, -limit, limit);
                if (a->dynamic) {
                    *a->velocity -= friction * a->inverse_mass;
                    *a->angular_velocity -= ra.cross(friction) * a->inverse_inertia;
                }
                if (b->dynamic) {
                    *b->velocity += friction * b->inverse_mass;
                    *b->angular_velocity += rb.cross(friction) * b->inverse_inertia;
                }
            }
        }
        ++stats.impacts;
        if (impact)
            impact(*fa, *fb, *first);
    }
    return stats;
}
} // namespace butter::physics2d
