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
    float normal_speed{}, angular_a{}, angular_b{}, separation{};
};
struct CcdStatistics {
    int impacts{};
    std::size_t sweeps{};
    std::size_t persistent_contacts{}; // Whole remaining trajectory certified on a contact plane.
    std::size_t bounds_tests{};
    // Overlap tests against stationary geometry, which is enumerated through its
    // own index rather than the mover sweep-and-prune. Reported separately
    // because it is the count that has to stay local when a scene adds walls.
    std::size_t stationary_tests{};
    bool limited{};
    float remaining_time{}; // Largest locally clamped interval, not discarded
                            // world time.
    float advanced_time{};
    std::size_t budget_exhaustions{}, non_convergences{}, zero_time_repeats{}, candidates{};
    std::vector<CcdDiagnostic> diagnostics;

    // Fold one window's counters into a running total. The continuous pass runs
    // once per integration window, so with sub-stepping enabled, publishing the
    // last window on its own would report a fraction of the frame's work as if
    // it were all of it. `advanced_time` sums because each window advanced its
    // own share of the frame; `limited` is an "any window needed the clamp"
    // flag, and the widest locally clamped interval is the one worth reporting.
    void accumulate(const CcdStatistics &window) {
        impacts += window.impacts;
        sweeps += window.sweeps;
        persistent_contacts += window.persistent_contacts;
        bounds_tests += window.bounds_tests;
        stationary_tests += window.stationary_tests;
        limited = limited || window.limited;
        remaining_time = std::max(remaining_time, window.remaining_time);
        advanced_time += window.advanced_time;
        budget_exhaustions += window.budget_exhaustions;
        non_convergences += window.non_convergences;
        zero_time_repeats += window.zero_time_repeats;
        candidates += window.candidates;
        diagnostics.insert(diagnostics.end(), window.diagnostics.begin(),
                           window.diagnostics.end());
    }
};

namespace ccd_detail {
// Conservative advancement needs a support function and a bounded projection,
// and a capsule has both in closed form: it is a segment plus a radius. It used
// to be missing from this list entirely, which meant a capsule was simply not
// swept -- a fast capsule passed straight through a thin floor while the
// discrete narrow phase happily generated its resting contact. Listing it here
// requires the three functions below to agree with `capsule_segment`, which
// they now do.
inline bool supported(const Shape &s) {
    return std::holds_alternative<Circle>(s) || std::holds_alternative<Box>(s) ||
           std::holds_alternative<Capsule>(s) ||
           (std::holds_alternative<Polygon>(s) && std::get<Polygon>(s).vertices.size() >= 3);
}
inline float radius(const Shape &s) {
    if (auto *c = std::get_if<Circle>(&s))
        return c->radius;
    if (auto *box = std::get_if<Box>(&s))
        return box->half_extents.length();
    if (auto *capsule = std::get_if<Capsule>(&s))
        return capsule->half_length + capsule->radius;
    float result = 0;
    const shape_detail::WorldVertices vertices(s, {});
    for (auto v : vertices.view())
        result = std::max(result, v.length());
    return result;
}
inline Vec2 support(const Shape &s, const Transform &t, Vec2 axis) {
    if (auto *c = std::get_if<Circle>(&s))
        return t.position + axis * c->radius;
    if (auto *capsule = std::get_if<Capsule>(&s))
        return capsule_support(*capsule, t, axis);
    const shape_detail::WorldVertices storage(s, t);
    const auto vertices = storage.view();
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
    // A capsule contributes its segment, not a tessellation. Building the
    // eighteen-gon here cost a heap allocation on every separation evaluation --
    // the innermost loop of conservative advancement -- and its facets are not
    // faces of the capsule anyway.
    std::array<Vec2, 2> capsule_a{}, capsule_b{};
    std::optional<shape_detail::WorldVertices> storage_a, storage_b;
    auto vertices_of = [](const Shape &shape, const Transform &t, auto &storage,
                          std::array<Vec2, 2> &capsule_storage) -> std::span<const Vec2> {
        if (const auto *capsule = std::get_if<Capsule>(&shape)) {
            const CapsuleSegment segment = capsule_segment(*capsule, t);
            capsule_storage = {segment.a, segment.b};
            return {capsule_storage.data(), capsule_storage.size()};
        }
        storage.emplace(shape, t);
        return storage->view();
    };
    const auto va = vertices_of(a, ta, storage_a, capsule_a);
    const auto vb = vertices_of(b, tb, storage_b, capsule_b);
    auto cached_support = [](const Shape &shape, const Transform &t,
                             std::span<const Vec2> vertices, Vec2 n) {
        if (auto *circle = std::get_if<Circle>(&shape))
            return t.position + n * circle->radius;
        // Exact, and deliberately *not* the vertex list: `vertices` for a capsule
        // is an eighteen-gon, and taking its support would put the separating
        // plane on a facet instead of on the round end the direction points at.
        if (auto *capsule = std::get_if<Capsule>(&shape))
            return capsule_support(*capsule, t, n);
        float best = -std::numeric_limits<float>::infinity();
        Vec2 point = t.position;
        for (auto v : vertices)
            if (v.dot(n) > best) {
                best = v.dot(n);
                point = v;
            }
        return point;
    };
    auto axis_test = [&](Vec2 axis) {
        if (axis.length_squared() < 1.0e-16f)
            return;
        axis = axis.normalized();
        for (Vec2 n : {axis, -axis}) {
            Vec2 pa = cached_support(a, ta, va, n), pb = cached_support(b, tb, vb, -n);
            float gap = (pb - pa).dot(n);
            if (gap > result.distance)
                result = {gap, n, (pa + pb) * 0.5f};
        }
    };
    auto faces = [&](std::span<const Vec2> vertices) {
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            Vec2 edge = vertices[(i + 1) % vertices.size()] - vertices[i];
            axis_test({-edge.y, edge.x});
        }
    };
    faces(va);
    faces(vb);
    // A capsule's rounded ends are not described by any face normal of either
    // shape: the separating plane for an end cap against a vertex runs through
    // the closest pair, not along an edge. `capsule_convex` reaches the same
    // answer by measuring distance instead, and these axes are what make the
    // separating-plane form agree with it, so the CCD gate and the discrete
    // contact cannot disagree about the same pair.
    auto capsule_axes = [&](const Shape &capsule_shape, const Transform &ct,
                            const Shape &other_shape, const Transform &ot,
                            std::span<const Vec2> other_vertices) {
        const auto *capsule = std::get_if<Capsule>(&capsule_shape);
        if (!capsule)
            return;
        const CapsuleSegment segment = capsule_segment(*capsule, ct);
        const Vec2 along = segment.b - segment.a;
        // The capsule's own flat sides.
        axis_test({-along.y, along.x});
        if (const auto *other_capsule = std::get_if<Capsule>(&other_shape)) {
            // Two rounded interiors: the perpendicular of the other segment is
            // the remaining face of the Minkowski difference.
            const Vec2 other_along = capsule_segment(*other_capsule, ot).b -
                                     capsule_segment(*other_capsule, ot).a;
            axis_test({-other_along.y, other_along.x});
            return;
        }
        if (std::holds_alternative<Circle>(other_shape)) {
            // A circle has no vertices, so the per-vertex loop below would leave
            // the direction from its centre to the barrel untested -- and that is
            // the one the support function maximises for a barrel resting against
            // a round side. Without it `separation` could report a gap smaller
            // than `capsule_circle` computes for the same pose, and the CCD gate
            // and the discrete contact would disagree.
            axis_test(ot.position - closest_on_segment(ot.position, segment.a, segment.b));
            return;
        }
        for (auto v : other_vertices)
            axis_test(v - closest_on_segment(v, segment.a, segment.b));
    };
    capsule_axes(a, ta, b, tb, vb);
    capsule_axes(b, tb, a, ta, va);
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
    auto feature = [&](const Shape &shape, const Transform &t, std::span<const Vec2> vertices,
                       Vec2 n) {
        Vec2 p = cached_support(shape, t, vertices, n);
        float lo = std::numeric_limits<float>::infinity(), hi = -lo;
        if (std::holds_alternative<Circle>(shape))
            return std::make_pair(p.dot(tangent), p.dot(tangent));
        // A capsule's support is `radius` outside its face list, so the plane
        // through `p` is not the plane its own endpoints sit on. Comparing
        // without the offset left *no* vertex on the support plane for a barrel
        // at rest -- the whole tangent span came back as infinity.
        const auto *capsule = std::get_if<Capsule>(&shape);
        const float offset = capsule ? capsule->radius : 0.0f;
        for (auto v : vertices)
            if (p.dot(n) - v.dot(n) - offset < 1.0e-5f) {
                lo = std::min(lo, v.dot(tangent));
                hi = std::max(hi, v.dot(tangent));
            }
        return std::make_pair(lo, hi);
    };
    const auto [alo, ahi] = feature(a, ta, va, result.normal);
    const auto [blo, bhi] = feature(b, tb, vb, -result.normal);
    const float along = (std::max(alo, blo) + std::min(ahi, bhi)) * 0.5f;
    result.point = result.normal * result.point.dot(result.normal) + tangent * along;
    return result;
}
// Bounds every projected vertex over the entire angular interval, including
// interior extrema. Endpoint-only checks would miss a rotating shape's return.
inline std::pair<double, double> projected_motion(const Shape &shape, const ShapeSweep &sweep,
                                                  Vec2 axis, Vec2 origin, Vec2 drift) {
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    const double start = sweep.start.angle, angle = double(sweep.end.angle) - start;
    const double linear = double(drift.x) * axis.x + double(drift.y) * axis.y;
    const double center = double(sweep.start.position.x - origin.x) * axis.x +
                          double(sweep.start.position.y - origin.y) * axis.y;
    auto vertex = [&](Vec2 q) {
        double x = double(axis.x) * q.x + double(axis.y) * q.y;
        double y = double(axis.y) * q.x - double(axis.x) * q.y;
        auto sample = [&](double t) {
            double theta = start + angle * t;
            double value = center + linear * t + x * std::cos(theta) + y * std::sin(theta);
            lo = std::min(lo, value);
            hi = std::max(hi, value);
        };
        sample(0);
        sample(1);
        constexpr double pi = 3.14159265358979323846;
        if (std::abs(angle) > 2 * pi) {
            // Constant-time conservative bound for arbitrarily many revolutions.
            double r = std::hypot(x, y);
            lo = std::min(lo, center + std::min(0.0, linear) - r);
            hi = std::max(hi, center + std::max(0.0, linear) + r);
            return;
        }
        double ca = angle * y, sa = -angle * x, r = std::hypot(ca, sa);
        if (r == 0 || std::abs(linear) > r)
            return;
        double phase = std::atan2(sa, ca), root = std::acos(std::clamp(-linear / r, -1.0, 1.0));
        double begin = std::min(start, start + angle), end = std::max(start, start + angle);
        for (double base : {phase - root, phase + root}) {
            double theta = base + std::ceil((begin - base) / (2 * pi)) * (2 * pi);
            for (int i = 0; i < 2 && theta <= end; ++i, theta += 2 * pi) {
                double t = (theta - start) / angle;
                if (t >= 0 && t <= 1)
                    sample(t);
            }
        }
    };
    if (auto *circle = std::get_if<Circle>(&shape)) {
        vertex(sweep.local.position);
        lo -= circle->radius;
        hi += circle->radius;
    } else if (auto *capsule = std::get_if<Capsule>(&shape)) {
        // Same shape as the circle case, for two points instead of one: bound the
        // segment over the whole sweep, then inflate by the radius. Tessellating
        // the capsule here instead would inflate a polygon that is inscribed in
        // the real one, so the certified interval would be too narrow and the
        // sweep could report a miss for a capsule that does touch.
        const Vec2 axis = rotate({0, capsule->half_length}, sweep.local.angle);
        vertex(sweep.local.position - axis);
        vertex(sweep.local.position + axis);
        lo -= capsule->radius;
        hi += capsule->radius;
    } else {
        const shape_detail::WorldVertices vertices(shape, sweep.local);
        for (auto q : vertices.view())
            vertex(q);
    }
    return {lo, hi};
}
inline bool separated_during_sweep(const Shape &a, const ShapeSweep &sa, const Shape &b,
                                   const ShapeSweep &sb, Vec2 normal, float tolerance) {
    // Remove common translation so a moving pair is not rejected merely because
    // its endpoints share a large displacement in world space.
    auto pa = projected_motion(a, sa, normal, sa.start.position, {});
    Vec2 relative = (sb.end.position - sb.start.position) - (sa.end.position - sa.start.position);
    auto pb = projected_motion(b, sb, normal, sa.start.position, relative);
    const double roundoff = 8 * std::numeric_limits<float>::epsilon() *
                            std::max({1.0, std::abs(pa.second), std::abs(pb.first)});
    return std::isfinite(pa.second) && std::isfinite(pb.first) &&
           pb.first - pa.second >= -double(tolerance) - roundoff;
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
    Vec2 local_center{}; // Centroid offset in body space; pure rotation pivots here.
    float inverse_mass{}, inverse_inertia{};
    bool dynamic{}, kinematic{}, bullet{};
    std::vector<CcdCollider> colliders;
    bool blocked{}; // Local conservative clamp, reset for each advance call.
    bool moves() const { return !blocked && (kinematic || (dynamic && (!sleeping || !*sleeping))); }
    Vec2 linear() const { return moves() ? *velocity : Vec2{}; }
    float angular() const { return moves() ? *angular_velocity : 0; }
    // Terminal transform after `duration`. A rigid body rotates about its
    // centroid, so the body origin is recovered from the advanced centroid.
    Transform advanced(float duration) const {
        Transform result = *transform;
        const float angle = transform->angle + angular() * duration;
        if (local_center.x == 0 && local_center.y == 0)
            result.position = transform->position + linear() * duration;
        else {
            const Vec2 center = transform->position + rotate(local_center, Rot(transform->angle));
            result.position = center + linear() * duration - rotate(local_center, Rot(angle));
        }
        result.angle = angle;
        return result;
    }
    void integrate(float duration) { *transform = advanced(duration); }
    ShapeSweep sweep(float dt, const CcdCollider &f) const { return {*transform, advanced(dt), f.local}; }
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

// Caller-owned scratch storage; never shared between worlds or nested advances.
// Contents are reset per call, while vector capacity survives across steps.
struct CcdWorkspace {
    struct Hit {
        std::size_t i, j, fi, fj;
        SweepHit hit;
    };
    struct Bounds {
        std::size_t i;
        AABB box;
    };
    std::vector<int> counts;
    std::map<std::array<std::size_t, 4>, int> zero_hits;
    std::vector<Bounds> bounds;
    std::vector<Hit> hits;
    std::vector<std::size_t> statics;
    // The split of the group into the bodies that can move and the geometry
    // that cannot. It is derived once per call rather than rediscovered inside
    // every advance iteration, and it is where a change in the split is
    // noticed: a body that was created, destroyed or flipped between static and
    // dynamic changes the count, which invalidates the cached index below
    // without the caller having to report it.
    std::vector<std::size_t> movers;
    // Bounds of the bodies that cannot move, sorted on one axis. A stationary
    // body's swept extent does not depend on how much of the step is left, and
    // in a scene of walls and terrain it does not change between steps either,
    // so it is built once and reused. Rebuilding it inside every advance
    // iteration was the entire reason a world with eight thousand static
    // obstacles spent most of a 300-frame run inside CCD.
    std::vector<Bounds> stationary;
    // The entries of `stationary` the mover sweep line currently touches.
    std::vector<std::size_t> stationary_active;
    bool stationary_valid{false};
    // Axis `stationary` is ordered on, or -1 when it is not ordered at all.
    // The mover sweep picks the axis with the wider spread, which can differ
    // from the stationary set's own, so the order is redone when it does.
    int stationary_axis{-1};
    std::size_t stationary_seen{0};

    // --- Swept-island partition ---
    // `active` is what actually moves this step; `sleepers` is the dynamics it
    // is not moving, which only a bullet or a kinematic body can still reach.
    // A sleeping body cannot be left out of the picture the way a wall can, so
    // it is exactly the question of who can still reach one that decides whether
    // the partition is taken at all.
    std::vector<std::size_t> active, sleepers;
    std::vector<Bounds> sleeper_bounds;
    std::vector<std::size_t> sleeper_line;
    // One swept bound per mover, over the whole step rather than the remaining
    // part of it. It only decides which movers can possibly interact, and that
    // answer must not shrink as impacts are resolved.
    std::vector<Bounds> island_bounds;
    // Disjoint set over `island_bounds`, plus the sweep-and-prune that decides
    // which movers are even worth unioning.
    std::vector<std::size_t> island_parent, island_order, island_line;
    // The partition, indexed by body so that the group scheduler can read it
    // without knowing how the islands were discovered, and laid out as
    // island -> its members so that neither the swept bounds nor the integration
    // has to walk the whole group once per island.
    std::vector<std::size_t> island_of, island_root, body_island, island_offset, island_cursor,
        island_members;
};

// Optional swept-island partition of a group's movers. Bodies carrying the same
// island number are solved on one timeline; different numbers never see each
// other at all, so the group is solved one island at a time instead of advancing
// every mover to the earliest impact and then rescanning all of them. An empty
// partition means one island, which is what a caller that does not partition
// gets, and it leaves the sweep exactly as it was.
//
// `vertical` is the ordering axis for the whole partition. Choosing it per
// island would let the stationary index be re-sorted once per island, which is
// the one thing that has to stay proportional to the world rather than to the
// number of islands.
struct CcdIslands {
    const std::size_t *of{};
    std::size_t islands{};
    bool vertical{};
};

inline CcdStatistics advance_continuous_group(std::vector<CcdMotion> &bodies, float dt,
                                              const CcdSettings &settings = {},
                                              const CcdFilter &filter = {},
                                              const CcdImpact &impact = {},
                                              CcdWorkspace *workspace = nullptr,
                                              bool stationary_changed = true,
                                              const CcdIslands &islands = {}) {
    CcdWorkspace local;
    auto &scratch = workspace ? *workspace : local;
    CcdStatistics stats;
    if (!std::isfinite(dt) || dt <= 0)
        return stats;
    const bool partitioned = islands.of != nullptr && islands.islands != 0;
    const std::size_t island_count = partitioned ? islands.islands : 1;
    // A body's island, or 0 when the caller did not partition.
    auto island_of_body = [&](std::size_t i) {
        return partitioned ? islands.of[i] : std::size_t(0);
    };
    // ... and once that split is known, each island's movers in one contiguous
    // run.
    auto &island_offset = scratch.island_offset;
    auto &island_members = scratch.island_members;
    auto &movers = scratch.movers;
    {
        // One pass over the group separates what can move from what cannot, and
        // notices a change in that split. Nothing moves a body between the two
        // lists inside an advance, so the list survives the whole call.
        movers.clear();
        std::size_t immovable = 0;
        for (std::size_t i = 0; i < bodies.size(); ++i)
            if (bodies[i].dynamic || bodies[i].kinematic)
                movers.push_back(i);
            else
                ++immovable;
        if (immovable != scratch.stationary_seen)
            scratch.stationary_valid = false;
    }
    {
        // Testing every mover against every island instead would cost a pass
        // over the group per island, which is the same rescans the partition
        // exists to remove.
        island_offset.assign(island_count + 1, 0);
        for (std::size_t index : movers)
            ++island_offset[island_of_body(index) + 1];
        for (std::size_t id = 0; id < island_count; ++id)
            island_offset[id + 1] += island_offset[id];
        island_members.resize(movers.size());
        auto &cursor = scratch.island_cursor;
        cursor.assign(island_offset.begin(), island_offset.end() - 1);
        for (std::size_t index : movers)
            island_members[cursor[island_of_body(index)]++] = index;
    }
    auto advance = [&](std::size_t island, float duration) {
        for (std::size_t k = island_offset[island]; k < island_offset[island + 1]; ++k)
            if (bodies[island_members[k]].moves())
                bodies[island_members[k]].integrate(duration);
    };
    for (auto &b : bodies)
        b.blocked = false;
    if (!settings.enabled) {
        advance(0, dt);
        stats.advanced_time = dt;
        return stats;
    }
    auto &counts = scratch.counts;
    counts.assign(bodies.size(), 0);
    auto &zero_hits = scratch.zero_hits;
    zero_hits.clear();
    auto &bounds = scratch.bounds;
    bounds.reserve(bodies.size());
    auto &stationary = scratch.stationary;
    if (stationary_changed)
        scratch.stationary_valid = false;
    auto &hits = scratch.hits;
    hits.reserve(bodies.size());
    auto sort_stationary = [&](bool vertical) {
        std::sort(stationary.begin(), stationary.end(), [vertical](const auto &a, const auto &b) {
            const float la = vertical ? a.box.min.y : a.box.min.x;
            const float lb = vertical ? b.box.min.y : b.box.min.x;
            return la == lb ? a.i < b.i : la < lb;
        });
        scratch.stationary_axis = vertical ? 1 : 0;
    };
    // Bodies that cannot move at all are indexed separately. Their bound is the
    // same in every advance iteration, and behind a wall or terrain it is the
    // same across steps, so rebuilding them on the mover sweep-and-prune was
    // pure waste: with eight thousand parked platforms the mover sweep spent
    // its time re-deriving the bounds of obstacles it could not touch.
    auto stationary_index = [&] {
        if (scratch.stationary_valid)
            return;
        stationary.clear();
        scratch.stationary_seen = 0;
        Vec2 low{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()};
        Vec2 high{-low.x, -low.y};
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            auto &x = bodies[i];
            if (x.dynamic || x.kinematic)
                continue;
            ++scratch.stationary_seen;
            bool initialized = false;
            AABB total{};
            for (auto &f : x.colliders) {
                auto box = compute_aabb(*f.shape, x.sweep(0.0f, f).at(0));
                const float reach = x.local_center.length();
                box.min -= Vec2{reach, reach};
                box.max += Vec2{reach, reach};
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
            if (!initialized)
                continue;
            stationary.push_back({i, total});
            Vec2 center = (total.min + total.max) * .5f;
            low.x = std::min(low.x, center.x);
            low.y = std::min(low.y, center.y);
            high.x = std::max(high.x, center.x);
            high.y = std::max(high.y, center.y);
        }
        // Ordered on the axis the stationary set itself spreads along, which is
        // stable for as long as the set is, so the sort is paid once.
        sort_stationary(high.y - low.y > high.x - low.x);
        scratch.stationary_valid = true;
    };
    // Solve one island at a time. Every mover belongs to exactly one island,
    // and no two islands can see each other, so each is carried to its own end
    // of step. Advancing the whole group to the earliest impact and then
    // rescanning every mover is what made an unrelated bullet expensive.
    for (std::size_t island = 0; island < island_count; ++island) {
        float remaining = dt;
        zero_hits.clear();

        while (remaining > 0) {
            hits.clear();
            bounds.clear();
            float earliest = 1;
            // Sweep-and-prune over the movers; their swept extent depends on how
            // much of the step is left, so this is rebuilt after every impulse.
            for (std::size_t m = island_offset[island]; m < island_offset[island + 1]; ++m) {
                const std::size_t index = island_members[m];
                auto &x = bodies[index];
                bool initialized = false;
                AABB total{};
                for (auto &f : x.colliders) {
                    auto sweep = x.sweep(remaining, f);
                    // Rotation pivots on the centroid, so the bound must also cover
                    // the lever from the body origin to the centroid.
                    const float reach = x.local_center.length();
                    auto box = std::abs(sweep.end.angle - sweep.start.angle) < 1e-8f
                                   ? compute_aabb(*f.shape, sweep.at(0))
                                   : ccd_detail::swept_bounds(*f.shape, sweep);
                    box.min -= Vec2{reach, reach};
                    box.max += Vec2{reach, reach};
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
                    bounds.push_back({index, total});
            }
            // Select the axis with the greater center spread. A tall stack should
            // not enumerate every pair merely because their X projections overlap.
            Vec2 low{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()};
            Vec2 high{-low.x, -low.y};
            for (auto &entry : bounds) {
                Vec2 center = (entry.box.min + entry.box.max) * .5f;
                low.x = std::min(low.x, center.x);
                low.y = std::min(low.y, center.y);
                high.x = std::max(high.x, center.x);
                high.y = std::max(high.y, center.y);
            }
            // The partition orders every island on one axis, so the stationary
            // index below is sorted once for the whole call rather than once per
            // island. That is the one cost that must stay proportional to the
            // world instead of to the number of islands.
            bool vertical = partitioned ? islands.vertical
                                        : high.y - low.y > high.x - low.x;
            auto lower = [&](const AABB &b) { return vertical ? b.min.y : b.min.x; };
            auto upper = [&](const AABB &b) { return vertical ? b.max.y : b.max.x; };
            std::sort(bounds.begin(), bounds.end(), [&](auto &a, auto &b) {
                return lower(a.box) == lower(b.box) ? a.i < b.i : lower(a.box) < lower(b.box);
            });
            // One body pair, from the bounds that let it through to the impulse it
            // produces. The two bodies can come from either list -- two movers, or
            // a mover against the immovable geometry -- so this takes indices
            // rather than sweep entries.
            auto test_pair = [&](std::size_t i, std::size_t j) {
                    auto &x = bodies[i];
                    auto &y = bodies[j];
                    if ((!x.dynamic && !y.dynamic) || (!x.moves() && !y.moves()) ||
                        (x.dynamic && y.dynamic && !x.bullet && !y.bullet))
                        return;
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
                            if (!hit && ccd_detail::supported(*fx.shape) &&
                                ccd_detail::supported(*fy.shape)) {
                                auto initial = ccd_detail::separation(*fx.shape, x.sweep(0, fx).at(0),
                                                                      *fy.shape, y.sweep(0, fy).at(0));
                                // A discrete resting overlap must not disable motion
                                // protection under pressure from the rest of a stack.
                                if (initial.distance < -settings.tolerance)
                                    hit = SweepHit{0, initial.normal, initial.point, true};
                            }
                            if (!hit)
                                continue;
                            auto start_separation = ccd_detail::separation(
                                *fx.shape, x.sweep(0, fx).at(0), *fy.shape, y.sweep(0, fy).at(0));
                            const float contact_depth =
                                hit->fraction == 0
                                    ? std::max(settings.tolerance, -start_separation.distance)
                                    : settings.tolerance;
                            {
                                if (ccd_detail::separated_during_sweep(
                                        *fx.shape, x.sweep(remaining, fx), *fy.shape,
                                        y.sweep(remaining, fy), hit->normal, contact_depth)) {
                                    ++stats.persistent_contacts;
                                    continue;
                                }
                            }
                            if (hit->fraction == 0 || !hit->converged) {
                                // The current corner may be separating while another
                                // corner returns later. Advance only a certified prefix
                                // instead of repeatedly resolving the current t=0 point.
                                auto sx = x.sweep(remaining, fx), sy = y.sweep(remaining, fy);
                                auto prefix = [](ShapeSweep sweep, float t) {
                                    sweep.end = {sweep.start.position +
                                                     (sweep.end.position - sweep.start.position) * t,
                                                 sweep.start.angle +
                                                     (sweep.end.angle - sweep.start.angle) * t};
                                    return sweep;
                                };
                                const float initial_gap =
                                    ccd_detail::separation(*fx.shape, sx.at(0), *fy.shape, sy.at(0))
                                        .distance;
                                // Re-entry uses an inner tolerance, leaving room for
                                // roundoff before the persistent-contact certificate.
                                const float prefix_tolerance =
                                    std::min(contact_depth,
                                             .5f * (contact_depth + std::max(0.0f, -initial_gap)));
                                float lo = 0, hi = 1;
                                for (int iteration = 0;
                                     iteration < std::min(24, std::max(1, settings.max_iterations));
                                     ++iteration) {
                                    float mid = (lo + hi) * .5f;
                                    if (ccd_detail::separated_during_sweep(
                                            *fx.shape, prefix(sx, mid), *fy.shape, prefix(sy, mid),
                                            hit->normal, prefix_tolerance))
                                        lo = mid;
                                    else
                                        hi = mid;
                                }
                                // Stay inside the proved interval without inventing an early
                                // collision a macroscopic distance from the surface.
                                float fraction = std::nextafter(lo, 0.0f);
                                if (fraction > 1e-6f) {
                                    auto next = ccd_detail::separation(*fx.shape, sx.at(fraction),
                                                                       *fy.shape, sy.at(fraction));
                                    if (next.distance <= 2 * settings.tolerance)
                                        hit = SweepHit{fraction, next.normal, next.point, true};
                                }
                            }
                            earliest = std::min(earliest, hit->fraction);
                            hits.push_back({i, j, fi, fj, *hit});
                        }
            };
            // Two movers. The list is ordered on the axis with the wider spread, so
            // a mover is only paired with the entries whose projection reaches it.
            for (std::size_t bi = 0; bi < bounds.size(); ++bi)
                for (std::size_t bj = bi + 1;
                     bj < bounds.size() && lower(bounds[bj].box) <= upper(bounds[bi].box); ++bj) {
                    ++stats.bounds_tests;
                    if (!bounds[bi].box.overlaps(bounds[bj].box))
                        continue;
                    test_pair(std::min(bounds[bi].i, bounds[bj].i),
                              std::max(bounds[bi].i, bounds[bj].i));
                }
            // And the movers against the immovable geometry. This is the pairing
            // that has to stay local: a scene of walls, platforms and terrain is
            // mostly obstacles no mover can reach, and the answer to "which ones can
            // this mover reach" must not be "all of them, every iteration". Both
            // lists are ordered on the same axis, so this is a merge that carries an
            // active set between movers instead of a scan per mover.
            stationary_index();
            if (scratch.stationary_axis != (vertical ? 1 : 0))
                sort_stationary(vertical);
            {
                auto &active = scratch.stationary_active;
                active.clear();
                std::size_t cursor = 0;
                for (auto &mover : bounds) {
                    const float lo = lower(mover.box), hi = upper(mover.box);
                    // The movers are visited in ascending lower edge, so an entry
                    // whose upper edge is already behind the sweep line will never
                    // be reached again.
                    for (std::size_t k = 0; k < active.size();) {
                        if (upper(stationary[active[k]].box) < lo) {
                            active[k] = active.back();
                            active.pop_back();
                        } else
                            ++k;
                    }
                    while (cursor < stationary.size() && lower(stationary[cursor].box) <= hi)
                        active.push_back(cursor++);
                    for (std::size_t slot : active) {
                        // This is the counter that must not follow the size of the
                        // world: it is the number of parked obstacles a mover is
                        // measured against, and a wall a hundred metres away may
                        // not appear in it.
                        ++stats.stationary_tests;
                        ++stats.bounds_tests;
                        if (!mover.box.overlaps(stationary[slot].box))
                            continue;
                        const std::size_t other = stationary[slot].i;
                        // The cached index is only trusted while the split it was
                        // built from still holds. Bodies cannot change list inside
                        // an advance, but the cache can predate a type flip, so an
                        // entry that is no longer immovable is dropped instead of
                        // being tested -- against a mover it would be a duplicate
                        // pair, and against a mover it never was one.
                        if (bodies[other].dynamic || bodies[other].kinematic)
                            continue;
                        test_pair(std::min(mover.i, other), std::max(mover.i, other));
                    }
                }
            }
            if (hits.empty()) {
                // No impact left in this island: give it the rest of the step
                // and move on. Other islands still have their own timeline.
                advance(island, remaining);
                break;
            }
            if (earliest > 0)
                zero_hits.clear();
            advance(island, remaining * earliest);
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
                    Vec2 ra = first->point - a->transform->position,
                         rb = first->point - b->transform->position;
                    Vec2 rv = b->linear() + Vec2{-b->angular() * rb.y, b->angular() * rb.x} -
                              a->linear() - Vec2{-a->angular() * ra.y, a->angular() * ra.x};
                    float gap = ccd_detail::separation(*fa->shape, a->sweep(0, *fa).at(0), *fb->shape,
                                                       b->sweep(0, *fb).at(0))
                                    .distance;
                    stats.diagnostics.push_back({reason, entry.i, entry.j, entry.fi, entry.fj,
                                                 dt - remaining, remaining, rv.dot(first->normal),
                                                 a->angular(), b->angular(), gap});
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
    }

    stats.advanced_time = dt;
    return stats;
}
// Two movers only have to be solved on one timeline when CCD looks at the pair
// and one of them can be pushed by it. A kinematic body is prescribed motion and
// absorbs nothing, so although everything it sweeps through has to be in its
// island, the box it knocks aside cannot turn around and change where the
// platform goes; ordinary dynamics are not even tested against each other, so
// only a bullet ties a dynamic body to the rest of the world.
//
// The movers are therefore partitioned on exactly that relation, and the group
// scheduler -- which already owns the stationary index and the sweep-and-prune
// -- solves one island at a time. Without the partition a scene containing a
// single bullet was forced onto one timeline in its entirety, even when the
// bullet was nowhere near anything: the coupled scheduler advances every mover
// to the earliest impact and then rescans all of them, so forty unrelated
// impacts cost forty rescans of every moving body.
inline CcdStatistics advance_continuous(std::vector<CcdMotion> &bodies, float dt,
                                        const CcdSettings &settings = {},
                                        const CcdFilter &filter = {}, const CcdImpact &impact = {},
                                        CcdWorkspace *workspace = nullptr,
                                        bool stationary_changed = true) {
    CcdWorkspace local;
    auto &scratch = workspace ? *workspace : local;
    CcdStatistics total;
    if (!std::isfinite(dt) || dt <= 0)
        return total;
    auto &statics = scratch.statics;
    statics.clear();
    auto &sleepers = scratch.sleepers;
    sleepers.clear();
    auto &active = scratch.active;
    active.clear();
    bool hazard = false;       // A mover that can reach geometry which is not.
    bool sleeper_bullet = false;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        auto &b = bodies[i];
        b.blocked = false;
        if (!b.dynamic && !b.kinematic) {
            statics.push_back(i);
            continue;
        }
        if (b.moves())
            active.push_back(i);
        else if (b.dynamic) {
            sleepers.push_back(i);
            sleeper_bullet |= b.bullet;
            continue;
        }
        hazard |= b.bullet || b.kinematic;
    }
    // Every island is a group of its own, so the stationary index -- which holds
    // group indices -- is rebuilt for each one. Past a handful of immovable
    // bodies that stops being a trade worth making: with walls and terrain the
    // group scheduler builds that index once and reuses it, and that is the
    // right tool there.
    if (!settings.enabled || statics.size() > 16)
        return advance_continuous_group(bodies, dt, settings, filter, impact, &scratch,
                                        stationary_changed);
    if (active.empty()) {
        total.advanced_time = dt;
        return total;
    }
    // Swept bound of every mover over the whole step. A mover with no geometry
    // gets an inverted box: it cannot interact with anything, but it still has
    // to be carried by a group so that its motion is integrated.
    auto &bounds = scratch.island_bounds;
    bounds.clear();
    for (std::size_t i : active) {
        auto &x = bodies[i];
        AABB total_box{};
        bool initialized = false;
        for (auto &f : x.colliders) {
            auto sweep = x.sweep(dt, f);
            // Rotation pivots on the centroid, so the bound must also cover the
            // lever from the body origin to the centroid.
            const float reach = x.local_center.length();
            const bool turning = std::abs(sweep.end.angle - sweep.start.angle) >= 1e-8f;
            auto part = turning ? ccd_detail::swept_bounds(*f.shape, sweep)
                                : compute_aabb(*f.shape, sweep.at(0));
            part.min -= Vec2{reach, reach};
            part.max += Vec2{reach, reach};
            if (!turning) {
                auto end = compute_aabb(*f.shape, sweep.at(1));
                part.min.x = std::min(part.min.x, end.min.x);
                part.min.y = std::min(part.min.y, end.min.y);
                part.max.x = std::max(part.max.x, end.max.x);
                part.max.y = std::max(part.max.y, end.max.y);
            }
            if (!initialized) {
                total_box = part;
                initialized = true;
            } else {
                total_box.min.x = std::min(total_box.min.x, part.min.x);
                total_box.min.y = std::min(total_box.min.y, part.min.y);
                total_box.max.x = std::max(total_box.max.x, part.max.x);
                total_box.max.y = std::max(total_box.max.y, part.max.y);
            }
        }
        if (!initialized)
            total_box = {{std::numeric_limits<float>::infinity(),
                          std::numeric_limits<float>::infinity()},
                         {-std::numeric_limits<float>::infinity(),
                          -std::numeric_limits<float>::infinity()}};
        bounds.push_back({i, total_box});
    }
    const std::size_t movers = bounds.size();
    auto &parent = scratch.island_parent;
    parent.resize(movers);
    for (std::size_t k = 0; k < movers; ++k)
        parent[k] = k;
    auto find = [&](std::size_t k) {
        while (parent[k] != k) {
            parent[k] = parent[parent[k]];
            k = parent[k];
        }
        return k;
    };
    // Which movers are worth unioning at all. Sweep and prune rather than the
    // full pair product: a cloud of free-falling debris is a mover list of
    // thousands whose bounds are mutually disjoint, and it must stay linear.
    auto &order = scratch.island_order;
    order.resize(movers);
    for (std::size_t k = 0; k < movers; ++k)
        order[k] = k;
    Vec2 low{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()};
    Vec2 high{-low.x, -low.y};
    for (auto &entry : bounds) {
        if (!(entry.box.min.x <= entry.box.max.x)) // Inverted: no geometry.
            continue;
        Vec2 center = (entry.box.min + entry.box.max) * .5f;
        low.x = std::min(low.x, center.x);
        low.y = std::min(low.y, center.y);
        high.x = std::max(high.x, center.x);
        high.y = std::max(high.y, center.y);
    }
    const bool vertical = high.y - low.y > high.x - low.x;
    auto lower_edge = [&](const AABB &b) { return vertical ? b.min.y : b.min.x; };
    auto upper_edge = [&](const AABB &b) { return vertical ? b.max.y : b.max.x; };
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return lower_edge(bounds[a].box) == lower_edge(bounds[b].box)
                   ? a < b
                   : lower_edge(bounds[a].box) < lower_edge(bounds[b].box);
    });
    // A sleeping body can still be reached -- and the impulse that wakes it will
    // throw it somewhere this partition does not describe -- so the partition is
    // only taken while no mover can touch one. Answering that is one sweep over
    // bounds that are fixed for the step, rather than the per-iteration pass the
    // group scheduler makes over the same bodies, so a world that carries a
    // mover which can reach nothing new still gets its islands.
    if (!sleepers.empty() && (hazard || sleeper_bullet)) {
        auto &sleep_bounds = scratch.sleeper_bounds;
        sleep_bounds.clear();
        for (std::size_t i : sleepers) {
            auto &x = bodies[i];
            AABB box{};
            bool initialized = false;
            for (auto &f : x.colliders) {
                auto part = compute_aabb(*f.shape, x.sweep(0.0f, f).at(0));
                if (!initialized) {
                    box = part;
                    initialized = true;
                } else {
                    box.min.x = std::min(box.min.x, part.min.x);
                    box.min.y = std::min(box.min.y, part.min.y);
                    box.max.x = std::max(box.max.x, part.max.x);
                    box.max.y = std::max(box.max.y, part.max.y);
                }
            }
            if (initialized)
                sleep_bounds.push_back({i, box});
        }
        std::sort(sleep_bounds.begin(), sleep_bounds.end(), [&](auto &a, auto &b) {
            return lower_edge(a.box) == lower_edge(b.box) ? a.i < b.i
                                                          : lower_edge(a.box) < lower_edge(b.box);
        });
        auto &sleep_line = scratch.sleeper_line;
        sleep_line.clear();
        std::size_t cursor = 0;
        bool reachable = false;
        for (std::size_t oi = 0; oi < movers && !reachable; ++oi) {
            const std::size_t k = order[oi];
            const float begin = lower_edge(bounds[k].box), end = upper_edge(bounds[k].box);
            for (std::size_t s = 0; s < sleep_line.size();) {
                if (upper_edge(sleep_bounds[sleep_line[s]].box) < begin) {
                    sleep_line[s] = sleep_line.back();
                    sleep_line.pop_back();
                } else
                    ++s;
            }
            while (cursor < sleep_bounds.size() && lower_edge(sleep_bounds[cursor].box) <= end)
                sleep_line.push_back(cursor++);
            auto &mover = bodies[bounds[k].i];
            for (std::size_t s : sleep_line) {
                auto &sleeping = bodies[sleep_bounds[s].i];
                // Both dynamics, neither of them a bullet: the pair test skips
                // it, so a settled stack does not cost anyone its partition.
                if (mover.dynamic && !mover.bullet && !sleeping.bullet)
                    continue;
                if (bounds[k].box.overlaps(sleep_bounds[s].box)) {
                    reachable = true;
                    break;
                }
            }
        }
        if (reachable)
            return advance_continuous_group(bodies, dt, settings, filter, impact, &scratch,
                                            stationary_changed);
    }
    auto &line = scratch.island_line;
    line.clear();
    for (std::size_t oi = 0; oi < movers; ++oi) {
        const std::size_t k = order[oi];
        const float begin = lower_edge(bounds[k].box);
        // The sweep line only moves forward, so an entry whose upper edge is
        // already behind it will never be reached again.
        for (std::size_t s = 0; s < line.size();) {
            if (upper_edge(bounds[line[s]].box) < begin) {
                line[s] = line.back();
                line.pop_back();
            } else
                ++s;
        }
        auto &x = bodies[bounds[k].i];
        for (std::size_t m : line) {
            auto &y = bodies[bounds[m].i];
            // The relation is exactly the one the pair test uses, so that every
            // pair CCD will look at ends up inside a single island: two
            // dynamics only meet when one of them is a bullet, and anything a
            // kinematic body sweeps through is a pair regardless.
            if (x.dynamic && y.dynamic && !x.bullet && !y.bullet)
                continue;
            if (!x.dynamic && !y.dynamic)
                continue;
            if (!bounds[k].box.overlaps(bounds[m].box))
                continue;
            const std::size_t ra = find(k), rb = find(m);
            if (ra != rb)
                parent[ra] = rb;
        }
        line.push_back(k);
    }
    // Label each mover with its island. Ascending order means an island is named
    // after its lowest member, so the numbering does not depend on the order the
    // unions happened to be made in.
    auto &island_of = scratch.island_of;
    island_of.assign(movers, movers);
    auto &root_id = scratch.island_root;
    root_id.assign(movers, movers);
    std::size_t islands = 0;
    for (std::size_t k = 0; k < movers; ++k) {
        const std::size_t r = find(k);
        if (root_id[r] == movers)
            root_id[r] = islands++;
        island_of[k] = root_id[r];
    }
    // Hand the partition to the group scheduler instead of slicing the world
    // into one group per island. The stationary index is built from the bodies
    // that cannot move, and rebuilding it once per island was the entire reason
    // the partition used to be limited to a handful of walls; inside one call it
    // is built once and every island reads it.
    auto &body_island = scratch.body_island;
    body_island.assign(bodies.size(), 0);
    for (std::size_t k = 0; k < movers; ++k)
        body_island[bounds[k].i] = island_of[k];
    total = advance_continuous_group(bodies, dt, settings, filter, impact, &scratch,
                                     stationary_changed,
                                     CcdIslands{body_island.data(), islands, vertical});
    total.advanced_time = dt;
    return total;
}
} // namespace butter::physics2d
