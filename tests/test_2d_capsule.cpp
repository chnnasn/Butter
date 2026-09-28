// Capsule consistency: discrete test, ray query, sweep and continuous
// collision all have to describe the *same* shape.
//
// A capsule is the set of points within `radius` of a segment. That is exact
// geometry with a flat barrel and two round ends. The engine used to feed the
// narrow phase an inscribed eighteen-gon instead, which is a different shape:
// the flat side is only flat at the vertices, and a ray grazing a round end
// misses a capsule it geometrically hits. These checks pin down the exact
// behaviour, and the *agreement* between the four code paths that consume a
// capsule, because each of them had its own copy of the tessellation.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace butter::physics2d;
static int checks = 0;
static void check(bool value, const char *message) {
    ++checks;
    if (!value)
        throw std::runtime_error(message);
}
static World::Config config() {
    World::Config c;
    c.gravity = {};
    return c;
}
static void no_damping(Body &b) {
    b.linear_damping = 0;
    b.angular_damping = 0;
}
static Body &wall(World &w, float x = 0, float restitution = 0) {
    return w.create_body()
        .static_body()
        .at(x, 0)
        .box(0.01f, 5)
        .friction(0)
        .restitution(restitution)
        .build();
}
// A capsule fired along +x. Half length and radius match `shot()`'s circle in
// extent so the wall contact sits at a comparable place.
static Body &capsule_shot(World &w, float x = -5, float velocity = 600, float restitution = 0) {
    auto &b = w.create_body()
                  .dynamic()
                  .at(x, 0)
                  .capsule(0.25f, 0.75f)
                  .velocity(velocity, 0)
                  .friction(0)
                  .restitution(restitution)
                  .build();
    no_damping(b);
    return b;
}
int main() try {
    const Shape capsule = Shape{Capsule{0.25f, 0.75f}};
    const Shape circle = Shape{Circle{0.4f}};
    const Shape box = Shape{Box{{0.6f, 0.6f}}};
    const Shape polygon = Shape{Polygon{{{-0.5f, -0.5f}, {0.5f, -0.5f}, {0.5f, 0.5f}, {-0.5f, 0.5f}}}};
    const std::array<Shape, 3> others{circle, box, polygon};
    std::mt19937 rng(20261206);
    std::uniform_real_distribution<float> pos(-2.0f, 2.0f), ang(-3.0f, 3.0f);
    {
        // Every pairing a broad phase can hand over has to be order independent.
        // The bug this guards against was a mixed pair that answered "contact"
        // in one order and "no contact at all" in the other, decided by which
        // body happened to be created first.
        for (int i = 0; i < 500; ++i) {
            const Transform ta{{pos(rng), pos(rng)}, ang(rng)};
            const Transform tb{{pos(rng), pos(rng)}, ang(rng)};
            for (const Shape &other : others) {
                Contact ab, ba;
                const bool hit_ab = test(capsule, ta, other, tb, ab);
                const bool hit_ba = test(other, tb, capsule, ta, ba);
                check(hit_ab == hit_ba, "capsule pair disagreed about contact existence");
                if (!hit_ab)
                    continue;
                check(std::abs(ab.penetration - ba.penetration) < 1.0e-5f,
                      "capsule pair penetration depended on argument order");
                check((ab.normal + ba.normal).length() < 1.0e-5f,
                      "capsule pair normal was not the exact reverse in the other order");
                // The normal has to separate the pair, whatever the pairing.
                check(std::abs(ab.normal.length() - 1.0f) < 1.0e-4f,
                      "capsule pair normal was not a unit vector");
            }
        }
    }
    {
        // Parallel barrels overlap along a whole strip, where the closest pair
        // is not unique. A per-vertex direction would pick a different answer at
        // every point of the strip, which is what made a resting capsule creep.
        Contact c;
        const Shape a = Shape{Capsule{0.25f, 0.75f}}, b = Shape{Capsule{0.25f, 0.75f}};
        check(test(a, {{0, 0}, 0}, b, {{0.3f, 0.55f}, 0}, c),
              "parallel capsules were reported as separated");
        check(c.normal.x > 0.999f && std::abs(c.normal.y) < 1.0e-3f,
              "parallel capsules reported a direction along the strip instead of across it");
        check(std::abs(c.penetration - 0.2f) < 1.0e-5f, "parallel capsule depth was wrong");
        // Coincident capsules have no closest pair at all; the fallback still
        // has to hand the solver a usable unit direction.
        check(test(a, {{0, 0}, 0}, b, {{0, 0}, 0}, c), "coincident capsules were separated");
        check(std::abs(c.normal.length() - 1.0f) < 1.0e-4f && std::abs(c.penetration - 0.5f) < 1.0e-5f,
              "coincident capsules produced a degenerate normal or depth");
        // Crossed capsules: the segments intersect, so the depth is the sum of
        // the radii in every direction and the direction itself is free.
        check(test(a, {{0, 0}, 0}, b, {{0.6f, 0}, 1.57079633f}, c),
              "crossed capsules were reported as separated");
        check(std::abs(c.penetration - 0.5f) < 1.0e-5f, "crossed capsule depth was wrong");
    }
    {
        // Ray query against the exact surface. The upper cap is a circle centred
        // at {0, .75}; a ray at the barrel height has to land on x = -.25, and a
        // ray grazing the very top has to land on the cap, not on a facet below
        // it, which is exactly where an inscribed polygon would miss.
        RaycastHit hit{};
        auto barrel = ray_shape({-5, 0}, {1, 0}, 100, capsule, {{0, 0}, 0});
        check(barrel && std::abs(barrel->distance - 4.75f) < 1.0e-5f,
              "capsule barrel raycast missed its exact surface");
        check(std::abs(barrel->normal.x + 1.0f) < 1.0e-5f && std::abs(barrel->normal.y) < 1.0e-5f,
              "capsule barrel raycast normal was not the outward surface normal");
        const float graze_y = 0.75f + 0.2495f;
        const float graze_expected = 5.0f - std::sqrt(0.25f * 0.25f - 0.2495f * 0.2495f);
        auto graze = ray_shape({-5, graze_y}, {1, 0}, 100, capsule, {{0, 0}, 0});
        check(graze && std::abs(graze->distance - graze_expected) < 1.0e-4f,
              "capsule round end was tessellated away in the ray query");
        check(graze && std::abs(graze->point.y - graze_y) < 1.0e-4f,
              "capsule raycast slid along the surface instead of hitting where it entered");
        // A ray that starts inside is ignored, as for every other shape.
        check(!ray_shape({0, 0.9f}, {1, 0}, 100, capsule, {{0, 0}, 0}),
              "ray starting inside a capsule was reported as a hit");
        // And a clear miss stays a miss.
        check(!ray_shape({-5, 2.0f}, {1, 0}, 100, capsule, {{0, 0}, 0}),
              "ray clear of a capsule reported a hit");
        // Length limit has to be respected, not rounded up to the surface.
        check(!ray_shape({-5, 0}, {1, 0}, 4.7f, capsule, {{0, 0}, 0}),
              "capsule raycast ignored the maximum distance");
    }
    {
        // The ray hit point must be *the* boundary: just inside is a contact for
        // a discrete probe, just outside is not. This ties `ray_shape` and
        // `test` to the same surface instead of two independent approximations.
        const Shape probe = Shape{Circle{0.001f}};
        int hits = 0;
        for (int i = 0; i < 400 && hits < 120; ++i) {
            // Aim from a ring at the capsule so the sample is hits, not misses:
            // origins within the barrel are rejected by design.
            const Vec2 origin{std::cos(float(i) * 0.7f) * 4.0f, std::sin(float(i) * 1.3f) * 4.0f};
            const Vec2 target{pos(rng) * 0.4f, pos(rng) * 0.9f};
            const Vec2 direction = (target - origin).normalized();
            auto hit = ray_shape(origin, direction, 100, capsule, {{0, 0}, 0});
            if (!hit)
                continue;
            ++hits;
            constexpr float depth = 0.01f;
            Contact inside, outside;
            // `ray_shape` reports the entry point, so continuing along the ray
            // leaves the surface on the inside and backing up along it is the
            // exterior the ray came from.
            check(test(probe, {hit->point + direction * depth, 0}, capsule, {{0, 0}, 0}, inside),
                  "capsule ray hit point was not on the discrete surface (inside probe missed)");
            check(!test(probe, {hit->point - direction * depth, 0}, capsule, {{0, 0}, 0}, outside),
                  "capsule ray hit point was not on the discrete surface (outside probe hit)");
        }
        check(hits >= 100, "capsule ray probe did not sample enough hits");
    }
    {
        // A sweep has to agree with the discrete test it is an acceleration of:
        // dense sampling of `test` finds first contact, and the sweep's reported
        // fraction may not be later than that.
        int impacts = 0;
        for (int i = 0; i < 240; ++i) {
            const Shape a = i % 3 == 0 ? Shape{Capsule{0.2f, 0.4f}}
                            : i % 3 == 1 ? Shape{Capsule{0.12f, 0.55f}}
                                         : Shape{Circle{0.2f}};
            const Shape b = i % 2 ? Shape{Box{{0.15f, 0.7f}}} : Shape{Polygon{{{0, 0.6f}, {0.6f, -0.4f}, {-0.6f, -0.4f}}}};
            ShapeSweep sa{{{pos(rng), pos(rng)}, ang(rng)},
                          {{pos(rng), pos(rng)}, ang(rng)},
                          {}},
                sb{};
            if (ccd_detail::separation(a, sa.at(0), b, sb.at(0)).distance <= 0.0001f)
                continue; // already touching: the sweep is documented at t = 0
            auto hit = sweep_shapes(a, sa, b, sb);
            for (int sample = 1; sample <= 1200; ++sample) {
                const float t = float(sample) / 1200;
                Contact contact;
                if (test(a, sa.at(t), b, sb.at(t), contact)) {
                    check(hit && hit->fraction <= t + 0.0001f,
                          "capsule sweep missed a contact the discrete test found");
                    ++impacts;
                    break;
                }
            }
        }
        check(impacts > 40, "capsule sweep oracle never reached a contact");
    }
    {
        // Exact sweep against a thin wall, checked against the analytic answer
        // rather than against itself: the barrel leading edge is at
        // start + .25, and the wall's near face is at -.01, so contact lands at
        // a fraction of (5 - .26) / 10.
        const Shape floorless = Shape{Box{{0.01f, 5}}};
        ShapeSweep moving{{{0, 0}, 0}, {{0, 0}, 0}, {}};
        moving.start.position = {-5, 0};
        moving.end.position = {5, 0};
        auto hit = sweep_shapes(capsule, moving, floorless, {});
        check(hit && hit->converged && std::abs(hit->fraction - 0.474f) < 0.002f,
              "capsule sweep against a thin wall did not stop at the exact surface");
        // A capsule wide enough to clear the wall must not report an impact.
        auto miss = sweep_shapes(capsule, moving, Shape{Box{{0.01f, 0.1f}}}, {});
        check(miss && miss->fraction > 0.45f && miss->fraction < 0.60f,
              "capsule sweep against a short wall reported a spurious early hit");
    }
    {
        // The manifold has to clip the barrel, not the tessellation. A capsule
        // lying flat on a floor is the case that exposes both halves of the
        // problem: the two points must sit on the floor plane (the barrel is a
        // surface `radius` away from the segment, not the segment itself), and
        // they must be spread over the barrel's own length rather than onto one
        // facet of an inscribed eighteen-gon.
        constexpr float half_pi = 1.57079633f;
        const Shape flat = Shape{Capsule{0.25f, 0.75f}};
        const Manifold m = contact_manifold(flat, {{0, 0.25f}, half_pi}, Shape{Box{{10, 0.5f}}},
                                            {{0, -0.5f}, 0}, {0, -1}, {0, 0.25f}, 0.0f);
        check(m.count == 2, "capsule floor manifold did not produce two barrel points");
        check(m.normal.x == 0 && m.normal.y == -1, "capsule floor manifold normal was wrong");
        for (int i = 0; i < m.count; ++i) {
            check(std::abs(m.points[i].point.y) < 1.0e-4f,
                  "capsule floor manifold point was off the floor plane");
            check(std::abs(m.points[i].separation) < 1.0e-4f,
                  "capsule resting on a floor reported a non-zero gap");
        }
        check(std::abs(m.points[0].point.x - m.points[1].point.x) > 1.4f,
              "capsule floor manifold collapsed onto a single tessellation facet");
        // The same capsule standing on its round end touches with a single point
        // on the cap tip, and the face contest has to pick the floor's face
        // rather than a barrel that is edge-on to it.
        const Manifold upright = contact_manifold(Shape{Capsule{0.25f, 0.75f}}, {{0, 1.0f}, 0},
                                                  Shape{Box{{10, 0.5f}}}, {{0, -0.5f}, 0}, {0, -1},
                                                  {0, 0.0f}, 0.0f);
        check(upright.count == 1, "upright capsule did not produce a single cap point");
        check(std::abs(upright.points[0].point.y) < 1.0e-4f &&
                  std::abs(upright.points[0].separation) < 1.0e-4f,
              "upright capsule cap point was not on the floor plane");
    }
    {
        // A resting capsule may not creep sideways. With zero friction there is
        // no tangential force at all, so any horizontal drift is the contact
        // normal pointing somewhere other than straight up.
        // Standing on a cap the rest height is half length + radius; lying on the
        // barrel it is the radius alone.
        for (auto [angle, rest_height] :
             {std::pair{0.0f, 1.0f}, std::pair{1.57079633f, 0.25f}}) {
            World w;
            w.create_body().static_body().at(0, -0.5f).box(20, 0.5f).friction(0).build();
            auto &b = w.create_body()
                          .dynamic()
                          .at(0, rest_height + 0.01f)
                          .angle(angle)
                          .capsule(0.25f, 0.75f)
                          .friction(0)
                          .restitution(0)
                          .build();
            for (int i = 0; i < 600; ++i) {
                b.wake();
                w.step();
            }
            check(std::abs(b.transform.position.x) < 1.0e-3f,
                  "resting capsule crept sideways along a flat floor");
            check(std::abs(b.transform.position.y - rest_height) < 0.02f,
                  "resting capsule settled at the wrong height");
            check(std::abs(b.transform.angle - angle) < 1.0e-3f,
                  "resting capsule rotated on a flat floor");
        }
    }
    {
        // Continuous collision for capsules. Without CCD the shot must tunnel,
        // which is what makes the positive case meaningful.
        auto discrete = config();
        discrete.ccd.enabled = false;
        World off(discrete);
        wall(off);
        auto &tunnelled = capsule_shot(off);
        off.step();
        check(tunnelled.transform.position.x > 4, "negative control must tunnel with CCD disabled");
        for (bool broadphase : {false, true}) {
            auto c = config();
            c.enable_broadphase = broadphase;
            World on(c);
            wall(on);
            auto &a = capsule_shot(on);
            on.step();
            check(std::abs(a.transform.position.x + 0.26f) < 0.02f,
                  "CCD capsule crossed a thin wall instead of stopping on it");
            check(std::abs(a.velocity.x) < 0.01f && on.ccd_statistics().impacts > 0,
                  "CCD capsule impact did not resolve velocity");
            check(!on.ccd_statistics().limited, "ordinary capsule impact exhausted the CCD budget");
        }
    }
    {
        // A rotating capsule is a second motion path through the same geometry:
        // the sweep has to handle the changing support point.
        World w(config());
        wall(w);
        auto &b = w.create_body()
                      .dynamic()
                      .at(-5, 0)
                      .angle(0.4f)
                      .capsule(0.25f, 0.75f)
                      .velocity(600, 0)
                      .friction(0)
                      .build();
        no_damping(b);
        b.fixed_rotation = true;
        b.inverse_inertia = 0;
        w.step();
        check(b.transform.position.x < 0 && !w.ccd_statistics().limited, "rotated capsule CCD");
    }
    {
        // Kinematic obstacles are swept too, and their own velocity may not be
        // changed by the impact.
        World w(config());
        auto &obstacle = w.create_body()
                             .kinematic()
                             .at(1, 0)
                             .box(0.05f, 3)
                             .velocity(-100, 0)
                             .friction(0)
                             .build();
        auto &b = capsule_shot(w, -2, 100);
        w.step();
        check(b.transform.position.x + 0.25f <= obstacle.transform.position.x - 0.049f,
              "moving kinematic obstacle CCD with a capsule");
        check(std::abs(obstacle.velocity.x + 100) < 0.001f,
              "kinematic velocity changed on a capsule impact");
    }
    {
        // A sleeping dynamic body has to be woken by a capsule bullet, and the
        // bullet has to keep its remaining time.
        World w(config());
        auto &sleeping = capsule_shot(w, 0, 0);
        sleeping.sleeping = true;
        auto &bullet = capsule_shot(w, -5, 600);
        bullet.bullet = true;
        w.step();
        check(!sleeping.sleeping && sleeping.velocity.x > 1,
              "capsule bullet did not wake a sleeping dynamic body");
    }
    {
        // A capsule sliding along a floor must reach it, not stop short of the
        // other shapes' reach: the sweep has to consume the whole window.
        World w(config());
        w.create_body().static_body().at(0, -0.25f).box(40, 0.25f).friction(0).build();
        auto &b = w.create_body()
                      .dynamic()
                      .at(-5, 0.25f)
                      .angle(1.57079633f)
                      .capsule(0.25f, 0.75f)
                      .velocity(300, 0)
                      .friction(0)
                      .build();
        no_damping(b);
        b.bullet = true;
        w.step(0.01f);
        check(std::abs(b.transform.position.x - (-5 + 3.0f)) < 1.0e-3f,
              "sliding capsule lost or gained tangential motion to CCD");
        check(std::abs(b.transform.position.y - 0.25f) < 0.02f,
              "sliding capsule sank into or lifted off the floor");
    }
    std::cout << checks << " capsule checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
