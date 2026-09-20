#include <butter/physics2d/butter2d.h>
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
static Body &shot(World &w, float x = -5, float velocity = 600, float restitution = 0) {
    auto &b = w.create_body()
                  .dynamic()
                  .at(x, 0)
                  .circle(0.05f)
                  .velocity(velocity, 0)
                  .friction(0)
                  .restitution(restitution)
                  .build();
    no_damping(b);
    return b;
}
int main() try {
    {
        auto discrete = config();
        discrete.ccd.enabled = false;
        World off(discrete);
        wall(off);
        auto &b = shot(off);
        off.step();
        check(b.transform.position.x > 4, "negative control must tunnel with CCD disabled");
        for (bool broadphase : {false, true}) {
            auto c = config();
            c.enable_broadphase = broadphase;
            World on(c);
            wall(on);
            auto &a = shot(on);
            on.step();
            check(a.transform.position.x <= -0.059f && a.transform.position.x > -0.07f,
                  "CCD circle crossed thin wall");
            check(std::abs(a.velocity.x) < 0.01f && on.ccd_statistics().impacts > 0,
                  "CCD must resolve velocity");
            check(!on.ccd_statistics().limited, "ordinary impact exhausted CCD budget");
        }
    }
    {
        World w(config());
        wall(w);
        auto &b = w.create_body()
                      .dynamic()
                      .at(-5, 2)
                      .box(0.1f, 0.1f)
                      .angle(0.4f)
                      .velocity(600, 0)
                      .friction(0)
                      .build();
        no_damping(b);
        b.fixed_rotation = true;
        b.inverse_inertia = 0;
        w.step();
        check(b.transform.position.x < 0 && !w.ccd_statistics().limited, "rotated box CCD");
        World p(config());
        wall(p);
        auto &triangle = p.create_body()
                             .dynamic()
                             .at(-5, 0)
                             .polygon({{-0.1f, -0.1f}, {0.1f, 0}, {-0.1f, 0.1f}})
                             .velocity(600, 0)
                             .build();
        no_damping(triangle);
        triangle.fixed_rotation = true;
        triangle.inverse_inertia = 0;
        p.step();
        check(triangle.transform.position.x < 0, "convex polygon CCD");
    }
    {
        World w(config());
        auto &obstacle = w.create_body()
                             .kinematic()
                             .at(1, 0)
                             .box(0.05f, 3)
                             .velocity(-100, 0)
                             .friction(0)
                             .build();
        auto &b = shot(w, -2, 100);
        w.step();
        check(b.transform.position.x + 0.05f <= obstacle.transform.position.x - 0.049f,
              "moving kinematic obstacle CCD");
        check(std::abs(obstacle.velocity.x + 100) < 0.001f, "kinematic velocity changed on impact");
    }
    {
        World w(config());
        auto &obstacle = w.create_body()
                             .kinematic()
                             .at(-5, 0)
                             .box(0.01f, 1)
                             .velocity(600, 0)
                             .friction(0)
                             .build();
        auto &body = shot(w, 0, 0);
        body.sleeping = true;
        w.step();
        check(!body.sleeping && body.transform.position.x >= obstacle.transform.position.x + 0.059f,
              "default CCD failed to wake a dynamic body hit by a moving kinematic "
              "wall");
    }
    for (bool bullet : {false, true}) {
        World w(config());
        auto &a = shot(w, -5, 600, 1);
        auto &b = shot(w, 5, -600, 1);
        a.bullet = bullet;
        w.step();
        if (bullet)
            check(a.transform.position.x < b.transform.position.x && a.velocity.x < 0 &&
                      b.velocity.x > 0,
                  "dynamic bullet pair did not bounce");
        else
            check(a.transform.position.x > b.transform.position.x,
                  "ordinary dynamic pairs unexpectedly use CCD");
    }
    {
        World w(config());
        wall(w, 0, 1);
        auto &b = shot(w, -1, 120, 1);
        w.step(1.0f / 30);
        check(b.transform.position.x < -2.9f && b.velocity.x < -119,
              "remaining time after bounce was discarded");
        World ricochet(config());
        wall(ricochet, -1, 1);
        wall(ricochet, 1, 1);
        auto &ball = shot(ricochet, 0, 600, 1);
        ricochet.step(0.01f);
        check(ricochet.ccd_statistics().impacts >= 3 && std::abs(ball.transform.position.x) < 0.95f,
              "multiple impacts in one step");
        check(!ricochet.ccd_statistics().limited, "ricochet unexpectedly exhausted budget");
        auto limited = config();
        limited.ccd.max_impacts = 1;
        World bounded(limited);
        wall(bounded, -1, 1);
        wall(bounded, 1, 1);
        auto &fast = shot(bounded, 0, 600, 1);
        bounded.step(0.01f);
        check(bounded.ccd_statistics().limited && bounded.ccd_statistics().remaining_time > 0 &&
                  std::abs(fast.transform.position.x) < 0.95f,
              "budget fallback must not advance unswept time");
    }
    {
        World w(config());
        wall(w).trigger = true;
        auto &b = shot(w);
        w.step();
        check(b.transform.position.x > 4 && w.ccd_statistics().impacts == 0,
              "sensor blocked a sweep");
        World filtered(config());
        wall(filtered).collision_mask = 0;
        auto &f = shot(filtered);
        filtered.step();
        check(f.transform.position.x > 4, "mask filtered obstacle blocked CCD");
        World custom(config());
        wall(custom);
        auto &c = shot(custom);
        custom.contact_filter = [](const Fixture &, const Fixture &) { return false; };
        custom.step();
        check(c.transform.position.x > 4, "custom filter ignored by CCD");
    }
    {
        const Shape box = Box{{2, 0.02f}}, circle = Circle{0.1f};
        ShapeSweep spin{{{0, 0}, 0}, {{0, 0}, 3.14159265f}, {}};
        ShapeSweep stationary{{{0, 1.5f}, 0}, {{0, 1.5f}, 0}, {}};
        auto hit = sweep_shapes(box, spin, circle, stationary);
        check(hit && hit->fraction > 0 && hit->fraction < 1 && hit->converged,
              "pure rotation swept past circle");
        ShapeSweep offset = spin;
        offset.local.position = {2, 0};
        stationary.start.position = {0, 2};
        stationary.end = stationary.start;
        hit = sweep_shapes(circle, offset, circle, stationary);
        check(hit && hit->fraction < 0.5f, "rotating local collider offset was not swept");
        World rotating(config());
        auto &bar = rotating.create_body()
                        .dynamic()
                        .box(2, 0.02f)
                        .angular_velocity(3.14159265f * 60)
                        .build();
        no_damping(bar);
        rotating.create_body().static_body().at(0, 1.5f).circle(0.1f).build();
        rotating.step();
        check(rotating.ccd_statistics().impacts > 0, "world rotation CCD not invoked");
    }
    {
        World w(config());
        auto &staticBody = w.create_empty_body();
        staticBody.type = BodyType::Static;
        staticBody.inverse_mass = 0;
        staticBody.inverse_inertia = 0;
        Fixture obstacle;
        obstacle.shape = Box{{0.01f, 5}};
        w.add_fixture(staticBody, obstacle);
        auto &body = w.create_empty_body();
        body.transform.position = {-5, 0};
        body.velocity = {600, 0};
        no_damping(body);
        Fixture collider;
        collider.shape = Circle{0.05f};
        collider.local.position = {2, 0};
        collider.material.friction = 0;
        w.add_fixture(body, collider);
        int enters = 0, exits = 0;
        w.on_contact = [&](Fixture &, Fixture &, bool enter) { enter ? ++enters : ++exits; };
        w.step();
        check(body.transform.position.x < -2.05f && enters == 1 && exits == 0,
              "compound offset CCD and contact entry");
        w.step();
        check(enters == 1 && exits == 0, "resting CCD contact chatter");
        body.transform.position = {-5, 0};
        w.step();
        check(exits == 1, "CCD contact exit after teleport");
        body.transform.position = {-5, 0};
        body.velocity = {600, 0};
        body.fixtures.front()->material.restitution = 1;
        w.step();
        check(enters == 2 && exits == 2, "transient CCD bounce must emit enter and exit");
    }
    {
        World w(config());
        auto &sleeping = shot(w, 0, 0);
        sleeping.sleeping = true;
        auto &bullet = shot(w, -5, 600);
        bullet.bullet = true;
        w.step();
        check(!sleeping.sleeping && sleeping.velocity.x > 1,
              "bullet did not wake sleeping dynamic body");
        auto gravity = config();
        gravity.gravity = {0, -10};
        World g(gravity);
        wall(g, 0, 1);
        auto &b = shot(g, -1, 120, 1);
        g.step(1.0f / 30);
        check(std::abs(b.velocity.y + 10.0f / 30) < 1e-5f,
              "gravity applied more than once during CCD");
    }
    {
        World w(config());
        wall(w, 2);
        wall(w, 0);
        auto &ball = shot(w);
        ball.velocity.y = 120;
        w.step();
        check(ball.transform.position.x < 0 && std::abs(ball.transform.position.y - 2) < 1e-4f,
              "earliest obstacle or tangential remaining motion was lost");
        check(std::abs(ball.angular_velocity) < 1e-4f,
              "circle impact used an incorrect off-center witness");
        CcdSettings budget;
        budget.max_iterations = 1;
        ShapeSweep spin{{{0, 0}, 0}, {{0, 0}, 3.14159265f}, {}},
            target{{{0, 1.5f}, 0}, {{0, 1.5f}, 0}, {}};
        auto limited =
            sweep_shapes(Shape{Box{{2, 0.02f}}}, spin, Shape{Circle{0.1f}}, target, budget);
        check(limited && !limited->converged,
              "iteration budget silently missed a rotational sweep");
    }
    {
        // Deterministic dense temporal oracle catches missed
        // rotational/translational contacts.
        std::mt19937 rng(9127);
        std::uniform_real_distribution<float> position(-3, 3), angle(-3.14f, 3.14f);
        for (int i = 0; i < 150; ++i) {
            Shape a = i % 2 ? Shape{Circle{0.17f}} : Shape{Box{{0.3f, 0.06f}}},
                  b = Box{{0.05f, 0.8f}};
            ShapeSweep sa{{{position(rng), position(rng)}, angle(rng)},
                          {{position(rng), position(rng)}, angle(rng)},
                          {}},
                sb{};
            if (ccd_detail::separation(a, sa.at(0), b, sb.at(0)).distance <= 0.0001f)
                continue;
            auto hit = sweep_shapes(a, sa, b, sb);
            for (int sample = 1; sample <= 1200; ++sample) {
                float t = float(sample) / 1200;
                Contact contact;
                if (test(a, sa.at(t), b, sb.at(t), contact)) {
                    check(hit && hit->fraction <= t + 0.0001f,
                          "dense oracle found missed CCD contact");
                    break;
                }
            }
        }
        ShapeSweep miss{{{-2, 1}, 0}, {{2, 1}, 0}, {}}, fixed{};
        check(!sweep_shapes(Shape{Circle{0.1f}}, miss, Shape{Box{{0.1f, 0.1f}}}, fixed),
              "clear miss reported a collision");
    }
    {
        World w(config());
        for (int i = 0; i < 1000; ++i)
            wall(w, 100 + float(i));
        wall(w);
        auto &body = shot(w);
        w.step();
        check(body.transform.position.x < 0 && w.ccd_statistics().sweeps < 5,
              "distant static obstacles were not rejected before narrow-phase CCD");
    }
    {
        // A failing local ricochet must not discard an unrelated body's time.
        auto c = config();
        c.ccd.max_impacts = 1;
        World w(c);
        wall(w, -1, 1);
        wall(w, 1, 1);
        auto &ball = shot(w, 0, 600, 1);
        auto &free = w.create_body().at(0, 20).circle(.1f).velocity(3, 0).build();
        no_damping(free);
        w.step(.01f);
        check(std::abs(free.transform.position.x - .03f) < 1e-6f,
              "local budget failure stopped unrelated motion");
        check(std::abs(ball.transform.position.x) < .95f, "local clamp tunneled");
        auto &stats = w.ccd_statistics();
        check(stats.budget_exhaustions == 1 && stats.non_convergences == 0 &&
                  !stats.diagnostics.empty(),
              "missing budget diagnostic");
        auto d = stats.diagnostics.front();
        check(d.reason == CcdFailure::ImpactBudget && d.remaining_time > 0 &&
                  std::abs(d.advanced_time + d.remaining_time - .01f) < 1e-6f,
              "incorrect diagnostic times");
        check(std::abs(stats.advanced_time - .01f) < 1e-6f,
              "world time was discarded on local failure");
    }
    {
        auto c = config();
        c.ccd.max_iterations = 1;
        World w(c);
        auto &bar = w.create_body().box(2, .02f).angular_velocity(3.14159265f * 60).build();
        no_damping(bar);
        w.create_body().static_body().at(0, 1.5f).circle(.1f).build();
        auto &free = w.create_body().at(20, 20).circle(.1f).velocity(3, 0).build();
        no_damping(free);
        w.step();
        check(w.ccd_statistics().non_convergences > 0, "missing non-convergence diagnostic");
        check(std::abs(free.transform.position.x - 20.05f) < 1e-5f,
              "non-convergence blocked unrelated body");
    }
    {
        // Joint projection is a second motion path, and needs its own sweep.
        auto c = config();
        World w(c);
        w.create_body().static_body().at(0, -.005f).box(2, .005f).build();
        auto &small = w.create_body().at(0, .08f).box(.02f, .02f).build();
        auto &anchor = w.create_body().static_body().at(0, -1).circle(.01f).build();
        w.add_distance_joint(small, anchor, 0);
        w.step();
        check(small.transform.position.y >= .019f, "joint position correction crossed thin floor");
        check(w.step_statistics().position_clamps > 0, "projection CCD guard was not exercised");
    }
    {
        World w(config());
        w.create_body().static_body().at(0, -.005f).box(2, .005f).build();
        auto &small = w.create_body().at(0, .05f).box(.02f, .02f).build();
        auto &upper = w.create_body().at(0, .5f).box(.5f, .5f).mass(1000).build();
        small.fixed_rotation = upper.fixed_rotation = true;
        bool traced = false;
        w.on_contact_diagnostic = [&](const World::ContactDiagnostic &d) {
            if (!d.after_integration && d.a == &small) {
                traced = true;
                check(d.after_a.position.y >= .019f, "contact projection trace crossed floor");
            }
        };
        w.step();
        check(traced && small.transform.position.y >= .019f,
              "contact position correction crossed floor outside original pairs");
        check(small.angular_velocity == 0 && upper.angular_velocity == 0,
              "fixed rotation changed during contact solving");
    }
    {
        // Sleep-cache snapshots must notice direct public geometry edits, even
        // when the new shapes' AABBs still overlap.
        World w;
        auto &ground = w.create_empty_body();
        ground.type = BodyType::Static;
        ground.inverse_mass = ground.inverse_inertia = 0;
        ground.transform.position = {0, -.5f};
        Fixture floor;
        floor.shape = Box{{1, .5f}};
        auto &f = w.add_fixture(ground, floor);
        auto &body = w.create_empty_body();
        body.transform.position = {.9f, .5f};
        Fixture circle;
        circle.shape = Circle{.5f};
        w.add_fixture(body, circle);
        int exits = 0;
        w.on_contact = [&](Fixture &, Fixture &, bool enter) {
            if (!enter)
                ++exits;
        };
        for (int i = 0; i < 120; ++i)
            w.step();
        check(body.sleeping, "sleep-cache test did not settle");
        f.shape = Circle{.6f};
        w.step();
        check(exits == 1, "sleep cache ignored direct shape mutation");
        f.shape = Box{{1, .5f}};
        w.step();
        f.local.position.x = -.6f;
        w.step();
        check(exits == 2, "sleep cache ignored fixture offset mutation");
        f.local.position.x = 0;
        w.step();
        ground.transform.position.x = -.6f;
        w.step();
        check(exits == 3, "sleep cache ignored public body transform mutation");
    }
    std::cout << checks << " CCD checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
