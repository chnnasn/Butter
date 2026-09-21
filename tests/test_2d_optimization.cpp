#include <butter/physics2d/butter2d.h>
#include <iostream>
#include <stdexcept>
#include <tuple>
using namespace butter::physics2d;
static void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
int main() try {
    // Stack storage and the unbounded fallback must preserve the public vertex
    // conversion, including empty/degenerate polygons and capsule fallback.
    {
        std::vector<Shape> shapes{Box{{.7f, .3f}}, Capsule{.2f, .6f}, Circle{.3f}};
        for (int n : {0, 1, 3, 16, 17, 64}) {
            Polygon p;
            for (int i = 0; i < n; ++i) {
                float angle = 6.2831853f * i / n;
                p.vertices.push_back({std::cos(angle), std::sin(angle)});
            }
            shapes.push_back(p);
        }
        for (const auto &shape : shapes) {
            Transform transform{{2.3f, -4.1f}, .731f};
            const auto expected = world_vertices(shape, transform);
            const shape_detail::WorldVertices storage(shape, transform);
            const auto actual = storage.view();
            require(actual.size() == expected.size() &&
                        std::equal(actual.begin(), actual.end(), expected.begin()),
                    "temporary vertex storage changed geometry");
        }
    }
    {
        World::Config cfg;
        cfg.gravity = {};
        cfg.ccd.enabled = false;
        World w(cfg);
        auto &a = w.create_empty_body();
        auto &b = w.create_empty_body();
        b.transform.position = {2.001f, 0};
        Fixture f;
        f.shape = Circle{1};
        w.add_fixture(a, f);
        w.add_fixture(b, f);
        int enters = 0;
        w.on_contact = [&](Fixture &, Fixture &, bool enter) { enters += enter; };
        w.step();
        w.step();
        require(enters == 0, "proximity constraint became a discrete contact event");
    }
    // Reused broadphase buffers must match brute force after public edits,
    // changes to proxy size/type, body insertion and body removal.
    World::Config config;
    config.gravity = {};
    World indexed(config);
    config.enable_broadphase = false;
    World oracle(config);
    std::vector<Body *> a, b;
    for (int i = 0; i < 100; ++i) {
        for (auto pair : {std::make_pair(&indexed, &a), std::make_pair(&oracle, &b)}) {
            auto &body = pair.first->create_body()
                             .at(float(i % 10), float(i / 10))
                             .box(.6f, .6f)
                             .trigger()
                             .build();
            body.user_data = i;
            pair.second->push_back(&body);
        }
    }
    using Event = std::tuple<std::uint64_t, std::uint64_t, bool>;
    std::vector<Event> x, y;
    indexed.on_trigger = [&](Body &p, Body &q, bool enter) {
        x.emplace_back(p.user_data, q.user_data, enter);
    };
    oracle.on_trigger = [&](Body &p, Body &q, bool enter) {
        y.emplace_back(p.user_data, q.user_data, enter);
    };
    for (int frame = 0; frame < 12; ++frame) {
        if (frame == 2 || frame == 6)
            for (int i = 0; i < 100; ++i) {
                a[i]->transform.position.x += float(i % 3) * .2f;
                b[i]->transform = a[i]->transform;
                a[i]->transform.angle = b[i]->transform.angle = float(i % 4) * .2f;
            }
        if (frame == 4) {
            a[0]->shape = b[0]->shape = Box{{100, .1f}};
        }
        if (frame == 5) {
            a[1]->type = b[1]->type = BodyType::Static;
        }
        if (frame == 7) {
            a[2]->collision_mask = b[2]->collision_mask = 0;
        }
        x.clear();
        y.clear();
        indexed.step();
        oracle.step();
        require(indexed.broadphase_candidate_count() == oracle.broadphase_candidate_count(),
                "broadphase missed a pair after mutation");
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        require(x == y, "cached event transitions differ from brute force");
    }
    indexed.destroy_body(*a.back());
    oracle.destroy_body(*b.back());
    indexed.create_body().at(4, 4).circle(.4f).trigger().build();
    oracle.create_body().at(4, 4).circle(.4f).trigger().build();
    indexed.step();
    oracle.step();
    require(indexed.broadphase_candidate_count() == oracle.broadphase_candidate_count(),
            "replaced body reused stale broadphase indices");

    // Resting groups retain contact state, but do no iterative solver work.
    // Force, teleport, geometry, local offset and filtering changes must wake
    // supported bodies and invalidate geometry/contact reuse.
    for (int mutation = 0; mutation < 8; ++mutation) {
        World w;
        auto &floor = w.create_empty_body();
        floor.type = BodyType::Static;
        floor.inverse_mass = floor.inverse_inertia = 0;
        Fixture f;
        f.shape = Box{{10, .5f}};
        floor.transform.position = {0, -.5f};
        auto &fixture = w.add_fixture(floor, f);
        auto &box = w.create_body().at(0, .5f).box(.5f, .5f).restitution(0).build();
        std::size_t reused = 0;
        for (int i = 0; i < 180; ++i) {
            w.step();
            reused += w.step_statistics().cached_manifolds;
        }
        require(box.sleeping && reused > 0, "resting geometry did not reuse its manifold");
        w.step();
        require(w.step_statistics().active_constraints == 0 &&
                    w.step_statistics().position_corrections == 0 &&
                    w.step_statistics().stationary_steps == 1,
                "sleeping group still solved constraints");
        switch (mutation) {
        case 0:
            box.force = {20, 0};
            break;
        case 1:
            floor.transform.position.y = -10;
            break;
        case 2:
            fixture.shape = Circle{.01f};
            break;
        case 3:
            fixture.local.position.y = -10;
            break;
        case 4:
            fixture.collision_mask = 0;
            break;
        case 5:
            w.contact_filter = [](const Fixture &, const Fixture &) { return false; };
            break;
        case 6:
            w.destroy_fixture(fixture);
            break;
        case 7:
            w.destroy_body(floor);
            break;
        }
        w.step();
        require(!box.sleeping, "public edit or force failed to wake supported body");
        require(w.step_statistics().stationary_steps == 0,
                "woken body incorrectly used stationary world path");
        if (mutation == 0)
            require(box.velocity.x > 0, "force applied to sleeping body was lost");
        else
            require(box.velocity.y < 0, "invalidated support prevented gravity");
    }

    {
        World w;
        auto &body = w.create_empty_body();
        body.type = BodyType::Kinematic;
        w.step();
        require(w.step_statistics().stationary_steps == 1,
                "stationary kinematic prevented resting path");
        body.velocity = {60, 0};
        w.step();
        require(std::abs(body.transform.position.x - 1) < 1e-5f &&
                    w.step_statistics().stationary_steps == 0,
                "moving kinematic skipped integration");
        body.velocity = {};
        body.angular_velocity = 60;
        w.step();
        require(std::abs(body.transform.angle - 1) < 1e-5f,
                "rotating kinematic skipped integration");
        require(std::abs(w.simulation_time() - 3.0 / 60) < 1e-7,
                "resting path discarded simulation time");
    }

    // A tall sparse bullet layout exercises the coupled CCD path. Its sweep
    // index should query Y intervals rather than all overlapping X intervals.
    Shape shape = Circle{.1f};
    std::vector<Transform> transforms(200);
    std::vector<Vec2> velocities(200, {1, 0});
    std::vector<float> angular(200);
    std::vector<CcdMotion> motions(200);
    for (int i = 0; i < 200; ++i) {
        transforms[i].position = {0, float(i) * 3};
        auto &m = motions[i];
        m.transform = &transforms[i];
        m.velocity = &velocities[i];
        m.angular_velocity = &angular[i];
        m.dynamic = m.bullet = true;
        m.inverse_mass = 1;
        m.colliders.push_back({&shape});
    }
    CcdWorkspace workspace;
    for (int step = 0; step < 3; ++step) {
        auto stats = advance_continuous(motions, 1.f / 60, {}, {}, {}, &workspace);
        require(stats.bounds_tests < 200 && stats.sweeps == 0,
                "CCD index enumerated unrelated vertical pairs");
        require(!stats.limited, "reused CCD workspace retained a limit");
    }
    std::cout << "2D optimization regressions passed\n";
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
