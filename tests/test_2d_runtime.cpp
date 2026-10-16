#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace butter::physics2d;
static void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}

static void queries() {
    World world;
    for (int i = 0; i < 100; ++i)
        world.create_body().static_body().at(5 + float(i) * 3, 0).box(.5f, .5f).build();
    const auto before = world.step_statistics().query_sync_scans;
    world.query_batch([&](const auto &batch) {
        for (int i = 0; i < 100; ++i) {
            auto hit = batch.raycast({0, 0}, {1, 0}, 1000);
            check(hit && hit->body_index == 0 && hit->distance == 4.5f, "batch ray missed nearest body");
            check(batch.query_aabb({{4, -1}, {6, 1}}).size() == 1, "batch AABB disagrees");
        }
        bool locked = false;
        try { world.create_body().build(); } catch (const std::logic_error &) { locked = true; }
        check(locked, "query batch allowed structural mutation");
    });
    check(world.step_statistics().query_sync_scans == before + 1, "batch scanned world per query");
    auto bodies = world.query_aabb({{4, -1}, {6, 1}});
    bodies[0]->transform.position.x = 2;
    check(world.raycast({0, 0}, {1, 0})->distance == 1.5f, "ordinary query lost direct-write compatibility");
    try { world.query_batch([](const auto &) { throw std::runtime_error("test"); }); }
    catch (const std::runtime_error &) {}
    check(!world.locked(), "throwing batch left world locked");

    DynamicTree tree;
    for (unsigned i = 0; i < 100; ++i)
        tree.create_proxy({{5 + float(i) * 3, -1}, {6 + float(i) * 3, 1}}, i, 0);
    int visited = 0;
    tree.raycast_clipped({0, 0}, {1, 0}, 1000, [&](unsigned, float distance) {
        ++visited;
        return distance;
    });
    check(visited == 1, "closest ray hit did not prune distant tree branches");
}

static void filters() {
    World::Config cfg;
    cfg.gravity = {};
    World world(cfg);
    auto &a = world.create_empty_body();
    auto &b = world.create_empty_body();
    b.transform.position.x = .8f;
    Fixture definition;
    definition.shape = Circle{.5f};
    world.add_fixture(a, definition);
    world.add_fixture(b, definition);
    auto &joint = world.add_distance_joint(a, b, .8f);
    joint.collide_connected = false;
    int enters = 0, exits = 0;
    world.on_contact = [&](Fixture &, Fixture &, bool enter) { enter ? ++enters : ++exits; };
    world.step();
    check(enters == 0, "joint exclusion ignored");
    joint.collide_connected = true;
    world.step();
    check(enters == 1, "public joint filter edit not observed");
    joint.collide_connected = false;
    world.step();
    check(exits == 1, "joint filter edit did not remove contact");
}

static void solver_paths() {
    World::Config base, fast;
    base.cache_solver_islands = base.buffered_velocities = false;
    base.ccd.tight_fixture_bounds = false;
    fast.detailed_profiling = true;
    int dispatched = 0;
    fast.parallel_for = [&](std::size_t count, const std::function<void(std::size_t)> &work) {
        ++dispatched;
        std::vector<std::thread> workers;
        for (std::size_t t = 0; t < std::min<std::size_t>(count, 3); ++t)
            workers.emplace_back([&, t] { for (std::size_t i = t; i < count; i += 3) work(i); });
        for (auto &worker : workers) worker.join();
    };
    World reference(base), optimized(fast);
    auto setup = [](World &world) {
        world.create_body().static_body().at(0, -.5f).box(100, .5f).build();
        std::vector<Body *> bodies;
        for (int column = 0; column < 4; ++column)
            for (int row = 0; row < 5; ++row)
                bodies.push_back(&world.create_body().at(float(column) * 5, .5f + row * 1.01f)
                    .box(.5f, .5f).build());
        return bodies;
    };
    const auto a = setup(reference), b = setup(optimized);
    std::size_t reuses = 0;
    for (int frame = 0; frame < 180; ++frame) {
        for (auto *body : a) body->wake();
        for (auto *body : b) body->wake();
        reference.step(); optimized.step();
        reuses += optimized.step_statistics().solver_island_reuses;
        for (std::size_t i = 0; i < a.size(); ++i) {
            check((a[i]->transform.position - b[i]->transform.position).length() < 1e-5f,
                  "optimized solver changed position");
            check((a[i]->velocity - b[i]->velocity).length() < 1e-5f,
                  "optimized solver changed velocity");
            check(std::abs(a[i]->transform.angle - b[i]->transform.angle) < 1e-5f &&
                  std::abs(a[i]->angular_velocity - b[i]->angular_velocity) < 1e-5f,
                  "optimized solver changed angular motion");
        }
    }
    check(dispatched > 0 && reuses > 0, "optimized paths were not exercised");
    check(reference.step_statistics().narrow_ms == 0, "disabled detailed clock still ran");
}

static void swept_envelopes() {
    for (const Shape &shape : {Shape{Box{{2, .1f}}}, Shape{Circle{.5f}}, Shape{Capsule{.2f, 1}}})
        for (float angle : {0.0f, .03f, 1.0f, 7.0f, -7.0f}) {
            ShapeSweep sweep{{{1, 2}, .2f}, {{-2, 3}, .2f + angle}, {{3, -1}, .6f}};
            const auto envelope = ccd_detail::tight_swept_bounds(shape, sweep, .0001f);
            for (int sample = 0; sample <= 100; ++sample)
                check(contains(envelope, compute_aabb(shape, sweep.at(float(sample) / 100))),
                      "CCD envelope failed to contain a rotated offset fixture");
        }
}

static void soft_contacts() {
    float previous = 0;
    for (float mass : {1.0f, 100.0f}) {
        World::Config cfg;
        cfg.contact_hertz = 20;
        cfg.substeps = 4;
        World world(cfg);
        world.create_body().static_body().at(0, -.5f).box(10, .5f).build();
        auto &b = world.create_body().at(0, 2).circle(.5f).mass(mass).build();
        for (int frame = 0; frame < 360; ++frame) {
            b.wake();
            world.step();
            check(std::isfinite(b.transform.position.y) && b.transform.position.y >= .45f,
                  "soft contact failed to support load");
        }
        check(std::abs(b.transform.position.y - .5f) < .02f && std::abs(b.velocity.y) < .05f,
              "soft contact failed to settle");
        if (previous != 0)
            check(std::abs(previous - b.transform.position.y) < .001f,
                  "softness changed with absolute body mass");
        previous = b.transform.position.y;
    }
}

static void simd_contacts() {
    World::Config cfg;
    cfg.simd_contacts = true;
    World scalar, packets(cfg);
    std::vector<Body *> a, b;
    auto populate = [](World &world, std::vector<Body *> &bodies) {
        world.create_body().static_body().at(0, -.5f).box(100, .5f).build();
        for (int i = 0; i < 12; ++i)
            bodies.push_back(&world.create_body().at(float(i) * 3, .5f).circle(.5f)
                .velocity(.1f, 0).friction(.5f).build());
    };
    populate(scalar, a); populate(packets, b);
    std::size_t count = 0;
    for (int frame = 0; frame < 120; ++frame) {
        for (auto *p : a) p->wake();
        for (auto *p : b) p->wake();
        scalar.step(); packets.step();
        count += packets.step_statistics().simd_packets;
        for (std::size_t i=0; i<a.size(); ++i) {
            check((a[i]->transform.position-b[i]->transform.position).length() < 1e-5f,
                  "SIMD packet changed independent body position");
            check((a[i]->velocity-b[i]->velocity).length() < 1e-5f &&
                  std::abs(a[i]->angular_velocity-b[i]->angular_velocity) < 1e-5f,
                  "SIMD packet changed contact response");
        }
    }
    check(count > 0, "SIMD regression never dispatched a packet");
}

static void translated_contacts() {
    World::Config cfg;
    cfg.solver_iterations = 2;
    World reference(cfg);
    cfg.translated_box_contacts = true;
    World optimized(cfg);
    auto setup = [](World &w) {
        w.create_body().static_body().at(0, -.5f).box(20, .5f).build();
        std::vector<Body *> bodies;
        for (int i = 0; i < 8; ++i) {
            auto &b = w.create_body().at(0, .5f + i).box(.5f, .5f).build();
            b.fixed_rotation = true;
            bodies.push_back(&b);
        }
        return bodies;
    };
    const auto a = setup(reference), b = setup(optimized);
    std::size_t updates = 0;
    for (int frame = 0; frame < 180; ++frame) {
        for (auto *p : a) p->wake();
        for (auto *p : b) p->wake();
        if (frame == 100) {
            a.back()->transform.angle = b.back()->transform.angle = .02f;
        }
        reference.step(); optimized.step();
        updates += optimized.step_statistics().translated_manifolds;
        for (std::size_t i = 0; i < a.size(); ++i)
            check((a[i]->transform.position - b[i]->transform.position).length() < 2e-4f,
                  "local-anchor translation changed supported box motion");
    }
    check(updates > 0, "translated-contact regression did not reuse anchors");
}

int main() try {
    queries(); filters(); solver_paths(); swept_envelopes(); soft_contacts(); simd_contacts(); translated_contacts();
    std::cout << "runtime optimization regressions passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
