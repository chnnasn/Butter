#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <cstdint>
using namespace butter::physics2d;

static int checks = 0;
static void check(bool ok, const char *message) {
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}

// A tree must answer exactly what a brute force scan answers. The world keeps
// its own tree, so this exercises it through the public API.
static void parity_against_bruteforce() {
    std::mt19937 rng(4242);
    std::uniform_real_distribution<float> position(-14, 14), size(0.2f, 1.4f);
    struct Spec {
        Vec2 position;
        float angle;
        Vec2 half_extents;
    };
    // Build the geometry once and replay it into both worlds: relying on the
    // evaluation order of repeated rng draws would make the two worlds differ.
    std::vector<Spec> spec;
    for (int i = 0; i < 240; ++i) {
        const float x = position(rng), y = position(rng), angle = position(rng) * 0.2f;
        const float hx = size(rng), hy = size(rng) * 0.7f;
        spec.push_back({{x, y}, angle, {hx, hy}});
    }
    auto build = [&](World &world, std::vector<Body *> &bodies) {
        for (std::size_t i = 0; i < spec.size(); ++i) {
            auto &body = world.create_body()
                             .at(spec[i].position.x, spec[i].position.y)
                             .angle(spec[i].angle)
                             .box(spec[i].half_extents.x, spec[i].half_extents.y)
                             .trigger()
                             .build();
            body.user_data = std::uint64_t(i);
            bodies.push_back(&body);
        }
    };
    World::Config indexed_config;
    indexed_config.gravity = {};
    indexed_config.ccd.enabled = false;
    World::Config oracle_config = indexed_config;
    oracle_config.enable_broadphase = false;
    World indexed(indexed_config), oracle(oracle_config);
    std::vector<Body *> a, b;
    build(indexed, a);
    build(oracle, b);
    for (std::size_t i = 0; i < a.size(); ++i)
        check(a[i]->shape.index() == b[i]->shape.index() &&
                  std::get<Box>(a[i]->shape).half_extents ==
                      std::get<Box>(b[i]->shape).half_extents &&
                  a[i]->transform.position == b[i]->transform.position,
              "parity worlds must start identical");
    for (int frame = 0; frame < 30; ++frame) {
        if (frame % 3 == 0) {
            std::uniform_real_distribution<float> step(-1.5f, 1.5f);
            for (std::size_t i = 0; i < a.size(); ++i) {
                Vec2 delta{step(rng), step(rng)};
                float turn = step(rng) * 0.1f;
                // Tiny sub-margin motion must not hide a pair change either.
                Vec2 nudge{};
                if (i % 7 == 0)
                    nudge = {step(rng) * 0.0005f, step(rng) * 0.0005f};
                for (Body *body : {a[i], b[i]}) {
                    body->transform.position += delta + nudge;
                    body->transform.angle += turn;
                }
            }
        }
        if (frame == 5) {
            a[3]->shape = Box{{6, 0.05f}};
            b[3]->shape = Box{{6, 0.05f}};
            a[9]->type = BodyType::Static;
            b[9]->type = BodyType::Static;
        }
        indexed.step();
        oracle.step();
        check(indexed.broadphase_candidate_count() == oracle.broadphase_candidate_count(),
              "dynamic tree candidate count differs from brute force");
        // Ray and AABB queries must agree with a linear scan as well.
        const AABB area{{-3.0f, -2.0f}, {4.0f, 5.0f}};
        const auto fast = indexed.query_aabb(area);
        const auto slow = oracle.query_aabb(area);
        // The two worlds own different Body objects, so compare identity, not
        // pointers: user_data carries the spec index in both worlds.
        std::vector<std::int64_t> fast_ids, slow_ids;
        for (auto *body : fast)
            fast_ids.push_back(std::int64_t(body->user_data));
        for (auto *body : slow)
            slow_ids.push_back(std::int64_t(body->user_data));
        std::sort(fast_ids.begin(), fast_ids.end());
        std::sort(slow_ids.begin(), slow_ids.end());
        if (fast_ids != slow_ids) {
            std::cerr << "frame " << frame << " query fast=" << fast_ids.size()
                      << " slow=" << slow_ids.size() << '\n';
            std::vector<std::int64_t> missing, extra;
            std::set_difference(slow_ids.begin(), slow_ids.end(), fast_ids.begin(), fast_ids.end(),
                                std::back_inserter(missing));
            std::set_difference(fast_ids.begin(), fast_ids.end(), slow_ids.begin(), slow_ids.end(),
                                std::back_inserter(extra));
            for (auto id : missing) {
                const Body *body = a[std::size_t(id)];
                const AABB box = compute_aabb(body->shape, body->transform);
                std::cerr << "  missing id=" << id << " box=" << box.min.x << ',' << box.min.y
                          << " .. " << box.max.x << ',' << box.max.y << '\n';
            }
            for (auto id : extra)
                std::cerr << "  extra id=" << id << '\n';
        }
        check(fast_ids == slow_ids, "tree AABB query differs from brute force");
        const auto fast_hit = indexed.raycast({-40, 0.5f}, {1, 0}, 100);
        const auto slow_hit = oracle.raycast({-40, 0.5f}, {1, 0}, 100);
        check(fast_hit.has_value() == slow_hit.has_value(),
              "tree raycast disagreed about a hit");
        if (fast_hit && slow_hit) {
            check(std::abs(fast_hit->distance - slow_hit->distance) < 1.0e-4f,
                  "tree raycast found a different distance");
            check(fast_hit->body_index == slow_hit->body_index,
                  "tree raycast found a different body");
        }
        const auto up_hit = indexed.raycast({1.25f, -50}, {0, 1}, 200);
        const auto up_slow = oracle.raycast({1.25f, -50}, {0, 1}, 200);
        check(up_hit.has_value() == up_slow.has_value(), "vertical raycast parity");
        if (up_hit && up_slow)
            check(std::abs(up_hit->distance - up_slow->distance) < 1.0e-4f,
                  "vertical raycast distance parity");
    }
}

// The tree must survive a stream of inserts, removals and teleports.
static void structural_churn() {
    World::Config config;
    config.gravity = {};
    World world(config);
    std::vector<Body *> bodies;
    for (int i = 0; i < 120; ++i)
        bodies.push_back(&world.create_body().at(float(i % 12) * 3, float(i / 12) * 3).box(1, 1).trigger().build());
    for (int frame = 0; frame < 40; ++frame) {
        if (frame % 4 == 0 && bodies.size() > 20) {
            world.destroy_body(*bodies.back());
            bodies.pop_back();
        }
        if (frame % 5 == 0)
            bodies.push_back(&world.create_body()
                                  .at(float(frame) * 1.5f - 30, 3)
                                  .circle(0.8f)
                                  .trigger()
                                  .build());
        if (frame % 6 == 0 && !bodies.empty()) {
            bodies.front()->transform.position = {float(frame) - 20, -6};
            bodies.front()->transform.angle = float(frame) * 0.3f;
        }
        world.step();
    }
    // Everything the tree reports must be findable by a brute force query.
    const AABB area{{-40, -40}, {40, 40}};
    check(world.query_aabb(area).size() == world.body_count(),
          "tree lost a body's proxy during churn");
    check(world.raycast({0, -1000}, {0, 1}, 2000).has_value(),
          "raycast lost every body during churn");
}

// Directly exercise the tree: structure, fat boxes and pruning behaviour.
static void tree_behaviour() {
    DynamicTree tree;
    std::vector<std::uint32_t> ids;
    for (int i = 0; i < 200; ++i) {
        AABB box{{float(i) * 2.0f, 0}, {float(i) * 2.0f + 1.0f, 1}};
        ids.push_back(tree.create_proxy(box, std::uint32_t(i), 0.05f));
    }
    check(tree.validate(), "tree structure invalid after inserts");
    check(tree.size() == 200, "tree lost proxies");
    int hits = 0;
    tree.query(AABB{{100, -1}, {100.5f, 2}}, 0.05f, [&](std::uint32_t) { ++hits; });
    check(hits == 1, "query must find exactly the overlapping proxy");
    // A move that stays inside the fat box must not refile the proxy.
    check(!tree.move_proxy(ids[0], AABB{{0, 0}, {1, 1}}, {0, 0}, 0.05f),
          "sub-margin move must not refile the proxy");
    // A large move must refile it and be visible at the destination only.
    const AABB destination{{1000, 1000}, {1001, 1001}};
    check(tree.move_proxy(ids[0], destination, {1000, 1000}, 0.05f), "large move refiles");
    check(tree.validate(), "tree structure invalid after move");
    hits = 0;
    tree.query(AABB{{-1, -1}, {1.5f, 2}}, 0.05f, [&](std::uint32_t) { ++hits; });
    check(hits == 0, "moved proxy still found at its old location");
    hits = 0;
    tree.query(AABB{{999, 999}, {1002, 1002}}, 0.05f, [&](std::uint32_t payload) {
        ++hits;
        check(payload == 0, "moved proxy kept its payload");
    });
    check(hits == 1, "moved proxy not found at its new location");
    for (std::size_t i = 0; i < ids.size(); i += 2)
        tree.destroy_proxy(ids[i]);
    check(tree.validate(), "tree structure invalid after removals");
    check(tree.size() == 100, "tree size wrong after removals");
    int remaining = 0;
    tree.query(AABB{{-1000, -1000}, {1000, 1000}}, 0.05f, [&](std::uint32_t) { ++remaining; });
    check(remaining == 100, "removed proxies still found");
    // Ray traversal must find every crossed proxy and visit the nearest first.
    DynamicTree rays;
    constexpr int ray_count = 60;
    for (int i = 0; i < ray_count; ++i)
        rays.create_proxy(AABB{{float(i) * 3 - 30, -1}, {float(i) * 3 - 29, 1}},
                          std::uint32_t(i), 0.01f);
    float first_entry = -1, previous_entry = -1;
    int visited = 0;
    // Range must cover the far end of the line of boxes (last entry ~247).
    rays.raycast({-100, 0}, {1, 0}, 1000, [&](std::uint32_t payload, float entry) {
        if (first_entry < 0) {
            first_entry = entry;
            check(payload == 0, "ray did not visit the nearest proxy first");
        }
        // Boxes are disjoint and sorted along x, so entries must not go backwards.
        check(entry >= previous_entry - 1.0e-6f, "ray visited proxies out of order");
        previous_entry = entry;
        ++visited;
        return true;
    });
    check(visited == ray_count, "ray missed a proxy");
    check(std::abs(first_entry - 70.0f) < 0.05f, "ray entry distance wrong");
}

int main() try {
    parity_against_bruteforce();
    structural_churn();
    tree_behaviour();
    std::cout << checks << " broadphase checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
