#include <butter/butter.h>

#include <iostream>

using namespace butter;
using namespace butter::math;

static int failures = 0;

static void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

int main() {
    World::Config config;
    config.gravity = {0, 0, 0};
    config.enable_sleeping = false;
    config.broadphase_cell_size = 1.0f;
    config.broadphase_max_cells_per_body = 1024;

    World broadphase_world(config);
    broadphase_world.create_body().static_body().at(0, 0, 0)
        .box(0.5f, 0.5f, 0.5f).build();
    broadphase_world.create_body().dynamic().at(0.8f, 0, 0)
        .box(0.5f, 0.5f, 0.5f).build();
    broadphase_world.create_body().dynamic().at(100, 100, 100)
        .sphere(1.0f).build();
    broadphase_world.create_body().dynamic().at(-100, -100, -100)
        .sphere(1.0f).build();

    int broadphase_collisions = 0;
    broadphase_world.on_collision = [&](const CollisionEvent&) {
        ++broadphase_collisions;
    };
    broadphase_world.step();
    check(broadphase_collisions > 0, "nearby bodies reach narrow phase");
    check(broadphase_world.broadphase_candidate_count() < 6,
          "spatial grid culls distant body pairs");

    const auto queried = broadphase_world.query_aabb({{-1, -1, -1}, {2, 1, 1}});
    check(queried.size() == 2, "query_aabb returns only nearby bodies");

    // A large proxy (typical ground plane) bypasses cell insertion but still
    // participates in exact broad-phase filtering.
    auto large_config = config;
    large_config.broadphase_max_cells_per_body = 16;
    World large_world(large_config);
    large_world.create_body().static_body().at(0, -1, 0)
        .box(100, 0.5f, 100).build();
    auto& ball = large_world.create_body().dynamic().at(0, 0, 0)
        .sphere(0.5f).build();
    int large_collisions = 0;
    large_world.on_collision = [&](const CollisionEvent&) {
        ++large_collisions;
    };
    large_world.step();
    check(large_collisions > 0, "large proxy still collides");
    check(ball.position.y >= -0.01f, "large proxy contact is solved");

    if (failures == 0) {
        std::cout << "test_broadphase: OK\n";
        return 0;
    }
    std::cerr << "test_broadphase: " << failures << " failure(s)\n";
    return 1;
}
