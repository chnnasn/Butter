#include <butter/butter.h>

#include <cmath>
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
    {
        World world;
        world.gravity = {0, -10, 0};

        auto& body = world.create_body()
            .dynamic()
            .at(0, 0, 0)
            .mass(1.0f)
            .build();

        world.step(1.0f);
        check(std::abs(body.velocity.y + 10.0f) < 0.02f,
              "gravity accelerates body for one second");
    }

    {
        World world;
        world.create_body()
            .static_body()
            .at(0, 0, 0)
            .box(1, 1, 1)
            .build();

        auto hit = world.raycast({0, 5, 0}, {0, -1, 0}, 10.0f);
        check(hit.has_value() && hit->hit, "raycast hits box");
        if (hit) {
            check(std::abs(hit->point.y - 1.0f) < 1.0e-4f, "raycast hit point");
            check(std::abs(hit->normal.y - 1.0f) < 1.0e-4f, "raycast hit normal");
        }
    }

    {
        World world;
        world.create_body()
            .static_body()
            .at(0, 0, 0)
            .box(10, 0.5f, 10)
            .build();
        auto& ball = world.create_body()
            .dynamic()
            .at(0, 5, 0)
            .sphere(0.5f)
            .build();

        for (int i = 0; i < 240; ++i) {
            world.step(1.0f / 60.0f);
        }
        check(ball.position.y > 0.0f && ball.position.y < 2.0f,
              "ball rests near the ground after simulation");
    }

    if (failures == 0) {
        std::cout << "test_world: OK\n";
        return 0;
    }
    std::cerr << "test_world: " << failures << " failure(s)\n";
    return 1;
}
