#include <butter/butter.h>

#include <iostream>

using namespace butter;
using namespace butter::math;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};

    world.create_body()
        .static_body()
        .at(0, 0, 0)
        .box(10, 0.5f, 10)
        .build();

    auto& ball = world.create_body()
        .dynamic()
        .at(0, 3, 0)
        .sphere(0.5f)
        .bounciness(0.5f)
        .build();

    auto hit = world.raycast()
        .from(0, 10, 0)
        .toward(0, -1, 0)
        .max_distance(20)
        .first();

    if (hit) {
        std::cout << "射线命中: " << hit->point.y << "\n";
    }

    for (int i = 0; i < 180; ++i) {
        world.step(1.0f / 60.0f);
    }

    std::cout << "球最终高度: " << ball.position.y << "\n";
    return 0;
}
