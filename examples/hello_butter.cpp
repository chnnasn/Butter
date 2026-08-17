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
        .box(50, 0.5f, 50)
        .friction(0.8f)
        .build();

    auto& ball = world.create_body()
        .dynamic()
        .at(0, 10, 0)
        .sphere(0.5f)
        .bounciness(0.7f)
        .build();

    world.on_collision = [](const CollisionEvent& e) {
        std::cout << "碰撞！冲量: " << e.impulse << "\n";
    };

    for (int i = 0; i < 120; ++i) {
        world.step(1.0f / 60.0f);
        std::cout << "球的位置: " << ball.position.y << "\n";
    }

    return 0;
}
