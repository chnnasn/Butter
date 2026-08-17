#include <butter/butter.h>

#include <chrono>
#include <iostream>

using namespace butter;
using namespace butter::math;
using Clock = std::chrono::steady_clock;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};
    world.create_body()
        .static_body()
        .at(0, 0, 0)
        .box(20, 0.5f, 20)
        .build();

    for (int i = 0; i < 100; ++i) {
        world.create_body()
            .dynamic()
            .at(static_cast<float>(i % 10) - 4.5f,
                2.0f + static_cast<float>(i / 10) * 0.8f,
                static_cast<float>((i / 10) % 2))
            .sphere(0.3f)
            .mass(1.0f)
            .build();
    }

    constexpr int steps = 600;
    const auto start = Clock::now();
    for (int i = 0; i < steps; ++i) {
        world.step(1.0f / 60.0f);
    }
    const auto end = Clock::now();

    const double seconds = std::chrono::duration<double>(end - start).count();
    std::cout << steps << " steps for " << world.dynamic_body_count()
              << " bodies in " << seconds << " s\n";
    return 0;
}
