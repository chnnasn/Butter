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
        .box(20, 0.5f, 20)
        .friction(0.9f)
        .build();

    for (int i = 0; i < 5; ++i) {
        world.create_body()
            .dynamic()
            .at(static_cast<float>(i - 2), 2.0f + static_cast<float>(i) * 1.2f, 0)
            .box(0.5f, 0.5f, 0.5f)
            .friction(0.6f)
            .bounciness(0.1f)
            .build();
    }

    for (int step = 0; step < 300; ++step) {
        world.step(1.0f / 60.0f);
        if (step % 30 == 0) {
            std::cout << "step " << step << ", dynamic bodies: "
                      << world.dynamic_body_count() << "\n";
        }
    }

    return 0;
}
