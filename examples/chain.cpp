#include <butter/butter.h>

#include <iostream>

using namespace butter;
using namespace butter::math;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};

    auto& anchor = world.create_body()
        .static_body()
        .at(0, 5, 0)
        .build();

    std::vector<Body*> links;
    links.push_back(&anchor);

    for (int i = 1; i <= 6; ++i) {
        auto& link = world.create_body()
            .dynamic()
            .at(0, 5.0f - static_cast<float>(i), 0)
            .mass(1.0f)
            .sphere(0.2f)
            .build();
        links.push_back(&link);
    }

    for (std::size_t i = 1; i < links.size(); ++i) {
        const Vec3 anchor_point = (*links[i - 1]).position;
        world.add_joint<HingeJoint>(*links[i - 1], *links[i], anchor_point);
    }

    for (int step = 0; step < 240; ++step) {
        world.step(1.0f / 60.0f);
        if (step % 60 == 0) {
            std::cout << "link 6 y = " << links[6]->position.y << "\n";
        }
    }

    return 0;
}
