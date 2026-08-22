#include <butter/physics2d/butter2d.h>
#include <iostream>
using namespace butter::physics2d;
int main() {
    World world;
    world.create_body().static_body().at(0, -1).box(20, 1).friction(0.8f).build();
    for (int i = 0; i < 8; ++i) world.create_body().dynamic().at(float(i % 4) - 1.5f, 1.0f + float(i / 4) * 1.1f).box(0.45f, 0.45f).friction(0.7f).build();
    for (int frame = 0; frame < 180; ++frame) world.step();
    std::cout << "2D playground simulated " << world.body_count() << " bodies\n";
}
