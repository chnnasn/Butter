#include <butter/physics2d/butter2d.h>
#include <chrono>
#include <iostream>
int main() {
    butter::physics2d::World world;
    world.create_body().static_body().at(0, -1).box(100, 1).build();
    for (int i = 0; i < 100; ++i) world.create_body().dynamic().at(float(i % 10), 1.0f + float(i / 10)).circle(0.25f).build();
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 600; ++i) world.step();
    const auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    std::cout << "2D 100 bodies / 600 steps: " << elapsed << " ms\n";
}
