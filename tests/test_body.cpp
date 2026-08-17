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
    World world;

    auto& body = world.create_body()
        .dynamic()
        .at(1, 2, 3)
        .mass(5.0f)
        .box(1, 1, 1)
        .friction(0.4f)
        .build();

    check(body.position == Vec3{1, 2, 3}, "builder sets position");
    check(std::abs(body.mass() - 5.0f) < 1.0e-6f, "builder sets mass");
    check(body.is_dynamic(), "builder creates dynamic body");
    check(body.colliders().size() == 1, "builder adds one collider");
    check(body.colliders()[0]->type() == Collider::Type::Box, "builder adds box collider");
    check(std::abs(body.colliders()[0]->material.friction - 0.4f) < 1.0e-6f,
          "builder sets material friction");

    body.apply_impulse({0, 10, 0});
    check(body.velocity.y == 2.0f, "impulse changes velocity by inverse mass");

    auto& static_body = world.create_body()
        .static_body()
        .at(0, 0, 0)
        .box(10, 0.5f, 10)
        .build();
    check(static_body.is_static(), "builder creates static body");
    check(static_body.inverse_mass() == 0.0f, "static body has zero inverse mass");

    if (failures == 0) {
        std::cout << "test_body: OK\n";
        return 0;
    }
    std::cerr << "test_body: " << failures << " failure(s)\n";
    return 1;
}
