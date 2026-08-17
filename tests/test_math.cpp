#include <butter/butter.h>

#include <cmath>
#include <iostream>

using namespace butter::math;

static int failures = 0;

static void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

int main() {
    const Vec3 a{1, 2, 3};
    const Vec3 b{4, 5, 6};

    check(a + b == Vec3{5, 7, 9}, "Vec3 addition");
    check(a - b == Vec3{-3, -3, -3}, "Vec3 subtraction");
    check(a * 2.0f == Vec3{2, 4, 6}, "Vec3 scalar multiplication");
    check(2.0f * a == Vec3{2, 4, 6}, "scalar multiplication friend");
    check(a.dot(b) == 32.0f, "Vec3 dot product");
    check(a.cross(b) == Vec3{-3, 6, -3}, "Vec3 cross product");
    check(std::abs(Vec3{3, 4, 0}.length() - 5.0f) < 1.0e-6f, "Vec3 length");

    Vec3 chain = Vec3::zero().set_x(1).set_y(2).set_z(3);
    check(chain == Vec3{1, 2, 3}, "Vec3 chained setters");

    auto [x, y, z] = chain;
    check(x == 1.0f && y == 2.0f && z == 3.0f, "Vec3 structured binding");

    const Vec2 v2{3, 4};
    check(std::abs(v2.length() - 5.0f) < 1.0e-6f, "Vec2 length");

    const Quat yaw = Quat::from_euler(0, 90_deg, 0);
    const Vec3 rotated = yaw.rotate(Vec3{1, 0, 0});
    check(std::abs(rotated.z + 1.0f) < 1.0e-5f, "Quat rotates X toward -Z");

    const AABB unit{{0, 0, 0}, {1, 1, 1}};
    check(unit.contains({0.5f, 0.5f, 0.5f}), "AABB contains point");
    check(!unit.contains({2, 0, 0}), "AABB rejects outside point");
    check(unit.overlaps({{0.5f, 0.5f, 0.5f}, {1.5f, 1.5f, 1.5f}}), "AABB overlaps");
    check(unit.center() == Vec3{0.5f, 0.5f, 0.5f}, "AABB center");

    if (failures == 0) {
        std::cout << "test_math: OK\n";
        return 0;
    }
    std::cerr << "test_math: " << failures << " failure(s)\n";
    return 1;
}
