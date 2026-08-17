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
    ContactInfo contact;

    {
        const SphereCollider a(0.5f);
        const SphereCollider b(0.5f);
        const Transform ta{{0, 0, 0}};
        const Transform tb{{0, 0.8f, 0}};
        check(Collider::test(a, ta, b, tb, contact), "sphere-sphere reports hit");
        check(std::abs(contact.penetration - 0.2f) < 1.0e-6f,
              "sphere-sphere penetration depth");
    }

    {
        const SphereCollider sphere(0.5f);
        const BoxCollider box(1, 1, 1);
        const Transform sphere_transform{{0, 1.2f, 0}};
        const Transform box_transform{{0, 0, 0}};
        check(Collider::test(sphere, sphere_transform, box, box_transform, contact),
              "sphere-box reports hit");
        check(contact.normal.y < -0.9f, "sphere-box normal points from sphere to box");
    }

    {
        const BoxCollider a(1, 1, 1);
        const BoxCollider b(1, 1, 1);
        const Transform ta{{0, 0, 0}};
        const Transform tb{{1.5f, 0, 0}};
        check(Collider::test(a, ta, b, tb, contact), "box-box reports hit");
        check(std::abs(contact.normal.x - 1.0f) < 1.0e-6f,
              "box-box normal points from A to B");
    }

    {
        const SphereCollider sphere(0.5f);
        const BoxCollider box(1, 1, 1);
        const Transform sphere_transform{{0, 3.0f, 0}};
        const Transform box_transform{{0, 0, 0}};
        check(!Collider::test(sphere, sphere_transform, box, box_transform, contact),
              "separated shapes do not collide");
    }

    if (failures == 0) {
        std::cout << "test_collision: OK\n";
        return 0;
    }
    std::cerr << "test_collision: " << failures << " failure(s)\n";
    return 1;
}
