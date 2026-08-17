#include <butter/butter.h>

#include <chrono>
#include <iostream>

using namespace butter;
using namespace butter::math;
using Clock = std::chrono::steady_clock;

int main() {
    const SphereCollider sphere(0.5f);
    const BoxCollider box(1, 1, 1);
    const Transform sphere_transform{{0, 1.2f, 0}};
    const Transform box_transform{{0, 0, 0}};
    ContactInfo contact;

    constexpr int iterations = 1'000'000;
    const auto start = Clock::now();
    for (int i = 0; i < iterations; ++i) {
        contact.reset();
        Collider::test(sphere, sphere_transform, box, box_transform, contact);
    }
    const auto end = Clock::now();

    const double seconds = std::chrono::duration<double>(end - start).count();
    std::cout << iterations << " sphere-box tests in " << seconds << " s ("
              << static_cast<double>(iterations) / seconds << " tests/s)\n";
    return 0;
}
