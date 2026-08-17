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
    {
        World world;
        world.gravity = {0, 0, 0};

        auto& body = world.create_body()
            .dynamic()
            .at(0, 0, 0)
            .mass(1.0f)
            .box(0.5f, 0.5f, 0.5f)
            .build();

        body.apply_impulse_at_point({0, 1, 0}, {0.5f, 0, 0});

        check(body.velocity.y > 0.0f, "apply_impulse_at_point changes linear velocity");
        check(body.angular_velocity.z > 0.0f, "apply_impulse_at_point changes angular velocity");
    }

    {
        World world;
        world.gravity = {0, 0, 0};
        bool collided = false;
        world.on_collision = [&](const CollisionEvent&) { collided = true; };

        auto& a = world.create_body().dynamic().at(0, 0, 0).sphere(0.5f).build();
        auto& b = world.create_body().dynamic().at(0, 0.5f, 0).sphere(0.5f).build();

        a.collision_group = 1;
        a.collision_mask = 0;
        b.collision_group = 2;
        b.collision_mask = 0;

        for (int i = 0; i < 2; ++i) world.step(1.0f / 60.0f);
        check(!collided, "collision filter blocks non-matching groups");

        collided = false;
        a.collision_mask = 0xFFFFFFFF;
        b.collision_mask = 0xFFFFFFFF;
        for (int i = 0; i < 2; ++i) world.step(1.0f / 60.0f);
        check(collided, "collision filter allows matching groups");
    }

    {
        World world;
        world.gravity = {0, 0, 0};

        auto& body = world.create_body().dynamic().at(0, 0, 0).mass(1.0f).build();
        body.velocity = {10, 0, 0};
        body.linear_damping = 0.5f;

        world.step(1.0f / 60.0f);
        check(body.velocity.x < 10.0f, "linear damping reduces velocity");
    }

    if (failures == 0) {
        std::cout << "test_features: OK\n";
        return 0;
    }
    std::cerr << "test_features: " << failures << " failure(s)\n";
    return 1;
}
