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
    // A stationary overlap should produce one enter, not a fresh event every
    // fixed step. Destroying the other body must close the pair with exactly
    // one exit before the body is removed.
    {
        World world;
        world.gravity = {0, 0, 0};

        int enters = 0;
        int stays = 0;
        int exits = 0;
        Body* trigger_body = nullptr;
        Body* other_body = nullptr;
        world.on_trigger = [&](const TriggerEvent& event) {
            check(&event.trigger == trigger_body, "trigger callback identifies trigger body");
            check(&event.other == other_body, "trigger callback identifies other body");
            if (event.is_enter) ++enters;
            if (event.is_stay) ++stays;
            if (event.is_exit) ++exits;
        };

        trigger_body = &world.create_body()
            .static_body()
            .at(0, 0, 0)
            .sphere(1.0f)
            .trigger()
            .build();
        other_body = &world.create_body()
            .dynamic()
            .at(0, 0, 0)
            .sphere(0.5f)
            .build();

        world.step();
        world.step();
        check(enters == 1, "overlap emits one enter event");
        check(stays == 0, "persistent overlap does not emit repeated events");
        check(exits == 0, "persistent overlap has no exit yet");

        world.destroy(*other_body);
        world.step();
        check(exits == 1, "destroyed body emits one trigger exit");
    }

    // Moving an overlapping body out of a trigger emits one exit, and moving
    // it back in starts a fresh enter.  Repeated separated steps stay quiet.
    {
        World world;
        world.gravity = {0, 0, 0};
        int enters = 0;
        int exits = 0;
        Body* other = nullptr;
        world.on_trigger = [&](const TriggerEvent& event) {
            if (event.is_enter) ++enters;
            if (event.is_exit) ++exits;
        };
        world.create_body().static_body().sphere(1.0f).trigger().build();
        other = &world.create_body().dynamic().sphere(0.5f).build();

        world.step();
        other->position = {5, 0, 0};
        world.step();
        world.step();
        check(enters == 1 && exits == 1,
              "separation emits exactly one trigger exit");

        other->position = {0, 0, 0};
        world.step();
        check(enters == 2 && exits == 1,
              "re-entry emits a fresh trigger enter");
    }

    // If only the second body owns the trigger collider, event direction must
    // still point at that body rather than following pair iteration order.
    {
        World world;
        world.gravity = {0, 0, 0};
        Body* trigger_body = nullptr;
        Body* other_body = nullptr;
        bool correctly_oriented = false;
        world.on_trigger = [&](const TriggerEvent& event) {
            correctly_oriented = (&event.trigger == trigger_body &&
                                  &event.other == other_body && event.is_enter);
        };

        other_body = &world.create_body()
            .static_body()
            .at(0, 0, 0)
            .sphere(1.0f)
            .build();
        trigger_body = &world.create_body()
            .dynamic()
            .at(0, 0, 0)
            .sphere(0.5f)
            .trigger()
            .build();

        world.step();
        check(correctly_oriented, "trigger event direction follows trigger material");
    }

    // Static trigger volumes should also observe static overlaps; the solid
    // static-static broad-phase cull must not suppress trigger pairs.
    {
        World world;
        world.gravity = {0, 0, 0};
        int directional_enters = 0;
        world.on_trigger = [&](const TriggerEvent& event) {
            if (event.is_enter) ++directional_enters;
        };
        world.create_body()
            .static_body()
            .at(0, 0, 0)
            .sphere(1.0f)
            .trigger()
            .build();
        world.create_body()
            .static_body()
            .at(0, 0, 0)
            .sphere(0.5f)
            .trigger()
            .build();
        world.step();
        check(directional_enters == 2,
              "both static trigger bodies receive directional enter events");
    }

    // PBD rebuilds contacts several times per fixed step. Trigger state is
    // still edge-triggered rather than one enter per projection iteration.
    {
        World world;
        world.gravity = {0, 0, 0};
        world.config.solver_mode = World::SolverMode::PBD;
        world.config.position_iterations = 8;
        int enters = 0;
        int exits = 0;
        world.on_trigger = [&](const TriggerEvent& event) {
            if (event.is_enter) ++enters;
            if (event.is_exit) ++exits;
        };
        world.create_body().static_body().sphere(1.0f).trigger().build();
        auto& body = world.create_body().dynamic().sphere(0.5f).build();
        world.step();
        world.step();
        check(enters == 1 && exits == 0,
              "PBD contact rebuilds do not duplicate trigger enter");
        body.position = {5, 0, 0};
        world.step();
        check(exits == 1, "PBD trigger separation emits one exit");
    }

    if (failures == 0) {
        std::cout << "test_triggers: OK\n";
        return 0;
    }
    std::cerr << "test_triggers: " << failures << " failure(s)\n";
    return 1;
}
