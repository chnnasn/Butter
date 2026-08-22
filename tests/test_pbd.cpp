#include <butter/butter.h>

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace butter;
using namespace butter::math;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};
    world.config.solver_mode = World::SolverMode::PBD;
    world.config.fixed_timestep = 1.0f / 120.0f;
    world.config.position_iterations = 8;
    world.config.pbd_compliance = 0.0f;

    world.create_body().static_body().at(0, 0, 0).box(10, 0.5f, 10).build();
    auto& box = world.create_body().dynamic().at(0, 5, 0)
                    .mass(1.0f).box(0.5f, 0.5f, 0.5f)
                    .friction(0.8f).bounciness(0.0f).build();

    for (int i = 0; i < 360; ++i) {
        world.step(1.0f / 120.0f);
    }

    // Ground top is y=.5 and box half-height is .5. A small tolerance allows
    // for the contact slop used by the position solver.
    assert(std::abs(box.position.y - 1.0f) < 0.08f);
    assert(std::abs(box.velocity.y) < 0.25f);

    World stack_world;
    stack_world.gravity = {0, -9.81f, 0};
    stack_world.config.solver_mode = World::SolverMode::PBD;
    stack_world.config.fixed_timestep = 1.0f / 120.0f;
    stack_world.config.position_iterations = 8;
    stack_world.config.velocity_iterations = 6;
    stack_world.config.enable_sleeping = false;
    stack_world.create_body().static_body().at(0, 0, 0).box(10, .5f, 10).build();
    for (int layer = 0; layer < 3; ++layer) {
        for (int x = 0; x < 3; ++x) {
            for (int z = 0; z < 3; ++z) {
                stack_world.create_body().dynamic()
                    .at((x - 1) * 1.02f, 2.0f + layer * 1.02f, (z - 1) * 1.02f)
                    .mass(2.0f).box(.5f, .5f, .5f).friction(.85f).build();
            }
        }
    }
    float max_speed = 0.0f;
    for (int i = 0; i < 600; ++i) {
        stack_world.step(1.0f / 120.0f);
        for (auto& body : stack_world.bodies()) {
            if (body.velocity.length() > max_speed) {
                max_speed = body.velocity.length();
            }
        }
    }
    assert(max_speed < 8.0f);

    // Convex support mapping participates in the same PBD contact path.
    World convex_world;
    convex_world.gravity = {0, -9.81f, 0};
    convex_world.config.solver_mode = World::SolverMode::PBD;
    convex_world.config.fixed_timestep = 1.0f / 120.0f;
    convex_world.config.position_iterations = 8;
    convex_world.config.velocity_iterations = 4;
    convex_world.create_body().static_body().at(0, 0, 0)
        .box(10, 0.5f, 10).build();
    const std::vector<Vec3> tetra{{-0.45f, -0.45f, -0.45f},
                                  {0.45f, -0.45f, -0.45f},
                                  {-0.45f, 0.45f, -0.45f},
                                  {-0.45f, -0.45f, 0.45f}};
    ConvexCollider convex_probe(tetra);
    BoxCollider ground_probe(10.0f, 0.5f, 10.0f);
    ContactInfo probe_contact;
    assert(Collider::test(convex_probe, Transform{{0, 0.9f, 0}},
                          ground_probe, Transform::identity(), probe_contact));
    assert(probe_contact.normal.y < 0.0f && probe_contact.penetration > 0.0f);
    probe_contact.reset();
    assert(Collider::test(ground_probe, Transform::identity(),
                          convex_probe, Transform{{0, 0.9f, 0}}, probe_contact));
    assert(probe_contact.normal.y > 0.0f && probe_contact.penetration > 0.0f);
    auto& convex = convex_world.create_body().dynamic().at(0, 4, 0)
                       .mass(1.0f).convex(tetra).friction(0.8f).build();
    int convex_collisions = 0;
    convex_world.on_collision = [&](const CollisionEvent&) { ++convex_collisions; };
    for (int i = 0; i < 360; ++i) convex_world.step(1.0f / 120.0f);
    assert(convex_collisions > 0);
    assert(convex.position.y > 0.8f && convex.position.y < 2.0f);
    assert(std::abs(convex.velocity.y) < 0.5f);

    World mesh_world;
    mesh_world.gravity = {0, -9.81f, 0};
    mesh_world.config.solver_mode = World::SolverMode::PBD;
    mesh_world.config.fixed_timestep = 1.0f / 120.0f;
    mesh_world.config.position_iterations = 6;
    mesh_world.config.velocity_iterations = 4;
    const std::vector<Vec3> mesh_vertices{{-5, 0, -5}, {5, 0, -5},
                                           {5, 0, 5}, {-5, 0, 5}};
    const std::vector<MeshCollider::Triangle> mesh_triangles{{0, 1, 2},
                                                              {0, 2, 3}};
    mesh_world.create_body().static_body().mesh(mesh_vertices, mesh_triangles).build();
    auto& mesh_ball = mesh_world.create_body().dynamic().at(0, 3, 0)
                          .sphere(0.5f).friction(0.8f).build();
    auto& mesh_convex = mesh_world.create_body().dynamic().at(2, 3, 0)
                            .convex(tetra).friction(0.8f).build();
    for (int i = 0; i < 300; ++i) mesh_world.step(1.0f / 120.0f);
    assert(mesh_ball.position.y > 0.4f && mesh_ball.position.y < 0.7f);
    assert(std::abs(mesh_ball.velocity.y) < 0.5f);
    assert(mesh_convex.position.y > 0.4f && mesh_convex.position.y < 2.0f);
    assert(std::abs(mesh_convex.velocity.y) < 0.5f);
    return 0;
}
