# Butter

> **Butter** - Make Physics Smooth Again!
> A modern, fluent, intuitive C++20 header-only physics engine.

[中文](README.zh-CN.md)

## Introduction

Butter aims for a C#/TypeScript-like fluent API while keeping C++ performance and zero-overhead abstraction.

## Quick Start

```cpp
#include <butter/butter.h>
#include <iostream>

using namespace butter;
using namespace butter::math;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};

    world.create_body()
        .static_body()
        .at(0, 0, 0)
        .box(50, 0.5f, 50)
        .friction(0.8f)
        .build();

    auto& ball = world.create_body()
        .dynamic()
        .at(0, 10, 0)
        .sphere(0.5f)
        .bounciness(0.7f);

    for (int i = 0; i < 120; ++i) {
        world.step(1.0f / 60.0f);
        std::cout << ball.position.y << "\n";
    }
}
```

## Build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Tests

Tests use lightweight assertions and do not depend on third-party frameworks.

| Test | Coverage |
| --- | --- |
| `test_math` | Vectors, matrices, quaternions, AABB |
| `test_body` | Builder, mass, impulse, colliders |
| `test_collision` | Sphere-sphere, sphere-box, box-box, separation |
| `test_world` | Gravity, raycast, resting ball |
| `test_features` | Impulse-at-point, collision filtering, damping |
| `test_pbd` | PBD gravity, contact projection, friction and stacks |
| `test_triggers` | Trigger enter/exit, direction and destruction cleanup |
| `test_shapes` | Capsule, convex and indexed-triangle narrow phase |
| `test_broadphase` | Spatial-hash broadphase, AABB query and large proxies |

Run all tests:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

## Modules

- `butter/math` - Vectors, matrices, quaternions, transforms, AABB
- `butter/core` - World, body, builder, material, events
- `butter/shapes` - Sphere, box, capsule, convex, indexed mesh, collision detection
- `butter/physics2d` - 2D circle, box, polygon, gravity and impulse contacts
- `butter/physics3d` - Explicit 3D facade over the existing Butter API
- `butter/constraints` - Distance, spring, hinge joints
- `butter/query` - Raycast and overlap queries

## Broadphase, PBD, triggers and collision shapes

World enables a dynamic spatial-hash broadphase by default. Aggregate body
AABBs are inserted into cells, then exact AABB and narrow-phase tests remove
false positives. Very large ground-like proxies automatically use a bounded
large-proxy fallback. Tune it per scene:

```cpp
World::Config config;
config.enable_broadphase = true;
config.broadphase_cell_size = 2.0f;
config.broadphase_fat_margin = 0.05f;
World world(config);

auto& hull = world.create_body().dynamic()
    .convex({{-1, 0, -1}, {1, 0, -1}, {0, 1, 0}, {0, 0, 1}})
    .build();

std::vector<Vec3> vertices{{-10, 0, -10}, {10, 0, -10},
                           {10, 0, 10}, {-10, 0, 10}};
std::vector<MeshCollider::Triangle> triangles{{0, 1, 2}, {0, 2, 3}};
world.create_body().static_body().mesh(std::move(vertices), std::move(triangles)).build();
```

`ConvexCollider` treats its point set as a convex hull through support mapping
(GJK/EPA). `MeshCollider` uses indexed triangle narrow phase for sphere, capsule,
box, convex and mesh pairs. Trigger callbacks are transition-based:
`event.is_enter` is emitted on entry and `event.is_exit` on exit; persistent
overlaps do not generate a callback every frame.

`broadphase_candidate_count()` reports the number of body pairs emitted by
the latest broadphase pass. `broadphase_max_cells_per_body` bounds how many
grid cells one body may occupy, so a large floor cannot flood the hash table.

PBD and trigger lifecycle can be configured explicitly:

```cpp
World::Config config;
config.solver_mode = World::SolverMode::PBD;
config.position_iterations = 8;
World world(config);

world.on_trigger = [](const TriggerEvent& event) {
    if (event.is_enter) std::cout << "enter\n";
    if (event.is_exit)  std::cout << "exit\n";
};

world.create_body().static_body().at(0, 1, 0)
    .sphere(3.0f).trigger().build();
```

Capsule tests use a finite segment plus radius (including capsule-box and
capsule-mesh cases). Convex data is solved with GJK/EPA. Mesh collision uses
indexed triangles with closest-point, SAT and GJK tests; invalid or degenerate
triangles are ignored. Collider offsets and rotated OBB AABBs are evaluated in
body-local transforms. Persistent trigger overlaps remain quiet after enter.

The spatial hash is currently a body-level broadphase. A mesh's narrow phase
walks its indexed triangles; an internal mesh BVH is a planned optimization
for very large meshes.

## Verification

```bash
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The current implementation passes all 9 test binaries. On Windows, inspect
the PBD crate stack and explosion demo with:

```powershell
cmake --build build --config Release --target exploding_crates
.\build\examples\Release\exploding_crates.exe
```

The demo settles a gravity-driven stack before an automatic blast. Press
`Space` to repeat the blast, `R` to reset, and `P` to pause.

## 2D module

The 2D API lives beside the 3D API under `butter::physics2d`, so it can be
adopted without changing existing 3D code:

```cpp
#include <butter/physics2d/butter2d.h>
using namespace butter::physics2d;

World world;
world.create_body().static_body().at(0, -1).box(20, 1).build();
auto& ball = world.create_body().dynamic().at(0, 5).circle(0.5f).build();
for (int i = 0; i < 120; ++i) world.step();
```

The 2D module includes circle, box, capsule, convex polygon and indexed
triangle mesh narrow phase, spatial-hash broadphase, PBD/impulse projection,
friction, angular motion, sleeping, collision filtering, enter/exit triggers,
raycast/AABB queries, distance constraints and a `physics2d_playground`
example. `bench_2d` measures a 100-body stack. CCD and a dedicated internal
mesh BVH remain future optimizations.

## Examples

```bash
cmake -S . -B build -DBUTTER_BUILD_EXAMPLES=ON
cmake --build build
./build/examples/hello_butter
./build/examples/falling_boxes
./build/examples/chain
./build/examples/playground
```

The GLFW-based exploding crate stack has its own documentation:
[examples/README.md](examples/README.md).

```bash
cmake -S . -B build -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --target exploding_crates
./build/examples/exploding_crates
```

On Windows, `E:/Github/glfw` is auto-detected when present; otherwise pass
`-DBUTTER_GLFW_DIR=<your GLFW source path>`. GLFW is linked as an external
dependency and is never vendored into Butter.

## Roadmap

### Done

- C++20 header-only core
- Math, rigid body, builder, world, events
- Sphere, box, capsule collision; box-box SAT
- Impulse solver, friction, position correction
- Contact angular velocity so spinning bodies are slowed by friction
- Collision filtering, linear/angular damping, sleeping
- PBD solver mode, dynamic spatial-hash broadphase, trigger enter/exit state
- Capsule, convex (GJK/EPA) and indexed triangle mesh narrow-phase collision
- Explosion + dynamic fracture demo

### Near-term

- CCD to eliminate high-speed tunneling
- Better stacking stability and contact quality
- More joints: slider, fixed, motor

### Mid-term

- Multithreaded solver
- Basic soft body, cloth, spring bones
- Vehicle and character controllers

### Long-term

- Fluids and particles
- Soft-body FEM
- Visual debugger and editor

## License

This project is licensed under the [MIT License](LICENSE).
