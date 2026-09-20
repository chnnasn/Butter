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
example. `bench_2d [count] [boxes]` checks landing correctness before reporting active and sleeping costs separately. A dedicated internal mesh BVH remains a future optimization. The 2D module now
includes conservative swept CCD for circles, boxes and convex polygons.

### Engine integration with explicit fixtures

`World::create_empty_body()` creates a dynamic body without the builder's default
shape. Use `add_fixture(body, definition)` for stable, independently removable
fixtures with local transforms, materials, trigger flags, and category/mask bits.
Do not add explicit fixtures to a builder-created single-shape body.

`contact_filter` adds an optional per-fixture filter. `on_contact(a, b, enter)`
reports transitions between explicit fixtures, including solid contacts and
triggers. Sleeping bodies retain contacts; fixture/body destruction emits exits
before releasing fixtures. Callbacks run with `locked() == true`: queue world
mutations until the callback or `step()` returns. Existing builder bodies retain
the `on_trigger` API. Mixed builder/explicit pairs participate in collision solving,
but do not publish transition callbacks.

`destroy_fixture`, `destroy_joint`, and `destroy_body` release individual objects;
body destruction also removes connected joints. `Body::user_data` and
`Fixture::user_data` are explicit 64-bit values. Engines may use `force`, `torque`,
`fixed_rotation`, and `BodyType::Kinematic`; engines that edit mass/type directly
must keep inverse mass/inertia consistent. Fixture density is metadata: mass and
inertia for explicit compound bodies are supplied by the engine integration.

Distance joints support rotated local anchors and `collide_connected`; soft
joints use `spring_stiffness` and `damping`. Deep penetration correction is bounded
by `Config::max_position_correction`. Ray queries test circles and convex polygons
(including rotated boxes) rather than their AABBs, and ignore origins inside them.
See the 2D CCD section below for supported shapes and remaining boundaries.

### 2D continuous collision detection (CCD)

2D worlds enable CCD by default for dynamic bodies against static and kinematic
bodies. Set `body.bullet = true` (or `.bullet()` on a builder) to additionally
sweep dynamic/dynamic pairs when either body is a bullet. `config.ccd.enabled =
false` restores discrete stepping. The 3D solver is unchanged.

CCD uses swept bounds and separating-plane conservative advancement for circles,
boxes and convex polygons. Both bodies' translation and unwrapped angular travel
are included, including rotating local fixture offsets. It advances to the
earliest time of impact, applies restitution/friction impulses, then sweeps the
remaining time again. Forces and damping are integrated once per outer step.
Fixture masks, the custom contact filter, and non-colliding connected joints are
respected. Explicit fixtures emit contact transitions for CCD impacts, including
an enter/exit pair for a bounce that separates within the same step.

`config.ccd.max_impacts` (default 32) bounds impacts **per moving body**;
`max_iterations` (64) and `tolerance` (0.0001 world units) bound conservative
advancement. Contacts at the same TOI are processed together. Budget exhaustion,
nonconvergence, or repeated zero-time hits conservatively clamp only the involved
motions for the rest of the step; unrelated bodies continue through CCD.
`world.ccd_statistics()` separates these failure reasons and records body/collider
indices plus elapsed/local remaining time. `remaining_time` now means the largest
locally clamped interval, not discarded world time. Monitor these diagnostics:
a local clamp is a conservative fallback, not equivalent-quality free motion.

Persistent contacts use clipped two-point polygon manifolds, feature/anchor
matching, accumulated normal/friction impulses, and warm starting. Velocity and
position solving are separate. Contact/joint position corrections are checked
against static and kinematic obstacles, including thin walls absent from the
original contact candidates. Connected dynamic bodies wake and sleep together;
a shared static floor does not join otherwise independent groups.

Capsules and meshes still use discrete collision detection. Sensors retain
endpoint overlap semantics and do not block CCD motion. Authored teleports are
not swept, and initial penetrations require discrete recovery. Polygons must be
convex and nondegenerate, and rotations are unwrapped (a full turn is `2*pi`, not
zero). Position guards cover the same convex shapes as CCD.

See [TomCat zero-time CCD follow-up](docs/ccd-zero-time-followup.md) for the original
fixture-based stack reproduction and independent static-environment timelines.

See [2D correctness and validation](docs/physics2d-correctness.md) for first-error
reproductions, diagnostic callbacks, stress-test commands, timing boundaries,
and remaining limitations. Historical speed ratios are not a comparison target:
verify physical time, penetration, sleep, and matching activity first.

`test_2d_ccd` includes a tunneling negative control, thin-wall circle/box/polygon
impacts, moving kinematics, bullet pairs, angular/offset sweeps, multiple rebounds,
filters, sleeping bodies, contact transitions, budget fallback, and a deterministic
dense-time collision oracle. Its checks stay enabled in Release builds.

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

- Extend CCD to 3D, capsules and meshes
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
