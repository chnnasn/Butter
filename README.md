# Butter

> **Butter** - Make Physics Smooth Again!
> A modern, fluent, intuitive C++20 header-only physics engine.

[中文](README.zh-CN.md)

## Introduction

Butter aims for a C#/TypeScript-like fluent API while keeping C++ performance and zero-overhead abstraction.

Both 3D and 2D APIs are available. The quick start below uses 3D; the latest
optimization and [benchmark results](#2d-benchmark-report-2026-09-21) concern 2D.

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
        .bounciness(0.7f)
        .build();

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
| `test_2d` | 2D bodies, collisions, queries and constraints |
| `test_2d_ccd` | Thin walls, angular/offset sweeps, compound position corrections, 4–64 vertex polygons and failure isolation |
| `test_2d_engine` | Explicit fixtures, engine integration and contact lifecycle |
| `test_3d_namespace` | Explicit 3D namespace compatibility |
| `test_2d_stability` | Large landings, long resting stacks and connected wake/sleep |
| `test_2d_tomcat_stack` | Original TomCat box parameters, time advancement and floor checks |
| `test_2d_tomcat_circles` | 1,000-circle TomCat regression (runs `test_2d_tomcat_stack 1000 circles`) |
| `test_2d_optimization` | Broadphase oracle, cache invalidation, sleeping-world fast path, kinematic motion and vertex/CCD buffers |

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

Revision `a8e4336` passes all 17 CTest cases (16 binaries) in
both Release and Debug. These are repository regressions,
separate from the external benchmark below. On Windows, inspect
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
triangle mesh narrow phase, spatial-grid broadphase, PBD/impulse projection,
friction, angular motion, sleeping, collision filtering, enter/exit triggers,
raycast/AABB queries, distance constraints and a `physics2d_playground`
example. `bench_2d [count] [circles|boxes] [awake]` checks landing correctness and
reports active and sleeping costs separately; `awake` explicitly keeps bodies
active for comparisons. This internal landing benchmark differs from the TomCat
stack benchmark below. A dedicated internal mesh BVH remains a future optimization. The 2D module now
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

## 2D benchmark report (2026-09-21)

### Latest: `a8e4336` against `9baad05`

The latest changes reduce CCD geometry allocations and reject distant obstacles
before rebuilding position-correction bounds. Potential hits still run the full
CCD guard. Earlier changes precompute velocity constraints, reuse contact geometry
and skip redundant CCD/detection in fully sleeping worlds. Solver iterations,
sleep thresholds and CCD budgets were not relaxed; the solver remains serial.

Environment and adapter: [TomCat Engine — `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter).
These local tests use an isolated copy of that adapter with Butter headers; they
do not update the engine branch's pinned dependency. Release, single thread,
1/60 s, 60 warm-up steps and 300 measured steps; six alternating old/new rounds,
discard the first and take the median of five. Timings are milliseconds per step.

| Scene | Bodies | `9baad05` | `a8e4336` |
| --- | ---: | ---: | ---: |
| Separated motion | 100 | 0.019314 | 0.018888 |
| Separated motion | 500 | 0.099302 | 0.098543 |
| Separated motion | 1,000 | 0.212034 | 0.210657 |
| Falling circles | 100 | 0.026110 | 0.027397 |
| Falling circles | 500 | 0.356339 | 0.338405 |
| Falling circles | 1,000 | 1.139450 | 1.033050 |
| Box stacks | 100 | 0.068556 | 0.059709 |
| Box stacks | 500 | 5.405940 | 4.081260 |
| Box stacks | 1,000 | 12.166300 | 9.203220 |

500/1,000 boxes take about **25%/24% less time**, and 1,000 circles about **9% less**.
Separated motion is comparable. The 100-circle scene is about 1.29 microseconds
slower (5%); this is not a universal speedup. These are Butter-to-Butter results,
not new Box2D ratios. A mean below 16.67 ms does not guarantee every step meets 60 Hz.

All nine scenes match the old version's per-step state hashes and active-body-step
counts. All timed runs advance five measured seconds with zero CCD limits, budget
exhaustion, nonconvergence, zero-time repeats, nonfinite positions or detected
below-ground centers. The stricter native 500/1,000-box tests check floor geometry
throughout 60 seconds: all bodies first sleep at **9.57/12.17 seconds** and remain
asleep. Their awake state at six seconds is not evidence of a permanent sleep failure.
These results validate the supplied scenes, not arbitrary physics workloads.

`world.step_statistics()` exposes CCD, detection, velocity, position, cache and
sleeping times, activity counts, `stationary_steps` and `projection_fast_rejections`.
In a separate 1,000-box diagnostic, position solving fell from about 4.03 to
1.49 ms/step with unchanged correction/candidate counts. Contact detection remains
a major cost; compare activity and trajectories before further performance claims.

See [implementation, invalidation rules and validation](docs/serial-optimization.md)
and [108 raw samples](docs/benchmarks/projection-followup-20260921.csv).
The raw file's `baseline9` means `9baad05`; `retest` means the code committed as `a8e4336`.
Reproduction commands are in the [examples and diagnostics guide](examples/README.md#2d-regressions-and-benchmarks).

### External Box2D reference: `52f15a5`

The supplied later TomCat report compared Butter `52f15a5` with Box2D 2.4.1 on an
i7-14650HX using Release, one thread and the same warm-up/measurement schedule:

| Scene (1,000 bodies) | Box2D | Butter `52f15a5` | Butter / Box2D |
| --- | ---: | ---: | ---: |
| Separated motion | 0.2118 ms | 0.2217 ms | 1.05× |
| Falling circles | 0.5340 ms | 1.4868 ms | 2.78× |
| Box stacks | 3.0651 ms | 13.6228 ms | 4.44× |

Separated motion was effectively on par; dense contacts remained several times
slower. The adapter used updated Butter in an isolated directory because the
branch still pinned an older dependency at the time. Different engine trajectories
and sleep schedules make these same-initial-scene costs, not equal-quality throughput.
Do not divide the latest local timings by these older Box2D timings to claim a new ratio.

<details>
<summary>Historical report: ea8ef90 and the first serial optimization</summary>

**Serial optimization follow-up:** a same-session rerun against `ea8ef90` measured
1,000-box stepping at **55.152 → 13.968 ms**, circles at **6.888 → 1.470 ms**, and
separated motion at **1.351 → 0.221 ms**. The 500/1,000-box scenes settle at about
9.57/12.17 seconds and stay asleep through 60 seconds. See the
[bilingual implementation and validation report](docs/serial-optimization.md)
for cache rules, matched-activity allocation diagnostics, raw samples and limits.
The report below preserves the earlier `ea8ef90` versus Box2D measurement.

Test environment and integration: [TomCat Engine — `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter).
This run tested **Butter `ea8ef90`** through that branch's physics adapter against
**Box2D 2.4.1** from TomCat's `main` branch. **At the time, `dev_butter` pinned an older
Butter revision**: these results used updated headers in an isolated test
directory, not an updated dependency in the engine branch.

Hardware and build: Intel Core i7-14650HX, Windows x64, MSVC 19.50.35724,
C++20 Release (`/O2 /DNDEBUG`), single-threaded workloads. Each scene uses a
1/60-second timestep, 60 warmup steps and 300 measured steps, with CCD and sleeping
enabled. Box2D uses `Step(8,3)`; the Butter adapter maps this to 8 iterations.
Six runs alternate library order; the first is discarded and the median of the
remaining five is reported. No concurrent tests or compilation ran during timing;
CPU affinity and power plans were not controlled. Timings exclude creation,
rendering, scripts and engine synchronization.

### Step time for 1,000 dynamic bodies

| Scene | Box2D | Butter | Butter / Box2D time |
| --- | ---: | ---: | ---: |
| Separated motion, no collisions | 0.206 ms | 1.108 ms | 5.38× |
| Falling circles | 0.537 ms | 5.867 ms | 10.92× |
| Box stacks | 2.967 ms | 51.237 ms | 17.27× |

The main observed improvement is in box stacking. Circle landings and separated
motion remain in roughly the same relative performance range.

| Box count | Previous run: Butter / Box2D time | This run: Butter / Box2D time |
| --- | ---: | ---: |
| 100 | 974.42× | 35.80× |
| 500 | 134.74× | 28.43× |
| 1,000 | 81.37× | 17.27× |

These cross-session ratios are observations, not controlled version speedup
factors. The engines follow different trajectories and sleep schedules, so the
results compare actual costs from the same initial scenes, not equal-quality
solver throughput.

### Correctness and stability checks

Separate diagnostics covered falling circles and box stacks at 100, 500 and
1,000 bodies (six scenes):

| Check | Result |
| --- | --- |
| Simulation time during 300 measured steps | Full 5 seconds in all six scenes |
| CCD-limited steps | 0 |
| Impact-budget exhaustion | 0 |
| CCD nonconvergence | 0 |
| Repeated zero-time collisions | 0 |
| Detected ground crossings within the platform | 0 |
| Nonfinite positions | 0 |

All circle scenes and the 100-box stack slept. The 500- and 1,000-box stacks
remained awake at the end. The ground check tests body centers within
`abs(x) < 99`; it is not a full geometry-level penetration proof. The previously
reproduced world stall, ground crossing and zero-time repeat failures did not
recur in these cases. This limited validation does not replace comprehensive
physics correctness testing.

At about **51 ms per step**, this historical `ea8ef90` result exceeded the
**16.67 ms** budget for 60 Hz. The newer measurements and sleep investigation
above supersede that performance assessment.

The supplied source, raw data and reproduction instructions are in the local
TomCat checkout at `.scratch/comparison-20260921-r2/README.md` (not published by
the branch link above). See also the repository's
[CCD follow-up](docs/ccd-zero-time-followup.md) for implementation diagnostics
and regression coverage; its local timing run is separate from this report.

</details>

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
- 2D convex CCD, persistent contact manifolds, warm starting and connected sleeping
- Original TomCat circle/box regressions and activity-aware timing diagnostics
- Cached velocity geometry, small CCD vertex buffers and sleeping-world fast path
- Conservative position-correction prefilter covering compound fixture offsets

### Near-term

- Extend CCD to 3D, capsules and meshes
- Better stacking stability and contact quality
- Reduce dense-contact detection and cache overhead with matched-activity benchmarks
- Preserve thin-wall CCD and long-term stability while optimizing serial work
- More joints: slider, fixed, motor

### Mid-term

- Multithreaded solver after serial correctness and profiling
- Basic soft body, cloth, spring bones
- Vehicle and character controllers

### Long-term

- Fluids and particles
- Soft-body FEM
- Visual debugger and editor

## License

This project is licensed under the [MIT License](LICENSE).
