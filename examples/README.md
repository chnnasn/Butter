# Butter examples and diagnostics

[中文](README.zh-CN.md)

Run the commands below from the repository root. Butter's 2D and 3D examples
use separate APIs; the latest 2D optimizations do not change the 3D solver.

| Target | Purpose | Extra dependency |
| --- | --- | --- |
| `hello_butter` | Minimal 3D simulation | None |
| `falling_boxes` | 3D falling bodies | None |
| `chain` | 3D joint chain | None |
| `playground` | 3D console example | None |
| `physics2d_playground` | Console 2D box simulation | None |
| `exploding_crates` | Interactive 3D explosion and fracture | GLFW and OpenGL |

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_EXAMPLES=ON
cmake --build build --config Release
```

Run `./build/examples/physics2d_playground` with a single-configuration generator,
or `./build/examples/Release/physics2d_playground.exe` with Visual Studio.
`physics2d_playground` is a console example, not an interactive renderer.

## Exploding crate stack

A GLFW-rendered 3D physics demo built with Butter.

## What It Shows

- Crates fall from the air under gravity.
- They collide, tumble, and settle into a stack through real physics.
- The demo triggers a blast automatically after the initial three-second
  settling period; press `Space` to trigger it earlier.
- Crates that receive a strong enough explosion impulse fracture into smaller debris.
- Fragments and surviving crates fly outward along different trajectories, then slow down and rest.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_EXAMPLES=ON -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --config Release --target exploding_crates
```

## Run

```bash
./build/examples/exploding_crates
```

On Windows, run:

```bash
build/examples/Release/exploding_crates.exe
```

## Controls

| Input | Action |
| --- | --- |
| `Space` | Detonate |
| `P` | Pause / resume |
| `R` | Reset the scene |
| Mouse drag | Rotate camera |
| Mouse wheel | Zoom |
| `Left` / `Right` | Slow down / speed up simulation |
| `Esc` | Quit |

## Physics Details

- Sustained blast force field pushes bodies over several frames.
- Off-center impulses give crates tumbling rotation.
- Contact velocity includes angular velocity, so spinning fragments are slowed by friction.
- Linear and angular damping let debris lose energy and settle.

## 2D benchmark results

This crate demo uses the 3D API. The separate TomCat 2D adapter benchmark,
including `a8e4336` results, historical Box2D references and correctness limits, is documented
in the [main README](../README.md#2d-benchmark-report-2026-09-21).
Its timing numbers do not describe this demo.

The latest local six-round comparison against `9baad05` measured 1,000 boxes at
**12.166 → 9.203 ms/step** (about 24% less) and 1,000 circles at
**1.139 → 1.033 ms/step** (about 9% less). The 100-circle scene increased by about
1.29 microseconds. Per-step state hashes and activity counts match across nine
scenes; these are not new Box2D ratios. Environment/adapter:
[TomCat Engine `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter).
Testing used isolated updated headers, not a dependency update to that branch.
See [bilingual validation and raw data](../docs/serial-optimization.md).

## 2D regressions and benchmarks

Enable tests and the internal benchmark explicitly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_TESTS=ON -DBUTTER_BUILD_BENCHMARKS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Visual Studio / Windows commands:

```powershell
.\build\tests\Release\test_2d_ccd.exe
.\build\tests\Release\test_2d_tomcat_stack.exe 1000 circles 360
.\build\tests\Release\test_2d_tomcat_stack.exe 500 boxes 3600
.\build\tests\Release\test_2d_tomcat_stack.exe 1000 boxes 3600
.\build\benchmarks\Release\bench_2d.exe 1000 circles
.\build\benchmarks\Release\bench_2d.exe 1000 boxes
.\build\benchmarks\Release\bench_2d.exe 1000 boxes awake
```

For single-configuration builds, omit `Release/` and `.exe`. The stack executable
takes body count, `circles`/`boxes`, then total steps (including the first 60 warm-up
steps). The 3,600-step runs observe 60 seconds, checking floor geometry, time and
CCD failures. At `a8e4336`, 500/1,000 boxes sleep at about 9.57/12.17 seconds and
stay asleep through the end. Thin-wall tests also cover rotation, compound offsets
and 4/16/17/64-vertex polygons. Release and Debug both pass all 17 CTest cases.

`bench_2d` is a separate internal landing workload, not a reproduction of the
TomCat adapter table. It reports active/sleeping phases, active-body steps, stage
times and allocations. Its `awake` mode deliberately wakes bodies every step;
compare that mode only with equivalent activity. The adapter tables exclude
rendering and scene creation. Do not use visual smoothness or FPS as a replacement
for physics time, floor checks, CCD diagnostics and sleeping validation.

Back to the main README: [../README.md](../README.md).
