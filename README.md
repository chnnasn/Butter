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

Run all tests:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

## Modules

- `butter/math` - Vectors, matrices, quaternions, transforms, AABB
- `butter/core` - World, body, builder, material, events
- `butter/shapes` - Sphere, box, capsule, convex, collision detection
- `butter/constraints` - Distance, spring, hinge joints
- `butter/query` - Raycast and overlap queries

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
- Explosion + dynamic fracture demo

### Near-term

- CCD to eliminate high-speed tunneling
- Faster broadphase (BVH / grid)
- Better stacking stability and contact quality
- More joints: slider, fixed, motor

### Mid-term

- Convex hull and mesh collision
- Multithreaded solver
- Basic soft body, cloth, spring bones
- Vehicle and character controllers

### Long-term

- Fluids and particles
- Soft-body FEM
- Visual debugger and editor

## License

This project is licensed under the [MIT License](LICENSE).
