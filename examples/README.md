# Exploding Crate Stack

[中文](README.zh-CN.md)

A GLFW-rendered physics demo built with Butter.

## What It Shows

- Crates fall from the air under gravity.
- They collide, tumble, and settle into a stack through real physics.
- The demo triggers a blast automatically after the initial three-second
  settling period; press `Space` to trigger it earlier.
- Crates that receive a strong enough explosion impulse fracture into smaller debris.
- Fragments and surviving crates fly outward along different trajectories, then slow down and rest.

## Build

```bash
cmake -S . -B build -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --target exploding_crates
```

## Run

```bash
./build/examples/exploding_crates
```

On Windows, run:

```bash
build/examples/Debug/exploding_crates.exe
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

Back to the main README: [../README.md](../README.md).
