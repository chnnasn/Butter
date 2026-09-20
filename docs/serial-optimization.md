# Serial 2D solver optimization

[中文](serial-optimization.zh-CN.md)

This follow-up builds on `ea8ef90`. It preserves the CCD impact budget, geometric
tolerance, velocity sleep thresholds and solver iteration count. No solver
threads were introduced. The historical TomCat report in the README remains a
separate measurement of `ea8ef90`.

## What kept the tall stacks awake?

Extending the original 1,000-box scene to 60 seconds showed that the old version
already settled before 15 seconds. At 5 seconds its maximum linear speed was
about 13.2 and at 10 seconds about 0.97. The six-second benchmark ended while
the tall stack was still moving; increasing sleep thresholds would hide motion.
The optimized version also settles before 15 seconds and remains asleep through
60 seconds. The independent time probe now advances exactly one X unit per
step, avoiding accumulated floating-point position error in long observations.

## Changes and invalidation rules

- Position correction queries cached static/kinematic AABBs with a bound covering
  the entire translation and angular interval. Rotation padding is bounded by
  `radius * min(abs(angle), 2)`. Distant obstacles are rejected before SAT or
  conservative advancement; nearby candidates still run the original sweep.
  Obstacle bounds are refreshed after integration, before position solving.
- The 2D grid uses retained flat cell and pair buffers instead of rebuilding
  hash containers and sets. Candidate pairs are sorted and deduplicated. Bounds
  and body types are checked on refresh; unchanged bounds reuse the pair list.
  The event pass reuses post-integration candidates because the intervening
  velocity solve does not move bodies.
- Contact geometry reuses exact detection-transform, shape and fixture-offset
  snapshots. A changed transform recomputes geometry; there is no approximate
  rotating-manifold shortcut. Friction and velocity targets are rebuilt from
  current inputs. Persistent event sets retain entries instead of reallocating
  every contact each step. Discrete proximity constraints remain distinct from
  actual contact events when CCD is disabled.
- Cached contacts detect edits to public transforms, geometry and offsets at
  step boundaries. Filtering is reevaluated. Removing support fixtures/bodies or
  joints wakes affected connected bodies; external force or velocity wakes a
  sleeping body. Destruction clears contact references before freeing objects.
  Mutations inside callbacks must still be deferred until `step()` returns.
- Connected-group topology and its buffers are reused while bodies and contact/
  joint edges are unchanged. Sleeping groups retain contacts but are omitted
  from iterative velocity and position work. A shared floor still does not join
  independent groups, and sleep thresholds were not increased.
- CCD retains caller-owned workspace and selects the sweep axis using spatial
  spread, avoiding quadratic X-interval enumeration in tall sparse layouts.
  Bounds are rebuilt after impulses because trajectories change. Existing
  bullet/kinematic coupling and the static-environment partition limit remain;
  this is not a general new island TOI scheduler.
- Polygon AABB construction no longer allocates a world-vertex array.

## Diagnostics and reproduction

`World::step_statistics()` adds cache and sleep/group timings, active constraint
count, manifold reuse, position correction/candidate/sweep counts, awake bodies,
moving/settling groups, maximum speed, angular speed, penetration and correction.
Sleep time includes wake propagation; cache maintenance is excluded from the
detection stage. `CcdStatistics::bounds_tests` counts spatial interval candidates.
Timing and counters describe the whole world, including diagnostic probes.

```powershell
cmake --build build --config Release --parallel 4
ctest --test-dir build -C Release --output-on-failure
build/tests/Release/test_2d_tomcat_stack.exe 500 boxes 3600
build/tests/Release/test_2d_tomcat_stack.exe 1000 boxes 3600
build/tests/Release/test_2d_tomcat_stack.exe 1000 circles
build/benchmarks/Release/bench_2d.exe 1000 boxes
build/benchmarks/Release/bench_2d.exe 1000 boxes awake
```

The last benchmark explicitly wakes bodies each frame. Compare it only with the
same forced-awake workload. Both benchmark modes report active-body steps and
ordinary C++ allocation counts; allocation timing is included in stage timings,
not an additional additive stage. The TomCat regression checks geometry against
the ground every frame, real probe motion, simulation time and all CCD failure
counters. Its optional fourth argument is total frames, including 60 warmup steps.

Release and Debug CTest each pass 17/17. New regression coverage compares grid
candidates and trigger events with brute force after edits and replacement,
checks cache invalidation and wake behavior, preserves discrete event semantics,
and checks the coupled CCD index on 200 vertically separated bullet bodies.
The existing thin-wall, angular re-entry and long-resting-stack tests also pass.

## Measurement boundaries

The adapter is from [TomCat `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter).
An isolated harness overlays either `ea8ef90` or the current headers; it does not
update the engine branch or its dependency pointer. Old and new binaries run
sequentially with alternating order, six rounds, discard round zero, 60 warmup
and 300 measured steps, 1/60 s, MSVC Release `/O2`, 1,000 bodies. Compilation and
tests finish before timing. No new Box2D run or comparison ratio is claimed.

Trajectories can differ even when the final awake counts agree. The timing table
is a same-scene cost comparison, not proof of identical trajectories or a
guarantee that every step on every machine fits 16.67 ms. Forced-awake allocation
results provide a separate workload with explicitly matched activity.

## Results

Six alternating rounds, median after discarding round zero, milliseconds per step:

| 1,000-body scene | `ea8ef90` measured in this run | Serial follow-up | Final awake bodies (both) |
| --- | ---: | ---: | ---: |
| Separated motion | 1.350690 | 0.220543 | 1,000 |
| Falling circles | 6.888180 | 1.469910 | 0 |
| Box stacks | 55.151900 | 13.968000 | 1,000 |

The five retained box samples range from 13.848 to 14.0705 ms. The mean step
time is below 16.67 ms in this workload; no tail-latency guarantee is implied.
All 36 timed runs report 5 measured simulation seconds, zero CCD limits and zero
nonfinite/below-ground centers. Box trajectory checksums differ between versions
(1881.63 vs 1914.06), so these are not identical-trajectory speedup claims.
[Raw six-round samples](benchmarks/serial-20260921.csv) include the discarded round;
`baseline` denotes `ea8ef90`, and `retest` denotes this serial follow-up.

The stricter native 500/1,000-box tests check body geometry against the floor
throughout 3,600 steps. They first all sleep at 9.56667/12.1667 seconds and remain
asleep through 60 seconds, with zero CCD limits, nonconvergence, repeated zero
hits or budget exhaustion. These are checks of the supplied scenes, not a proof
of correctness for arbitrary geometry or extreme inputs.

A separate 300-step native profile of 1,000 boxes reports 300,000 active-body
steps and 12.3505 ms/step: CCD 0.7275, detection 4.9264, velocity 1.7252, position
4.1093, cache 0.5225 and sleeping/group work 0.3327 ms. Other step bookkeeping
accounts for the remainder. Only 280,993 obstacle candidates were checked for
4,159,957 point corrections. A separate instrumented adapter run spends about
1.461 ms/step inside contact-filter/contact-event callbacks; that is an inclusive
subset of step time and incurs timer overhead, not an extra additive stage.

Ordinary allocation diagnostic, 1,000 single-layer boxes over 360 steps:

| Mode / phase | Activity matched in both versions | Old allocations | New allocations | Old / new total ms |
| --- | --- | ---: | ---: | ---: |
| Natural sleep: active | 49 frames, 49,000 active-body steps | 1,538,150 | 25,244 | 270.611 / 84.111 |
| Natural sleep: resting | 311 frames, 0 active-body steps | 8,096,540 | 622 | 1,087.820 / 226.242 |
| Forced awake | 360 frames, 360,000 active-body steps | 13,366,690 | 25,866 | 2,375.860 / 869.397 |

These allocation-instrumented runs are separate diagnostics, not six-round
timing medians. Source harness outputs are retained locally under
`build/tomcat-retest/serial-comparison.txt`, `build/tomcat-serial-*-long.txt`,
`build/tomcat-serial-stages.txt`, `build/tomcat-serial-adapter.txt` and
`build/alloc-{baseline,serial}-{natural,awake}.txt`. The adapter callback profile
uses scope timers in a copied adapter; the primary timing binaries do not.
