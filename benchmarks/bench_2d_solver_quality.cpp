// Solver-quality comparison harness.
//
// "Is the solver good enough" is not a question a single millisecond count can
// answer. This file measures the situations that actually expose a solver's
// limits -- a tall stack, a large mass ratio, a long joint chain, a moving
// platform and a fast spin -- and reports what each (iteration, sub-step)
// budget buys, so the next millisecond can be spent where the measurements say
// it belongs rather than where it feels like it should go.
//
// Two rules keep the comparison honest, because a benchmark that can be gamed
// is worse than no benchmark:
//
//   * Accuracy is read out of the bodies' own transforms -- penetration under
//     load, anchor separation, relative slide, residual spin -- never from an
//     engine counter that a solver is free to define away.
//   * Every dynamic body is woken on every frame. Sleeping is a legitimate
//     engine feature but it is not a solver improvement, and letting a
//     configuration win by freezing its stack earlier would measure the sleep
//     threshold rather than the solve. `bench_2d` already exposes the same idea
//     as an explicit `forced_awake` mode; this file simply always uses it.
//
// Nothing here relaxes an error tolerance, lowers a CCD guard or shortens a
// simulation to look better: every configuration runs the same number of frames
// over the same scene, and the only thing that changes is how the work inside a
// frame is arranged.
//
// Usage:
//   bench_2d_solver_quality [frames]
// With no argument it runs the full iteration x sub-step sweep.

#include <butter/physics2d/butter2d.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

using namespace butter::physics2d;
using Clock = std::chrono::steady_clock;

namespace {

constexpr float kDt = 1.0f / 60.0f;

// Two solver arrangements are not (iterations, sub-steps), and they are worth
// the same kind of measurement: how compliant a contact is, and whether a
// joint is solved inside the sweep that corrects the contacts it shares a body
// with. Both are read by `config()`, which every scenario goes through, so the
// sweep can vary them without threading another argument through five scenes.
float g_stiffness = 0;
bool g_interleave = false;
float g_relaxation = 0.2f;

struct Measurement {
    double milliseconds{};
    double primary{};   // The scenario's headline accuracy number. Lower is better.
    double secondary{}; // A second, independent symptom of the same failure.
};

World::Config config(int iterations, int substeps) {
    World::Config cfg;
    cfg.solver_iterations = iterations;
    cfg.substeps = substeps;
    cfg.contact_stiffness = g_stiffness;
    cfg.interleave_joints = g_interleave;
    cfg.position_relaxation = g_relaxation;
    return cfg;
}

// Waking every dynamic body each frame also resets `sleep_counter`, so nothing
// can fall asleep for the duration of a measurement.
void force_awake(const std::vector<Body *> &bodies) {
    for (auto *body : bodies)
        if (body->type == BodyType::Dynamic)
            body->wake();
}

double elapsed_ms(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Separation of a hinge's two anchors in world space. A revolute joint has to
// hold this at zero; anything else is the joint being pulled apart.
double anchor_gap(const HingeJoint &joint) {
    const Vec2 a =
        joint.a->transform.position + rotate(joint.anchor_a, Rot(joint.a->transform.angle));
    const Vec2 b =
        joint.b->transform.position + rotate(joint.anchor_b, Rot(joint.b->transform.angle));
    return double((b - a).length());
}

// The floor is deliberately enormous. A scene whose bodies can slide off the
// end of the world measures how long ago that happened, not how well the solver
// did: every "penetration" number here would otherwise be dominated by bodies
// in free fall. It is a *static* floor either way, so a body that reaches the
// edge has already failed for reasons this file should be reporting on.
//
// Its top surface is at y = 0.
void add_ground(World &world) {
    world.create_body()
        .static_body()
        .at(0, -0.05f)
        .box(400.0f, 0.05f)
        .friction(0.6f)
        .restitution(0)
        .build();
}

// Worst penetration of any body through the floor, in millimetres. This is the
// number the position pass exists to keep at zero, and it is read straight off
// the body's own AABB.
double worst_floor_penetration(const std::vector<Body *> &bodies) {
    double worst = 0;
    for (auto *body : bodies)
        worst = std::max(worst, -double(compute_aabb(body->shape, body->transform).min.y));
    return worst * 1000.0;
}

// ---------------------------------------------------------------------------
// 1. Tall stack
//
// Twenty unit boxes in a column. The bottom contact carries twenty bodies, so a
// contact that is not solved to convergence shows up as the column settling
// below its true resting height. The headline number is how far the top of the
// column ends up below where a perfectly rigid column would hold it, which is
// the accumulated penetration of all twenty contacts at once.
// ---------------------------------------------------------------------------
Measurement tall_stack(int iterations, int substeps, int frames) {
    constexpr int layers = 20;
    constexpr float half = 0.5f;
    constexpr float pitch = 2 * half + 0.01f;

    World world(config(iterations, substeps));
    add_ground(world);
    std::vector<Body *> boxes;
    for (int i = 0; i < layers; ++i)
        boxes.push_back(&world.create_body()
                             .at(0, half + float(i) * pitch)
                             .box(half, half)
                             .friction(0.6f)
                             .restitution(0)
                             .build());

    auto start = Clock::now();
    for (int frame = 0; frame < frames; ++frame) {
        force_awake(boxes);
        world.step(kDt);
    }
    const double milliseconds = elapsed_ms(start);

    // Compression of the column under its own weight, in millimetres: how far
    // the top box ended up from where a rigid in-contact column would put it.
    // Measured against the contact height, not the drop height -- the boxes
    // start with a 0.01 gap per level, and charging the solver for the settling
    // those gaps owe would report 190 mm of "penetration" in a column that is
    // perfect.
    //
    // The deviation is signed-and-measured in both directions on purpose. A
    // one-sided `max(0, ideal - top)` would report 0 mm for a column that never
    // finished settling and is *standing on its start gaps*, which is not
    // accuracy: it is the solver failing to fall. Leaving 190 mm of its own
    // drop-height unaccounted for has to cost the configuration a number as
    // large as sinking 190 mm does.
    const double ideal_top = double(half + (layers - 1) * (2 * half));
    const double compression = std::abs(ideal_top - double(boxes.back()->transform.position.y));
    // Lateral spread: the symptom that distinguishes "sinks a little" from
    // "collapses".
    double drift = 0;
    for (auto *box : boxes)
        drift = std::max(drift, std::abs(double(box->transform.position.x)));
    return {milliseconds, compression * 1000.0, drift * 1000.0};
}

// ---------------------------------------------------------------------------
// 2. Large mass ratio
//
// Ten light boxes with a crate fifty times heavier on top. The light bodies
// have almost no inertia to push back with, so an under-converged contact lets
// the heavy one sink through them: the classic failure of solving each contact
// as an isolated unit instead of solving the coupled stack.
// ---------------------------------------------------------------------------
Measurement mass_ratio(int iterations, int substeps, int frames) {
    constexpr int light_count = 10;
    constexpr float half = 0.5f;
    constexpr float pitch = 2 * half + 0.01f;

    World world(config(iterations, substeps));
    add_ground(world);
    std::vector<Body *> bodies;
    for (int i = 0; i < light_count; ++i)
        bodies.push_back(&world.create_body()
                             .at(0, half + float(i) * pitch)
                             .box(half, half)
                             .mass(1.0f)
                             .friction(0.6f)
                             .restitution(0)
                             .build());
    auto &heavy = world.create_body()
                      .at(0, half + float(light_count) * pitch)
                      .box(half, half)
                      .mass(50.0f)
                      .friction(0.6f)
                      .restitution(0)
                      .build();

    auto start = Clock::now();
    for (int frame = 0; frame < frames; ++frame) {
        force_awake(bodies);
        heavy.wake();
        world.step(kDt);
    }
    const double milliseconds = elapsed_ms(start);
    // Steel does not compress: a crate this heavy should load the stack without
    // the stack getting any shorter. Contact height again, not drop height, and
    // two-sided for the same reason as `tall_stack`: the crate starts with a gap
    // above the stack, and a solver that lets it stop there has not earned a
    // better number than one that lets it sink.
    const double ideal_top = double(half + light_count * (2 * half));
    const double drop = std::abs(ideal_top - double(heavy.transform.position.y));
    return {milliseconds, drop * 1000.0, worst_floor_penetration(bodies)};
}

// ---------------------------------------------------------------------------
// 3. Long joint chain
//
// Forty links pinned end to end, hanging from a static anchor with no contacts
// at all. This isolates the joint solver: the chain's own weight has to be
// carried through forty coupled constraints, so anything the velocity pass
// fails to converge the position pass has to make up for. The headline number
// is the largest anchor separation once the swing has died down.
// ---------------------------------------------------------------------------
Measurement joint_chain(int iterations, int substeps, int frames) {
    constexpr int links = 40;
    constexpr float link_length = 0.5f;
    constexpr float top = 6.0f;

    World world(config(iterations, substeps));
    auto &anchor = world.create_body().static_body().at(0, top).box(0.1f, 0.1f).build();
    std::vector<Body *> chain{&anchor};
    for (int i = 1; i <= links; ++i) {
        auto &link = world.create_body()
                         .at(0, top - link_length * float(i))
                         .box(0.04f, link_length * 0.5f)
                         .build();
        link.angular_damping = 0.05f;
        chain.push_back(&link);
    }
    std::vector<HingeJoint *> joints;
    for (int i = 1; i <= links; ++i)
        joints.push_back(&world.add_hinge_joint_at(*chain[std::size_t(i - 1)],
                                                   *chain[std::size_t(i)],
                                                   {0, top - link_length * float(i) +
                                                           link_length * 0.5f}));

    double worst_gap = 0, steady_gap = 0;
    auto start = Clock::now();
    for (int frame = 0; frame < frames; ++frame) {
        force_awake(chain);
        world.step(kDt);
        double gap = 0;
        for (const auto *joint : joints)
            gap = std::max(gap, anchor_gap(*joint));
        worst_gap = std::max(worst_gap, gap);
        if (frame + 120 >= frames)
            steady_gap = std::max(steady_gap, gap);
    }
    const double milliseconds = elapsed_ms(start);
    return {milliseconds, steady_gap * 1000.0, worst_gap * 1000.0};
}

// ---------------------------------------------------------------------------
// 4. Moving platform
//
// A kinematic platform sweeping back and forth under a row of boxes. The
// platform's motion is imposed rather than solved, so every frame presents the
// contact solver with a fresh boundary condition: the boxes have to be carried
// without sinking into the platform and without creeping across it.
// ---------------------------------------------------------------------------
Measurement moving_platform(int iterations, int substeps, int frames) {
    constexpr int count = 8;
    constexpr float half = 0.5f;
    constexpr float platform_top = 0.1f;
    constexpr float platform_half_width = 8.0f;

    World world(config(iterations, substeps));
    // The ground is only a safety net: a box that reaches it has already been
    // dropped by the platform, which the sink metric reports anyway.
    world.create_body()
        .static_body()
        .at(0, -4.05f)
        .box(400.0f, 0.05f)
        .friction(0.6f)
        .restitution(0)
        .build();
    auto &platform = world.create_body()
                         .kinematic()
                         .at(0, 0)
                         .box(platform_half_width, platform_top)
                         .friction(0.9f)
                         .restitution(0)
                         .build();
    std::vector<Body *> boxes;
    std::vector<float> initial_offset;
    for (int i = 0; i < count; ++i) {
        boxes.push_back(&world.create_body()
                            .at(-float(count) * 0.55f + float(i) * 1.1f, platform_top + half)
                            .box(half, half)
                            .friction(0.9f)
                            .restitution(0)
                            .build());
        initial_offset.push_back(boxes.back()->transform.position.x - platform.transform.position.x);
    }

    double worst_sink = 0, worst_slide = 0;
    auto start = Clock::now();
    for (int frame = 0; frame < frames; ++frame) {
        // Triangular sweep at one metre per second, reversing every two seconds.
        const double phase = std::fmod(double(frame) * kDt, 4.0);
        platform.velocity.x = phase < 2.0 ? 1.0f : -1.0f;
        platform.wake();
        force_awake(boxes);
        world.step(kDt);
        for (std::size_t i = 0; i < boxes.size(); ++i) {
            worst_sink = std::max(worst_sink, double(platform_top) + double(half) -
                                                  double(boxes[i]->transform.position.y));
            // How far the box has crept along the platform relative to where it
            // was placed. A box that is properly carried keeps this near zero.
            const double offset =
                double(boxes[i]->transform.position.x - platform.transform.position.x);
            worst_slide = std::max(worst_slide, std::abs(offset - double(initial_offset[i])));
        }
    }
    const double milliseconds = elapsed_ms(start);
    return {milliseconds, worst_sink * 1000.0, worst_slide * 1000.0};
}

// ---------------------------------------------------------------------------
// 5. Fast spin
//
// A box dropped while spinning at 40 rad/s. High angular velocity is the
// hardest case for a position pass that applies one large correction: the
// contact normal swings several degrees per frame, so a correction computed
// once and applied in full is already stale when it lands. The box also carries
// far more kinetic energy than the landing needs, so penetration depth during
// the impact is the direct symptom, and the spin it is still holding at the end
// says whether the impact was resolved or merely survived.
// ---------------------------------------------------------------------------
Measurement fast_spin(int iterations, int substeps, int frames) {
    World world(config(iterations, substeps));
    add_ground(world);
    auto &spinner = world.create_body()
                        .at(0, 1.2f)
                        .box(0.5f, 0.5f)
                        .friction(0.6f)
                        .restitution(0)
                        .angular_velocity(40.0f)
                        .build();

    double worst_sink = 0;
    auto start = Clock::now();
    for (int frame = 0; frame < frames; ++frame) {
        spinner.wake();
        world.step(kDt);
        worst_sink = std::max(worst_sink, -double(compute_aabb(spinner.shape, spinner.transform).min.y));
    }
    const double milliseconds = elapsed_ms(start);
    // How much spin the landing left behind. A solver that lets the box grind
    // along the floor keeps feeding it; one that resolves the impact leaves it
    // nearly still.
    const double residual = std::abs(double(spinner.angular_velocity));
    return {milliseconds, worst_sink * 1000.0, residual * 1000.0};
}

struct Scenario {
    const char *name;
    Measurement (*run)(int, int, int);
    const char *primary;
    const char *secondary;
};

const int kIterations[] = {2, 4, 8, 16};
const int kSubsteps[] = {1, 2, 4};
constexpr int kIterationCount = int(sizeof(kIterations) / sizeof(kIterations[0]));
constexpr int kSubstepCount = int(sizeof(kSubsteps) / sizeof(kSubsteps[0]));

}  // namespace

int main(int argc, char **argv) try {
    const int frames = argc > 1 ? std::atoi(argv[1]) : 600;
    if (frames < 60 || frames > 20000)
        throw std::runtime_error("frames must be 60..20000");

    const Scenario scenarios[] = {
        {"tall_stack", tall_stack, "compression_mm", "drift_mm"},
        {"mass_ratio", mass_ratio, "heavy_drop_mm", "floor_penetration_mm"},
        {"joint_chain", joint_chain, "steady_gap_mm", "worst_gap_mm"},
        {"moving_platform", moving_platform, "sink_mm", "slide_mm"},
        {"fast_spin", fast_spin, "penetration_mm", "residual_spin_mrad_s"},
    };

    std::printf("frames=%d  (all bodies force-awake every frame; nothing else is relaxed)\n",
                frames);
    std::printf("Each cell is measured over the full %d frames of the same scene.\n", frames);

    for (const auto &scenario : scenarios) {
        Measurement table[kSubstepCount][kIterationCount];
        double best = 1e30;
        int best_iterations = 0, best_substeps = 0;
        for (int s = 0; s < kSubstepCount; ++s)
            for (int i = 0; i < kIterationCount; ++i) {
                table[s][i] = scenario.run(kIterations[i], kSubsteps[s], frames);
                if (table[s][i].primary < best) {
                    best = table[s][i].primary;
                    best_iterations = kIterations[i];
                    best_substeps = kSubsteps[s];
                }
            }
        // The cheapest arrangement within one percent of the best accuracy the
        // scene can reach. That is what the engine's defaults should be chosen
        // from: accuracy nobody can measure is just cost.
        const double tolerance = best + std::max(1.0e-6, std::abs(best) * 0.01);
        double cheapest_ms = 1e30;
        int cheapest_iterations = 0, cheapest_substeps = 0;
        for (int s = 0; s < kSubstepCount; ++s)
            for (int i = 0; i < kIterationCount; ++i)
                if (table[s][i].primary <= tolerance && table[s][i].milliseconds < cheapest_ms) {
                    cheapest_ms = table[s][i].milliseconds;
                    cheapest_iterations = kIterations[i];
                    cheapest_substeps = kSubsteps[s];
                }

        std::printf("\nscenario=%s  primary=%s (lower is better)  secondary=%s\n", scenario.name,
                    scenario.primary, scenario.secondary);
        std::printf("  %6s", "sub\\it");
        for (int i = 0; i < kIterationCount; ++i)
            std::printf(" %22d", kIterations[i]);
        std::printf("\n");
        std::printf("  %6s", "");
        for (int i = 0; i < kIterationCount; ++i)
            std::printf(" %10s %11s", "primary", "ms");
        std::printf("\n");
        for (int s = 0; s < kSubstepCount; ++s) {
            std::printf("  %6d", kSubsteps[s]);
            for (int i = 0; i < kIterationCount; ++i)
                std::printf(" %10.4f %11.2f", table[s][i].primary, table[s][i].milliseconds);
            std::printf("\n");
        }
        std::printf("  %6s", "second");
        for (int i = 0; i < kIterationCount; ++i)
            std::printf(" %22.4f", table[0][i].secondary);
        std::printf("   (sub=1)\n");
        std::printf("  best primary=%.4f at iter=%d sub=%d | within 1%%, cheapest = "
                    "iter=%d sub=%d at %.2f ms\n",
                    best, best_iterations, best_substeps, cheapest_iterations, cheapest_substeps,
                    cheapest_ms);
        // The arrangements that are not (iterations, sub-steps). Comparing each
        // one's *best* budget against the rigid baseline's best budget is not a
        // comparison of arrangements -- it is a comparison of budgets, and on a
        // scene that settles it hands the win to whichever arrangement happened
        // to be measured at the luckier budget. So each arrangement is measured
        // twice and both readings are printed:
        //
        //   at the baseline's budget   same scene, same frames, same work per
        //                              step; only the arrangement differs, which
        //                              is the only honest way to price a lever
        //   best over all budgets      what the arrangement can reach if it is
        //                              also given more iterations or sub-steps
        //
        // The second line is what says whether the arrangement can reach an
        // accuracy the rigid solve cannot reach at any budget. Neither line is
        // "how much error does it hide": every configuration runs the same
        // frames over the same scene, and nothing is relaxed to make a number
        // look better -- the accuracy is read off the bodies' own transforms.
        std::printf("  arrangement sweep, all measured at iter=%d sub=%d:\n", best_iterations,
                    best_substeps);
        const float stiffnesses[] = {0.0f, 200.0f, 1000.0f};
        const float relaxations[] = {0.2f, 0.5f, 1.0f};
        double baseline_primary = 1e30, baseline_ms = 1e30;
        for (float relaxation : relaxations)
            for (float stiffness : stiffnesses)
                for (int interleave = 0; interleave < 2; ++interleave) {
                    g_stiffness = stiffness;
                    g_interleave = interleave != 0;
                    g_relaxation = relaxation;
                    const Measurement same =
                        scenario.run(best_iterations, best_substeps, frames);
                    if (baseline_primary > 1e29) {
                        baseline_primary = same.primary;
                        baseline_ms = same.milliseconds;
                    }
                    double reach_ms = 1e30, reach_primary = 1e30;
                    int reach_iterations = 0, reach_substeps = 0;
                    for (int s = 0; s < kSubstepCount; ++s)
                        for (int i = 0; i < kIterationCount; ++i) {
                            const Measurement m = scenario.run(kIterations[i], kSubsteps[s], frames);
                            if (m.primary < reach_primary) {
                                reach_primary = m.primary;
                                reach_ms = m.milliseconds;
                                reach_iterations = kIterations[i];
                                reach_substeps = kSubsteps[s];
                            }
                        }
                    // A baseline of zero has no ratio to report: dividing by a
                    // clamped epsilon turns a change from nothing into a
                    // billion-percent figure. Report the absolute change in the
                    // metric instead when there is nothing to divide by.
                    char delta[32];
                    if (baseline_primary > 1.0e-9)
                        std::snprintf(delta, sizeof(delta), "%+8.1f%%",
                                      100.0 * (same.primary - baseline_primary) /
                                          baseline_primary);
                    else
                        std::snprintf(delta, sizeof(delta), "%+8.4f abs",
                                      same.primary - baseline_primary);
                    std::printf("    relax=%.1f stiff=%-5.0f inter=%d | at that budget %10.4f "
                                "(%s) %8.2f ms (%+6.1f%%) | best %10.4f at iter=%2d sub=%d "
                                "for %8.2f ms\n",
                                double(relaxation), double(stiffness), interleave, same.primary,
                                delta, same.milliseconds,
                                100.0 * (same.milliseconds - baseline_ms) / baseline_ms,
                                reach_primary, reach_iterations, reach_substeps, reach_ms);
                }
        g_stiffness = 0;
        g_interleave = false;
        g_relaxation = 0.2f;
    }
    return 0;
} catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s\n", e.what());
    return 1;
}
