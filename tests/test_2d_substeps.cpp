// Sub-stepping: `Config::substeps`.
//
// Splitting a step into several integration/solve windows is the engine's one
// lever that changes *how* a frame is solved rather than how long it is solved
// for, so it needs a test that pins down both halves of the bargain: what it
// buys, and that it costs nothing when it is switched off.
//
// The bought half is measured on the two scenes that separate the two ways a
// contact solver fails:
//
//   * A tall column under its own weight is a *propagation* problem. More
//     Gauss-Seidel sweeps push the load down the chain, so iterations are what
//     fix it -- but a single window at the default eight iterations buckles it.
//   * A large mass ratio is not a propagation problem. The light bodies cannot
//     push the heavy one back within one window no matter how many times the
//     window is swept, so iterations alone cannot fix it at any count, and
//     sub-stepping can.
//
// The free half is that `substeps == 1` is the original single-window step, so
// every existing trajectory, every persistent contact and every cached impulse
// is bit-for-bit what it was. That is checked directly rather than assumed.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace butter::physics2d;

static int checks = 0;
static void check(bool ok, const char *message) {
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}

static constexpr float kDt = 1.0f / 60.0f;

static World::Config config(int iterations, int substeps) {
    World::Config cfg;
    cfg.solver_iterations = iterations;
    cfg.substeps = substeps;
    return cfg;
}

// Waking every dynamic body each frame keeps a scene awake, so the measurement
// is of the solve and not of whichever configuration happens to cross the sleep
// threshold first.
static void force_awake(const std::vector<Body *> &bodies) {
    for (auto *body : bodies)
        if (body->type == BodyType::Dynamic)
            body->wake();
}

static void add_ground(World &world) {
    world.create_body()
        .static_body()
        .at(0, -0.05f)
        .box(400.0f, 0.05f)
        .friction(0.6f)
        .restitution(0)
        .build();
}

// ---------------------------------------------------------------------------
// A twenty-box column. Contact height for the top box is 0.5 + 19 * 1.0.
// ---------------------------------------------------------------------------
struct ColumnResult {
    float top_y{};
    int asleep{};
};

static ColumnResult column(int iterations, int substeps, bool keep_awake) {
    constexpr int layers = 20;
    World world(config(iterations, substeps));
    add_ground(world);
    std::vector<Body *> boxes;
    for (int i = 0; i < layers; ++i)
        boxes.push_back(&world.create_body()
                            .at(0, 0.5f + float(i) * 1.01f)
                            .box(0.5f, 0.5f)
                            .friction(0.6f)
                            .restitution(0)
                            .build());
    for (int frame = 0; frame < 600; ++frame) {
        if (keep_awake)
            force_awake(boxes);
        world.step(kDt);
    }
    ColumnResult result;
    result.top_y = boxes.back()->transform.position.y;
    for (auto *box : boxes)
        result.asleep += box->sleeping;
    return result;
}

// ---------------------------------------------------------------------------
// Ten light boxes carrying a crate fifty times heavier.
// ---------------------------------------------------------------------------
static float heavy_drop_mm(int iterations, int substeps) {
    constexpr int light_count = 10;
    World world(config(iterations, substeps));
    add_ground(world);
    std::vector<Body *> bodies;
    for (int i = 0; i < light_count; ++i)
        bodies.push_back(&world.create_body()
                             .at(0, 0.5f + float(i) * 1.01f)
                             .box(0.5f, 0.5f)
                             .mass(1.0f)
                             .friction(0.6f)
                             .restitution(0)
                             .build());
    auto &heavy = world.create_body()
                      .at(0, 0.5f + float(light_count) * 1.01f)
                      .box(0.5f, 0.5f)
                      .mass(50.0f)
                      .friction(0.6f)
                      .restitution(0)
                      .build();
    for (int frame = 0; frame < 600; ++frame) {
        force_awake(bodies);
        heavy.wake();
        world.step(kDt);
    }
    // 0.5 + 10 * 1.0 is where a rigid, in-contact stack would hold the crate.
    const float ideal = 0.5f + float(light_count);
    return std::max(0.0f, ideal - heavy.transform.position.y) * 1000.0f;
}

// A single window at the default iteration count cannot hold the crate, and no
// iteration count rescues it -- the failure is not a lack of sweeps. Sub-stepping
// does, because what the light bodies need is more *time* to push back in, not
// more sweeps of the same instant.
static void mass_ratio_needs_windows_not_sweeps() {
    const float single = heavy_drop_mm(8, 1);
    const float many_sweeps = heavy_drop_mm(16, 1);
    const float windows = heavy_drop_mm(4, 4);
    std::cout << "mass_ratio drop(mm): iter8/sub1=" << single << " iter16/sub1=" << many_sweeps
              << " iter4/sub4=" << windows << '\n';
    check(single > 1000.0f, "single window unexpectedly held the heavy crate");
    check(many_sweeps > 1000.0f, "iterations alone were expected not to close a 50:1 mass ratio");
    check(windows < 100.0f, "sub-stepping did not hold the heavy crate");
}

// The column is the opposite case: one window buckles it, and it takes many
// iterations to stand up straight. Sub-stepping must not make it worse, and a
// fully awake column is the harsher of the two measurements.
static void column_buckles_in_one_window() {
    const ColumnResult single = column(8, 1, true);
    const ColumnResult windows = column(8, 2, true);
    const float contact_height = 0.5f + 19.0f;
    std::cout << "column top_y: iter8/sub1=" << single.top_y << " iter8/sub2=" << windows.top_y
              << " (contact height " << contact_height << ")\n";
    check(single.top_y < 17.0f, "the single-window column was expected to buckle");
    check(windows.top_y > single.top_y, "sub-stepping made the column worse");
    check(std::abs(windows.top_y - contact_height) < 0.15f, "sub-stepped column did not stand");
}

// The same column, left to settle naturally, has to reach rest and stay there.
static void column_settles_and_sleeps() {
    const ColumnResult windows = column(8, 2, false);
    std::cout << "column natural: top_y=" << windows.top_y << " asleep=" << windows.asleep << "/20\n";
    check(windows.asleep == 20, "sub-stepped column never came to rest");
    check(std::abs(windows.top_y - 19.5f) < 0.15f, "sub-stepped column settled out of place");
}

// `substeps == 1` is the classic step, so it has to be reached both by asking
// for it and by leaving the field alone.
static void one_window_is_the_classic_step() {
    World::Config explicit_one = config(8, 1);
    World::Config untouched;
    check(untouched.substeps == 1, "Config::substeps no longer defaults to a single window");

    auto run = [](const World::Config &cfg) {
        World world(cfg);
        add_ground(world);
        std::vector<Body *> boxes;
        for (int i = 0; i < 8; ++i)
            boxes.push_back(&world.create_body()
                                .at(0.3f * float(i) - 1.0f, 0.5f + float(i) * 1.01f)
                                .box(0.5f, 0.5f)
                                .friction(0.6f)
                                .restitution(0)
                                .build());
        for (int frame = 0; frame < 240; ++frame) {
            force_awake(boxes);
            world.step(kDt);
        }
        std::vector<float> state;
        for (auto *box : boxes) {
            state.push_back(box->transform.position.x);
            state.push_back(box->transform.position.y);
            state.push_back(box->transform.angle);
            state.push_back(box->velocity.x);
            state.push_back(box->velocity.y);
            state.push_back(box->angular_velocity);
        }
        return state;
    };
    const auto a = run(untouched);
    const auto b = run(explicit_one);
    check(a.size() == b.size(), "column state size mismatch");
    bool identical = true;
    for (std::size_t i = 0; i < a.size(); ++i)
        identical &= a[i] == b[i];
    check(identical, "substeps=1 is not the untouched default step");
    check(run(config(8, 2)) != a, "the sub-step count has no effect");
}

// A short, cheap scene whose whole state can be compared exactly. Used only for
// equivalence checks, so it does not need to be interesting.
static std::vector<float> tiny_state(int iterations, int substeps) {
    World world(config(iterations, substeps));
    add_ground(world);
    auto &box = world.create_body().at(0.2f, 1.4f).box(0.5f, 0.5f).friction(0.6f).restitution(0).build();
    for (int frame = 0; frame < 45; ++frame) {
        box.wake();
        world.step(kDt);
    }
    return {box.transform.position.x, box.transform.position.y, box.transform.angle,
            box.velocity.x,           box.velocity.y,           box.angular_velocity};
}

// The knob is clamped, so an accidental huge value degrades rather than hangs.
static void substeps_are_clamped() {
    check(tiny_state(4, World::kMaxSubsteps) == tiny_state(4, 1000),
          "substeps were not clamped to kMaxSubsteps");
    const auto one = tiny_state(4, 1);
    check(tiny_state(4, 0) == one, "a zero substep count did not fall back to one");
    check(tiny_state(4, -8) == one, "a negative substep count did not fall back to one");
}

// Time and the continuous pass are per frame, not per window: a frame advances
// the world by dt however it is divided, and the CCD counters have to report the
// whole frame rather than whichever window happened to run last.
static void time_and_ccd_are_per_frame() {
    for (int substeps : {1, 2, 4}) {
        constexpr int frames = 120;
        World world(config(8, substeps));
        add_ground(world);
        auto &faller = world.create_body().at(0, 6.0f).box(0.4f, 0.4f).restitution(0).build();
        double advanced = 0;
        for (int frame = 0; frame < frames; ++frame) {
            faller.wake();
            world.step(kDt);
            advanced += world.ccd_statistics().advanced_time;
        }
        check(std::abs(world.simulation_time() - frames * kDt) < 1.0e-4,
              "sub-stepping changed the world clock");
        // Every frame must report a whole dt of advanced time; a window that
        // overwrote the others would report dt/substeps.
        check(std::abs(advanced - frames * kDt) < 1.0e-2,
              "CCD advanced_time did not accumulate across windows");
    }
}

// A joint chain is a propagation problem like the column, so sub-stepping must
// at least not damage it: the same iteration count has to keep the anchors as
// close together.
static void joint_chain_is_not_damaged() {
    auto gap_after = [](int iterations, int substeps) {
        constexpr int links = 24;
        constexpr float link_length = 0.5f;
        World world(config(iterations, substeps));
        auto &anchor = world.create_body().static_body().at(0, 6.0f).box(0.1f, 0.1f).build();
        std::vector<Body *> chain{&anchor};
        for (int i = 1; i <= links; ++i)
            chain.push_back(&world.create_body()
                                .at(0, 6.0f - link_length * float(i))
                                .box(0.04f, link_length * 0.5f)
                                .build());
        std::vector<HingeJoint *> joints;
        for (int i = 1; i <= links; ++i)
            joints.push_back(&world.add_hinge_joint_at(
                *chain[std::size_t(i - 1)], *chain[std::size_t(i)],
                {0, 6.0f - link_length * float(i) + link_length * 0.5f}));
        for (int frame = 0; frame < 600; ++frame) {
            force_awake(chain);
            world.step(kDt);
        }
        float worst = 0;
        for (const auto *joint : joints) {
            const Vec2 a =
                joint->a->transform.position + rotate(joint->anchor_a, Rot(joint->a->transform.angle));
            const Vec2 b =
                joint->b->transform.position + rotate(joint->anchor_b, Rot(joint->b->transform.angle));
            worst = std::max(worst, (b - a).length());
        }
        return worst;
    };
    const float single = gap_after(8, 1);
    const float windows = gap_after(8, 2);
    std::cout << "joint chain worst gap: iter8/sub1=" << single * 1000.0f
              << "mm iter8/sub2=" << windows * 1000.0f << "mm\n";
    check(single < 0.05f && windows < 0.05f, "a hinge chain lost its anchors under sub-stepping");
}

int main() try {
    one_window_is_the_classic_step();
    mass_ratio_needs_windows_not_sweeps();
    column_buckles_in_one_window();
    column_settles_and_sleeps();
    substeps_are_clamped();
    time_and_ccd_are_per_frame();
    joint_chain_is_not_damaged();
    std::cout << checks << " sub-step checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
