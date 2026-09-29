// Position-correction obstacles at scale.
//
// The guard that stops a solver correction from pushing a body through static
// geometry used to walk a sorted list of every obstacle in the world. The
// acceptance scene for it was "a box pile and a floor", which is exactly the
// scene where a linear scan is already cheap. The scene that matters is a pile
// plus thousands of walls, platforms and terrain blocks.
//
// Two things are asserted here, and they are different in kind.
//
//   * The amount of guard work is *exactly* the same whether the world contains
//     no parked platform or eight thousand: the same candidates, the same
//     sweeps, the same fast rejections. That is a bit-for-bit claim and it is
//     the index's whole contract.
//   * The result is the same too, but only up to a stated tolerance, and the
//     reason is worth writing down. A twelve-box tower creeps sideways as it
//     settles and that creep is chaotic; once the world holds sixteen or more
//     static bodies the first step already differs in the last bits of one
//     angular velocity, and 300 frames later that has grown to a few
//     centimetres. The counters are identical throughout, CCD-off makes the
//     difference vanish entirely, and the number of parked platforms is
//     irrelevant once there are sixteen of them -- so this is a pre-existing
//     floating-point ordering wart of the solver, not something the index
//     scales with. Asserting bit equality on a chaotic tower would only be
//     asserting a coincidence, so the tower is checked against the physical
//     quantities it is supposed to preserve, and bit equality is asserted
//     where it is genuinely stable: between two obstacle counts.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using namespace butter::physics2d;
static int checks = 0;
static void check(bool value, const char *message) {
    ++checks;
    if (!value)
        throw std::runtime_error(message);
}
static World::Config config() {
    World::Config c;
    c.gravity = {};
    return c;
}
namespace {
constexpr int kFrames = 300;
constexpr int kTower = 12;
// Far enough that no correction envelope the tower can produce ever reaches a
// parked platform, with tens of metres of margin.
constexpr float kParked = 60.0f;
struct Scene {
    std::size_t candidates{}, sweeps{}, clamps{}, corrections{}, fast_rejections{};
    double total_ms{};
    // Where the frame went, so the claim below is checkable rather than
    // asserted: the guard's own phase, the continuous pass, the two per-body
    // compatibility walks, and the event scan.
    double ccd_ms{}, cache_ms{}, position_ms{}, detection_ms{};
    std::vector<Vec2> positions;
    std::vector<float> angles;
    int asleep{};
};
// A `kTower`-box column on a floor, plus `obstacles` static platforms parked
// `kParked` units to the side. `force_awake` keeps the guard running on every
// frame instead of letting the column fall asleep and skip the work entirely.
std::unique_ptr<World> column(int obstacles, bool force_awake, Scene &out) {
    auto world = std::make_unique<World>();
    auto &w = *world;
    w.create_body().static_body().at(0, -0.5f).box(40, 0.5f).friction(0.6f).build();
    std::vector<Body *> stack;
    for (int i = 0; i < kTower; ++i)
        stack.push_back(&w.create_body()
                             .dynamic()
                             .at(0, 0.5f + float(i) * 0.995f)
                             .box(0.5f, 0.5f)
                             .friction(0.6f)
                             .build());
    int side = 1;
    while (side * side < obstacles)
        ++side;
    for (int i = 0; i < obstacles; ++i)
        w.create_body()
            .static_body()
            .at(kParked + float(i % side) * 6.0f, -20.0f + float(i / side) * 6.0f)
            .box(2.0f, 0.2f)
            .build();
    const auto start = std::chrono::steady_clock::now();
    for (int frame = 0; frame < kFrames; ++frame) {
        if (force_awake)
            for (Body *b : stack)
                b->wake();
        w.step();
        const auto &s = w.step_statistics();
        out.candidates += s.projection_candidates;
        out.sweeps += s.projection_sweeps;
        out.clamps += s.position_clamps;
        out.corrections += s.position_corrections;
        out.fast_rejections += s.projection_fast_rejections;
        out.ccd_ms += s.ccd_ms;
        out.cache_ms += s.cache_ms;
        out.position_ms += s.position_ms;
        out.detection_ms += s.detection_ms;
    }
    out.total_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
            .count();
    for (Body *b : stack) {
        out.positions.push_back(b->transform.position);
        out.angles.push_back(b->transform.angle);
        out.asleep += b->sleeping ? 1 : 0;
    }
    return world;
}
// Every quantity that says *how much work* the position guard did. These are the
// index's contract, and they have to match exactly.
bool same_guard_work(const Scene &a, const Scene &b) {
    return a.candidates == b.candidates && a.sweeps == b.sweeps && a.clamps == b.clamps &&
           a.corrections == b.corrections && a.fast_rejections == b.fast_rejections;
}
} // namespace
int main() try {
    {
        // Settling tower, sleeping naturally. The guard work must be identical
        // for every obstacle count, and two large counts must agree with each
        // other bit for bit -- the count of parked platforms cannot matter.
        Scene base{};
        column(0, false, base);
        check(base.candidates > 0, "the settling tower never exercised the projection guard");
        check(base.corrections > 0, "the settling tower produced no position corrections");
        check(base.fast_rejections > 0, "the settling tower never used the cheap rejection");
        Scene two_hundred{};
        column(200, false, two_hundred);
        Scene eight_thousand{};
        column(8000, false, eight_thousand);
        check(same_guard_work(base, two_hundred), "200 parked platforms changed the guard's work");
        check(same_guard_work(base, eight_thousand),
              "8000 parked platforms changed the guard's work");
        check(two_hundred.positions == eight_thousand.positions &&
                  two_hundred.angles == eight_thousand.angles,
              "the number of parked platforms changed the tower");
        check(eight_thousand.clamps == base.clamps, "parked platforms changed the clamp count");
        check(eight_thousand.asleep == base.asleep, "parked platforms changed the sleep count");
        // Physical outcome: the tower does not topple, does not sink, and its
        // lateral creep stays inside a fifth of a box.
        for (int i = 0; i < kTower; ++i) {
            check(std::abs(eight_thousand.positions[i].y - base.positions[i].y) < 1.0e-3f,
                  "parked platforms changed a contact height in the tower");
            check(std::abs(eight_thousand.positions[i].x - base.positions[i].x) < 0.1f,
                  "parked platforms changed the tower's lateral creep");
            check(std::abs(eight_thousand.angles[i]) < 0.02f, "the tower toppled");
        }
        std::cout << "settling tower: candidates " << base.candidates << " sweeps " << base.sweeps
                  << " fast " << base.fast_rejections << " base_ms " << base.total_ms << '\n';
    }
    {
        // The same scene while the tower is forced to stay active, so the guard
        // really runs on every frame. A regression back to a linear scan would
        // show up as the guard work and the frame cost growing with the obstacle
        // count. The baseline for *that* comparison is a second large count
        // rather than the empty world: keeping the tower awake makes its chaotic
        // creep diverge from the empty-world run, which perturbs the guard
        // counts for reasons that have nothing to do with the index.
        Scene two_thousand{};
        column(2000, true, two_thousand);
        Scene eight_thousand{};
        column(8000, true, eight_thousand);
        check(two_thousand.candidates > 200,
              "the active tower never exercised the projection guard");
        check(same_guard_work(two_thousand, eight_thousand),
              "active tower: 8000 parked platforms did 4x the guard work of 2000");
        check(two_thousand.positions == eight_thousand.positions &&
                  two_thousand.angles == eight_thousand.angles,
              "active tower: the number of parked platforms changed the trajectory");
        Scene empty{};
        column(0, true, empty);
        auto phases = [](const char *label, const Scene &s) {
            std::cout << "  " << label << ": total " << s.total_ms << " ccd " << s.ccd_ms << " guard "
                      << s.position_ms << " compat_walks " << s.cache_ms << " events "
                      << s.detection_ms << '\n';
        };
        phases("empty       ", empty);
        phases("2000 parked ", two_thousand);
        phases("8000 parked ", eight_thousand);
        // What one parked obstacle may cost per frame. The guard's own work is
        // exactly scale-free -- the counters above are the whole proof of that --
        // but a world of eight thousand bodies still pays to *find out* that
        // none of them changed, and that is deliberate: transforms are public
        // fields, so a step compares them against the snapshot it filed instead
        // of trusting a version counter. What this budget rejects is a wall
        // costing a *query*: a bounding box against every mover, or a sweep, or
        // a projection candidate. Tens of nanoseconds buys a pose comparison and
        // nothing more, so a regression to a linear scan inside the guard lands
        // two orders of magnitude above it.
        //
        // "Tens of nanoseconds" is a claim about optimised code. An unoptimised
        // build does the same work -- the same comparisons, the same counters --
        // through real calls instead of inlined ones, and measures about seven
        // times as much: 560-600 ns per obstacle-frame against 80-110 ns here.
        // The budget is therefore scaled by configuration rather than fixed. It
        // still rejects exactly what it is for: the class it guards against sits
        // two orders of magnitude above *both* figures.
#ifdef NDEBUG
        constexpr double kMarginalNsPerObstacleFrame = 300.0;
#else
        constexpr double kMarginalNsPerObstacleFrame = 2000.0;
#endif
        const double marginal_ns =
            (eight_thousand.total_ms - empty.total_ms) * 1.0e6 / (8000.0 * kFrames);
        std::cout << "active tower:   candidates " << two_thousand.candidates << " sweeps "
                  << two_thousand.sweeps << " fast " << two_thousand.fast_rejections << " empty_ms "
                  << empty.total_ms << " worst_ms " << eight_thousand.total_ms
                  << " marginal_ns_per_obstacle_frame " << marginal_ns << '\n';
        check(marginal_ns < kMarginalNsPerObstacleFrame,
              "a parked obstacle costs more than a pose comparison per frame");
    }
    {
        // Obstacles placed *inside* the tower's reach must still be respected.
        // The index may not become so eager to prune that it loses a real
        // obstacle: without the wall the tower slides right, with it the tower
        // has to stop short of the wall's face.
        auto run = [&](bool walled) {
            World w;
            w.create_body().static_body().at(0, -0.5f).box(40, 0.5f).friction(0.6f).build();
            if (walled)
                w.create_body().static_body().at(1.2f, 20.0f).box(0.1f, 40.0f).friction(0).build();
            std::vector<Body *> stack;
            for (int i = 0; i < kTower; ++i)
                stack.push_back(&w.create_body()
                                    .dynamic()
                                    .at(0, 0.5f + float(i) * 0.98f)
                                    .box(0.5f, 0.5f)
                                    .friction(0)
                                    .velocity(1.0f, 0)
                                    .build());
            for (int frame = 0; frame < 240; ++frame) {
                for (Body *b : stack)
                    b->wake();
                w.step();
            }
            float furthest = -1.0e9f;
            for (Body *b : stack)
                furthest = std::max(furthest, b->transform.position.x);
            return furthest;
        };
        const float free_slide = run(false);
        const float stopped = run(true);
        check(free_slide > 2.0f, "the unblocked tower never travelled far enough to matter");
        // The wall's left face is at 1.2 - 0.1 = 1.1, and a box is a unit wide.
        check(stopped < 0.9f, "the tower crossed a wall that stood inside its reach");
        check(stopped < free_slide - 1.0f, "the wall inside the tower's reach had no effect");
    }
    {
        // A kinematic obstacle has to be filed at the pose the position solver
        // actually sees. Integration moves it *before* the solve, so if the index
        // were refreshed earlier, a correction would still be gated by bounds the
        // obstacle has already left.
        World w(config());
        auto &wall = w.create_body().kinematic().at(0, 0).box(0.01f, 1.0f).friction(0).build();
        auto &small = w.create_body().at(-0.04f, 0).box(0.02f, 0.02f).friction(0).build();
        small.fixed_rotation = true;
        auto &anchor = w.create_body().static_body().at(5, 0).circle(0.01f).build();
        w.add_distance_joint(small, anchor, 5.0f);
        // The joint pulls `small` towards the wall's face, so the projection
        // guard has to stop it there.
        w.step();
        check(w.step_statistics().position_clamps > 0,
              "a stationary kinematic obstacle did not gate a projection");
        const float gated = small.transform.position.x;
        check(gated <= -0.0299f, "the correction crossed the kinematic obstacle's face");
        // Now carry the wall away in a single step. Its old bounds must be gone
        // by the time the position solver runs, so the same correction is free.
        wall.velocity = {600, 0};
        small.transform.position = {-0.04f, 0};
        w.step();
        check(wall.transform.position.x > 5, "the kinematic obstacle did not move");
        check(w.step_statistics().position_clamps == 0,
              "expired kinematic obstacle bounds were reused by the position guard");
        check(small.transform.position.x > gated,
              "the correction was still held back by the obstacle's old position");
    }
    std::cout << checks << " obstacle checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
