// Contacts as maintained objects, not rebuilt ones.
//
// A step detects contacts twice: once before integration and once after it. The
// second pass exists to refresh the manifolds, the anchors and the restitution
// targets of the pairs that *moved*. Rebuilding the entire contact set for that
// re-ran the narrow phase, the filter and the fixture pairing for every pair in
// the world, and moved every constraint in the array -- which is the order the
// solver visits them in.
//
// The scene where that matters is a large world that is mostly at rest with a
// handful of bodies awake. This file asserts two things about it, and they are
// different in kind:
//
//   * The result is *bit for bit* the same as rebuilding everything. The
//     in-place refresh is an implementation choice, so the test runs the same
//     scene both ways and compares every position and angle. Anything the
//     refresh cannot express in place -- a contact appearing or separating, a
//     manifold that vanishes, an event membership that flips -- makes it defer
//     to the rebuild by itself, which is why the two runs can agree exactly.
//
//   * The work follows the movers. The number of candidate pairs the contact
//     build examines is the awake bodies' neighbourhood, not the world, so a
//     hundred sleeping boxes cost nothing to leave alone. That is asserted on
//     the candidate-pair counter, which is deterministic, rather than on a
//     wall-clock ratio that a shared machine would make noisy.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace butter::physics2d;
static int checks = 0;
static void check(bool value, const char *message) {
    ++checks;
    if (!value)
        throw std::runtime_error(message);
}
namespace {
constexpr int kSettle = 480;   // frames allowed for the pile to fall asleep
constexpr int kFrames = 240;   // frames measured afterwards
struct Run {
    std::vector<Vec2> positions;
    std::vector<float> angles;
    std::size_t candidate_pairs{}, fixture_pairs{}, narrow_tests{}, refreshed{}, unchanged{},
        rebuilt{}, asleep{};
    double ms{};
};
// A floor with a 12 x 8 pile on the left and a small tower on the right. The
// pile settles and falls asleep and is never touched again; the tower is woken
// every frame, so its bodies are re-detected on every step while the world they
// belong to stands still. That is the case the refresh exists for: the pile's
// contacts are live, correct and untouched, and re-deriving them every frame is
// pure waste.
Run run(bool incremental) {
    World::Config c;
    c.incremental_contacts = incremental;
    World w(c);
    w.create_body().static_body().at(0, -0.5f).box(30, 0.5f).friction(0.6f).build();
    std::vector<Body *> pile;
    for (int row = 0; row < 12; ++row)
        for (int column = 0; column < 8; ++column)
            pile.push_back(&w.create_body()
                                .dynamic()
                                .at(-3.75f + 0.51f * float(column), 0.5f + 1.002f * float(row))
                                .box(0.25f, 0.5f)
                                .friction(0.6f)
                                .build());
    std::vector<Body *> tower;
    for (int row = 0; row < 4; ++row)
        tower.push_back(&w.create_body()
                            .dynamic()
                            .at(14.0f, 0.5f + 1.002f * float(row))
                            .box(0.5f, 0.5f)
                            .friction(0.6f)
                            .build());
    for (int frame = 0; frame < kSettle; ++frame)
        w.step();
    Run out{};
    const auto start = std::chrono::steady_clock::now();
    for (int frame = 0; frame < kFrames; ++frame) {
        // Waking a resting tower is enough: its boxes shift by fractions of a
        // micron every frame, which is a move the broad phase has to report but
        // not a pair that appears or disappears.
        for (Body *b : tower)
            b->wake();
        w.step();
        const auto &s = w.step_statistics();
        out.candidate_pairs += s.candidate_pairs;
        out.fixture_pairs += s.fixture_pairs;
        out.narrow_tests += s.narrow_tests;
        out.refreshed += s.refreshed_contact_steps;
        out.unchanged += s.unchanged_contact_steps;
        out.rebuilt += (s.unchanged_contact_steps == 0 && s.refreshed_contact_steps == 0) ? 1 : 0;
    }
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                 .count();
    for (Body *b : pile) {
        out.positions.push_back(b->transform.position);
        out.angles.push_back(b->transform.angle);
        out.asleep += b->sleeping ? 1 : 0;
    }
    return out;
}
} // namespace
int main() try {
    const Run incremental = run(true);
    const Run rebuilt = run(false);
    // The scene has to be the one the refresh is for: a sleeping pile with a
    // handful of waking bodies, and the incremental run has to have taken the
    // incremental path, or this file proves nothing.
    check(incremental.asleep == 96, "the pile did not settle");
    check(incremental.refreshed > kFrames / 2,
          "a sleeping pile with an awake tower never used the in-place contact refresh");
    check(rebuilt.refreshed == 0, "the refresh ran while it was switched off");
    // Both ways must agree exactly. A rebuilt contact set and a refreshed one
    // differ only in how the same contacts were arrived at, so any difference
    // at all -- however small -- is a bug in the refresh.
    check(incremental.positions == rebuilt.positions,
          "the in-place contact refresh changed the trajectory");
    check(incremental.angles == rebuilt.angles, "the in-place contact refresh changed a rotation");
    // Work follows the movers: four awake boxes must not make the contact build
    // walk the whole pile. The rebuilt run examines every pair every frame.
    std::cout << "refresh: candidate_pairs " << incremental.candidate_pairs << " fixture_pairs "
              << incremental.fixture_pairs << " narrow_tests " << incremental.narrow_tests
              << " refreshed_frames " << incremental.refreshed << " unchanged_frames "
              << incremental.unchanged << " rebuilt_frames " << incremental.rebuilt << " ms "
              << incremental.ms << '\n';
    std::cout << "rebuild: candidate_pairs " << rebuilt.candidate_pairs << " fixture_pairs "
              << rebuilt.fixture_pairs << " narrow_tests " << rebuilt.narrow_tests << " ms "
              << rebuilt.ms << '\n';
    check(rebuilt.candidate_pairs > 0, "the rebuilt run examined no candidate pair");
    check(incremental.candidate_pairs * 8 < rebuilt.candidate_pairs,
          "four awake boxes still made the contact build walk the whole world");
    std::cout << checks << " contact refresh checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
