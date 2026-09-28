// The solver arrangements that are not (iterations, sub-steps).
//
// Three arrangements are not a budget and are therefore not covered by the
// sub-step sweep in test_2d_substeps or by the benchmark's table: how much of a
// contact's measured overlap one position sweep removes, whether a compliant
// contact relaxes its own impulse, and when an island's joints are solved
// relative to its contacts.
//
// A benchmark can show a lever moves a number; only a test can say what the
// lever *is*. This file pins the mechanisms:
//
//   * A fully relaxed position pass is a sound stack, a stack that is left to
//     settle actually settles, and the default relaxation is what makes a
//     single-sub-step stack collapse -- which is the whole reason the knob
//     exists, stated as a measurement rather than as an opinion.
//   * Relaxation is monotone: less overlap left per sweep is never more sink.
//   * A fully relaxed correction still respects the world -- a dropped box
//     lands on the floor and a body moving at 50 m/s does not pass through a
//     wall. Cranking the correction up may not become a way to tunnel.
//   * Interleaving a joint into the contact sweep holds a chain at least as
//     well as the separate sweep, at the same budget.
//
// Every number is read out of the bodies' own transforms, never out of a
// counter a solver is free to define away, and every comparison is at a fixed
// budget: nothing here is relaxed, shortened or slept early to make a number
// look better.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
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
// A `rows`-box column on a floor, plus the worst penetration of any box through
// the floor's top surface. The world is owned here because `boxes` points into
// it: a helper that built the column in a local world and returned the body
// pointers would hand back dangling ones.
struct Column {
    std::unique_ptr<World> world;
    std::vector<Body *> boxes;
    float worst_floor_penetration{};
};
void settle_column(World &world, int rows, float spacing, Column &out) {
    world.create_body().static_body().at(0, -0.05f).box(40.0f, 0.05f).friction(0.6f).build();
    for (int i = 0; i < rows; ++i)
        out.boxes.push_back(&world.create_body()
                                .dynamic()
                                .at(0, 0.5f + spacing * float(i))
                                .box(0.5f, 0.5f)
                                .friction(0.6f)
                                .restitution(0)
                                .build());
    for (int frame = 0; frame < 600; ++frame) {
        for (Body *b : out.boxes)
            b->wake();
        world.step();
    }
    for (Body *b : out.boxes)
        out.worst_floor_penetration =
            std::max(out.worst_floor_penetration, 0.5f - b->transform.position.y);
}
// Compression of the column under its own weight, in millimetres: how far the
// top box ended up from where a rigid in-contact column would hold it. Measured
// against the contact height, so the 0.01 gap per level the boxes settle through
// is not charged to the solver -- and two-sided, so a column that never finished
// settling and is still standing on its start gaps reports that failure instead
// of reporting zero error.
float compression_mm(const Column &c) {
    const float ideal = 0.5f + float(c.boxes.size() - 1) * 1.0f;
    return std::abs(ideal - c.boxes.back()->transform.position.y) * 1000.0f;
}
Column column_at(int rows, float spacing, int iterations, int substeps, float relaxation) {
    World::Config cfg;
    cfg.solver_iterations = iterations;
    cfg.substeps = substeps;
    cfg.position_relaxation = relaxation;
    Column out{};
    out.world = std::make_unique<World>(cfg);
    settle_column(*out.world, rows, spacing, out);
    return out;
}
} // namespace
int main() try {
    {
        // The mechanism: the position sweep runs `solver_iterations` times a
        // step and removes this fraction of the overlap it measures, so a
        // resting stack settles where that removal and the sink gravity adds
        // each step balance. Less left per sweep is therefore never more sink --
        // measured as a chain, not as one comparison, because a single pair of
        // budgets can only ever show one point of it.
        const float relaxations[] = {0.2f, 0.4f, 0.6f, 0.8f, 1.0f};
        float previous = 1.0e30f;
        for (float relaxation : relaxations) {
            const Column c = column_at(20, 1.01f, 8, 4, relaxation);
            const float compression = compression_mm(c);
            std::cout << "relaxation " << relaxation << ": compression " << compression
                      << " mm (floor " << c.worst_floor_penetration * 1000.0f << " mm)\n";
            check(c.worst_floor_penetration < 0.02f, "a relaxed column sank through the floor");
            check(std::isfinite(compression), "a relaxed column diverged");
            check(compression <= previous,
                  "relaxing the position sweep further put the column lower, not higher");
            previous = compression;
        }
        check(previous < 100.0f, "a fully relaxed position sweep collapsed the column");
    }
    {
        // What the lever is worth. The default 0.2 relaxation leaves four fifths
        // of every measured overlap in place, which a single-sub-step stack
        // cannot absorb: twenty contacts each keeping a slice of the step's sink
        // is a column that walks itself apart. Fully relaxed, the same stack at
        // the same budget holds its contact height, and holds it as well as the
        // default arrangement does with four times the sub-steps.
        const Column default_single = column_at(20, 1.01f, 8, 1, 0.2f);
        const Column full_single = column_at(20, 1.01f, 8, 1, 1.0f);
        const Column default_four = column_at(20, 1.01f, 16, 4, 0.2f);
        const Column full_four = column_at(20, 1.01f, 16, 4, 1.0f);
        std::cout << "single sub-step: default " << compression_mm(default_single)
                  << " mm, fully relaxed " << compression_mm(full_single) << " mm\n";
        std::cout << "four sub-steps:  default " << compression_mm(default_four)
                  << " mm, fully relaxed " << compression_mm(full_four) << " mm\n";
        check(full_single.worst_floor_penetration < 0.02f,
              "a fully relaxed single sub-step let the column sink through the floor");
        check(compression_mm(full_single) < 100.0f,
              "a fully relaxed single sub-step collapsed the column");
        // The same budget, one knob apart: the default relaxation is the reason
        // a single sub-step cannot hold a stack at all.
        check(compression_mm(full_single) * 20.0f < compression_mm(default_single),
              "the default position relaxation is not what breaks a single-sub-step stack");
        // The same accuracy, a quarter of the sub-steps and half the iterations.
        check(compression_mm(full_single) < compression_mm(default_four),
              "a fully relaxed single sub-step did not reach what the default needs four for");
        // And at an adequate budget more relaxation is still less sink.
        check(compression_mm(full_four) < compression_mm(default_four),
              "full relaxation regressed at four sub-steps");
    }
    {
        // A bigger correction may not become a way to over-correct a body into
        // something else or to tunnel through it. Both scenes run with the
        // correction at its maximum and the default CCD guard in place.
        World::Config cfg;
        cfg.position_relaxation = 1.0f;
        World world(cfg);
        world.create_body().static_body().at(0, -0.05f).box(40.0f, 0.05f).friction(0.6f).build();
        world.create_body().static_body().at(4, 0.5f).box(0.5f, 0.5f).friction(0.6f).build();
        auto &faller = world.create_body()
                           .dynamic()
                           .at(0, 3.0f)
                           .box(0.5f, 0.5f)
                           .friction(0.6f)
                           .restitution(0)
                           .build();
        auto &slider = world.create_body()
                           .dynamic()
                           .at(-20.0f, 0.5f)
                           .box(0.5f, 0.5f)
                           .friction(0.2f)
                           .restitution(0)
                           .velocity(50.0f, 0)
                           .build();
        for (int frame = 0; frame < 600; ++frame) {
            faller.wake();
            slider.wake();
            world.step();
        }
        std::cout << "over-correction: faller y " << faller.transform.position.y << ", slider x "
                  << slider.transform.position.x << '\n';
        check(std::abs(faller.transform.position.y - 0.5f) < 0.01f,
              "full relaxation did not let a dropped box settle on the floor");
        // The wall's near face is at x = 3.5 and the slider's half width is 0.5,
        // so a centre past 4.0 has been pushed into the wall.
        check(slider.transform.position.x < 4.05f,
              "a fully relaxed correction pushed a body into a static obstacle");
        check(slider.transform.position.x > 2.0f,
              "a fully relaxed correction threw a body far off a static obstacle");
    }
    {
        // Interleaving changes *when* a joint sees a contact correction, so the
        // joint has to hold at least as well as it does in the separate sweep.
        // The links are placed face to face, so this is an island holding both
        // joints and contacts, and the gap is read off the joint's own anchors:
        // anything else would be measuring the scene rather than the solve.
        auto chain = [](bool interleave) {
            constexpr float link_length = 0.5f, top = 6.0f;
            World::Config cfg;
            cfg.gravity = {0, -9.81f};
            cfg.interleave_joints = interleave;
            World w(cfg);
            auto &anchor = w.create_body().static_body().at(0, top).box(0.1f, 0.1f).build();
            std::vector<Body *> links{&anchor};
            std::vector<HingeJoint *> joints;
            for (int i = 1; i <= 24; ++i) {
                auto &link = w.create_body()
                                 .dynamic()
                                 .at(0, top - link_length * float(i))
                                 .box(0.04f, link_length * 0.5f)
                                 .friction(0.4f)
                                 .build();
                links.push_back(&link);
                // The boundary between the two links, in world space, so both
                // anchors land on the same point and the joint starts satisfied.
                joints.push_back(&w.add_hinge_joint_at(*links[std::size_t(i - 1)], link,
                                                       {0, top - link_length * float(i) +
                                                               link_length * 0.5f}));
            }
            // The gap is sampled every frame, not just at the end. Once the
            // swing has died down the chain settles to a gap at the float
            // resolution of a six-metre world -- about 7e-7 -- so a settled-only
            // reading cannot tell two arrangements apart. The worst gap during
            // the swing is where a joint that is being corrected a sweep late
            // shows up.
            float worst = 0, settled = 0;
            for (int frame = 0; frame < 600; ++frame) {
                for (Body *link : links)
                    link->wake();
                w.step();
                settled = 0;
                for (const HingeJoint *joint : joints) {
                    const Vec2 a = joint->a->transform.position +
                                   rotate(joint->anchor_a, Rot(joint->a->transform.angle));
                    const Vec2 b = joint->b->transform.position +
                                   rotate(joint->anchor_b, Rot(joint->b->transform.angle));
                    settled = std::max(settled, (b - a).length());
                }
                worst = std::max(worst, settled);
            }
            return std::make_pair(worst, settled);
        };
        const auto separate = chain(false);
        const auto interleaved = chain(true);
        std::cout << "chain: separate gap " << separate.first << " m worst, " << separate.second
                  << " m settled\n";
        std::cout << "chain: interleaved gap " << interleaved.first << " m worst, "
                  << interleaved.second << " m settled\n";
        check(separate.first < 0.05f, "the separate joint sweep failed to hold the chain");
        check(interleaved.first < 0.05f, "interleaving the joint sweep failed to hold the chain");
        check(interleaved.first <= separate.first * 1.05f,
              "interleaving the joint sweep stretched the chain");
    }
    std::cout << checks << " solver option checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
