// Revolute (hinge) joint behaviour.
//
// The hinge used to be a zero-length DistanceJoint, and DistanceJoint::solve
// bails out as soon as the two anchors are closer than 1e-6: the "joint" was
// therefore a no-op the moment it was created, and a chain simply fell apart.
// These checks pin down the properties a real revolute joint must have:
// the anchors stay coincident along both axes, relative rotation stays free,
// the motor reaches its target speed, and the limits hold.
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

// Distance between the two anchor points in world space. A working hinge keeps
// this at zero in both the radial and the tangential direction.
static float anchor_gap(const HingeJoint &joint) {
    const Vec2 a =
        joint.a->transform.position + rotate(joint.anchor_a, Rot(joint.a->transform.angle));
    const Vec2 b =
        joint.b->transform.position + rotate(joint.anchor_b, Rot(joint.b->transform.angle));
    return (b - a).length();
}

static World::Config config(Vec2 gravity) {
    World::Config c;
    c.gravity = gravity;
    c.ccd.enabled = false;
    return c;
}

// A pendulum pinned at one end must swing, and the pin must not drift.
static void pendulum() {
    World world(config({0, -9.81f}));
    auto &pin = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &rod = world.create_body().at(1, 0).box(0.5f, 0.1f).build();
    rod.linear_damping = 0;
    rod.angular_damping = 0;
    auto &joint = world.add_hinge_joint_at(pin, rod, {0, 0});

    float max_gap = 0, min_angle = 0, max_angle = -10, min_radius = 10, max_radius = 0;
    // A frictionless pendulum released from horizontal swings down to -pi and
    // back up to horizontal. Twelve seconds covers several periods.
    float first_peak = 0, last_peak = 0;
    for (int step = 0; step < 720; ++step) {
        world.step(kDt);
        const float rate = std::abs(rod.angular_velocity);
        if (step < 120)
            first_peak = std::max(first_peak, rate);
        if (step >= 600)
            last_peak = std::max(last_peak, rate);
        max_gap = std::max(max_gap, anchor_gap(joint));
        min_angle = std::min(min_angle, joint.angle());
        max_angle = std::max(max_angle, joint.angle());
        const float radius = rod.transform.position.length();
        min_radius = std::min(min_radius, radius);
        max_radius = std::max(max_radius, radius);
    }
    check(max_gap < 0.01f, "hinge anchor drifted");
    // The rod hangs from a pin one unit away, so its centre never leaves that circle.
    check(min_radius > 0.995f && max_radius < 1.005f, "hinge did not hold the rod at its radius");
    // Rotation must stay completely free: a horizontal rod swings down through
    // vertical and up to horizontal on the far side, i.e. all the way to -pi.
    check(min_angle < -2.9f, "hinge blocked rotation");
    // ... while never swinging past its release height: that would mean the
    // solver created energy. It may only lose energy.
    check(min_angle > -3.3f && max_angle < 0.05f, "hinge let the rod swing past its release height");
    check(max_angle > -0.2f, "pendulum lost its swing after one period");
    check(last_peak <= first_peak * 1.01f, "hinge injected energy into the swing");

    // Damping has to be strong: angular damping scales the body's own spin, and
    // the pendulum's effective inertia about the pin is I_c + m*r^2, which here
    // is 12.5x I_c. A weak coefficient barely touches the swing.
    rod.angular_damping = 30.0f;
    for (int step = 0; step < 2400; ++step)
        world.step(kDt);
    check(std::abs(rod.transform.angle - (-1.5708f)) < 0.15f, "pendulum did not settle hanging");
    check(std::abs(rod.transform.position.x) < 0.02f, "settled pendulum is not vertical");
    check(std::abs(rod.transform.position.length() - 1.0f) < 0.02f, "settled pendulum left its radius");
    check(joint.angle() == rod.transform.angle - pin.transform.angle, "hinge angle convenience");
}

// A chain of hinges has to stay straight and connected; the old no-op joint
// let every link fall independently.
static void chain() {
    World world(config({0, -9.81f}));
    auto &anchor = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    std::vector<Body *> links{&anchor};
    const int link_count = 8;
    for (int i = 1; i <= link_count; ++i)
        links.push_back(
            &world.create_body().at(0, -0.5f * float(i)).box(0.25f, 0.05f).build());

    std::vector<HingeJoint *> joints;
    for (int i = 1; i <= link_count; ++i) {
        const float pin_y = -0.5f * float(i) + 0.25f;
        joints.push_back(&world.add_hinge_joint_at(*links[std::size_t(i - 1)], *links[std::size_t(i)],
                                                   {0, pin_y}));
    }

    float max_gap = 0, max_drift = 0;
    for (int step = 0; step < 600; ++step) {
        world.step(kDt);
        for (const auto *joint : joints)
            max_gap = std::max(max_gap, anchor_gap(*joint));
        // A straight hanging chain never picks up sideways motion.
        for (int i = 1; i <= link_count; ++i)
            max_drift = std::max(max_drift, std::abs(links[std::size_t(i)]->transform.position.x));
    }
    check(max_gap < 0.02f, "chain hinge stretched");
    check(max_drift < 0.05f, "chain sagged sideways");
    // Every link still hangs below the one above it.
    for (int i = 1; i <= link_count; ++i)
        check(links[std::size_t(i)]->transform.position.y < links[std::size_t(i - 1)]->transform.position.y,
              "chain link order broke");
}

// Warm starting and the 2x2 block have to make a heavily loaded hinge rigid.
// A radial-only joint lets a load like this slide along the tangent.
static void loaded_hinge() {
    World world(config({0, -9.81f}));
    auto &anchor = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &heavy = world.create_body().at(2, 0).box(1.0f, 0.2f).mass(40).build();
    auto &joint = world.add_hinge_joint_at(anchor, heavy, {0, 0});
    float max_gap = 0, min_angle = 0, min_arm = 10, max_arm = 0;
    for (int step = 0; step < 600; ++step) {
        world.step(kDt);
        max_gap = std::max(max_gap, anchor_gap(joint));
        min_angle = std::min(min_angle, joint.angle());
        const float arm = heavy.transform.position.length();
        min_arm = std::min(min_arm, arm);
        max_arm = std::max(max_arm, arm);
    }
    check(max_gap < 0.02f, "loaded hinge stretched");
    check(min_arm > 1.98f && max_arm < 2.02f, "loaded hinge arm changed");
    check(min_angle < -1.45f, "loaded hinge never swung past vertical");
}

// The motor drives the relative angle at the requested speed.
static void motor() {
    World world(config({0, 0}));
    auto &pin = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &rod = world.create_body().at(1, 0).box(0.5f, 0.1f).build();
    rod.linear_damping = 0;
    rod.angular_damping = 0;
    auto &joint = world.add_hinge_joint_at(pin, rod, {0, 0});
    joint.motor(2.0f, 500.0f);
    check(joint.max_motor_torque == 500.0f, "motor setter");
    for (int step = 0; step < 120; ++step)
        world.step(kDt);
    // Two seconds at 2 rad/s.
    check(std::abs(joint.angle() - 4.0f) < 0.2f, "motor did not reach its speed");
    check(std::abs(rod.angular_velocity - 2.0f) < 0.05f, "motor settled at the wrong rate");
    joint.motor(-3.0f, 500.0f);
    for (int step = 0; step < 60; ++step)
        world.step(kDt);
    check(std::abs(rod.angular_velocity + 3.0f) < 0.05f, "motor did not reverse");
    // A torque cap of zero is the documented way to switch the motor off.
    joint.motor(0, 0);
    const float coasting = rod.angular_velocity;
    world.step(kDt);
    // A stale motor impulse would spend itself as an immediate kick.
    check(std::abs(rod.angular_velocity - coasting) < 0.02f, "disabled motor kicked the body");
    for (int step = 0; step < 29; ++step)
        world.step(kDt);
    // Afterwards it merely coasts: only the solver's own drift may remain.
    check(std::abs(rod.angular_velocity - coasting) < std::abs(coasting) * 0.08f,
          "disabled motor still pulled");
    check(rod.angular_velocity * coasting > 0, "disabled motor reversed the spin");
}

// Limits hold the relative angle inside the requested window.
static void limits() {
    World world(config({0, 0}));
    auto &pin = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &rod = world.create_body().at(1, 0).box(0.5f, 0.1f).build();
    rod.linear_damping = 0;
    rod.angular_damping = 0;
    auto &joint = world.add_hinge_joint_at(pin, rod, {0, 0});
    joint.motor(5.0f, 1000.0f).limit(-0.4f, 0.4f);
    check(joint.enable_limit && joint.lower_angle == -0.4f && joint.upper_angle == 0.4f,
          "limit setter");

    float max_angle = -10, min_angle = 10;
    for (int step = 0; step < 240; ++step) {
        world.step(kDt);
        max_angle = std::max(max_angle, joint.angle());
        min_angle = std::min(min_angle, joint.angle());
    }
    check(max_angle <= 0.42f, "upper limit leaked");
    check(max_angle >= 0.38f, "upper limit never reached");
    check(min_angle >= -0.02f, "rod moved backwards off the upper limit");

    joint.motor(-5.0f, 1000.0f);
    for (int step = 0; step < 240; ++step) {
        world.step(kDt);
        min_angle = std::min(min_angle, joint.angle());
        max_angle = std::max(max_angle, joint.angle());
    }
    check(min_angle >= -0.42f, "lower limit leaked");
    check(min_angle <= -0.38f, "lower limit never reached");

    // A symmetric limit about an offset reference pose still reads relative angles.
    joint.no_limit();
    check(!joint.enable_limit, "no_limit");
    for (int step = 0; step < 240; ++step)
        world.step(kDt);
    check(joint.angle() < -1.0f, "limit stayed active after no_limit");
}

// A locked hinge (lower == upper) behaves like a weld.
static void locked_hinge() {
    World world(config({0, -9.81f}));
    auto &pin = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &rod = world.create_body().at(1, 0).box(0.5f, 0.1f).build();
    auto &joint = world.add_hinge_joint_at(pin, rod, {0, 0});
    joint.limit(0.0f, 0.0f);
    float max_gap = 0, max_tilt = 0;
    for (int step = 0; step < 300; ++step) {
        world.step(kDt);
        max_gap = std::max(max_gap, anchor_gap(joint));
        max_tilt = std::max(max_tilt, std::abs(joint.angle()));
    }
    check(max_gap < 0.02f, "locked hinge anchor drifted");
    check(max_tilt < 0.05f, "locked hinge rotated");
    check(std::abs(rod.transform.position.x - 1.0f) < 0.03f, "locked hinge translated");
}

// Distance and spring joints must keep working exactly as before.
static void distance_regression() {
    World world(config({0, -9.81f}));
    auto &anchor = world.create_body().static_body().at(0, 0).box(0.1f, 0.1f).build();
    auto &bob = world.create_body().at(0, -2).box(0.1f, 0.1f).build();
    bob.angular_damping = 0;
    world.add_distance_joint(anchor, bob, 2.0f);
    float max_error = 0;
    for (int step = 0; step < 120; ++step) {
        world.step(kDt);
        max_error = std::max(max_error, std::abs((bob.transform.position - anchor.transform.position).length() - 2.0f));
    }
    check(max_error < 0.02f, "distance joint no longer holds its length");

    auto &spring = world.create_body().at(2, 0).box(0.1f, 0.1f).build();
    auto &link = world.add_spring_joint(anchor, spring, 1.0f, 0.9f);
    check(link.spring_stiffness == 0.9f, "spring joint lost its stiffness");
}

int main() try {
    pendulum();
    chain();
    loaded_hinge();
    motor();
    limits();
    locked_hinge();
    distance_regression();
    std::cout << checks << " hinge checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
