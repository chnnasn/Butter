#include <butter/physics2d/butter2d.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace butter::physics2d;
static void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
static void falling(int count, bool boxes, bool stack, int layers = 5, bool must_sleep = true) {
    World::Config cfg;
    cfg.solver_iterations = 8;
    World w(cfg);
    w.create_body()
        .static_body()
        .at(0, -0.05f)
        .box(2000, 0.05f)
        .friction(.6f)
        .restitution(0)
        .build();
    std::vector<Body *> bodies;
    for (int i = 0; i < count; ++i) {
        int columns = stack ? count / layers : count;
        int column = i % columns, row = i / columns;
        auto builder = w.create_body()
                           .at(float(column) * 1.1f - float(columns) * .55f, 1.0f + row * 1.05f)
                           .friction(.6f)
                           .restitution(0);
        if (boxes)
            builder.box(.5f, .5f);
        else
            builder.circle(.5f);
        bodies.push_back(&builder.build());
    }
    auto &clock = w.create_body().kinematic().at(-2200, 0).velocity(3, 0).build();
    int limits = 0;
    std::size_t warm = 0;
    float min_bottom = 100;
    int frames = stack ? 600 : 360;
    for (int f = 0; f < frames; ++f) {
        w.step();
        limits += w.ccd_statistics().limited;
        warm += w.step_statistics().warm_started_points;
        for (auto *b : bodies) {
            auto bounds = compute_aabb(b->shape, b->transform);
            min_bottom = std::min(min_bottom, bounds.min.y);
            if (bounds.min.y < -.025f || !std::isfinite(b->transform.position.y)) {
                std::cerr << "first error count=" << count << " frame=" << f
                          << " y=" << b->transform.position.y << " bottom=" << bounds.min.y
                          << " vy=" << b->velocity.y << " angle=" << b->transform.angle << '\n';
                throw std::runtime_error("body penetrated floor");
            }
        }
    }
    int asleep = 0;
    float speed = 0;
    for (auto *b : bodies) {
        asleep += b->sleeping;
        speed = std::max(speed, b->velocity.length());
    }
    std::cout << (boxes ? "boxes" : "circles") << " count=" << count << " stack=" << stack
              << " sleeping=" << asleep << " limits=" << limits << " bottom=" << min_bottom
              << " speed=" << speed << " warm=" << warm << std::endl;
    require(std::abs(w.simulation_time() - frames / 60.0) < 1e-4, "world clock lost time");
    require(std::abs((clock.transform.position.x + 2200) - frames * .05f) < .08f,
            "unrelated motion lost time");
    if (must_sleep)
        require(asleep == count, "landings failed to sleep");
    if (!stack)
        require(limits == 0, "simple landings required CCD fallback");
}
static void low_stack() {
    World::Config cfg;
    cfg.solver_iterations = 8;
    World w(cfg);
    w.create_body().static_body().at(0, -.05f).box(10, .05f).friction(.6f).build();
    std::vector<Body *> boxes;
    for (int i = 0; i < 5; ++i)
        boxes.push_back(&w.create_body()
                             .at(0, .5f + i * 1.01f)
                             .box(.5f, .5f)
                             .friction(.6f)
                             .restitution(0)
                             .build());
    for (int i = 0; i < 900; ++i)
        w.step();
    std::vector<Transform> settled;
    int asleep = 0;
    for (auto b : boxes) {
        settled.push_back(b->transform);
        asleep += b->sleeping;
    }
    std::cout << "low stack sleeping=" << asleep << " base=" << boxes[0]->transform.position.y
              << std::endl;
    require(asleep == 5, "low stack failed to sleep");
    for (int i = 0; i < 1800; ++i)
        w.step();
    for (int i = 0; i < 5; ++i) {
        require((boxes[i]->transform.position - settled[i].position).length() < .001f,
                "resting stack drifted");
        require(std::abs(boxes[i]->transform.angle - settled[i].angle) < .001f,
                "resting stack jittered");
        require(boxes[i]->velocity.length() < .001f && std::abs(boxes[i]->angular_velocity) < .001f,
                "sleep residual velocity");
        require(boxes[i]->transform.position.y > .48f + i * .98f, "low stack sank or collapsed");
    }
    boxes.back()->wake();
    boxes.back()->velocity.x = 2;
    w.step();
    for (auto b : boxes)
        require(!b->sleeping, "wake did not propagate through island");
}
int main(int argc, char **) try {
    if (argc > 1) {
        falling(500, true, true, 25, false);
        falling(1000, true, true, 50, false);
        return 0;
    }
    for (int n : {100, 500, 1000})
        falling(n, false, false);
    for (int n : {500, 1000})
        falling(n, true, false);
    low_stack();
    for (int n : {500, 1000})
        falling(n, true, true);
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
