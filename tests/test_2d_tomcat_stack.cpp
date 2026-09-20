#include <butter/physics2d/butter2d.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace butter::physics2d;
static void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
// Geometry, mass/inertia, damping and layout match TomCat dev_butter's adapter
// and .scratch/comparison-20260921/diagnostic.cpp; no renderer or adapter timing.
int main(int argc, char **argv) try {
    int count = argc > 1 ? std::stoi(argv[1]) : 100;
    World w;
    bool circles = argc > 2 && std::string(argv[2]) == "circles";
    std::vector<Body *> bodies;
    auto add = [&](float x, float y, bool ground) {
        auto &b = w.create_empty_body();
        b.type = ground ? BodyType::Static : BodyType::Dynamic;
        b.transform.position = {x, y};
        b.linear_damping = b.angular_damping = 0;
        float hx = ground ? 100.f : .25f, hy = ground ? .5f : .25f;
        Polygon p{{{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}}};
        Fixture f;
        f.shape = p;
        if (circles && !ground)
            f.shape = Circle{.25f};
        f.material.friction = .3f;
        f.material.restitution = 0;
        f.material.density = ground ? 0 : 1;
        w.add_fixture(b, f);
        float area = 0, moment = 0;
        for (int i = 0; i < 4; ++i) {
            auto a = p.vertices[i], v = p.vertices[(i + 1) % 4];
            float cross = a.cross(v);
            area += cross;
            moment += cross * (a.dot(a) + a.dot(v) + v.dot(v));
        }
        b.mass = ground ? 0 : std::abs(area) * .5f;
        b.inverse_mass = ground ? 0 : 1 / b.mass;
        b.inertia = ground ? 0 : std::abs(moment) / 12;
        b.inverse_inertia = ground ? 0 : 1 / b.inertia;
        if (circles && !ground) {
            b.mass = 3.14159265358979323846f * .25f * .25f;
            b.inverse_mass = 1 / b.mass;
            b.inertia = .5f * b.mass * .25f * .25f;
            b.inverse_inertia = 1 / b.inertia;
        }
        return &b;
    };
    add(0, -.5f, true);
    int cols = circles ? 50 : 10;
    for (int i = 0; i < count; ++i)
        bodies.push_back(add(float(i % cols) * .65f - 16, .3f + float(i / cols) * .65f, false));
    auto &probe = w.create_body().at(1000, 1000).circle(.1f).velocity(1, 0).build();
    probe.linear_damping = probe.angular_damping = 0;
    std::size_t zero = 0, nonconv = 0, budget = 0;
    int limits = 0;
    double ccd = 0, detect = 0, velocity = 0, position = 0;
    for (int frame = 0; frame < 360; ++frame) {
        w.step();
        for (auto b : bodies) {
            auto bounds = compute_aabb(b->fixtures[0]->shape, b->transform);
            require(std::isfinite(b->transform.position.y), "nonfinite stack state");
            if (std::abs(b->transform.position.x) < 99 && bounds.min.y < -.002f) {
                std::cerr << "first_error frame=" << frame << " bottom=" << bounds.min.y
                          << " x=" << b->transform.position.x << " y=" << b->transform.position.y
                          << " vy=" << b->velocity.y << std::endl;
                require(false, "stack crossed floor");
            }
        }
        if (frame >= 60) {
            auto &s = w.ccd_statistics();
            zero += s.zero_time_repeats;
            nonconv += s.non_convergences;
            budget += s.budget_exhaustions;
            limits += s.limited;
            auto &t = w.step_statistics();
            ccd += t.ccd_ms;
            detect += t.detection_ms;
            velocity += t.velocity_ms;
            position += t.position_ms;
        }
    }
    int asleep = 0;
    for (auto b : bodies)
        asleep += b->sleeping;
    std::cout << "shape=" << (circles ? "circles" : "boxes") << " count=" << count
              << " sleeping=" << asleep << " limited=" << limits << " zero=" << zero
              << " nonconvergence=" << nonconv << " budget=" << budget << " ccd_ms=" << ccd / 300
              << " detection_ms=" << detect / 300 << " velocity_ms=" << velocity / 300
              << " position_ms=" << position / 300 << std::endl;
    require(std::abs(w.simulation_time() - 6) < 1e-5, "world time lost");
    require(std::abs(probe.transform.position.x - 1006) < .003f, "independent motion time lost");
    require(zero == 0 && nonconv == 0 && budget == 0 && limits == 0,
            "persistent contacts still trigger CCD protection");
    if (count == 100 || circles)
        require(asleep == count && limits == 0, "supported stack failed to settle");
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
