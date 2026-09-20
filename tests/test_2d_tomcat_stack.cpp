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
    int frames = argc > 3 ? std::stoi(argv[3]) : 360;
    require(count > 0 && count <= 10000 && frames > 60,
            "count must be 1..10000 and frames must exceed warmup");
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
    auto &probe = w.create_body().at(1000, 100000).circle(.1f).velocity(60, 0).build();
    probe.linear_damping = probe.angular_damping = 0;
    std::size_t zero = 0, nonconv = 0, budget = 0;
    int limits = 0;
    double ccd = 0, detect = 0, velocity = 0, position = 0;
    double total = 0, sleeping_ms = 0, cache_ms = 0, active_ms = 0, resting_ms = 0;
    std::size_t active_body_steps = 0, projection_candidates = 0, corrections = 0, cached = 0;
    int active_frames = 0, resting_frames = 0, first_sleep = -1;
    for (int frame = 0; frame < frames; ++frame) {
        int active = 0;
        for (auto b : bodies)
            active += !b->sleeping;
        auto begin = std::chrono::steady_clock::now();
        w.step();
        double elapsed =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin)
                .count();
        bool all_sleep = true;
        for (auto b : bodies)
            all_sleep &= b->sleeping;
        if (all_sleep && first_sleep < 0)
            first_sleep = frame + 1;
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
            total += elapsed;
            sleeping_ms += t.sleeping_ms;
            cache_ms += t.cache_ms;
            active_body_steps += active;
            if (active) {
                ++active_frames;
                active_ms += elapsed;
            } else {
                ++resting_frames;
                resting_ms += elapsed;
            }
            projection_candidates += t.projection_candidates;
            corrections += t.position_corrections;
            cached += t.cached_manifolds + t.sleeping_contacts;
        }
        if (frames > 360 && (frame + 1) % 300 == 0) {
            float speed = 0, angular = 0;
            int awake = 0;
            for (auto b : bodies) {
                speed = std::max(speed, b->velocity.length());
                angular = std::max(angular, std::abs(b->angular_velocity));
                awake += !b->sleeping;
            }
            const auto &s = w.step_statistics();
            std::cout << "seconds=" << (frame + 1) / 60 << " awake=" << awake << " speed=" << speed
                      << " angular=" << angular << " penetration=" << s.max_penetration
                      << " correction=" << s.max_correction << " moving_groups=" << s.moving_groups
                      << " settling_groups=" << s.settling_groups << std::endl;
        }
    }
    int asleep = 0;
    for (auto b : bodies)
        asleep += b->sleeping;
    std::cout << "shape=" << (circles ? "circles" : "boxes") << " count=" << count
              << " sleeping=" << asleep << " limited=" << limits << " zero=" << zero
              << " nonconvergence=" << nonconv << " budget=" << budget
              << " ccd_ms=" << ccd / (frames - 60) << " detection_ms=" << detect / (frames - 60)
              << " velocity_ms=" << velocity / (frames - 60)
              << " position_ms=" << position / (frames - 60) << std::endl;
    std::cout << "total_ms=" << total / (frames - 60) << " cache_ms=" << cache_ms / (frames - 60)
              << " sleep_ms=" << sleeping_ms / (frames - 60) << " active_frames=" << active_frames
              << " active_body_steps=" << active_body_steps << " active_ms=" << active_ms
              << " resting_frames=" << resting_frames << " resting_ms=" << resting_ms
              << " first_sleep_seconds=" << first_sleep / 60.0 << " corrections=" << corrections
              << " projection_candidates=" << projection_candidates << " cached_contacts=" << cached
              << std::endl;
    require(std::abs(w.simulation_time() - frames / 60.0) < 1e-4, "world time lost");
    require(std::abs(probe.transform.position.x - (1000.0f + frames)) < .001f,
            "independent motion time lost");
    require(zero == 0 && nonconv == 0 && budget == 0 && limits == 0,
            "persistent contacts still trigger CCD protection");
    if (count == 100 || circles)
        require(asleep == count && limits == 0, "supported stack failed to settle");
    if (frames >= 1800)
        require(asleep == count, "long stack observation did not settle");
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
