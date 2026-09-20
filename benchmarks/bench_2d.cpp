#include <butter/physics2d/butter2d.h>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

// Counts ordinary C++ allocations made on this serial benchmark thread.
// Allocation timing is included in total/stage timings, not an additive stage.
static std::size_t allocation_calls = 0, allocation_bytes = 0;
static double allocation_ms = 0;
void *operator new(std::size_t bytes) {
    auto start = std::chrono::steady_clock::now();
    void *result = std::malloc(bytes ? bytes : 1);
    allocation_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    if (!result)
        throw std::bad_alloc();
    ++allocation_calls;
    allocation_bytes += bytes;
    return result;
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void *operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
using namespace butter::physics2d;
struct Totals {
    int frames{};
    std::size_t active_body_steps{}, allocations{}, bytes{};
    double total{}, ccd{}, detection{}, velocity{}, position{}, allocation{};
};
int main(int argc, char **argv) try {
    int count = argc > 1 ? std::stoi(argv[1]) : 100;
    if (count < 1 || count > 10000)
        throw std::invalid_argument("body count must be 1..10000");
    bool boxes = argc > 2 && std::string(argv[2]) == "boxes";
    World w;
    w.create_body().static_body().at(0, -.05f).box(float(count), .05f).build();
    std::vector<Body *> bodies;
    for (int i = 0; i < count; ++i) {
        auto b = w.create_body().at((i - count * .5f) * 1.1f, 1).restitution(0);
        if (boxes)
            b.box(.5f, .5f);
        else
            b.circle(.5f);
        bodies.push_back(&b.build());
    }
    auto &probe = w.create_body().kinematic().at(-float(count) - 10, 0).velocity(0, 1).build();
    Totals totals[2];
    int limited = 0;
    for (int frame = 0; frame < 360; ++frame) {
        int active = 0;
        for (auto *b : bodies)
            active += !b->sleeping;
        auto calls = allocation_calls, bytes = allocation_bytes;
        double ams = allocation_ms;
        auto start = std::chrono::steady_clock::now();
        w.step();
        double elapsed =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        auto &t = totals[active ? 0 : 1];
        auto &s = w.step_statistics();
        ++t.frames;
        t.active_body_steps += active;
        t.total += elapsed;
        t.ccd += s.ccd_ms;
        t.detection += s.detection_ms;
        t.velocity += s.velocity_ms;
        t.position += s.position_ms;
        t.allocations += allocation_calls - calls;
        t.bytes += allocation_bytes - bytes;
        t.allocation += allocation_ms - ams;
        limited += w.ccd_statistics().limited;
        for (auto *b : bodies)
            if (compute_aabb(b->shape, b->transform).min.y < -.025f)
                throw std::runtime_error("correctness gate: floor penetration");
        if (std::abs(probe.transform.position.y - (frame + 1) / 60.0f) > .001f)
            throw std::runtime_error("correctness gate: lost independent motion time");
    }
    for (auto *b : bodies)
        if (!b->sleeping)
            throw std::runtime_error("correctness gate: failed to sleep");
    if (std::abs(w.simulation_time() - 6) > 1e-5)
        throw std::runtime_error("correctness gate: lost simulation time");
    if (limited)
        throw std::runtime_error("correctness gate: simple landings required CCD clamps");
    std::cout << "PASS count=" << count << " shape=" << (boxes ? "boxes" : "circles")
              << " simulated_seconds=" << w.simulation_time() << " sleeping=" << count << '\n';
    for (int i = 0; i < 2; ++i) {
        auto &t = totals[i];
        std::cout << (i ? "sleeping" : "active") << " frames=" << t.frames
                  << " active_body_steps=" << t.active_body_steps << " total_ms=" << t.total
                  << " ccd_ms=" << t.ccd << " detection_ms=" << t.detection
                  << " velocity_ms=" << t.velocity << " position_ms=" << t.position
                  << " ordinary_allocations=" << t.allocations << " allocated_bytes=" << t.bytes
                  << " allocation_ms_inclusive=" << t.allocation << '\n';
    }
    std::cout << "TomCat adapter: unavailable in this repository. Compare only "
                 "matching activity counts and correctness gates; no historical "
                 "speedup ratio.\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
