#include <butter/physics2d/butter2d.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace butter::physics2d;

static int checks = 0;
static void check(bool ok, const char *message) {
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}
static bool near(float a, float b, float tolerance = 1.0e-4f) {
    return std::abs(a - b) <= tolerance * std::max(1.0f, std::max(std::abs(a), std::abs(b)));
}

int main() try {
    constexpr float pi = 3.14159265358979f;
    // Shape-level mass properties: closed forms against the analytic values.
    {
        const float r = 0.75f, density = 2.5f;
        auto data = mass_data(Shape{Circle{r}}, density);
        check(near(data.mass, density * pi * r * r), "circle mass");
        check(near(data.inertia, 0.5f * data.mass * r * r), "circle inertia");
        check(data.center.length() == 0, "circle centroid");
    }
    {
        const float hx = 0.4f, hy = 1.3f, density = 3.0f;
        auto data = mass_data(Shape{Box{{hx, hy}}}, density);
        check(near(data.mass, density * 4 * hx * hy), "box mass");
        check(near(data.inertia, data.mass * (hx * hx + hy * hy) / 3.0f), "box inertia");
    }
    {
        // A square outline must reproduce the analytic box properties, and its
        // centroid must be the middle of the outline, not the vertex average.
        const float h = 1.2f, density = 1.7f;
        Polygon square{{{-h, -h}, {h, -h}, {h, h}, {-h, h}}};
        auto data = mass_data(Shape{square}, density);
        auto box = mass_data(Shape{Box{{h, h}}}, density);
        check(near(data.mass, box.mass), "polygon mass matches box");
        check(near(data.inertia, box.inertia), "polygon inertia matches box");
        check(data.center.length() < 1.0e-6f, "polygon centroid at origin");
    }
    {
        // An off-centre outline must report its real centroid and inertia about
        // that centroid, not about the body origin.
        const float density = 1.0f;
        Polygon free{{{0, 0}, {2, 0}, {2, 2}, {0, 2}}};
        auto data = mass_data(Shape{free}, density);
        check(near(data.mass, 4.0f), "offset polygon mass");
        check(near(data.center.x, 1.0f) && near(data.center.y, 1.0f), "offset polygon centroid");
        check(near(data.inertia, data.mass * (4.0f + 4.0f) / 12.0f), "offset polygon inertia");
    }
    {
        const float r = 0.3f, hl = 1.1f, density = 1.0f;
        auto data = mass_data(Shape{Capsule{r, hl}}, density);
        const float rect = density * 4 * r * hl, cap = density * pi * r * r;
        check(near(data.mass, rect + cap), "capsule mass");
        check(near(data.inertia, rect * (r * r + hl * hl) / 3.0f + 0.5f * cap * r * r),
              "capsule inertia");
    }
    {
        // Degenerate outlines carry no mass instead of inventing one.
        Polygon line{{{0, 0}, {1, 0}, {2, 0}}};
        check(mass_data(Shape{line}, 1.0f).mass == 0, "degenerate outline mass");
        Polygon empty{};
        check(mass_data(Shape{empty}, 1.0f).mass == 0, "empty outline mass");
    }
    // Builder bodies derive inertia (and the centroid) from their shape while
    // keeping the requested mass.
    {
        World w;
        auto &ball = w.create_body().circle(0.5f).build();
        check(near(ball.mass, 1.0f), "builder keeps requested mass");
        check(near(ball.inertia, 0.5f * 1.0f * 0.25f), "builder circle inertia");
        auto &cube = w.create_body().box(0.5f, 0.5f).mass(4.0f).build();
        check(near(cube.mass, 4.0f), "builder explicit mass");
        check(near(cube.inertia, 4.0f * (0.25f + 0.25f) / 3.0f), "builder box inertia");
        auto &dense = w.create_body().box(0.5f, 0.5f).density(3.0f).build();
        check(near(dense.mass, 3.0f), "builder density mass");
        auto &off = w.create_body().polygon({{0, 0}, {2, 0}, {2, 2}, {0, 2}}).build();
        check(near(off.local_center.x, 1.0f) && near(off.local_center.y, 1.0f),
              "builder polygon centroid");
        check(off.center_of_mass().x == off.transform.position.x + 1.0f,
              "builder centroid in world space");
    }
    // Compound fixture bodies recompute mass/centroid/inertia on every change.
    {
        World w;
        auto &body = w.create_empty_body();
        Fixture f;
        f.shape = Box{{1, 1}};
        f.material.density = 2.0f;
        w.add_fixture(body, f);
        check(near(body.mass, 8.0f), "fixture mass");
        check(near(body.inertia, 8.0f * (1 + 1) / 3.0f), "fixture inertia");
        check(body.local_center.length() == 0, "centred fixture centroid");
        // Shift the second fixture and the combined centroid must move with it.
        f.local.position = {5, 0};
        w.add_fixture(body, f);
        check(near(body.mass, 16.0f), "compound mass");
        check(near(body.local_center.x, 2.5f), "compound centroid");
        check(near(body.inertia, 2 * (8.0f * (1 + 1) / 3.0f) + 2 * 8.0f * 2.5f * 2.5f),
              "compound inertia about the combined centroid");
        // Density zero geometry must not contribute mass.
        f.material.density = 0;
        auto &extra = w.add_fixture(body, f);
        check(near(body.mass, 16.0f), "zero-density fixture ignored");
        w.destroy_fixture(extra);
        check(near(body.mass, 16.0f), "destroy recomputes mass");
        // Explicit mass data opts out of the automatic recompute.
        body.set_mass_data({7.0f, {0, 0}, 9.0f});
        f.material.density = 1.0f;
        w.add_fixture(body, f);
        check(near(body.mass, 7.0f) && near(body.inertia, 9.0f), "explicit mass data preserved");
    }
    // A body whose centroid is off the origin must rotate about the centroid.
    {
        World::Config cfg;
        cfg.gravity = {};
        World w(cfg);
        auto &body = w.create_empty_body();
        body.transform.position = {0, 0};
        body.angular_damping = 0;
        Fixture f;
        f.shape = Box{{0.5f, 0.5f}};
        f.local.position = {2, 0};
        w.add_fixture(body, f);
        check(near(body.local_center.x, 2.0f), "offset fixture centroid");
        const Vec2 center = body.center_of_mass();
        body.angular_velocity = 1.0f;
        for (int i = 0; i < 30; ++i)
            w.step(1.0f / 60.0f);
        const Vec2 moved = body.center_of_mass();
        check((moved - center).length() < 1.0e-4f, "rotation must not translate the centroid");
        check(std::abs(body.transform.angle - 0.5f) < 1.0e-4f, "rotation applied");
        check(std::abs(body.transform.position.x - (-2.0f)) > 0.5f,
              "body origin orbits the centroid");
    }
    std::cout << checks << " mass-property checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
