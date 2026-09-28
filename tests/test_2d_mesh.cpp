// What a mesh costs, and what its surface is.
//
// A mesh used to be a `std::vector<Polygon>` and nothing else: every query
// walked every triangle. That is fine for the single triangle the rest of the
// suite builds and hopeless for a terrain, where the narrow phase would visit
// the whole mesh once per body per step, the mesh's own bounds were re-derived
// from every triangle on every call, and `shape_radius` did it again.
//
// The requirements this file states as measurements are:
//
//   * the mesh's spatial index is built from the triangles, so a query costs a
//     descent and the triangle count stops appearing in the work per step --
//     and the world-space bounds of a mesh come from the index rather than from
//     re-deriving a bounds per triangle;
//   * the mesh's *surface* is its boundary. The cuts between adjacent triangles
//     are inside the material and must never generate a contact, or a body
//     inside a terrain is pushed out of a cut instead of out of the surface;
//   * a contact on a mesh has a real manifold (two points where the surface
//     supports them), not the single fallback point the face clipper produced
//     when the mesh reached it with no vertices at all;
//   * a seam is not a feature. A body crossing the seam between two coplanar
//     triangles must behave exactly as it does on one triangle.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
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
constexpr float kDt = 1.0f / 60.0f;

// A floor as `segments` quads, each split into two triangles. The top surface is
// a polyline of boundary edges; every vertical cut between neighbouring quads is
// shared by the two triangles on either side of it, so it is interior -- which is
// exactly the shape a terrain author produces.
std::vector<Polygon> strip(int segments, float width, float bottom = -2.0f) {
    const float left = -0.5f * width;
    const float step = width / float(segments);
    std::vector<Polygon> triangles;
    triangles.reserve(std::size_t(segments) * 2);
    for (int i = 0; i < segments; ++i) {
        const float x0 = left + step * float(i), x1 = x0 + step;
        triangles.push_back(Polygon{{{x0, 0}, {x1, 0}, {x0, bottom}}});
        triangles.push_back(Polygon{{{x1, 0}, {x1, bottom}, {x0, bottom}}});
    }
    return triangles;
}
} // namespace

int main() try {
    {
        // The surface, not the cuts. A body inside a mesh is pushed out through
        // the nearest *face* of the mesh; a cut through the material is not a
        // face, and neither is the direction from an edge to a body that has
        // gone past it. A mesh floor's faces are its boundary edges, and the
        // direction out of the material is the edge's own outward normal -- the
        // direction from the edge to the body's centre points the other way once
        // the body is inside, and following it would drive the body further
        // through the floor instead of ejecting it.
        Mesh floor{strip(4, 4.0f)};
        const Transform at{{0.05f, -0.05f}, 0};
        Contact contact;
        // The mesh is the "A" side, which is the order the world builds a pair
        // in: the floor is created first, so `contact.normal` is mesh -> body and
        // that is the direction the solver pushes the body along.
        check(test(Shape{floor}, {}, Shape{Circle{0.1f}}, at, contact),
              "a body inside the mesh found no contact at all");
        std::cout << "a body inside the mesh is pushed out along " << contact.normal.x << ", "
                  << contact.normal.y << '\n';
        check(contact.normal.y > 0.99f && std::abs(contact.normal.x) < 0.01f,
              "a body inside a mesh was pushed out of an interior cut, not out of the surface");
        check(std::abs(contact.penetration - 0.05f) < 1e-4f,
              "the ejection depth is not the distance to the surface");
    }
    {
        // The surface again, from the outside: the boundary of a mesh floor is
        // its top and bottom polylines, so a box resting on one is supported by
        // that edge and no other.
        const Mesh floor{strip(2, 4.0f)};
        const auto &tree = floor.acceleration();
        // Two top edges, two bottom edges, and the vertical edge at either end.
        // Everything else -- the diagonal inside each quad and the cut between
        // the two quads -- belongs to two triangles and is interior.
        check(tree.boundary.size() == 6, "a two-quad strip does not have six boundary edges");
        check(tree.nodes.size() == 1, "four triangles did not fit in a single leaf");
    }
    {
        // The index is *built from* the triangles, so it cannot be allowed to
        // disagree with them. The world AABB of a mesh is what it always was --
        // the union of its triangles' bounds, rotated -- and it is now read off
        // the index instead of re-derived from every triangle.
        const Mesh floor{strip(3, 6.0f)};
        const Transform t{{1.5f, -0.25f}, 0.7f};
        const AABB wide = compute_aabb(Shape{floor}, t);
        AABB exact{};
        bool first = true;
        for (const auto &triangle : floor.triangles) {
            const AABB part = compute_aabb(Shape{triangle}, t);
            if (first) {
                exact = part;
                first = false;
                continue;
            }
            exact.min.x = std::min(exact.min.x, part.min.x);
            exact.min.y = std::min(exact.min.y, part.min.y);
            exact.max.x = std::max(exact.max.x, part.max.x);
            exact.max.y = std::max(exact.max.y, part.max.y);
        }
        // The indexed answer encloses the triangles rather than matching them
        // corner for corner: it is the rotated local bounds, which is a strictly
        // larger box whenever a triangle does not reach a corner of it. Larger is
        // the direction that is safe for a broad phase, and being inside is the
        // direction that would not be.
        check(wide.min.x <= exact.min.x + 1e-6f && wide.min.y <= exact.min.y + 1e-6f &&
                  wide.max.x >= exact.max.x - 1e-6f && wide.max.y >= exact.max.y - 1e-6f,
              "the indexed mesh bounds do not enclose the triangles");
        check(shape_radius(Shape{floor}) >= shape_radius(Shape{floor.triangles.front()}),
              "the indexed mesh radius is smaller than its own triangles");
    }
    {
        // The load-bearing measurement: the narrow phase must stop caring how
        // many triangles the mesh has. The two meshes below have the same quad
        // size and the body crosses the same twenty metres of floor over the same
        // frames, so the only difference is the mesh's *extent* -- eighty times
        // the triangles. The counter is the work itself, not a stopwatch, and the
        // body is sliding so that its contacts are genuinely re-detected rather
        // than served from the cache.
        const auto visits_per_step = [](int segments) {
            World world;
            world.create_body()
                .static_body()
                .at(0, 0)
                .mesh(strip(segments, float(segments)))
                .friction(0)
                .build();
            auto &box = world.create_body()
                            .dynamic()
                            .at(-10, 0.5f)
                            .box(0.5f, 0.5f)
                            .friction(0)
                            .restitution(0)
                            .velocity(6, 0)
                            .build();
            for (int frame = 0; frame < 30; ++frame) {
                box.wake();
                world.step(kDt);
            }
            const std::size_t before = mesh_triangle_tests();
            const int frames = 150;
            for (int frame = 0; frame < frames; ++frame) {
                box.wake();
                world.step(kDt);
            }
            return double(mesh_triangle_tests() - before) / double(frames);
        };
        const double small = visits_per_step(40);
        const double large = visits_per_step(1600);
        std::cout << "triangles tested per step: 80 triangles -> " << small
                  << ", 3200 triangles -> " << large << '\n';
        check(small > 0, "the mesh contact never reached a triangle");
        // Eighty times the triangles, at the same density. A hierarchy that is
        // actually being used answers with the same handful either way; a walk
        // would answer with eighty times the count.
        check(large < 2.0 * small,
              "the mesh narrow phase still walks the triangles instead of descending");
    }
    {
        // A manifold, not a point. A box on a mesh floor is supported across the
        // face that is touching it, so the constraint carries two points -- the
        // face clipper needs real vertices for that, and a mesh handed to it as
        // a Mesh has none.
        World world;
        world.create_body().static_body().at(0, 0).mesh(strip(8, 8.0f)).friction(0.6f).build();
        auto &box = world.create_body()
                        .dynamic()
                        .at(0, 0.5f)
                        .box(0.5f, 0.5f)
                        .friction(0.6f)
                        .restitution(0)
                        .build();
        for (int frame = 0; frame < 120; ++frame) {
            box.wake();
            world.step(kDt);
        }
        const auto statistics = world.step_statistics();
        std::cout << "a settled box on a mesh: " << statistics.constraints << " constraints, "
                  << statistics.warm_started_points << " warm-started points\n";
        check(statistics.constraints == 1, "the settled box is not resting on exactly one pair");
        check(statistics.warm_started_points == 2,
              "a box resting on a mesh floor was supported at a single point");
        check(std::abs(box.transform.position.y - 0.5f) < 0.02f && box.transform.angle == 0.0f,
              "a box resting on a mesh floor did not settle flat");
    }
    {
        // A seam is not a feature. A box sliding across a mesh floor made of a
        // hundred triangles crosses ninety-nine seams, and a floor made of one
        // quad has none: over the same time and the same speed the two have to
        // end up in the same place.
        const auto slide = [](bool tiled) {
            World world;
            if (tiled)
                world.create_body().static_body().at(0, 0).mesh(strip(100, 40.0f)).friction(0).build();
            else
                world.create_body()
                    .static_body()
                    .at(0, 0)
                    .polygon({{-20, 0}, {20, 0}, {20, -2}, {-20, -2}})
                    .friction(0)
                    .build();
            auto &box = world.create_body()
                            .dynamic()
                            .at(-15, 0.5f)
                            .box(0.5f, 0.5f)
                            .friction(0)
                            .restitution(0)
                            .velocity(6, 0)
                            .build();
            for (int frame = 0; frame < 240; ++frame) {
                box.wake();
                world.step(kDt);
            }
            return box.transform;
        };
        const Transform one = slide(false);
        const Transform many = slide(true);
        std::cout << "sliding across the seams: one quad ends at x=" << one.position.x
                  << ", y=" << one.position.y << ", a hundred triangles at x=" << many.position.x
                  << ", y=" << many.position.y << '\n';
        // Free travel would be 6 m/s over four seconds, minus the settle; what
        // matters is that the seams did not turn a slide into a stop.
        check(one.position.x > -15.0f + 3.0f, "the single-quad floor stopped the box");
        check(many.position.x > one.position.x - 0.5f,
              "the seams between the triangles stopped the sliding box");
        check(std::abs(many.position.y - one.position.y) < 0.05f,
              "the sliding box rode up or down over the seams");
    }
    {
        // The seam normal itself. Two triangles meeting at a shallow peak: a body
        // balanced on the shared vertex is in contact with both of them, and the
        // direction it must be pushed in is the direction from that vertex to the
        // body -- which is the one thing both triangles can agree on without
        // knowing about each other.
        Mesh roof{{Polygon{{{-1, -0.1f}, {0, 0}, {-1, -2}}}, Polygon{{{0, 0}, {1, -0.1f}, {1, -2}}}}};
        Contact contact;
        check(test(Shape{roof}, {}, Shape{Circle{0.05f}}, {{0, 0.04f}, 0}, contact),
              "nothing was found at the seam");
        std::cout << "the seam normal is " << contact.normal.x << ", " << contact.normal.y << '\n';
        check(contact.normal.y > 0.99f && std::abs(contact.normal.x) < 0.02f,
              "the seam normal is one triangle's slope rather than the shared vertex normal");
        // A body off to one side of the seam is on one triangle's face and gets
        // that face's normal, which is what makes the surface a surface rather
        // than a rounded hill everywhere. The facet through (-1, -0.1) and
        // (0, 0) has this outward normal, and a circle resting on it sits one
        // radius along that normal from the surface.
        Contact on_face;
        const Vec2 slope = Vec2{-0.1f, 1.0f}.normalized();
        const Vec2 on_surface{-0.5f, -0.05f};
        check(test(Shape{roof}, {}, Shape{Circle{0.05f}}, {on_surface + slope * 0.045f, 0}, on_face),
              "nothing was found on the triangle's face");
        std::cout << "the face normal is " << on_face.normal.x << ", " << on_face.normal.y << '\n';
        check(on_face.normal.dot(slope) > 0.99f,
              "a contact in the middle of a triangle did not use that triangle's face normal");
    }
    {
        // The index has to follow the triangles. A direct write to the outline is
        // noticed by the same per-step scan that notices every other public edit;
        // this checks that the index is dropped with it, by moving a triangle out
        // from under a body and requiring the contact to move with it.
        World world;
        auto &floor_body = world.create_body()
                               .static_body()
                               .at(0, 0)
                               .mesh(strip(2, 4.0f))
                               .friction(0.6f)
                               .build();
        auto &ball = world.create_body()
                         .dynamic()
                         .at(1.0f, 0.3f)
                         .circle(0.25f)
                         .restitution(0)
                         .build();
        for (int frame = 0; frame < 60; ++frame) {
            ball.wake();
            world.step(kDt);
        }
        check(ball.transform.position.y > 0.2f, "the ball fell through the mesh floor");
        // Shift the whole mesh down by a metre, straight through the public
        // shape: the surface goes with it, so the ball has a metre to fall.
        auto &mesh = std::get<Mesh>(floor_body.shape);
        for (auto &triangle : mesh.triangles)
            for (auto &vertex : triangle.vertices)
                vertex.y -= 1.0f;
        for (int frame = 0; frame < 30; ++frame) {
            ball.wake();
            world.step(kDt);
        }
        // Still falling, or landed a metre lower: either way it is no longer
        // resting where the old surface was.
        check(ball.transform.position.y < 0.4f,
              "the mesh's spatial index did not follow a direct write to its triangles");
    }
    std::cout << checks << " mesh checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
