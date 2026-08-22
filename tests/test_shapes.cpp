#include <butter/butter.h>

#include <cassert>
#include <cmath>
#include <vector>

using namespace butter;
using namespace butter::math;

int main() {
    ContactInfo contact;

    CapsuleCollider capsule(0.25f, 2.0f);
    BoxCollider box(0.5f, 0.5f, 0.5f);
    assert(Collider::test(capsule, Transform{{0, 0, 0}},
                           box, Transform{{0.6f, 0, 0}}, contact));
    assert(contact.penetration > 0.0f);
    contact.reset();
    assert(!Collider::test(capsule, Transform{{0, 0, 0}},
                            box, Transform{{2.0f, 0, 0}}, contact));
    CapsuleCollider other_capsule(0.2f, 1.0f);
    contact.reset();
    assert(Collider::test(capsule, Transform{{0, 0, 0}},
                          other_capsule, Transform{{0.35f, 0.1f, 0}}, contact));
    contact.reset();
    assert(!Collider::test(capsule, Transform{{0, 0, 0}},
                           other_capsule, Transform{{2.0f, 0, 0}}, contact));

    const std::vector<Vec3> tetra{
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    ConvexCollider convex_a(tetra);
    ConvexCollider convex_b(tetra);
    assert(Collider::test(convex_a, Transform{{0, 0, 0}},
                          convex_b, Transform{{0.25f, 0.25f, 0.25f}}, contact));
    assert(contact.penetration > 0.01f);
    assert(contact.normal.x > 0.0f && contact.normal.y > 0.0f &&
           contact.normal.z > 0.0f);
    contact.reset();
    assert(!Collider::test(convex_a, Transform{{0, 0, 0}},
                            convex_b, Transform{{3, 3, 3}}, contact));
    contact.reset();
    assert(Collider::test(convex_a, Transform::identity(),
                          convex_b, Transform{{1, 0, 0}}, contact));
    contact.reset();
    assert(!Collider::test(convex_a, Transform::identity(),
                           convex_b, Transform{{1.001f, 0, 0}}, contact));

    SphereCollider convex_sphere(0.12f);
    contact.reset();
    assert(Collider::test(convex_a, Transform::identity(),
                          convex_sphere, Transform{{0.3f, 0.3f, 0.3f}}, contact));
    contact.reset();
    assert(!Collider::test(convex_a, Transform::identity(),
                           convex_sphere, Transform{{2, 2, 2}}, contact));
    contact.reset();
    // Coincident convex centers use a deterministic fallback search direction.
    assert(Collider::test(convex_a, Transform::identity(),
                          convex_b, Transform::identity(), contact));

    BoxCollider ground_box(10.0f, 0.5f, 10.0f);
    contact.reset();
    assert(Collider::test(convex_a, Transform{{0, 0.4f, 0}},
                          ground_box, Transform::identity(), contact));
    assert(contact.normal.y < 0.0f && contact.penetration > 0.0f);

    const std::vector<Vec3> vertices{{-2, 0, -2}, {2, 0, -2},
                                     {2, 0, 2}, {-2, 0, 2}};
    const std::vector<MeshCollider::Triangle> triangles{{0, 1, 2}, {0, 2, 3}};
    MeshCollider mesh(vertices, triangles);
    SphereCollider sphere(0.5f);
    contact.reset();
    assert(Collider::test(sphere, Transform{{0, 0.25f, 0}},
                          mesh, Transform::identity(), contact));
    contact.reset();
    assert(!Collider::test(sphere, Transform{{0, 2, 0}},
                           mesh, Transform::identity(), contact));

    contact.reset();
    assert(Collider::test(box, Transform{{0, 0.25f, 0}},
                          mesh, Transform::identity(), contact));
    assert(contact.penetration > 0.0f);
    contact.reset();
    assert(Collider::test(mesh, Transform::identity(),
                          box, Transform{{0, 0.25f, 0}}, contact));

    // The mesh narrow phase also handles a finite capsule segment, not just
    // its bounding sphere.  This capsule crosses the indexed triangle plane.
    CapsuleCollider mesh_capsule(0.15f, 1.0f);
    contact.reset();
    assert(Collider::test(mesh_capsule, Transform{{0, 0, 0}},
                          mesh, Transform::identity(), contact));
    contact.reset();
    assert(!Collider::test(mesh_capsule, Transform{{3, 0, 0}},
                           mesh, Transform::identity(), contact));

    // A long segment can run parallel to a triangle plane with both segment
    // endpoints outside the triangle's footprint. The interior projection is
    // still the closest pair and must be considered by capsule-mesh.
    CapsuleCollider parallel_capsule(0.15f, 6.0f);
    const Quat rotate_to_x = Quat::from_axis_angle({0, 0, 1},
                                                    1.57079632679f);
    contact.reset();
    assert(Collider::test(parallel_capsule,
                          Transform{{0, 0.1f, 0}, rotate_to_x},
                          MeshCollider({{-1, 0, -1}, {1, 0, -1},
                                        {0, 0, 1}},
                                       {{0, 1, 2}}),
                          Transform::identity(), contact));

    // Convex-vs-mesh uses support mapping against every triangle rather than
    // falling back to the mesh AABB.
    contact.reset();
    assert(Collider::test(convex_a, Transform{{0.2f, -0.1f, 0.2f}},
                          mesh, Transform::identity(), contact));
    assert(contact.penetration > 0.0f);
    contact.reset();
    assert(!Collider::test(convex_a, Transform{{0, 3, 0}},
                           mesh, Transform::identity(), contact));

    // Indexed mesh pairs use triangle SAT and preserve coplanar contact.
    const std::vector<Vec3> shifted_vertices{{-1.5f, 0, -1.5f},
                                              {2.5f, 0, -1.5f},
                                              {2.5f, 0, 2.5f},
                                              {-1.5f, 0, 2.5f}};
    MeshCollider shifted_mesh(shifted_vertices, triangles);
    contact.reset();
    assert(Collider::test(mesh, Transform::identity(),
                          shifted_mesh, Transform::identity(), contact));
    contact.reset();
    assert(!Collider::test(mesh, Transform{{0, 0, 10}},
                           shifted_mesh, Transform::identity(), contact));

    MeshCollider invalid_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}},
                              {{0, 1, 2}});
    contact.reset();
    assert(!Collider::test(sphere, Transform::identity(),
                           invalid_mesh, Transform::identity(), contact));

    ConvexCollider empty_convex;
    contact.reset();
    assert(!Collider::test(empty_convex, Transform::identity(),
                           sphere, Transform::identity(), contact));

    // A point cloud is interpreted as its actual convex hull. Check a cube
    // hull against spheres at several outside/inside positions so GJK does
    // not silently regress to an AABB overlap test.
    const std::vector<Vec3> cube_points{
        {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
        {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
    ConvexCollider cube_hull(cube_points);
    SphereCollider probe(0.2f);
    for (int ix = -3; ix <= 3; ++ix) {
        for (int iy = -2; iy <= 2; ++iy) {
            for (int iz = -3; iz <= 3; ++iz) {
                const Vec3 position{ix * 0.6f, iy * 0.6f, iz * 0.6f};
                const float dx = std::max(0.0f, std::abs(position.x) - 1.0f);
                const float dy = std::max(0.0f, std::abs(position.y) - 1.0f);
                const float dz = std::max(0.0f, std::abs(position.z) - 1.0f);
                const bool expected = std::sqrt(dx * dx + dy * dy + dz * dz) <= 0.20001f;
                contact.reset();
                const bool hit = Collider::test(cube_hull, Transform::identity(),
                                                probe, Transform{position}, contact);
                assert(hit == expected);
            }
        }
    }
    contact.reset();
    assert(!Collider::test(empty_convex, Transform{{0, 0, 10}},
                           sphere, Transform{{0, 0, 10}}, contact));

    // A 2-D convex hull lying on a mesh surface should still be treated as a
    // touching/overlapping pair; the GJK fallback must not require a 3-D
    // tetrahedral hull.
    ConvexCollider flat_convex({{-0.5f, 0, -0.5f}, {0.8f, 0, -0.5f},
                                {0.8f, 0, 0.8f}, {-0.5f, 0, 0.8f}});
    contact.reset();
    assert(Collider::test(flat_convex, Transform::identity(),
                          mesh, Transform::identity(), contact));
    contact.reset();
    assert(!Collider::test(flat_convex, Transform{{5, 0, 0}},
                           mesh, Transform::identity(), contact));
    contact.reset();
    assert(!Collider::test(flat_convex, Transform{{0, 0.01f, 0}},
                           mesh, Transform::identity(), contact));

    // The triangle proxy used by convex-mesh GJK carries its world centroid
    // in the transform, so translated meshes must not be tested against the
    // origin by the centre-direction fallback.
    contact.reset();
    assert(Collider::test(convex_a, Transform{{100.2f, -0.1f, 100.2f}},
                          mesh, Transform{{100, 0, 100}}, contact));

    // Convex point clouds need not be centred at the body origin.  The
    // contact normal follows the actual hull centres, not the body positions.
    ConvexCollider offset_hull({{10, 0, 0}, {11, 0, 0},
                                {10, 1, 0}, {10, 0, 1}});
    contact.reset();
    assert(Collider::test(offset_hull, Transform::identity(),
                          convex_a, Transform{{9.5f, 0, 0}}, contact));
    assert(contact.normal.x < -0.5f);

    ConvexCollider line_hull({{-1, 0, 0}, {1, 0, 0}, {1, 0, 0}});
    contact.reset();
    assert(Collider::test(line_hull, Transform::identity(),
                          sphere, Transform{{0, 0.4f, 0}}, contact));
    contact.reset();
    assert(!Collider::test(line_hull, Transform::identity(),
                           sphere, Transform{{0, 1.0f, 0}}, contact));

    // Collider offsets live in body-local space.  A rotated body must rotate
    // the offset translation for both narrow phase and broad-phase AABBs.
    SphereCollider offset_sphere(0.25f);
    offset_sphere.set_offset(Transform{{1, 0, 0}});
    const Quat quarter_turn = Quat::from_axis_angle({0, 0, 1}, 1.57079632679f);
    contact.reset();
    assert(Collider::test(offset_sphere, Transform{{0, 0, 0}, quarter_turn},
                           sphere, Transform{{0, 1, 0}}, contact));
    const AABB offset_bounds = offset_sphere.compute_aabb(
        Transform{{0, 0, 0}, quarter_turn});
    assert(offset_bounds.center().y > 0.9f &&
           std::abs(offset_bounds.center().x) < 1.0e-4f);

    BoxCollider rotated_box(1.0f, 2.0f, 3.0f);
    const AABB rotated_bounds = rotated_box.compute_aabb(
        Transform{{0, 0, 0}, quarter_turn});
    assert(std::abs(rotated_bounds.extents().x - 4.0f) < 1.0e-4f);
    assert(std::abs(rotated_bounds.extents().y - 2.0f) < 1.0e-4f);
    assert(std::abs(rotated_bounds.extents().z - 6.0f) < 1.0e-4f);

    // A sphere exactly on a mesh plane gets a deterministic A->B normal even
    // when the triangle winding is reversed.
    MeshCollider reversed_mesh({{-1, 0, -1}, {1, 0, -1}, {0, 0, 1}},
                               {{0, 2, 1}});
    contact.reset();
    assert(Collider::test(sphere, Transform{{0, 0, 0}},
                          reversed_mesh, Transform::identity(), contact));
    assert(contact.normal.y > 0.0f);
    return 0;
}
