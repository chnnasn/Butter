#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "butter/shapes/box.h"
#include "butter/shapes/capsule.h"
#include "butter/shapes/collider.h"
#include "butter/shapes/convex.h"
#include "butter/shapes/mesh.h"
#include "butter/shapes/sphere.h"

namespace butter {
namespace detail {

inline math::Vec3 collider_center(const Collider& collider, const math::Transform& transform) {
    // The offset is expressed in the rigid body's local frame.  Adding it
    // directly to the world position loses the body's rotation and makes an
    // offset collider drift to the wrong side when the body rotates.  Keep
    // this in lock-step with Transform::operator* and the point-cloud/mesh
    // paths below.
    return transform.transform_point(collider.offset.position);
}

inline math::Quat collider_rotation(const Collider& collider, const math::Transform& transform) {
    return (transform * collider.offset).rotation;
}

// A collider's body origin is not necessarily the geometric centre: convex
// point clouds and offset child colliders are both valid.  Use the world AABB
// centre when available for GJK's initial direction and normal orientation;
// fall back to the body/offset centre for empty or malformed shapes.
inline math::Vec3 collider_shape_center(const Collider& collider,
                                        const math::Transform& transform) {
    const math::AABB bounds = collider.compute_aabb(transform);
    return bounds.is_empty() ? collider_center(collider, transform)
                             : bounds.center();
}

inline math::Vec3 closest_point_on_segment(const math::Vec3& p,
                                           const math::Vec3& a,
                                           const math::Vec3& b) {
    const math::Vec3 ab = b - a;
    const float length_squared = ab.length_squared();
    if (length_squared < 1.0e-12f) return a;
    float t = (p - a).dot(ab) / length_squared;
    t = std::clamp(t, 0.0f, 1.0f);
    return a + ab * t;
}

inline void segment_endpoints(const CapsuleCollider& capsule,
                              const math::Transform& transform,
                              math::Vec3& a,
                              math::Vec3& b) {
    const math::Vec3 center = collider_center(capsule, transform);
    const math::Quat rotation = collider_rotation(capsule, transform);
    const float h = capsule.half_height();
    a = center + rotation.rotate({0, -h, 0});
    b = center + rotation.rotate({0, h, 0});
}

// Squared closest-point distance between two finite segments.  The helper is
// declared here because capsule-capsule narrow phase uses the same exact
// segment primitive as capsule-mesh below.
inline float closest_points_segment_segment(const math::Vec3& p1,
                                            const math::Vec3& q1,
                                            const math::Vec3& p2,
                                            const math::Vec3& q2,
                                            math::Vec3& c1,
                                            math::Vec3& c2);

inline bool test_sphere_sphere(const SphereCollider& a, const math::Transform& ta,
                               const SphereCollider& b, const math::Transform& tb,
                               ContactInfo& contact) {
    const math::Vec3 center_a = collider_center(a, ta);
    const math::Vec3 center_b = collider_center(b, tb);
    const float radius_sum = a.radius() + b.radius();

    math::Vec3 delta = center_b - center_a;
    const float distance_squared = delta.length_squared();
    if (distance_squared > radius_sum * radius_sum) return false;

    const float distance = std::sqrt(distance_squared);
    const math::Vec3 normal = distance > 1.0e-6f ? delta / distance : math::Vec3{0, 1, 0};
    contact.normal = normal;
    contact.penetration = radius_sum - distance;
    contact.point = center_a + normal * a.radius();
    contact.hit = true;
    return true;
}

inline bool test_sphere_box(const SphereCollider& sphere, const math::Transform& sphere_transform,
                            const BoxCollider& box, const math::Transform& box_transform,
                            ContactInfo& contact) {
    const math::Vec3 sphere_center = collider_center(sphere, sphere_transform);
    const math::Vec3 box_center = collider_center(box, box_transform);
    const math::Quat box_rotation = collider_rotation(box, box_transform);
    const math::Quat inv_rotation = box_rotation.conjugate();
    const math::Vec3 local_center = inv_rotation.rotate(sphere_center - box_center);
    const math::Vec3& h = box.half_extents();

    const math::Vec3 closest{
        std::clamp(local_center.x, -h.x, h.x),
        std::clamp(local_center.y, -h.y, h.y),
        std::clamp(local_center.z, -h.z, h.z),
    };

    const math::Vec3 delta = local_center - closest;
    const float distance_squared = delta.length_squared();
    const float radius = sphere.radius();

    math::Vec3 normal_local;
    float penetration = 0;

    if (distance_squared > 1.0e-12f) {
        const float distance = std::sqrt(distance_squared);
        if (distance > radius) return false;
        normal_local = -delta / distance;
        penetration = radius - distance;
        contact.point = box_center + box_rotation.rotate(closest);
    } else {
        // Sphere center is inside the box; push it out along the axis of least penetration.
        const math::Vec3 distances{
            h.x - std::abs(local_center.x),
            h.y - std::abs(local_center.y),
            h.z - std::abs(local_center.z),
        };
        if (distances.x < distances.y && distances.x < distances.z) {
            normal_local = {local_center.x >= 0 ? 1.0f : -1.0f, 0, 0};
            penetration = radius + distances.x;
        } else if (distances.y < distances.z) {
            normal_local = {0, local_center.y >= 0 ? 1.0f : -1.0f, 0};
            penetration = radius + distances.y;
        } else {
            normal_local = {0, 0, local_center.z >= 0 ? 1.0f : -1.0f};
            penetration = radius + distances.z;
        }
        contact.point = box_center + box_rotation.rotate(local_center);
    }

    // Normal points from sphere (A) toward box (B).
    contact.normal = box_rotation.rotate(normal_local);
    contact.penetration = penetration;
    contact.hit = true;
    return true;
}

inline bool test_sphere_capsule(const SphereCollider& sphere, const math::Transform& sphere_transform,
                                const CapsuleCollider& capsule, const math::Transform& capsule_transform,
                                ContactInfo& contact) {
    const math::Vec3 sphere_center = collider_center(sphere, sphere_transform);
    math::Vec3 segment_a;
    math::Vec3 segment_b;
    segment_endpoints(capsule, capsule_transform, segment_a, segment_b);

    const math::Vec3 closest = closest_point_on_segment(sphere_center, segment_a, segment_b);
    math::Vec3 delta = closest - sphere_center;
    const float distance_squared = delta.length_squared();
    const float radius_sum = sphere.radius() + capsule.radius();
    if (distance_squared > radius_sum * radius_sum) return false;

    const float distance = std::sqrt(distance_squared);
    math::Vec3 normal = distance > 1.0e-6f ? delta / distance : math::Vec3{0, 1, 0};
    // normal currently points from sphere center toward the capsule segment.
    contact.normal = normal;
    contact.penetration = radius_sum - distance;
    contact.point = sphere_center + normal * sphere.radius();
    contact.hit = true;
    return true;
}

inline bool test_capsule_capsule(const CapsuleCollider& a, const math::Transform& ta,
                                 const CapsuleCollider& b, const math::Transform& tb,
                                 ContactInfo& contact) {
    math::Vec3 a0, a1, b0, b1;
    segment_endpoints(a, ta, a0, a1);
    segment_endpoints(b, tb, b0, b1);

    math::Vec3 pa;
    math::Vec3 pb;
    closest_points_segment_segment(a0, a1, b0, b1, pa, pb);

    math::Vec3 delta = pb - pa;
    const float distance_squared = delta.length_squared();
    const float radius_sum = a.radius() + b.radius();
    if (distance_squared > radius_sum * radius_sum) return false;

    const float distance = std::sqrt(distance_squared);
    const math::Vec3 normal = distance > 1.0e-6f ? delta / distance : math::Vec3{0, 1, 0};
    contact.normal = normal;
    contact.penetration = radius_sum - distance;
    contact.point = pa + normal * a.radius();
    contact.hit = true;
    return true;
}

inline math::Vec3 closest_point_on_aabb(const math::Vec3& point,
                                        const math::Vec3& half_extents) {
    return {std::clamp(point.x, -half_extents.x, half_extents.x),
            std::clamp(point.y, -half_extents.y, half_extents.y),
            std::clamp(point.z, -half_extents.z, half_extents.z)};
}

// Exact enough for a convex box: distance from the capsule's centre segment
// to the OBB is minimized in the box's local space. The squared distance is a
// convex function of segment parameter, so a bounded ternary search converges
// without replacing the capsule by a sphere or an AABB.
inline bool test_capsule_box(const CapsuleCollider& capsule,
                             const math::Transform& capsule_transform,
                             const BoxCollider& box,
                             const math::Transform& box_transform,
                             ContactInfo& contact) {
    math::Vec3 segment_a, segment_b;
    segment_endpoints(capsule, capsule_transform, segment_a, segment_b);
    const math::Vec3 box_center = collider_center(box, box_transform);
    const math::Quat box_rotation = collider_rotation(box, box_transform);
    const math::Quat inverse_rotation = box_rotation.conjugate();
    const math::Vec3 local_a = inverse_rotation.rotate(segment_a - box_center);
    const math::Vec3 local_b = inverse_rotation.rotate(segment_b - box_center);
    const math::Vec3 half_extents = box.half_extents();

    auto evaluate = [&](float t, math::Vec3* out_segment,
                        math::Vec3* out_box) {
        const math::Vec3 point = local_a + (local_b - local_a) * t;
        const math::Vec3 closest = closest_point_on_aabb(point, half_extents);
        if (out_segment) *out_segment = point;
        if (out_box) *out_box = closest;
        return (point - closest).length_squared();
    };

    float lo = 0.0f;
    float hi = 1.0f;
    for (int i = 0; i < 28; ++i) {
        const float m1 = lo + (hi - lo) / 3.0f;
        const float m2 = hi - (hi - lo) / 3.0f;
        if (evaluate(m1, nullptr, nullptr) < evaluate(m2, nullptr, nullptr)) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    const float t = (lo + hi) * 0.5f;
    math::Vec3 local_segment;
    math::Vec3 local_box;
    const float distance_squared = evaluate(t, &local_segment, &local_box);
    const float radius = capsule.radius();
    if (distance_squared > radius * radius) return false;

    const float distance = std::sqrt(std::max(0.0f, distance_squared));
    math::Vec3 normal_local;
    if (distance > 1.0e-6f) {
        // Normal points from capsule (A) toward box (B).
        normal_local = (local_box - local_segment) / distance;
    } else {
        const math::Vec3 centre = (local_a + local_b) * 0.5f;
        const math::Vec3 d{half_extents.x - std::abs(centre.x),
                           half_extents.y - std::abs(centre.y),
                           half_extents.z - std::abs(centre.z)};
        if (d.x <= d.y && d.x <= d.z) {
            normal_local = {centre.x >= 0.0f ? 1.0f : -1.0f, 0, 0};
        } else if (d.y <= d.z) {
            normal_local = {0, centre.y >= 0.0f ? 1.0f : -1.0f, 0};
        } else {
            normal_local = {0, 0, centre.z >= 0.0f ? 1.0f : -1.0f};
        }
    }

    const math::Vec3 normal = box_rotation.rotate(normal_local).normalized();
    contact.normal = normal;
    contact.penetration = radius - distance;
    contact.point = box_center + box_rotation.rotate(local_segment + normal_local * radius);
    contact.hit = true;
    return true;
}

struct GjkSupportPoint {
    math::Vec3 minkowski{};
    math::Vec3 point_a{};
    math::Vec3 point_b{};
};

inline math::Vec3 support_point(const Collider& collider,
                                const math::Transform& transform,
                                const math::Vec3& direction) {
    const math::Vec3 unit_direction = direction.length_squared() > 1.0e-12f
        ? direction.normalized() : math::Vec3::up();
    const math::Vec3 center = collider_center(collider, transform);
    const math::Quat rotation = collider_rotation(collider, transform);
    const math::Vec3 local_direction = rotation.conjugate().rotate(unit_direction);

    switch (collider.type()) {
        case Collider::Type::Sphere: {
            const auto& sphere = static_cast<const SphereCollider&>(collider);
            return center + unit_direction * sphere.radius();
        }
        case Collider::Type::Box: {
            const auto& box = static_cast<const BoxCollider&>(collider);
            const math::Vec3 h = box.half_extents();
            return center + rotation.rotate({
                local_direction.x >= 0.0f ? h.x : -h.x,
                local_direction.y >= 0.0f ? h.y : -h.y,
                local_direction.z >= 0.0f ? h.z : -h.z});
        }
        case Collider::Type::Capsule: {
            const auto& capsule = static_cast<const CapsuleCollider&>(collider);
            const float sign = local_direction.y >= 0.0f ? 1.0f : -1.0f;
            const math::Vec3 endpoint = rotation.rotate({0, sign * capsule.half_height(), 0});
            return center + endpoint + unit_direction * capsule.radius();
        }
        case Collider::Type::Convex: {
            const auto& convex = static_cast<const ConvexCollider&>(collider);
            math::Vec3 best = center;
            float best_dot = -std::numeric_limits<float>::infinity();
            const math::Transform world = transform * collider.offset;
            for (const auto& point : convex.points()) {
                const math::Vec3 world_point = world.transform_point(point);
                const float value = world_point.dot(direction);
                if (value > best_dot) {
                    best_dot = value;
                    best = world_point;
                }
            }
            return best;
        }
        default:
            return center;
    }
}

inline GjkSupportPoint gjk_support(const Collider& a, const math::Transform& ta,
                                   const Collider& b, const math::Transform& tb,
                                   const math::Vec3& direction) {
    const math::Vec3 point_a = support_point(a, ta, direction);
    const math::Vec3 point_b = support_point(b, tb, -direction);
    return {point_a - point_b, point_a, point_b};
}

inline bool handle_gjk_simplex(std::vector<GjkSupportPoint>& simplex,
                               math::Vec3& direction) {
    const math::Vec3 origin{};
    const auto& a = simplex.back();
    const math::Vec3 ao = origin - a.minkowski;

    if (simplex.size() == 2) {
        const math::Vec3 ab = simplex[0].minkowski - a.minkowski;
        if (ab.dot(ao) > 0.0f) {
            direction = ab.cross(ao).cross(ab);
            if (direction.length_squared() < 1.0e-12f) {
                direction = math::Vec3{-ab.y, ab.x, 0.0f};
                if (direction.length_squared() < 1.0e-12f)
                    direction = math::Vec3{0.0f, -ab.z, ab.y};
            }
        } else {
            simplex = {a};
            direction = ao;
        }
        return false;
    }

    if (simplex.size() == 3) {
        const auto& b = simplex[1];
        const auto& c = simplex[0];
        const math::Vec3 ab = b.minkowski - a.minkowski;
        const math::Vec3 ac = c.minkowski - a.minkowski;
        const math::Vec3 abc = ab.cross(ac);

        if (abc.cross(ac).dot(ao) > 0.0f) {
            if (ac.dot(ao) > 0.0f) {
                simplex = {c, a};
                direction = ac.cross(ao).cross(ac);
            } else {
                simplex = {a};
                direction = ao;
            }
            return false;
        }
        if (ab.cross(abc).dot(ao) > 0.0f) {
            if (ab.dot(ao) > 0.0f) {
                simplex = {b, a};
                direction = ab.cross(ao).cross(ab);
            } else {
                simplex = {a};
                direction = ao;
            }
            return false;
        }
        if (abc.dot(ao) > 0.0f) {
            direction = abc;
        } else {
            simplex = {b, c, a};
            direction = -abc;
        }
        return false;
    }

    // Tetrahedron: if the origin lies outside any face adjacent to A, keep
    // that face and continue; otherwise the Minkowski volumes intersect.
    const auto& b = simplex[2];
    const auto& c = simplex[1];
    const auto& d = simplex[0];
    const math::Vec3 ab = b.minkowski - a.minkowski;
    const math::Vec3 ac = c.minkowski - a.minkowski;
    const math::Vec3 ad = d.minkowski - a.minkowski;
    const math::Vec3 ao3 = -a.minkowski;
    const math::Vec3 abc = ab.cross(ac);
    const math::Vec3 acd = ac.cross(ad);
    const math::Vec3 adb = ad.cross(ab);
    if (abc.dot(ao3) > 0.0f) {
        simplex = {c, b, a};
        direction = abc;
        return false;
    }
    if (acd.dot(ao3) > 0.0f) {
        simplex = {d, c, a};
        direction = acd;
        return false;
    }
    if (adb.dot(ao3) > 0.0f) {
        simplex = {b, d, a};
        direction = adb;
        return false;
    }
    return true;
}

struct EpaFace {
    int a{0};
    int b{0};
    int c{0};
    math::Vec3 normal{0, 1, 0};
    float distance{0.0f};
};

struct EpaResult {
    math::Vec3 normal{0, 1, 0};
    float penetration{0.0f};
};

// Build an outward-facing EPA face.  `opposite` is supplied for the initial
// tetrahedron so the face can be oriented away from the tetrahedron interior;
// newly-created horizon faces only need the origin-side distance check.
inline bool make_epa_face(const std::vector<GjkSupportPoint>& vertices,
                          int ia, int ib, int ic, int opposite,
                          EpaFace& face) {
    math::Vec3 normal =
        (vertices[ib].minkowski - vertices[ia].minkowski).cross(
            vertices[ic].minkowski - vertices[ia].minkowski);
    const float length_squared = normal.length_squared();
    if (length_squared <= 1.0e-14f) return false;
    normal /= std::sqrt(length_squared);

    if (opposite >= 0) {
        const math::Vec3 toward_opposite =
            vertices[opposite].minkowski - vertices[ia].minkowski;
        if (normal.dot(toward_opposite) > 0.0f) {
            std::swap(ib, ic);
            normal = -normal;
        }
    }

    float distance = normal.dot(vertices[ia].minkowski);
    if (distance < 0.0f) {
        std::swap(ib, ic);
        normal = -normal;
        distance = -distance;
    }
    face = {ia, ib, ic, normal, distance};
    return std::isfinite(distance);
}

inline bool epa_penetration(const Collider& a, const math::Transform& ta,
                            const Collider& b, const math::Transform& tb,
                            const std::vector<GjkSupportPoint>& simplex,
                            EpaResult& result) {
    if (simplex.size() < 4) return false;

    std::vector<GjkSupportPoint> vertices = simplex;
    std::vector<EpaFace> faces;
    faces.reserve(64);
    // Four faces of the tetrahedron returned by GJK.  The opposite index is
    // used to orient each face consistently before the first expansion.
    const int initial_faces[4][4] = {
        {0, 1, 2, 3}, {0, 3, 1, 2}, {0, 2, 3, 1}, {1, 3, 2, 0}};
    for (const auto& indices : initial_faces) {
        EpaFace face;
        if (make_epa_face(vertices, indices[0], indices[1], indices[2],
                          indices[3], face)) {
            faces.push_back(face);
        }
    }
    if (faces.size() < 4) return false;

    struct Edge {
        int a{0};
        int b{0};
    };
    auto add_horizon_edge = [](std::vector<Edge>& edges, int from, int to) {
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (it->a == to && it->b == from) {
                edges.erase(it);
                return;
            }
        }
        edges.push_back({from, to});
    };

    constexpr int max_iterations = 64;
    constexpr std::size_t max_vertices = 128;
    constexpr float absolute_tolerance = 1.0e-5f;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (faces.empty()) return false;

        std::size_t closest_index = 0;
        for (std::size_t i = 1; i < faces.size(); ++i) {
            if (faces[i].distance < faces[closest_index].distance)
                closest_index = i;
        }
        const EpaFace closest = faces[closest_index];
        const GjkSupportPoint support =
            gjk_support(a, ta, b, tb, closest.normal);
        const float support_distance = support.minkowski.dot(closest.normal);
        if (!std::isfinite(support_distance)) return false;

        const float tolerance = absolute_tolerance *
            std::max(1.0f, std::abs(support_distance));
        if (support_distance - closest.distance <= tolerance) {
            result.normal = closest.normal;
            result.penetration = std::max(0.0f, closest.distance);
            return true;
        }

        bool duplicate = false;
        for (const auto& vertex : vertices) {
            if ((vertex.minkowski - support.minkowski).length_squared() <
                absolute_tolerance * absolute_tolerance) {
                duplicate = true;
                break;
            }
        }
        if (duplicate || vertices.size() >= max_vertices) {
            result.normal = closest.normal;
            result.penetration = std::max(0.0f, closest.distance);
            return true;
        }

        const int new_index = static_cast<int>(vertices.size());
        vertices.push_back(support);

        std::vector<Edge> horizon;
        horizon.reserve(faces.size() * 3);
        for (std::size_t i = faces.size(); i-- > 0;) {
            const EpaFace& face = faces[i];
            const bool visible = face.normal.dot(
                support.minkowski - vertices[face.a].minkowski) > absolute_tolerance;
            if (!visible) continue;
            add_horizon_edge(horizon, face.a, face.b);
            add_horizon_edge(horizon, face.b, face.c);
            add_horizon_edge(horizon, face.c, face.a);
            faces.erase(faces.begin() + static_cast<std::ptrdiff_t>(i));
        }
        if (horizon.empty()) return false;

        for (const Edge& edge : horizon) {
            EpaFace face;
            if (make_epa_face(vertices, edge.a, edge.b, new_index, -1, face)) {
                faces.push_back(face);
            }
        }
    }
    return false;
}

inline float projected_overlap(const Collider& a, const math::Transform& ta,
                               const Collider& b, const math::Transform& tb,
                               const math::Vec3& raw_axis) {
    const float length_squared = raw_axis.length_squared();
    if (length_squared <= 1.0e-14f) return 0.0f;
    const math::Vec3 axis = raw_axis / std::sqrt(length_squared);
    const float a_min = support_point(a, ta, -axis).dot(axis);
    const float a_max = support_point(a, ta, axis).dot(axis);
    const float b_min = support_point(b, tb, -axis).dot(axis);
    const float b_max = support_point(b, tb, axis).dot(axis);
    return std::max(0.0f, std::min(a_max, b_max) - std::max(a_min, b_min));
}

inline math::Vec3 stable_support_contact(const Collider& a,
                                         const math::Transform& ta,
                                         const Collider& b,
                                         const math::Transform& tb,
                                         const math::Vec3& normal) {
    const math::Vec3 point_a = support_point(a, ta, normal);
    const math::Vec3 point_b = support_point(b, tb, -normal);
    // A support query on a large box returns an arbitrary corner when the
    // direction is exactly vertical.  Averaging those raw points would put a
    // ground contact several metres away from the body and create a huge
    // artificial rotational effective mass. Preserve the support planes but
    // anchor tangential coordinates on the smaller/finite-thickness shape
    // (usually the dynamic object); this also gives a useful contact for a
    // zero-thickness mesh triangle.
    const math::AABB bounds_a = a.compute_aabb(ta);
    const math::AABB bounds_b = b.compute_aabb(tb);
    const math::Vec3 centre_a = bounds_a.is_empty()
        ? collider_center(a, ta) : bounds_a.center();
    const math::Vec3 centre_b = bounds_b.is_empty()
        ? collider_center(b, tb) : bounds_b.center();
    math::Vec3 tangent_anchor = (centre_a + centre_b) * 0.5f;
    if (!bounds_a.is_empty() && !bounds_b.is_empty()) {
        const math::Vec3 extents_a = bounds_a.extents();
        const math::Vec3 extents_b = bounds_b.extents();
        const bool surface_a = std::min({extents_a.x, extents_a.y, extents_a.z}) < 1.0e-4f;
        const bool surface_b = std::min({extents_b.x, extents_b.y, extents_b.z}) < 1.0e-4f;
        const float volume_a = extents_a.x * extents_a.y * extents_a.z;
        const float volume_b = extents_b.x * extents_b.y * extents_b.z;
        if (surface_a && !surface_b) {
            tangent_anchor = centre_b;
        } else if (surface_b && !surface_a) {
            tangent_anchor = centre_a;
        } else if (volume_a > volume_b * 4.0f) {
            tangent_anchor = centre_b;
        } else if (volume_b > volume_a * 4.0f) {
            tangent_anchor = centre_a;
        }
    }
    const float plane = (point_a.dot(normal) + point_b.dot(normal)) * 0.5f;
    return tangent_anchor + normal * (plane - tangent_anchor.dot(normal));
}

inline bool test_convex_pair(const Collider& a, const math::Transform& ta,
                             const Collider& b, const math::Transform& tb,
                             ContactInfo& contact) {
    // An empty point cloud has no support point and therefore represents no
    // geometric volume.  Treat it as non-collidable instead of letting the
    // fallback support-at-origin accidentally overlap another shape.
    if ((a.type() == Collider::Type::Convex &&
         static_cast<const ConvexCollider&>(a).points().empty()) ||
        (b.type() == Collider::Type::Convex &&
         static_cast<const ConvexCollider&>(b).points().empty())) {
        return false;
    }
    math::Vec3 direction = collider_shape_center(b, tb) -
                           collider_shape_center(a, ta);
    if (direction.length_squared() < 1.0e-12f) direction = math::Vec3::up();

    const auto orient_contact_normal = [&](math::Vec3 normal) {
        const math::Vec3 centre_delta = collider_shape_center(b, tb) -
                                         collider_shape_center(a, ta);
        if (centre_delta.length_squared() > 1.0e-12f &&
            centre_delta.dot(normal) < 0.0f) {
            normal = -normal;
        }
        return normal;
    };

    const auto report_touching_contact = [&](math::Vec3 normal) {
        if (normal.length_squared() < 1.0e-12f) {
            normal = collider_shape_center(b, tb) -
                     collider_shape_center(a, ta);
        }
        if (normal.length_squared() < 1.0e-12f) normal = math::Vec3::up();
        normal = orient_contact_normal(normal);
        normal = normal.normalized();
        contact.normal = normal;
        contact.penetration = 0.0f;
        contact.point = stable_support_contact(a, ta, b, tb, normal);
        contact.hit = true;
        return true;
    };

    std::vector<GjkSupportPoint> simplex;
    simplex.reserve(4);
    simplex.push_back(gjk_support(a, ta, b, tb, direction));
    if (simplex.back().minkowski.dot(direction) < -1.0e-7f) return false;
    if (simplex.back().minkowski.length_squared() <= 1.0e-14f)
        return report_touching_contact(direction);
    direction = -simplex.back().minkowski;

    for (int iteration = 0; iteration < 32; ++iteration) {
        if (direction.length_squared() < 1.0e-12f)
            return report_touching_contact(direction);
        const GjkSupportPoint point = gjk_support(a, ta, b, tb, direction);
        if (point.minkowski.dot(direction) < -1.0e-7f) return false;
        if (point.minkowski.length_squared() <= 1.0e-14f)
            return report_touching_contact(direction);
        simplex.push_back(point);
        if (handle_gjk_simplex(simplex, direction)) {
            EpaResult epa;
            math::Vec3 normal;
            float penetration = 0.0f;
            if (epa_penetration(a, ta, b, tb, simplex, epa)) {
                normal = epa.normal;
                penetration = epa.penetration;
            } else {
                normal = direction.length_squared() > 1.0e-12f
                    ? direction.normalized()
                    : (collider_shape_center(b, tb) -
                       collider_shape_center(a, ta)).normalized();
            }
            normal = orient_contact_normal(normal);
            if (normal.length_squared() < 1.0e-12f) normal = math::Vec3::up();
            contact.normal = normal;
            // EPA normally gives the minimum translation depth.  A GJK
            // simplex can be coplanar for thin/near-degenerate hulls, in
            // which case EPA may converge to a zero-distance face; use the
            // exact interval overlap on that face as a stable fallback.
            const float axis_overlap = projected_overlap(a, ta, b, tb, normal);
            contact.penetration = penetration > 1.0e-4f
                ? penetration : axis_overlap;
            contact.point = stable_support_contact(a, ta, b, tb, normal);
            contact.hit = true;
            return true;
        }
    }
    return false;
}

inline math::Vec3 closest_point_on_triangle(const math::Vec3& p,
                                             const math::Vec3& a,
                                             const math::Vec3& b,
                                             const math::Vec3& c) {
    const math::Vec3 ab = b - a;
    const math::Vec3 ac = c - a;
    const math::Vec3 ap = p - a;
    const float d1 = ab.dot(ap);
    const float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f) return a;

    const math::Vec3 bp = p - b;
    const float d3 = ab.dot(bp);
    const float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3) return b;

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        const float v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const math::Vec3 cp = p - c;
    const float d5 = ab.dot(cp);
    const float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6) return c;

    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        const float w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        const math::Vec3 bc = c - b;
        const float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + bc * w;
    }

    const float denominator = 1.0f / (va + vb + vc);
    const float v = vb * denominator;
    const float w = vc * denominator;
    return a + ab * v + ac * w;
}

inline bool point_in_triangle(const math::Vec3& point,
                              const math::Vec3& a,
                              const math::Vec3& b,
                              const math::Vec3& c,
                              const math::Vec3& unit_normal) {
    constexpr float epsilon = 1.0e-6f;
    const float side_a = unit_normal.dot((b - a).cross(point - a));
    const float side_b = unit_normal.dot((c - b).cross(point - b));
    const float side_c = unit_normal.dot((a - c).cross(point - c));
    return side_a >= -epsilon && side_b >= -epsilon && side_c >= -epsilon;
}

inline bool segment_triangle_intersection(const math::Vec3& p0,
                                          const math::Vec3& p1,
                                          const math::Vec3& a,
                                          const math::Vec3& b,
                                          const math::Vec3& c,
                                          math::Vec3& out_point) {
    constexpr float epsilon = 1.0e-8f;
    const math::Vec3 direction = p1 - p0;
    const math::Vec3 normal = (b - a).cross(c - a);
    if (normal.length_squared() <= epsilon * epsilon) return false;

    const float denominator = normal.dot(direction);
    if (std::abs(denominator) <= epsilon) return false;
    const float t = normal.dot(a - p0) / denominator;
    if (t < -epsilon || t > 1.0f + epsilon) return false;

    const math::Vec3 point = p0 + direction * std::clamp(t, 0.0f, 1.0f);
    const math::Vec3 c0 = (b - a).cross(point - a);
    const math::Vec3 c1 = (c - b).cross(point - b);
    const math::Vec3 c2 = (a - c).cross(point - c);
    if (normal.dot(c0) < -epsilon || normal.dot(c1) < -epsilon ||
        normal.dot(c2) < -epsilon) return false;
    out_point = point;
    return true;
}

inline float closest_points_segment_segment(const math::Vec3& p1,
                                            const math::Vec3& q1,
                                            const math::Vec3& p2,
                                            const math::Vec3& q2,
                                            math::Vec3& c1,
                                            math::Vec3& c2) {
    const math::Vec3 d1 = q1 - p1;
    const math::Vec3 d2 = q2 - p2;
    const math::Vec3 r = p1 - p2;
    const float a = d1.dot(d1);
    const float e = d2.dot(d2);
    const float f = d2.dot(r);
    float s = 0.0f;
    float t = 0.0f;
    if (a <= 1.0e-12f && e <= 1.0e-12f) {
        c1 = p1;
        c2 = p2;
        return (c1 - c2).length_squared();
    }
    if (a <= 1.0e-12f) {
        s = 0.0f;
        t = std::clamp(f / e, 0.0f, 1.0f);
    } else {
        const float c = d1.dot(r);
        if (e <= 1.0e-12f) {
            t = 0.0f;
            s = std::clamp(-c / a, 0.0f, 1.0f);
        } else {
            const float b = d1.dot(d2);
            const float denominator = a * e - b * b;
            s = denominator != 0.0f ? std::clamp((b * f - c * e) / denominator,
                                                 0.0f, 1.0f) : 0.0f;
            t = (b * s + f) / e;
            if (t < 0.0f) {
                t = 0.0f;
                s = std::clamp(-c / a, 0.0f, 1.0f);
            } else if (t > 1.0f) {
                t = 1.0f;
                s = std::clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    return (c1 - c2).length_squared();
}

inline bool mesh_triangle(const MeshCollider& mesh, const math::Transform& transform,
                          const MeshCollider::Triangle& triangle,
                          math::Vec3& a, math::Vec3& b, math::Vec3& c) {
    if (triangle.a >= mesh.vertices().size() ||
        triangle.b >= mesh.vertices().size() ||
        triangle.c >= mesh.vertices().size()) return false;
    const math::Transform world = transform * mesh.offset;
    a = world.transform_point(mesh.vertices()[triangle.a]);
    b = world.transform_point(mesh.vertices()[triangle.b]);
    c = world.transform_point(mesh.vertices()[triangle.c]);
    // Repeated/collinear indices do not describe a surface triangle and
    // should not produce an infinite SAT penetration or a bogus contact.
    return (b - a).cross(c - a).length_squared() > 1.0e-14f;
}

inline bool test_sphere_mesh(const SphereCollider& sphere,
                             const math::Transform& sphere_transform,
                             const MeshCollider& mesh,
                             const math::Transform& mesh_transform,
                             ContactInfo& contact) {
    const math::Vec3 centre = collider_center(sphere, sphere_transform);
    float best_distance_squared = std::numeric_limits<float>::infinity();
    math::Vec3 best_normal{0, 1, 0};
    for (const auto& triangle : mesh.triangles()) {
        math::Vec3 a, b, c;
        if (!mesh_triangle(mesh, mesh_transform, triangle, a, b, c)) continue;
        const math::Vec3 point = closest_point_on_triangle(centre, a, b, c);
        const math::Vec3 delta = point - centre;
        const float distance_squared = delta.length_squared();
        if (distance_squared >= best_distance_squared) continue;
        best_distance_squared = distance_squared;
        const math::Vec3 face_normal = (b - a).cross(c - a).normalized();
        best_normal = distance_squared > 1.0e-12f
            ? delta.normalized() : (face_normal.length_squared() > 0.0f
                ? face_normal : math::Vec3::up());
        // A sphere centre exactly on a triangle plane has no separating
        // direction from the closest-point query.  Orient the fallback face
        // normal from the sphere toward the triangle so reversing the mesh
        // winding cannot make the solver push the sphere through the surface.
        if (distance_squared <= 1.0e-12f &&
            best_normal.dot(((a + b + c) / 3.0f) - centre) < 0.0f) {
            best_normal = -best_normal;
        }
    }
    if (best_distance_squared == std::numeric_limits<float>::infinity() ||
        best_distance_squared > sphere.radius() * sphere.radius()) return false;
    const float distance = std::sqrt(std::max(0.0f, best_distance_squared));
    contact.normal = best_normal;
    contact.penetration = sphere.radius() - distance;
    contact.point = centre + best_normal * sphere.radius();
    contact.hit = true;
    return true;
}

inline bool test_capsule_mesh(const CapsuleCollider& capsule,
                              const math::Transform& capsule_transform,
                              const MeshCollider& mesh,
                              const math::Transform& mesh_transform,
                              ContactInfo& contact) {
    math::Vec3 segment_a, segment_b;
    segment_endpoints(capsule, capsule_transform, segment_a, segment_b);
    float best_distance_squared = std::numeric_limits<float>::infinity();
    math::Vec3 best_segment{};
    math::Vec3 best_triangle{};
    math::Vec3 best_face_normal{0, 1, 0};
    math::Vec3 best_face_center{};
    for (const auto& triangle : mesh.triangles()) {
        math::Vec3 a, b, c;
        if (!mesh_triangle(mesh, mesh_transform, triangle, a, b, c)) continue;
        const math::Vec3 plane_normal = (b - a).cross(c - a);
        const float plane_length_squared = plane_normal.length_squared();
        const math::Vec3 face_normal = plane_length_squared > 1.0e-12f
            ? plane_normal / std::sqrt(plane_length_squared) : math::Vec3::up();
        const math::Vec3 face_center = (a + b + c) / 3.0f;
        math::Vec3 intersection{};
        if (plane_length_squared > 1.0e-12f &&
            segment_triangle_intersection(segment_a, segment_b, a, b, c,
                                           intersection)) {
            best_distance_squared = 0.0f;
            best_segment = intersection;
            best_triangle = intersection;
            best_face_normal = face_normal;
            best_face_center = face_center;
        }
        auto consider = [&](const math::Vec3& segment_point,
                            const math::Vec3& triangle_point) {
            const float distance_squared = (triangle_point - segment_point).length_squared();
            if (distance_squared < best_distance_squared) {
                best_distance_squared = distance_squared;
                best_segment = segment_point;
                best_triangle = triangle_point;
                best_face_normal = face_normal;
                best_face_center = face_center;
            }
        };

        // When the capsule segment runs parallel to a triangle plane, the
        // closest pair can be an interior point of the triangle and an
        // interior point of the segment. Endpoint-to-face and edge-to-segment
        // checks alone miss a long segment whose endpoints lie outside the
        // triangle. Project the segment onto the plane and test that 2D
        // segment against the triangle boundary/interior in this case.
        const float signed_distance_a = face_normal.dot(segment_a - a);
        const float signed_distance_b = face_normal.dot(segment_b - a);
        if (std::abs(signed_distance_b - signed_distance_a) <= 1.0e-6f) {
            const math::Vec3 projected_a = segment_a - face_normal * signed_distance_a;
            const math::Vec3 projected_b = segment_b - face_normal * signed_distance_b;
            if (point_in_triangle(projected_a, a, b, c, face_normal))
                consider(segment_a, projected_a);
            if (point_in_triangle(projected_b, a, b, c, face_normal))
                consider(segment_b, projected_b);

            const math::Vec3 triangle_vertices[3] = {a, b, c};
            const math::Vec3 triangle_edges[3] = {b, c, a};
            for (int edge = 0; edge < 3; ++edge) {
                math::Vec3 projected_segment_point;
                math::Vec3 edge_point;
                closest_points_segment_segment(projected_a, projected_b,
                                               triangle_vertices[edge],
                                               triangle_edges[edge],
                                               projected_segment_point, edge_point);
                const math::Vec3 projected_direction = projected_b - projected_a;
                const float projected_length_squared = projected_direction.length_squared();
                const float parameter = projected_length_squared > 1.0e-12f
                    ? std::clamp((projected_segment_point - projected_a).dot(
                                     projected_direction) / projected_length_squared,
                                 0.0f, 1.0f)
                    : 0.0f;
                consider(segment_a + (segment_b - segment_a) * parameter, edge_point);
            }
        }

        consider(segment_a, closest_point_on_triangle(segment_a, a, b, c));
        consider(segment_b, closest_point_on_triangle(segment_b, a, b, c));
        const math::Vec3 triangle_vertices[3] = {a, b, c};
        const math::Vec3 triangle_edges[3] = {b, c, a};
        for (int edge = 0; edge < 3; ++edge) {
            math::Vec3 segment_point;
            math::Vec3 edge_point;
            closest_points_segment_segment(segment_a, segment_b,
                                           triangle_vertices[edge], triangle_edges[edge],
                                           segment_point, edge_point);
            consider(segment_point, edge_point);
        }

    }
    if (best_distance_squared == std::numeric_limits<float>::infinity() ||
        best_distance_squared > capsule.radius() * capsule.radius()) return false;
    const float distance = std::sqrt(std::max(0.0f, best_distance_squared));
    const math::Vec3 delta = best_triangle - best_segment;
    math::Vec3 normal = distance > 1.0e-6f ? delta / distance : best_face_normal;
    if (distance <= 1.0e-6f &&
        normal.dot(best_face_center - best_segment) < 0.0f) {
        normal = -normal;
    }
    contact.normal = normal;
    contact.penetration = capsule.radius() - distance;
    contact.point = best_segment + normal * capsule.radius();
    contact.hit = true;
    return true;
}

inline bool test_box_mesh(const BoxCollider& box,
                          const math::Transform& box_transform,
                          const MeshCollider& mesh,
                          const math::Transform& mesh_transform,
                          ContactInfo& contact) {
    const math::Vec3 box_center = collider_center(box, box_transform);
    const math::Quat box_rotation = collider_rotation(box, box_transform);
    const math::Quat inverse_rotation = box_rotation.conjugate();
    const math::Vec3 h = box.half_extents();
    float best_overlap = std::numeric_limits<float>::infinity();
    math::Vec3 best_normal_local{0, 1, 0};
    math::Vec3 best_triangle_center{};
    bool hit = false;

    for (const auto& triangle : mesh.triangles()) {
        math::Vec3 wa, wb, wc;
        if (!mesh_triangle(mesh, mesh_transform, triangle, wa, wb, wc)) continue;
        const math::Vec3 a = inverse_rotation.rotate(wa - box_center);
        const math::Vec3 b = inverse_rotation.rotate(wb - box_center);
        const math::Vec3 c = inverse_rotation.rotate(wc - box_center);
        const std::array<math::Vec3, 3> edges{b - a, c - b, a - c};
        std::array<math::Vec3, 13> axes{};
        axes[0] = {1, 0, 0};
        axes[1] = {0, 1, 0};
        axes[2] = {0, 0, 1};
        int axis_count = 3;
        const math::Vec3 triangle_normal = (b - a).cross(c - a);
        if (triangle_normal.length_squared() > 1.0e-12f)
            axes[axis_count++] = triangle_normal.normalized();
        const math::Vec3 box_axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (const auto& edge : edges) {
            for (const auto& box_axis : box_axes) {
                const math::Vec3 axis = edge.cross(box_axis);
                if (axis.length_squared() > 1.0e-12f && axis_count < 13)
                    axes[axis_count++] = axis.normalized();
            }
        }

        float triangle_overlap = std::numeric_limits<float>::infinity();
        math::Vec3 triangle_axis{0, 1, 0};
        bool separated = false;
        for (int i = 0; i < axis_count; ++i) {
            const math::Vec3 axis = axes[i];
            const float box_radius = h.x * std::abs(axis.x) +
                                     h.y * std::abs(axis.y) +
                                     h.z * std::abs(axis.z);
            const float p0 = a.dot(axis);
            const float p1 = b.dot(axis);
            const float p2 = c.dot(axis);
            const float tri_min = std::min({p0, p1, p2});
            const float tri_max = std::max({p0, p1, p2});
            float overlap = std::min(box_radius, tri_max) -
                            std::max(-box_radius, tri_min);
            // Along the triangle's face normal the triangle has zero
            // thickness.  The interval length is therefore always zero even
            // when the box has sunk deeply through the surface; use the
            // box-radius minus plane-distance penetration for that axis.
            if (triangle_normal.length_squared() > 1.0e-12f &&
                std::abs(axis.dot(triangle_normal.normalized())) > 1.0f - 1.0e-4f &&
                std::abs(tri_max - tri_min) <= 1.0e-4f) {
                overlap = box_radius - std::abs((tri_min + tri_max) * 0.5f);
            }
            // A triangle is zero-thickness along its face normal. An overlap
            // of exactly zero therefore still means the solid box crosses
            // the triangle plane; only a negative separation rejects it.
            if (overlap < -1.0e-6f) {
                separated = true;
                break;
            }
            if (overlap < triangle_overlap) {
                triangle_overlap = overlap;
                triangle_axis = axis;
            }
        }
        if (separated || triangle_overlap >= best_overlap) continue;
        const math::Vec3 centre = (a + b + c) / 3.0f;
        if (centre.dot(triangle_axis) < 0.0f) triangle_axis = -triangle_axis;
        best_overlap = triangle_overlap;
        best_normal_local = triangle_axis;
        best_triangle_center = (wa + wb + wc) / 3.0f;
        hit = true;
    }
    if (!hit) return false;
    const math::Vec3 normal = box_rotation.rotate(best_normal_local).normalized();
    const float box_radius = h.x * std::abs(best_normal_local.x) +
                             h.y * std::abs(best_normal_local.y) +
                             h.z * std::abs(best_normal_local.z);
    contact.normal = normal;
    contact.penetration = best_overlap;
    contact.point = box_center + normal * box_radius;
    (void)best_triangle_center;
    contact.hit = true;
    return true;
}

inline bool test_triangle_triangle(const math::Vec3& a0, const math::Vec3& a1,
                                   const math::Vec3& a2, const math::Vec3& b0,
                                   const math::Vec3& b1, const math::Vec3& b2,
                                   ContactInfo& contact) {
    const math::Vec3 a_edges[3] = {a1 - a0, a2 - a1, a0 - a2};
    const math::Vec3 b_edges[3] = {b1 - b0, b2 - b1, b0 - b2};
    const math::Vec3 a_normal = a_edges[0].cross(a_edges[1]);
    const math::Vec3 b_normal = b_edges[0].cross(b_edges[1]);
    if (a_normal.length_squared() <= 1.0e-14f ||
        b_normal.length_squared() <= 1.0e-14f) {
        return false;
    }
    // In addition to the usual triangle normals and edge x edge axes, add
    // edge x face-normal axes.  The latter are essential for coplanar
    // triangles: all edge x edge axes can be zero while the triangles are
    // still separated in their common plane.
    std::array<math::Vec3, 17> axes{};
    int count = 0;
    if (a_normal.length_squared() > 1.0e-12f) axes[count++] = a_normal.normalized();
    if (b_normal.length_squared() > 1.0e-12f) axes[count++] = b_normal.normalized();
    for (const auto& ea : a_edges) {
        for (const auto& eb : b_edges) {
            const math::Vec3 axis = ea.cross(eb);
            if (axis.length_squared() > 1.0e-12f && count < static_cast<int>(axes.size()))
                axes[count++] = axis.normalized();
        }
    }
    if (a_normal.length_squared() > 1.0e-12f) {
        const math::Vec3 n = a_normal.normalized();
        for (const auto& edge : a_edges) {
            const math::Vec3 axis = edge.cross(n);
            if (axis.length_squared() > 1.0e-12f && count < static_cast<int>(axes.size()))
                axes[count++] = axis.normalized();
        }
    }
    if (b_normal.length_squared() > 1.0e-12f) {
        const math::Vec3 n = b_normal.normalized();
        for (const auto& edge : b_edges) {
            const math::Vec3 axis = edge.cross(n);
            if (axis.length_squared() > 1.0e-12f && count < static_cast<int>(axes.size()))
                axes[count++] = axis.normalized();
        }
    }
    const math::Vec3 centre_a = (a0 + a1 + a2) / 3.0f;
    const math::Vec3 centre_b = (b0 + b1 + b2) / 3.0f;
    float best_overlap = std::numeric_limits<float>::infinity();
    math::Vec3 best_axis{0, 1, 0};
    for (int i = 0; i < count; ++i) {
        const math::Vec3 axis = axes[i];
        const float pa[3] = {a0.dot(axis), a1.dot(axis), a2.dot(axis)};
        const float pb[3] = {b0.dot(axis), b1.dot(axis), b2.dot(axis)};
        const float amin = std::min({pa[0], pa[1], pa[2]});
        const float amax = std::max({pa[0], pa[1], pa[2]});
        const float bmin = std::min({pb[0], pb[1], pb[2]});
        const float bmax = std::max({pb[0], pb[1], pb[2]});
        const float overlap = std::min(amax, bmax) - std::max(amin, bmin);
        if (overlap < -1.0e-6f) return false;
        if (overlap < best_overlap) {
            best_overlap = overlap;
            best_axis = axis;
        }
    }
    if ((centre_b - centre_a).dot(best_axis) < 0.0f) best_axis = -best_axis;
    contact.normal = best_axis;
    contact.penetration = best_overlap;
    contact.point = (centre_a + centre_b) * 0.5f;
    contact.hit = true;
    return true;
}

inline bool test_mesh_mesh(const MeshCollider& a, const math::Transform& ta,
                           const MeshCollider& b, const math::Transform& tb,
                           ContactInfo& contact) {
    for (const auto& triangle_a : a.triangles()) {
        math::Vec3 a0, a1, a2;
        if (!mesh_triangle(a, ta, triangle_a, a0, a1, a2)) continue;
        for (const auto& triangle_b : b.triangles()) {
            math::Vec3 b0, b1, b2;
            if (!mesh_triangle(b, tb, triangle_b, b0, b1, b2)) continue;
            if (test_triangle_triangle(a0, a1, a2, b0, b1, b2, contact)) return true;
        }
    }
    return false;
}

inline bool test_convex_mesh(const Collider& convex,
                             const math::Transform& convex_transform,
                             const MeshCollider& mesh,
                             const math::Transform& mesh_transform,
                             ContactInfo& contact) {
    for (const auto& triangle : mesh.triangles()) {
        math::Vec3 a, b, c;
        if (!mesh_triangle(mesh, mesh_transform, triangle, a, b, c)) continue;
        // A triangle is itself a convex support shape. GJK handles the
        // convex-vs-triangle case without reducing the mesh to its AABB.
        // Keep the triangle in a local frame and carry its world-space
        // centroid in the transform. Feeding world coordinates with an
        // identity transform would make the GJK centre/normal fallback point
        // at the origin for meshes translated away from it.
        const math::Vec3 centre = (a + b + c) / 3.0f;
        ConvexCollider triangle_shape({a - centre, b - centre, c - centre});
        if (test_convex_pair(convex, convex_transform, triangle_shape,
                             math::Transform{centre}, contact)) return true;
    }
    return false;
}

inline void obb_axes(const math::Quat& rotation, math::Vec3 axes[3]) {
    const math::Mat3 m = rotation.to_mat3();
    axes[0] = {m.m[0][0], m.m[1][0], m.m[2][0]};
    axes[1] = {m.m[0][1], m.m[1][1], m.m[2][1]};
    axes[2] = {m.m[0][2], m.m[1][2], m.m[2][2]};
}

inline float project_obb_radius(const math::Vec3& axis,
                                const math::Vec3 axes[3],
                                const math::Vec3& half_extents) {
    return half_extents.x * std::abs(axis.dot(axes[0])) +
           half_extents.y * std::abs(axis.dot(axes[1])) +
           half_extents.z * std::abs(axis.dot(axes[2]));
}

inline bool test_box_box(const BoxCollider& a, const math::Transform& ta,
                         const BoxCollider& b, const math::Transform& tb,
                         ContactInfo& contact) {
    const math::Vec3 center_a = collider_center(a, ta);
    const math::Vec3 center_b = collider_center(b, tb);
    const math::Quat rot_a = collider_rotation(a, ta);
    const math::Quat rot_b = collider_rotation(b, tb);

    math::Vec3 axes_a[3];
    math::Vec3 axes_b[3];
    obb_axes(rot_a, axes_a);
    obb_axes(rot_b, axes_b);

    math::Vec3 axes[15];
    for (int i = 0; i < 3; ++i) {
        axes[i] = axes_a[i];
        axes[3 + i] = axes_b[i];
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            axes[6 + i * 3 + j] = axes_a[i].cross(axes_b[j]);
        }
    }

    const math::Vec3 delta = center_b - center_a;
    float min_overlap = std::numeric_limits<float>::infinity();
    math::Vec3 best_axis{0, 0, 0};

    for (const math::Vec3& raw_axis : axes) {
        math::Vec3 axis = raw_axis;
        const float length_squared = axis.length_squared();
        if (length_squared < 1.0e-12f) continue;
        axis = axis / std::sqrt(length_squared);

        const float radius_a = project_obb_radius(axis, axes_a, a.half_extents());
        const float radius_b = project_obb_radius(axis, axes_b, b.half_extents());
        const float center_distance = std::abs(delta.dot(axis));
        const float overlap = radius_a + radius_b - center_distance;
        if (overlap <= 0.0f) return false;

        if (overlap < min_overlap) {
            min_overlap = overlap;
            best_axis = axis;
        }
    }

    if (delta.dot(best_axis) < 0.0f) best_axis = -best_axis;
    contact.normal = best_axis;
    contact.penetration = min_overlap;
    // Use the support point on B along the contact normal. Averaging the two
    // centers (the old approximation) puts a ground contact halfway toward
    // the world origin in X/Z, creating a false lever arm and turning normal
    // support into a persistent rolling/sideways impulse for stacked boxes.
    // The support point remains on the actual contacting face and is stable
    // for static-vs-dynamic contacts.
    const float radius_b = project_obb_radius(best_axis, axes_b, b.half_extents());
    contact.point = center_b - best_axis * radius_b;
    contact.hit = true;
    return true;
}

} // namespace detail

inline bool Collider::test(const Collider& a, const math::Transform& ta,
                           const Collider& b, const math::Transform& tb,
                           ContactInfo& contact) {
    contact.reset();

    const Type type_a = a.type();
    const Type type_b = b.type();

    // An empty point set has no convex hull. Treat it as a non-collidable
    // shape rather than letting the GJK support fallback (the collider centre)
    // manufacture a false point-vs-shape contact.
    if ((type_a == Type::Convex &&
         static_cast<const ConvexCollider&>(a).points().empty()) ||
        (type_b == Type::Convex &&
         static_cast<const ConvexCollider&>(b).points().empty())) {
        return false;
    }

    if (type_a == Type::Sphere && type_b == Type::Sphere) {
        return detail::test_sphere_sphere(
            static_cast<const SphereCollider&>(a), ta,
            static_cast<const SphereCollider&>(b), tb, contact);
    }
    if (type_a == Type::Sphere && type_b == Type::Box) {
        return detail::test_sphere_box(
            static_cast<const SphereCollider&>(a), ta,
            static_cast<const BoxCollider&>(b), tb, contact);
    }
    if (type_a == Type::Box && type_b == Type::Sphere) {
        const bool hit = detail::test_sphere_box(
            static_cast<const SphereCollider&>(b), tb,
            static_cast<const BoxCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Sphere && type_b == Type::Capsule) {
        return detail::test_sphere_capsule(
            static_cast<const SphereCollider&>(a), ta,
            static_cast<const CapsuleCollider&>(b), tb, contact);
    }
    if (type_a == Type::Capsule && type_b == Type::Sphere) {
        const bool hit = detail::test_sphere_capsule(
            static_cast<const SphereCollider&>(b), tb,
            static_cast<const CapsuleCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Capsule && type_b == Type::Capsule) {
        return detail::test_capsule_capsule(
            static_cast<const CapsuleCollider&>(a), ta,
            static_cast<const CapsuleCollider&>(b), tb, contact);
    }
    if (type_a == Type::Capsule && type_b == Type::Box) {
        return detail::test_capsule_box(
            static_cast<const CapsuleCollider&>(a), ta,
            static_cast<const BoxCollider&>(b), tb, contact);
    }
    if (type_a == Type::Box && type_b == Type::Capsule) {
        const bool hit = detail::test_capsule_box(
            static_cast<const CapsuleCollider&>(b), tb,
            static_cast<const BoxCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Sphere && type_b == Type::Mesh) {
        return detail::test_sphere_mesh(
            static_cast<const SphereCollider&>(a), ta,
            static_cast<const MeshCollider&>(b), tb, contact);
    }
    if (type_a == Type::Mesh && type_b == Type::Sphere) {
        const bool hit = detail::test_sphere_mesh(
            static_cast<const SphereCollider&>(b), tb,
            static_cast<const MeshCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Capsule && type_b == Type::Mesh) {
        return detail::test_capsule_mesh(
            static_cast<const CapsuleCollider&>(a), ta,
            static_cast<const MeshCollider&>(b), tb, contact);
    }
    if (type_a == Type::Mesh && type_b == Type::Capsule) {
        const bool hit = detail::test_capsule_mesh(
            static_cast<const CapsuleCollider&>(b), tb,
            static_cast<const MeshCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Box && type_b == Type::Mesh) {
        return detail::test_box_mesh(
            static_cast<const BoxCollider&>(a), ta,
            static_cast<const MeshCollider&>(b), tb, contact);
    }
    if (type_a == Type::Mesh && type_b == Type::Box) {
        const bool hit = detail::test_box_mesh(
            static_cast<const BoxCollider&>(b), tb,
            static_cast<const MeshCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    if (type_a == Type::Mesh && type_b == Type::Mesh) {
        return detail::test_mesh_mesh(
            static_cast<const MeshCollider&>(a), ta,
            static_cast<const MeshCollider&>(b), tb, contact);
    }
    if (type_a == Type::Convex && type_b == Type::Mesh) {
        return detail::test_convex_mesh(
            a, ta, static_cast<const MeshCollider&>(b), tb, contact);
    }
    if (type_a == Type::Mesh && type_b == Type::Convex) {
        const bool hit = detail::test_convex_mesh(
            b, tb, static_cast<const MeshCollider&>(a), ta, contact);
        contact.normal = -contact.normal;
        return hit;
    }
    const bool convex_a = type_a == Type::Convex || type_a == Type::Box ||
                          type_a == Type::Sphere || type_a == Type::Capsule;
    const bool convex_b = type_b == Type::Convex || type_b == Type::Box ||
                          type_b == Type::Sphere || type_b == Type::Capsule;
    if (convex_a && convex_b &&
        (type_a == Type::Convex || type_b == Type::Convex)) {
        return detail::test_convex_pair(a, ta, b, tb, contact);
    }
    if (type_a == Type::Box && type_b == Type::Box) {
        return detail::test_box_box(
            static_cast<const BoxCollider&>(a), ta,
            static_cast<const BoxCollider&>(b), tb, contact);
    }

    // Fallback: AABB overlap check.
    const math::AABB aabb_a = a.compute_aabb(ta);
    const math::AABB aabb_b = b.compute_aabb(tb);
    if (!aabb_a.overlaps(aabb_b)) return false;
    contact.normal = math::Vec3::up();
    contact.penetration = 0.05f;
    contact.point = (aabb_a.center() + aabb_b.center()) * 0.5f;
    contact.hit = true;
    return true;
}

} // namespace butter
