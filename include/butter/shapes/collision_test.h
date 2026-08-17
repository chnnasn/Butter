#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "butter/shapes/box.h"
#include "butter/shapes/capsule.h"
#include "butter/shapes/collider.h"
#include "butter/shapes/convex.h"
#include "butter/shapes/sphere.h"

namespace butter {
namespace detail {

inline math::Vec3 collider_center(const Collider& collider, const math::Transform& transform) {
    return transform.position + collider.offset.position;
}

inline math::Quat collider_rotation(const Collider& collider, const math::Transform& transform) {
    return (transform.rotation * collider.offset.rotation).normalized();
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

    // A simple iterative closest-point approximation is fine for the reference solver.
    math::Vec3 pa = (a0 + a1) * 0.5f;
    math::Vec3 pb = (b0 + b1) * 0.5f;
    for (int i = 0; i < 8; ++i) {
        pa = closest_point_on_segment(pb, a0, a1);
        pb = closest_point_on_segment(pa, b0, b1);
    }

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
    contact.point = (center_a + center_b) * 0.5f;
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
