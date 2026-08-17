#pragma once

#include "butter/shapes/collider.h"
#include "butter/shapes/collision_test.h"

namespace butter {

// Convenience overlap query. Collider::test is defined by collision_test.h,
// which is included by butter.h.
inline bool overlap_test(const Collider& a, const math::Transform& ta,
                         const Collider& b, const math::Transform& tb) {
    ContactInfo contact;
    return Collider::test(a, ta, b, tb, contact);
}

} // namespace butter
