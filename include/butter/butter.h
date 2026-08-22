#pragma once

// 🧈 Butter - Make Physics Smooth Again!
// One header to include them all.

#include "butter/math/vec2.h"
#include "butter/math/vec3.h"
#include "butter/math/mat3.h"
#include "butter/math/mat4.h"
#include "butter/math/quat.h"
#include "butter/math/transform.h"
#include "butter/math/aabb.h"

#include "butter/core/material.h"
#include "butter/core/events.h"

#include "butter/shapes/collider.h"
#include "butter/shapes/sphere.h"
#include "butter/shapes/box.h"
#include "butter/shapes/capsule.h"
#include "butter/shapes/convex.h"
#include "butter/shapes/mesh.h"
#include "butter/shapes/collision_test.h"

#include "butter/core/body.h"
#include "butter/core/body_builder.h"

#include "butter/constraints/joint.h"
#include "butter/constraints/distance.h"
#include "butter/constraints/spring.h"
#include "butter/constraints/hinge.h"

#include "butter/query/raycast.h"
#include "butter/query/overlap.h"

#include "butter/core/world.h"

// Optional 2D module (kept in the physics2d namespace).
#include "butter/physics2d/butter2d.h"
