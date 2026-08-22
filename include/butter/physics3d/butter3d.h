#pragma once

// Explicit 3D facade. The original namespace butter API remains compatible.
#include "butter/butter.h"

namespace butter::physics3d {
using math::Vec3;
using math::Quat;
using math::Transform;
using math::AABB;
using Material = butter::Material;
using World = butter::World;
using Body = butter::Body;
using BodyBuilder = butter::BodyBuilder;
using BodyType = butter::BodyType;
using CollisionEvent = butter::CollisionEvent;
using TriggerEvent = butter::TriggerEvent;
using SphereCollider = butter::SphereCollider;
using BoxCollider = butter::BoxCollider;
using CapsuleCollider = butter::CapsuleCollider;
using ConvexCollider = butter::ConvexCollider;
using MeshCollider = butter::MeshCollider;
using ContactInfo = butter::ContactInfo;
}
