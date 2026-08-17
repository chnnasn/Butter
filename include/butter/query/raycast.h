#pragma once

#include <limits>
#include <optional>
#include <vector>

#include "butter/math/vec3.h"

namespace butter {

class Body;
class World;

struct RaycastHit {
    Body* body{nullptr};
    math::Vec3 point{};
    math::Vec3 normal{};
    float distance{std::numeric_limits<float>::infinity()};
    bool hit{false};

    explicit operator bool() const { return hit; }
};

class QueryBuilder {
public:
    explicit QueryBuilder(World& world) : world_(world) {}

    QueryBuilder& from(const math::Vec3& origin) {
        origin_ = origin;
        return *this;
    }
    QueryBuilder& from(float x, float y, float z) {
        origin_ = {x, y, z};
        return *this;
    }

    QueryBuilder& toward(const math::Vec3& direction) {
        direction_ = direction.normalized();
        return *this;
    }
    QueryBuilder& toward(float x, float y, float z) {
        direction_ = math::Vec3{x, y, z}.normalized();
        return *this;
    }

    QueryBuilder& max_distance(float distance) {
        max_distance_ = distance;
        return *this;
    }

    QueryBuilder& ignore(Body& body) {
        ignore_list_.push_back(&body);
        return *this;
    }

    QueryBuilder& only_dynamic() {
        dynamic_only_ = true;
        return *this;
    }

    // Defined in world.h after World is complete.
    std::optional<RaycastHit> first();
    std::vector<RaycastHit> all();

private:
    World& world_;
    math::Vec3 origin_{0, 0, 0};
    math::Vec3 direction_{0, 0, -1};
    float max_distance_{std::numeric_limits<float>::max()};
    std::vector<Body*> ignore_list_;
    bool dynamic_only_{false};

    friend class World;
};

} // namespace butter
