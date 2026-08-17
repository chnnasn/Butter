#pragma once

namespace butter {

struct Material {
    float friction{0.5f};
    float restitution{0.0f};
    float density{1.0f};
    bool is_trigger{false};

    constexpr Material() = default;
    constexpr Material(float friction, float restitution, float density = 1.0f)
        : friction(friction), restitution(restitution), density(density) {}

    Material& set_friction(float value) { friction = value; return *this; }
    Material& set_restitution(float value) { restitution = value; return *this; }
    Material& set_density(float value) { density = value; return *this; }
};

} // namespace butter
