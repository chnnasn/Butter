#pragma once
#include <algorithm>
#include <cmath>

namespace butter::physics2d {
// Mass-independent spring/damper coefficients, prepared once per substep.
// The impulse term is essential: merely modifying effective mass is not a
// compliant constraint. Frequency is limited relative to the substep rate.
struct ContactSoftness {
    float bias_rate{0}, mass_scale{1}, impulse_scale{0};
    static ContactSoftness make(float hertz, float damping_ratio, float dt) {
        if (!(hertz > 0) || !(dt > 0)) return {};
        const float omega = 6.28318530718f * std::min(hertz, .25f / dt);
        const float a = 2 * std::max(0.0f, damping_ratio) + dt * omega;
        const float b = dt * omega * a;
        return {omega / a, b / (1 + b), 1 / (1 + b)};
    }
};
}
