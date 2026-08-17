#pragma once

#include "vec3.h"

namespace butter::math {

struct Mat4 {
    // Column-major storage, matching OpenGL conventions. Default is identity.
    float m[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
    };

    constexpr Mat4() = default;

    static constexpr Mat4 identity() { return {}; }

    static constexpr Mat4 translation(const Vec3& t) {
        Mat4 result;
        result.m[3][0] = t.x;
        result.m[3][1] = t.y;
        result.m[3][2] = t.z;
        return result;
    }

    static constexpr Mat4 scale(const Vec3& s) {
        Mat4 result;
        result.m[0][0] = s.x;
        result.m[1][1] = s.y;
        result.m[2][2] = s.z;
        return result;
    }

    constexpr Mat4 operator*(const Mat4& other) const {
        Mat4 result;
        for (int col = 0; col < 4; ++col) {
            for (int row = 0; row < 4; ++row) {
                result.m[col][row] = m[0][row] * other.m[col][0] +
                                     m[1][row] * other.m[col][1] +
                                     m[2][row] * other.m[col][2] +
                                     m[3][row] * other.m[col][3];
            }
        }
        return result;
    }
};

} // namespace butter::math
