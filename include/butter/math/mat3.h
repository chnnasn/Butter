#pragma once

#include "vec3.h"

namespace butter::math {

struct Mat3 {
    // Row-major storage. The default constructor produces the identity matrix.
    float m[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f},
    };

    constexpr Mat3() = default;
    constexpr Mat3(float m00, float m01, float m02,
                   float m10, float m11, float m12,
                   float m20, float m21, float m22)
        : m{{m00, m01, m02}, {m10, m11, m12}, {m20, m21, m22}} {}

    constexpr Vec3 operator*(const Vec3& v) const {
        return {
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
        };
    }

    constexpr Mat3 operator*(const Mat3& other) const {
        Mat3 result;
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                result.m[row][col] = m[row][0] * other.m[0][col] +
                                     m[row][1] * other.m[1][col] +
                                     m[row][2] * other.m[2][col];
            }
        }
        return result;
    }

    constexpr Mat3 transposed() const {
        return {
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2],
        };
    }

    static constexpr Mat3 zero() {
        return {
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        };
    }

    static constexpr Mat3 identity() {
        return {};
    }

    static constexpr Mat3 diagonal(float a, float b, float c) {
        return {
            a, 0, 0,
            0, b, 0,
            0, 0, c,
        };
    }

    static constexpr Mat3 scale(float s) {
        return diagonal(s, s, s);
    }
};

} // namespace butter::math
