#pragma once

#include <algorithm>
#include <cmath>

#include "mat3.h"
#include "vec3.h"

namespace butter::math {

struct Quat {
    float w{1};
    float x{0};
    float y{0};
    float z{0};

    constexpr Quat() = default;
    constexpr Quat(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}

    static constexpr Quat identity() { return {}; }

    static Quat from_euler(float x_rad, float y_rad, float z_rad) {
        const float cx = std::cos(x_rad * 0.5f);
        const float sx = std::sin(x_rad * 0.5f);
        const float cy = std::cos(y_rad * 0.5f);
        const float sy = std::sin(y_rad * 0.5f);
        const float cz = std::cos(z_rad * 0.5f);
        const float sz = std::sin(z_rad * 0.5f);

        return {
            cx * cy * cz + sx * sy * sz,
            sx * cy * cz - cx * sy * sz,
            cx * sy * cz + sx * cy * sz,
            cx * cy * sz - sx * sy * cz,
        };
    }

    static Quat from_axis_angle(const Vec3& axis, float angle_rad) {
        const Vec3 n = axis.normalized();
        const float half = angle_rad * 0.5f;
        const float s = std::sin(half);
        return {std::cos(half), n.x * s, n.y * s, n.z * s};
    }

    constexpr Quat operator*(const Quat& q) const {
        return {
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w,
        };
    }

    constexpr Quat conjugate() const { return {w, -x, -y, -z}; }

    float length_squared() const { return w * w + x * x + y * y + z * z; }
    float length() const { return std::sqrt(length_squared()); }

    Quat normalized() const {
        const float len = length();
        return len > 1.0e-12f ? Quat{w / len, x / len, y / len, z / len} : Quat{};
    }

    Quat& normalize() { *this = normalized(); return *this; }

    Quat inverse() const {
        const float len2 = length_squared();
        const Quat c = conjugate();
        if (len2 < 1.0e-12f) return Quat{};
        return {c.w / len2, c.x / len2, c.y / len2, c.z / len2};
    }

    Vec3 rotate(const Vec3& v) const {
        return to_mat3() * v;
    }

    Mat3 to_mat3() const {
        const float xx = x * x;
        const float yy = y * y;
        const float zz = z * z;
        const float xy = x * y;
        const float xz = x * z;
        const float yz = y * z;
        const float wx = w * x;
        const float wy = w * y;
        const float wz = w * z;

        return {
            1.0f - 2.0f * (yy + zz), 2.0f * (xy - wz),       2.0f * (xz + wy),
            2.0f * (xy + wz),       1.0f - 2.0f * (xx + zz), 2.0f * (yz - wx),
            2.0f * (xz - wy),       2.0f * (yz + wx),       1.0f - 2.0f * (xx + yy),
        };
    }

    static Quat slerp(const Quat& a, const Quat& b, float t) {
        float dot = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
        Quat qb = b;
        if (dot < 0.0f) {
            dot = -dot;
            qb = {-b.w, -b.x, -b.y, -b.z};
        }
        dot = std::clamp(dot, -1.0f, 1.0f);
        const float theta = std::acos(dot);
        const float sin_theta = std::sin(theta);
        if (sin_theta < 1.0e-6f) {
            const Quat lerped{
                a.w + (qb.w - a.w) * t,
                a.x + (qb.x - a.x) * t,
                a.y + (qb.y - a.y) * t,
                a.z + (qb.z - a.z) * t,
            };
            return lerped.normalized();
        }
        const float wa = std::sin((1.0f - t) * theta) / sin_theta;
        const float wb = std::sin(t * theta) / sin_theta;
        return {
            a.w * wa + qb.w * wb,
            a.x * wa + qb.x * wb,
            a.y * wa + qb.y * wb,
            a.z * wa + qb.z * wb,
        };
    }
};

} // namespace butter::math
