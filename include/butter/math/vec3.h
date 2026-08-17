#pragma once

#include <cmath>
#include <cstddef>
#include <tuple>

namespace butter::math {

constexpr float pi = 3.14159265358979323846f;
constexpr float tau = pi * 2.0f;
constexpr float deg_to_rad(float deg) { return deg * (pi / 180.0f); }
constexpr float rad_to_deg(float rad) { return rad * (180.0f / pi); }

constexpr float operator""_deg(long double value) {
    return static_cast<float>(value) * (pi / 180.0f);
}
constexpr float operator""_deg(unsigned long long value) {
    return static_cast<float>(value) * (pi / 180.0f);
}

struct Vec3 {
    float x{0};
    float y{0};
    float z{0};

    constexpr Vec3() = default;
    constexpr Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    constexpr explicit Vec3(float scalar) : x(scalar), y(scalar), z(scalar) {}

    constexpr float& operator[](std::size_t i) {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    constexpr const float& operator[](std::size_t i) const {
        return i == 0 ? x : (i == 1 ? y : z);
    }

    constexpr Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    constexpr Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    constexpr Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(float s) const { return {x / s, y / s, z / s}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }

    constexpr Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    constexpr Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
    constexpr Vec3& operator/=(float s) { x /= s; y /= s; z /= s; return *this; }

    constexpr bool operator==(const Vec3& v) const { return x == v.x && y == v.y && z == v.z; }
    constexpr bool operator!=(const Vec3& v) const { return !(*this == v); }

    constexpr Vec3& set_x(float value) { x = value; return *this; }
    constexpr Vec3& set_y(float value) { y = value; return *this; }
    constexpr Vec3& set_z(float value) { z = value; return *this; }

    float length_squared() const { return x * x + y * y + z * z; }
    float length() const { return std::sqrt(length_squared()); }
    Vec3 normalized() const {
        const float len = length();
        return len > 1.0e-12f ? *this / len : Vec3{};
    }
    Vec3& normalize() { *this = normalized(); return *this; }

    constexpr float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    constexpr Vec3 cross(const Vec3& v) const {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }

    constexpr Vec3 component_min(const Vec3& v) const {
        return {x < v.x ? x : v.x, y < v.y ? y : v.y, z < v.z ? z : v.z};
    }
    constexpr Vec3 component_max(const Vec3& v) const {
        return {x > v.x ? x : v.x, y > v.y ? y : v.y, z > v.z ? z : v.z};
    }

    static constexpr Vec3 zero() { return {}; }
    static constexpr Vec3 one() { return {1, 1, 1}; }
    static constexpr Vec3 up() { return {0, 1, 0}; }
    static constexpr Vec3 down() { return {0, -1, 0}; }
    static constexpr Vec3 forward() { return {0, 0, 1}; }
    static constexpr Vec3 back() { return {0, 0, -1}; }
    static constexpr Vec3 left() { return {-1, 0, 0}; }
    static constexpr Vec3 right() { return {1, 0, 0}; }

    friend constexpr Vec3 operator*(float s, const Vec3& v) { return v * s; }
};

using Point3 = Vec3;
using Direction3 = Vec3;
using Force3 = Vec3;
using Velocity3 = Vec3;

} // namespace butter::math


