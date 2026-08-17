#pragma once

#include <cmath>
#include <cstddef>
#include <tuple>

namespace butter::math {

struct Vec2 {
    float x{0};
    float y{0};

    constexpr Vec2() = default;
    constexpr Vec2(float x, float y) : x(x), y(y) {}
    constexpr explicit Vec2(float scalar) : x(scalar), y(scalar) {}

    constexpr float& operator[](std::size_t i) { return i == 0 ? x : y; }
    constexpr const float& operator[](std::size_t i) const { return i == 0 ? x : y; }

    constexpr Vec2 operator+(const Vec2& v) const { return {x + v.x, y + v.y}; }
    constexpr Vec2 operator-(const Vec2& v) const { return {x - v.x, y - v.y}; }
    constexpr Vec2 operator*(float s) const { return {x * s, y * s}; }
    constexpr Vec2 operator/(float s) const { return {x / s, y / s}; }
    constexpr Vec2 operator-() const { return {-x, -y}; }

    constexpr Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    constexpr Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    constexpr Vec2& operator*=(float s) { x *= s; y *= s; return *this; }
    constexpr Vec2& operator/=(float s) { x /= s; y /= s; return *this; }

    constexpr bool operator==(const Vec2& v) const { return x == v.x && y == v.y; }
    constexpr bool operator!=(const Vec2& v) const { return !(*this == v); }

    constexpr Vec2& set_x(float value) { x = value; return *this; }
    constexpr Vec2& set_y(float value) { y = value; return *this; }

    float length_squared() const { return x * x + y * y; }
    float length() const { return std::sqrt(length_squared()); }
    Vec2 normalized() const {
        const float len = length();
        return len > 1.0e-12f ? *this / len : Vec2{};
    }
    Vec2& normalize() { *this = normalized(); return *this; }

    constexpr float dot(const Vec2& v) const { return x * v.x + y * v.y; }
    constexpr float cross(const Vec2& v) const { return x * v.y - y * v.x; }

    static constexpr Vec2 zero() { return {}; }
    static constexpr Vec2 one() { return {1, 1}; }
    static constexpr Vec2 up() { return {0, 1}; }
    static constexpr Vec2 down() { return {0, -1}; }
    static constexpr Vec2 left() { return {-1, 0}; }
    static constexpr Vec2 right() { return {1, 0}; }

    friend constexpr Vec2 operator*(float s, const Vec2& v) { return v * s; }
};

using Point2 = Vec2;
using Direction2 = Vec2;

} // namespace butter::math


