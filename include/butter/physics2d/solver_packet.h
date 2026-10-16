#pragma once
#include <algorithm>
#include <array>
#if defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
#endif

namespace butter::physics2d::solver_detail {
using Lanes = std::array<float, 4>;
// Four independent constraints. Callers must guarantee that no writable body
// is shared by lanes. Unsupported architectures retain the same scalar math.
inline Lanes normal_packet(const Lanes &previous, const Lanes &target, const Lanes &speed,
                           const Lanes &bias, const Lanes &mass, float scale, float impulse_scale) {
    Lanes result;
#if defined(__SSE2__) || defined(_M_X64)
    auto delta = _mm_div_ps(_mm_mul_ps(_mm_set1_ps(scale),
        _mm_sub_ps(_mm_sub_ps(_mm_loadu_ps(target.data()), _mm_loadu_ps(speed.data())),
                   _mm_loadu_ps(bias.data()))), _mm_loadu_ps(mass.data()));
    auto old = _mm_loadu_ps(previous.data());
    delta = _mm_sub_ps(delta, _mm_mul_ps(_mm_set1_ps(impulse_scale), old));
    _mm_storeu_ps(result.data(), _mm_max_ps(_mm_setzero_ps(), _mm_add_ps(old, delta)));
#else
    for (int i=0;i<4;++i) {
        const float delta = scale * (target[i]-speed[i]-bias[i])/mass[i] - impulse_scale*previous[i];
        result[i] = std::max(0.0f, previous[i] + delta);
    }
#endif
    return result;
}
inline Lanes friction_packet(const Lanes &previous, const Lanes &speed, const Lanes &mass,
                             const Lanes &limit) {
    Lanes result;
#if defined(__SSE2__) || defined(_M_X64)
    auto value = _mm_sub_ps(_mm_loadu_ps(previous.data()),
                           _mm_div_ps(_mm_loadu_ps(speed.data()), _mm_loadu_ps(mass.data())));
    auto upper = _mm_loadu_ps(limit.data());
    _mm_storeu_ps(result.data(), _mm_min_ps(upper, _mm_max_ps(_mm_sub_ps(_mm_setzero_ps(),upper),value)));
#else
    for (int i=0;i<4;++i) result[i]=std::clamp(previous[i]-speed[i]/mass[i],-limit[i],limit[i]);
#endif
    return result;
}
}
