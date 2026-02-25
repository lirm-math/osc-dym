// real_t.h
#pragma once

#include <cfloat>
#include <limits>
#include <string>

#if defined(REAL_4)
using real_t = float;
static constexpr real_t pi = 3.14159265358979323846f;
inline real_t operator""_r(const char *str, std::size_t)
{
    return std::stof(str);
}
#elif defined(REAL_10)
using real_t = long double;
static constexpr real_t pi = 3.14159265358979323846L;
inline real_t operator""_r(const char *str, std::size_t)
{
    return std::stold(str);
}
#elif defined(REAL_16)
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp> 
using real_t = boost::multiprecision::cpp_bin_float_quad;
static inline const real_t pi = boost::math::constants::pi<real_t>();
inline real_t operator""_r(const char *str, std::size_t)
{
    return real_t(str);
}
#elif defined(REAL_8)
// 默认 double
using real_t = double;
static constexpr real_t pi = 3.14159265358979323846;
inline real_t operator""_r(const char *str, std::size_t)
{
    return std::stod(str);
}
#else
// 默认 double
using real_t = double;
static constexpr real_t pi = 3.14159265358979323846;
inline real_t operator""_r(const char *str, std::size_t)
{
    return std::stod(str);
}
#endif

inline real_t operator""_r(unsigned long long val)
{
    return real_t(val);
}

inline real_t operator""_r(long double val)
{
    return real_t(val);
}