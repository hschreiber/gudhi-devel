/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_UINT128_H_
#define GUDHI_UINT128_H_
#include <gudhi/Debug_utils.h>
#include <cstdint>
#include <utility>
#include <boost/container_hash/hash.hpp>

/* What about boost::multiprecision::uint128_t?
 * It works on linux, mac, windows, etc.
 * On linux x86_64, it has low overhead since it delegates to unsigned __int128.
 * It is a bit slower than Fake_uint128, probably because of checks in << and >>, but I only noticed a 10%
 * overhead compared to unsigned __int128 on 1 testcase.
 * On windows, it is implemented as { size_t; uint64_t[2]; }, that's 50% overhead on storage. On the same
 * testcase, forcing linux-gcc to use the windows code led to 4x slowdown. Other testcases were not affected.
 * The pathological testcase is computing homology in all dimensions for the vertices of a regular 24-gon.
 */

// GUDHI_FORCE_FAKE_UINT128 is only used for tests
#if !defined __SIZEOF_INT128__ || defined GUDHI_FORCE_FAKE_UINT128
namespace Gudhi::numbers {
class Fake_uint128 {
  // Debug
#ifdef __SIZEOF_INT128__
  unsigned __int128 native() const { return ((unsigned __int128)high << 64) + low; }
#define GUDHI_VERIF(X) GUDHI_CHECK(res.native() == (X), "")
#else
#define GUDHI_VERIF(X)
#endif
 public:
  constexpr Fake_uint128() : high(0), low(0) {}
  constexpr Fake_uint128(std::uint64_t a) : high(0), low(a) {}
  // Arithmetic
  // (multiplication and division are not needed for now)
  friend Fake_uint128 operator+(Fake_uint128 a, Fake_uint128 b) {
    Fake_uint128 res;
    res.low = a.low + b.low;
    res.high = a.high + b.high + (res.low < a.low);
    GUDHI_VERIF(a.native() + b.native());
    return res;
  }
  friend Fake_uint128 operator-(Fake_uint128 a, Fake_uint128 b) {
    Fake_uint128 res;
    res.low = a.low - b.low;
    res.high = a.high - b.high - (res.low > a.low);
    GUDHI_VERIF(a.native() - b.native());
    return res;
  }
  friend Fake_uint128 operator<<(Fake_uint128 a, uint8_t b) {
    Fake_uint128 res;
    GUDHI_CHECK(b < 128, "");
    if (b >= 64) {
      res.low = 0;
      res.high = a.low << (b - 64);
    } else if (b == 0) {
      res = a;
    } else {
      res.low = a.low << b;
      res.high = a.high << b | a.low >> (64 - b);
    }
    GUDHI_VERIF(a.native() << b);
    return res;
  }
  friend Fake_uint128 operator>>(Fake_uint128 a, uint8_t b) {
    Fake_uint128 res;
    GUDHI_CHECK(b < 128, "");
    if (b >= 64) {
      res.high = 0;
      res.low = a.high >> (b - 64);
    } else if (b == 0) {
      res = a;
    } else {
      res.high = a.high >> b;
      res.low = a.low >> b | a.high << (64 - b);
    }
    GUDHI_VERIF(a.native() >> b);
    return res;
  }
  friend Fake_uint128 operator&(Fake_uint128 a, Fake_uint128 b) {
    Fake_uint128 res;
    res.low = a.low & b.low;
    res.high = a.high & b.high;
    GUDHI_VERIF(a.native() & b.native());
    return res;
  }
  friend Fake_uint128 operator|(Fake_uint128 a, Fake_uint128 b) {
    Fake_uint128 res;
    res.low = a.low | b.low;
    res.high = a.high | b.high;
    GUDHI_VERIF(a.native() | b.native());
    return res;
  }
  friend Fake_uint128 operator~(Fake_uint128 a) {
    Fake_uint128 res;
    res.low = ~a.low;
    res.high = ~a.high;
    return res;
  }
  // In-place arithmetic
  Fake_uint128& operator+=(Fake_uint128 a) {
    *this = *this + a;
    return *this;
  }
  Fake_uint128& operator-=(Fake_uint128 a) {
    *this = *this - a;
    return *this;
  }
  Fake_uint128& operator++() {
    if (++low == 0) ++high;
    return *this;
  }
  Fake_uint128& operator--() {
    if (low-- == 0) --high;
    return *this;
  }
  Fake_uint128& operator<<=(uint8_t a) {
    *this = *this << a;
    return *this;
  }
  Fake_uint128& operator>>=(uint8_t a) {
    *this = *this >> a;
    return *this;
  }
  Fake_uint128& operator&=(Fake_uint128 a) {
    *this = *this & a;
    return *this;
  }
  Fake_uint128& operator|=(Fake_uint128 a) {
    *this = *this | a;
    return *this;
  }
  // Comparisons
  friend bool operator==(Fake_uint128 a, Fake_uint128 b) { return a.low == b.low && a.high == b.high; }
  friend bool operator!=(Fake_uint128 a, Fake_uint128 b) { return a.low != b.low || a.high != b.high; }
  friend bool operator<(Fake_uint128 a, Fake_uint128 b) {
    return a.high < b.high || (a.high == b.high && a.low < b.low);
  }
  friend bool operator>(Fake_uint128 a, Fake_uint128 b) {
    return a.high > b.high || (a.high == b.high && a.low > b.low);
  }
  friend bool operator<=(Fake_uint128 a, Fake_uint128 b) {
    return a.high < b.high || (a.high == b.high && a.low <= b.low);
  }
  friend bool operator>=(Fake_uint128 a, Fake_uint128 b) {
    return a.high > b.high || (a.high == b.high && a.low >= b.low);
  }
  // Misc
  friend std::size_t hash_value(Fake_uint128 a) {
    typedef std::pair<std::uint64_t, std::uint64_t> P;
    return boost::hash_value(P(a.high, a.low));
  }
  template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
  explicit operator T() const {
    GUDHI_CHECK(high == 0 && low <= std::numeric_limits<T>::max(), "");
    return static_cast<T>(low);
  }

 private:
  std::uint64_t high, low;  // does the order matter?
#undef GUDHI_VERIF
};
typedef Fake_uint128 uint128_t;
}  // namespace Gudhi::numbers
template <>
class std::numeric_limits<Gudhi::numbers::Fake_uint128> {
 public:
  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr bool has_infinity = false;
  static constexpr bool is_modulo = true;
  static constexpr int digits = 128;
  static constexpr int radix = 2;
  // etc
};
#else
namespace Gudhi::numbers {
__extension__ using uint128_t = unsigned __int128;
}  // namespace Gudhi::numbers
#endif

// ------------------------------------------------------------------------------------------------------------
// Wide arithmetic on 64 bits words.
//
// Fake_uint128 has no multiplication nor division, and emulating them through its operators would be slow. The few
// primitives below are what is needed to do modular arithmetic with 64 bits moduli (see gudhi/arithmetic.h).
// They only use std::uint64_t in their interface, with a native implementation (unsigned __int128) when available,
// and a portable one otherwise (or when GUDHI_FORCE_FAKE_UINT128 is defined, which is used to test it).
// Compiler intrinsics (_umul128, __umulh, _udiv128 for MSVC x64) could be plugged in the portable branch.
namespace Gudhi::numbers {

#if !defined __SIZEOF_INT128__ || defined GUDHI_FORCE_FAKE_UINT128

/** @brief Full product of two 64 bits integers: returns the high 64 bits and stores the low 64 bits in @p low. */
constexpr std::uint64_t mul_wide(std::uint64_t a, std::uint64_t b, std::uint64_t& low) noexcept {
  const std::uint64_t mask = 0xFFFFFFFFu;
  const std::uint64_t a0 = a & mask, a1 = a >> 32, b0 = b & mask, b1 = b >> 32;
  const std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const std::uint64_t mid = (p00 >> 32) + (p01 & mask) + (p10 & mask);  // does not overflow
  low = (mid << 32) | (p00 & mask);
  return p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
}

/**
 * @brief Divides the 128 bits integer `hi * 2^64 + lo` by @p d. Stores the remainder in @p rem and returns the
 * quotient. Requires `hi < d` (so that the quotient fits in 64 bits). Meant for set-up code: shift-subtract.
 */
constexpr std::uint64_t div_wide(std::uint64_t hi, std::uint64_t lo, std::uint64_t d, std::uint64_t& rem) noexcept {
  std::uint64_t q = 0, r = hi;
  for (int i = 63; i >= 0; --i) {
    const bool carry = (r >> 63) != 0;  // r < d, so r * 2 + bit < 2 * d may overflow 64 bits
    r = (r << 1) | ((lo >> i) & 1u);
    q <<= 1;
    if (carry || r >= d) {
      r -= d;
      q |= 1u;
    }
  }
  rem = r;
  return q;
}

#else

constexpr std::uint64_t mul_wide(std::uint64_t a, std::uint64_t b, std::uint64_t& low) noexcept {
  const uint128_t p = static_cast<uint128_t>(a) * b;
  low = static_cast<std::uint64_t>(p);
  return static_cast<std::uint64_t>(p >> 64);
}

constexpr std::uint64_t div_wide(std::uint64_t hi, std::uint64_t lo, std::uint64_t d, std::uint64_t& rem) noexcept {
  const uint128_t n = (static_cast<uint128_t>(hi) << 64) | lo;
  const uint128_t q = n / d;
  rem = static_cast<std::uint64_t>(n - q * d);
  return static_cast<std::uint64_t>(q);
}

#endif

/** @brief High 64 bits of the product of two 64 bits integers. */
constexpr std::uint64_t mul_hi(std::uint64_t a, std::uint64_t b) noexcept {
  std::uint64_t low = 0;
  return mul_wide(a, b, low);
}

}  // namespace Gudhi::numbers

#endif  // GUDHI_UINT128_H_