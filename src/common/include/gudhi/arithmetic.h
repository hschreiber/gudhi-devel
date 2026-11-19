/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 *      - 2026/11 Hannah Schreiber: modular arithmetic helpers (reducers, extended Euclid, mul_mod, pow_mod) shared
 *        by the persistence fields, and exact primality test for large integers.
 */

#ifndef GUDHI_ARITHMETIC_H_
#define GUDHI_ARITHMETIC_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <gudhi/uint128.h>

namespace Gudhi {

namespace numbers {

/**
 * @brief Returns `(a * b) % m` for any `m > 0` and any `a`, `b`. Not meant to be used in a hot loop: see
 * @ref Modular_reducer for that.
 */
constexpr std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b, std::uint64_t m) noexcept {
  a %= m;
  b %= m;
  if (m <= 0xFFFFFFFFu) return a * b % m;  // a * b < 2^64
  std::uint64_t low = 0, rem = 0;
  const std::uint64_t high = mul_wide(a, b, low);  // high < m as a * b < m^2
  div_wide(high, low, m, rem);
  return rem;
}

/** @brief Returns `(base ^ exp) % m` for any `m > 0`. */
constexpr std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exp, std::uint64_t m) noexcept {
  std::uint64_t res = 1 % m;
  base %= m;
  while (exp > 0) {
    if ((exp & 1U) != 0U) res = mul_mod(res, base, m);
    base = mul_mod(base, base, m);
    exp >>= 1;
  }
  return res;
}

/**
 * @brief Extended Euclid algorithm on `(a, m)`, with `0 <= a < m`.
 *
 * Stores `g = gcd(a, m)` in @p g (`g = m` if `a == 0`) and returns the Bezout coefficient `t` in `[0, m)` such that
 * `t * a = g` modulo `m`. In particular, if `g == 1`, `t` is the inverse of `a` modulo `m`; and if `a == 0`,
 * `t = 0`.
 *
 * @tparam T Unsigned integer type of at most 64 bits.
 * @tparam Signed Signed integer type used for the coefficients, which have to fit in it: `m` has to be at most its
 * maximal value. The default always works for 32 bits `T` and for 64 bits `T` if `m < 2^63`; for smaller `m`,
 * `std::make_signed_t<T>` is faster (`m < 2^31` or `m < 2^63`).
 */
template <typename T, typename Signed = std::int64_t>
constexpr T extended_euclid(T a, T m, T& g) noexcept {
  static_assert(std::is_unsigned_v<T> && sizeof(T) <= 8, "Unsigned integer types of at most 64 bits only.");
  static_assert(std::is_signed_v<Signed>, "A signed integer type is needed for the coefficients.");
  Signed t0 = 0, t1 = 1;
  T r0 = m, r1 = a;
  while (r1 != 0) {
    const T q = r0 / r1;
    const T r2 = r0 - (q * r1);
    r0 = r1;
    r1 = r2;
    const auto t2 = static_cast<Signed>(t0 - (static_cast<Signed>(q) * t1));  // |t| never exceeds m
    t0 = t1;
    t1 = t2;
  }
  g = r0;
  return static_cast<T>(t0 < 0 ? t0 + static_cast<Signed>(m) : t0);
}

/** @brief Returns the inverse of `a` modulo `m`, if it exists (`gcd(a, m) == 1`), see @ref extended_euclid. */
template <typename T, typename Signed = std::int64_t>
constexpr T mod_inverse(T a, T m) noexcept {
  T g = 0;
  return extended_euclid<T, Signed>(a, m, g);
}

/**
 * @brief Fast `x % m` and modular multiplications for a fixed modulus `2 <= m < 2^32`, without hardware division
 * (Barrett reduction: one multiplication to estimate the quotient, one conditional correction).
 *
 * Also provides multiplication by a fixed element with a precomputed quotient (Shoup), which needs `m < 2^31`.
 */
class Modular_reducer_32 {
 public:
  using Element = std::uint32_t;

  constexpr Modular_reducer_32() noexcept : m_(0), inv_(0) {}
  /** @param m Modulus, `m >= 2`. */
  constexpr explicit Modular_reducer_32(Element m) noexcept : m_(m), inv_(UINT64_MAX / m) {}

  [[nodiscard]] constexpr Element modulus() const noexcept { return m_; }

  /** @brief `x % m` for any 64 bits integer `x`. */
  [[nodiscard]] constexpr Element reduce(std::uint64_t x) const noexcept {
    const std::uint64_t q = mul_hi(x, inv_);                          // floor(x / m) or one less
    const std::uint64_t r = x - (q * m_);                             // in [0, 2m)
    return static_cast<Element>(std::min<std::uint64_t>(r, r - m_));  // r - m wraps around if r < m
  }
  /** @brief `(a * b) % m`, with `a, b < m`. */
  [[nodiscard]] constexpr Element mul(Element a, Element b) const noexcept { return reduce(std::uint64_t(a) * b); }
  /** @brief `(a * b + c) % m`, with `a, b, c < m`. */
  [[nodiscard]] constexpr Element mul_add(Element a, Element b, Element c) const noexcept {
    return reduce((std::uint64_t(a) * b) + c);
  }

  /** @brief Pre-computation for @ref mul_shoup, `w < m`. */
  [[nodiscard]] constexpr Element shoup_quotient(Element w) const noexcept {
    return static_cast<Element>((std::uint64_t(w) << 32) / m_);
  }
  /** @brief `(x * w) % m` for any `x`, with `w < m` and `wq = shoup_quotient(w)`. Requires `m < 2^31`. */
  [[nodiscard]] constexpr Element mul_shoup(Element x, Element w, Element wq) const noexcept {
    const auto q = static_cast<Element>((std::uint64_t(x) * wq) >> 32);
    const auto r = ((x * w) - (q * m_));  // exact value in [0, 2m)
    return std::min<Element>(r, r - m_);
  }

 private:
  Element m_;
  std::uint64_t inv_;  // floor((2^64 - 1) / m)
};

/**
 * @brief Same as @ref Modular_reducer_32 for a modulus `2 <= m < 2^64`. Division by an invariant integer
 * (Möller & Granlund, "Improved division by invariant integers", 2011, Algorithm 4): a reduction costs two 64x64
 * multiplications and no division, whatever the compiler support for 128 bits integers is.
 */
class Modular_reducer_64 {
 public:
  using Element = std::uint64_t;

  constexpr Modular_reducer_64() noexcept : m_(0), shift_(0), d_(0), v_(0) {}
  /** @param m Modulus, `m >= 2`. */
  constexpr explicit Modular_reducer_64(Element m) noexcept : m_(m), shift_(0), d_(m), v_(0) {
    while ((d_ >> 63) == 0) {  // normalize: top bit of d_ is set
      d_ <<= 1;
      ++shift_;
    }
    std::uint64_t rem = 0;
    v_ = div_wide(~d_, ~std::uint64_t(0), d_, rem);  // floor((2^128 - 1) / d) - 2^64
  }

  [[nodiscard]] constexpr Element modulus() const noexcept { return m_; }

  /** @brief `(hi * 2^64 + lo) % m`, with `hi < m`. */
  [[nodiscard]] constexpr Element reduce_wide(std::uint64_t hi, std::uint64_t lo) const noexcept {
    // shifted dividend (u1, u0): u1 < d_ because hi < m
    const std::uint64_t u1 = (hi << shift_) | ((lo >> 1) >> (63 - shift_));
    const std::uint64_t u0 = lo << shift_;
    std::uint64_t q0 = 0;
    std::uint64_t q1 = mul_wide(v_, u1, q0);  // (q1, q0) = v * u1
    q0 += u0;                                 // (q1, q0) += (u1, u0)
    q1 += u1 + static_cast<std::uint64_t>(q0 < u0);
    ++q1;
    std::uint64_t r = u0 - (q1 * d_);  // all operations modulo 2^64 from here
    if (r > q0) r += d_;
    if (r >= d_) r -= d_;
    return r >> shift_;
  }
  /** @brief `x % m` for any 64 bits integer `x`. */
  [[nodiscard]] constexpr Element reduce(std::uint64_t x) const noexcept { return reduce_wide(0, x); }
  /** @brief `(a * b) % m`, with `a, b < m`. */
  [[nodiscard]] constexpr Element mul(Element a, Element b) const noexcept {
    std::uint64_t lo = 0;
    const std::uint64_t hi = mul_wide(a, b, lo);
    return reduce_wide(hi, lo);
  }
  /** @brief `(a * b + c) % m`, with `a, b, c < m`. */
  [[nodiscard]] constexpr Element mul_add(Element a, Element b, Element c) const noexcept {
    std::uint64_t lo = 0;
    std::uint64_t hi = mul_wide(a, b, lo);
    lo += c;
    hi += static_cast<std::uint64_t>(lo < c);  // a * b + c <= m * (m - 1) < m * 2^64, so hi < m still holds
    return reduce_wide(hi, lo);
  }

  /** @brief Pre-computation for @ref mul_shoup, `w < m`. */
  [[nodiscard]] constexpr Element shoup_quotient(Element w) const noexcept {
    std::uint64_t rem = 0;
    return div_wide(w, 0, m_, rem);  // floor(w * 2^64 / m)
  }
  /** @brief `(x * w) % m` for any `x`, with `w < m` and `wq = shoup_quotient(w)`. Requires `m < 2^63`. */
  [[nodiscard]] constexpr Element mul_shoup(Element x, Element w, Element wq) const noexcept {
    const Element q = mul_hi(x, wq);
    const Element r = (x * w) - (q * m_);  // exact value in [0, 2m)
    return std::min<Element>(r, r - m_);
  }

 private:
  Element m_;
  unsigned int shift_;
  std::uint64_t d_;  // m << shift_
  std::uint64_t v_;  // floor((2^128 - 1) / d_) - 2^64
};

/** @brief The reducer fitting the given unsigned integer type (32 or 64 bits). */
template <typename Unsigned_integer_type>
using Modular_reducer = std::conditional_t<sizeof(Unsigned_integer_type) == 4, Modular_reducer_32, Modular_reducer_64>;

/**
 * @brief Deterministic Miller-Rabin primality test, exact for every 64 bits integer. (The 12 first primes as bases
 * are enough up to 3.3 * 10^24; 2, 7 and 61 are enough up to 4 759 123 141.)
 */
constexpr bool is_prime_miller_rabin(std::uint64_t n) noexcept {
  constexpr std::array<std::uint64_t, 12> smallPrimes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  if (n < 2) return false;
  for (std::uint64_t q : smallPrimes) {
    if (n == q) return true;
    if (n % q == 0) return false;
  }
  std::uint64_t d = n - 1;
  unsigned int s = 0;
  while ((d & 1U) == 0) {
    d >>= 1;
    ++s;
  }
  constexpr std::array<std::uint64_t, 3> basesSmall = {2, 7, 61};
  const bool small = n < 4759123141ULL;
  const std::size_t nbBases = small ? 3 : 12;
  for (std::size_t i = 0; i < nbBases; ++i) {
    const std::uint64_t a = small ? basesSmall[i] : smallPrimes[i];
    std::uint64_t x = pow_mod(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool composite = true;
    for (unsigned int r = 1; r < s && composite; ++r) {
      x = mul_mod(x, x, n);
      if (x == n - 1) composite = false;
    }
    if (composite) return false;
  }
  return true;
}

}  // namespace numbers

// integers below 2^16 (and types wider than 64 bits) use the standard modulo.
// integral types of at most 64 bits above that use an exact Miller-Rabin test instead, so that it is also usable
// for large primes (e.g. the characteristic of a field of 64 bits elements).
template <typename Integer>
static constexpr bool is_prime(const Integer p) {
  if (p <= 1) return false;
  if (p <= 3) return true;
  if (p % 2 == 0 || p % 3 == 0) return false;

  if constexpr (std::is_integral_v<Integer> && sizeof(Integer) <= sizeof(std::uint64_t)) {
    if (p >= 65536) return numbers::is_prime_miller_rabin(static_cast<std::uint64_t>(p));
  }

  for (long i = 5; i * i <= p; i = i + 6)
    if (p % i == 0 || p % (i + 2) == 0) return false;

  return true;
}

}  // namespace Gudhi

#endif  // GUDHI_ARITHMETIC_H_