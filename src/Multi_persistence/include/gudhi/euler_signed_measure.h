/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux, Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file euler_signed_measure.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Contains the method @ref Gudhi::multi_persistence::compute_euler_signed_measure.
 */

#ifndef GUDHI_MP_EULER_SIGNED_H_
#define GUDHI_MP_EULER_SIGNED_H_

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtration_value.h>
#include <gudhi/Multi_persistence/tensor_helpers.h>
#include <gudhi/Multi_persistence/utils.h>

namespace Gudhi {
namespace multi_persistence {

namespace detail {

// Number of trailing zero bits of x. Precondition: x != 0.
constexpr unsigned int count_trailing_zero_u64(unsigned long long int x) noexcept {
  unsigned int n = 0;
  while ((x & 1ULL) == 0) {
    x >>= 1;
    ++n;
  }
  return n;
}

// Adds `value` to every element of the cone
template <typename T, class StoragePolicy>
inline void _add_cone(const Tensor<T> &tensor,
                      const multi_filtration::Multi_parameter_filtration_value<StoragePolicy> &basepoint, T value) {
  GUDHI_CHECK(tensor.rank() == basepoint.num_parameters(), std::invalid_argument("Invalid coordinates for cone."));
  GUDHI_CHECK(tensor.rank() != 0, std::logic_error("Number of parameters cannot be 0."));

  const std::size_t r = tensor.rank();
  T *data = tensor.data_handle();

  if (r == 0) {
    data[0] += value;
    return;
  }

  std::vector<std::ptrdiff_t> len(r), cnt(r, 0);
  std::ptrdiff_t offset = 0;
  for (std::size_t d = 0; d < r; ++d) {
    const auto b = static_cast<std::ptrdiff_t>(basepoint(0, d));
    len[d] = static_cast<std::ptrdiff_t>(tensor.extent(d)) - b;
    if (len[d] <= 0) return;  // empty cone
    offset += b * static_cast<std::ptrdiff_t>(tensor.stride(d));
  }

  const std::ptrdiff_t n = len[r - 1];

  for (;;) {
    T *p = data + offset;
    for (std::ptrdiff_t k = 0; k < n; ++k) p[k] += value;

    // Advance the odometer over the outer dimensions r-2 ... 0.
    std::ptrdiff_t d = static_cast<std::ptrdiff_t>(r) - 2;
    for (; d >= 0; --d) {
      const std::ptrdiff_t sd = static_cast<std::ptrdiff_t>(tensor.stride(d));
      offset += sd;
      if (++cnt[d] < len[d]) break;
      offset -= len[d] * sd;
      cnt[d] = 0;
    }
    if (d < 0) return;
  }
}

template <typename T, class StoragePolicy>
inline void _add_cone_boundary(const Tensor<T> &tensor,
                               const multi_filtration::Multi_parameter_filtration_value<StoragePolicy> &basepoint,
                               T value) {
  const std::size_t r = tensor.rank();
  GUDHI_CHECK(r < 64, std::invalid_argument("Number of parameters are limited to 63 for optimization purposes."));
  T *data = tensor.data_handle();

  std::vector<std::ptrdiff_t> delta(r);  // hi - lo, per dimension
  std::ptrdiff_t off = 0;                // start at the all-"lo" corner
  for (std::size_t d = 0; d < r; ++d) {
    const auto st = static_cast<std::ptrdiff_t>(tensor.stride(d));
    const auto lo = static_cast<std::ptrdiff_t>(basepoint(0, d));
    const auto hi = static_cast<std::ptrdiff_t>(tensor.extent(d)) - 1;
    off += lo * st;
    delta[d] = (hi - lo) * st;
  }

  // constructing a mask to only cover the boundary
  const unsigned long long int corners = 1ULL << r;  // requires r < 64
  unsigned long long int gray = 0;
  data[off] += value;

  for (unsigned long long int i = 1; i < corners; ++i) {
    const unsigned int d = count_trailing_zero_u64(i);  // bit flipped at step i
    gray ^= 1ULL << d;
    off += ((gray >> d) & 1ULL) ? delta[d] : -delta[d];
    data[off] += (i & 1ULL) ? -value : value;
  }
}

template <class FilteredComplex, typename T>
inline void _get_euler_surface(const Tensor<T> &out, const FilteredComplex &cpx, bool mobiusInversion, bool zeroPad) {
  std::vector<typename Tensor<T>::index_type> coordinates(mobiusInversion && !zeroPad ? out.rank() : 0);
  for (auto sh : cpx.complex_simplex_range()) {
    const auto &fil = cpx.filtration(sh);
    GUDHI_CHECK(fil.num_parameters() == out.rank(),
                std::logic_error("Tensor does not have the right shape for given complex."));
    int sign = 1 - 2 * (cpx.dimension(sh) % 2);
    if (mobiusInversion && zeroPad) {
      detail::_add_cone_boundary(out, fil, sign);
    } else if (mobiusInversion) {
      for (std::size_t p = 0; p < fil.num_parameters(); ++p) coordinates[p] = fil(0, p);
      out[coordinates] += sign;
    } else {
      detail::_add_cone(out, fil, sign);
    }
  }
}

}  // namespace detail

// only takes into account first generator
template <typename T = int, typename I = std::uint16_t, class FilteredComplex, typename U>
inline Signed_measure<I, T> compute_euler_signed_measure(const FilteredComplex &cpx, const std::vector<U> &shape,
                                                         bool zeroPad, bool mobiusInversion = true,
                                                         bool verbose = false) {
  std::size_t total = 1;
  for (U value : shape) {
    total *= static_cast<std::size_t>(value);
  }
  std::vector<T> container(total, 0);
  Tensor<T> tensor(container, shape);

  print_progress_wo_line_break(verbose, "Container shape : ");
  for (auto r : container.extents()) print_progress_wo_line_break(verbose, r, ", ");
  print_progress_wo_line_break(verbose, "\nContainer size : ", container.size(),
                               "\nComputing Euler Characteristic ...");

  detail::_get_euler_surface(tensor, cpx, mobiusInversion, zeroPad);

  print_progress(verbose, "Done.");
  print_progress_wo_line_break(verbose, "Sparsifying the measure ...");

  auto res = sparsify<I>(tensor);

  print_progress(verbose, "Done.");

  return res;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_EULER_SIGNED_H_