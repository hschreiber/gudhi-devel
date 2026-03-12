/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Complex_pcoh_interface.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::Complex_pcoh_interface class.
 */

#ifndef PY_COMPLEX_PCOH_INTERFACE_H_INCLUDED
#define PY_COMPLEX_PCOH_INTERFACE_H_INCLUDED

#include <cstddef>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>

// Must come BEFORE any boost/range headers
namespace std {
template <>
struct iterator_traits<nanobind::iterator> {
  using difference_type = std::ptrdiff_t;
  using value_type = nanobind::handle;
  using pointer = const nanobind::handle *;
  using reference = const nanobind::handle &;
  using iterator_category = std::input_iterator_tag;  // single-pass only
};
}  // namespace std

#include <boost/range/iterator_range_core.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gudhi/Debug_utils.h>

namespace Gudhi {

class Complex_pcoh_interface {
 public:
  using Complex = nanobind::object;
  using Simplex_key = std::uint32_t;    /**< Simplex_key type */
  using Simplex_handle = std::uint32_t; /**< Simplex_handle type. Has to go from 0 to num_simplices */
  using Filtration_value = double;      /**< Internal filtration value type */
  using Dimension = int;                /**< Internal dimension type */

  Complex_pcoh_interface(Complex &cpx)
      : cpx_(cpx),
        null_key_(-1),
        null_simplex_(nanobind::cast<Simplex_handle>(cpx_.attr("null_handle")())),
        num_simplices_(nanobind::cast<std::size_t>(cpx_.attr("num_cells")())),
        filtration_f_(cpx_.attr("filtration")),
        dimension_f_(cpx_.attr("dimension")),
        endpoints_f_(cpx_.attr("endpoints")),
        filtration_simplex_range_f_(cpx_.attr("filtration_cell_range")),
        boundary_simplex_range_f_(cpx_.attr("boundary_cell_range")),
        handleToKey_(num_simplices_, null_key_),
        keyToHandle_(num_simplices_, null_simplex_) {}

  friend void swap(Complex_pcoh_interface &be1, Complex_pcoh_interface &be2) noexcept { std::swap(be1.cpx_, be2.cpx_); }

  [[nodiscard]] std::size_t num_simplices() const { return num_simplices_; }

  [[nodiscard]] Filtration_value filtration(Simplex_handle sh) const {
    return nanobind::cast<Filtration_value>(filtration_f_(sh));
  }

  [[nodiscard]] Dimension dimension() const { return nanobind::cast<Dimension>(dimension_f_()); }

  [[nodiscard]] Dimension dimension(Simplex_handle sh) const { return nanobind::cast<Dimension>(dimension_f_(sh)); }

  void assign_key(Simplex_handle sh, Simplex_key key) {
    if (key != null_key_)
      keyToHandle_[key] = sh;
    else
      keyToHandle_[handleToKey_[sh]] = null_simplex_;
    handleToKey_[sh] = key;
  }

  [[nodiscard]] Simplex_key key(Simplex_handle sh) const {
    if (sh == null_simplex_) return null_key_;
    return handleToKey_[sh];
  }

  [[nodiscard]] Simplex_key null_key() const { return null_key_; }

  [[nodiscard]] Simplex_handle simplex(Simplex_key key) const {
    if (key == null_key_) return null_simplex_;
    return keyToHandle_[key];
  }

  [[nodiscard]] Simplex_handle null_simplex() const { return null_simplex_; }

  // only used in update_cohomology_groups_edge, so not used without optimizations
  [[nodiscard]] std::pair<Simplex_handle, Simplex_handle> endpoints(Simplex_handle sh) const {
    return nanobind::cast<std::pair<Simplex_handle, Simplex_handle>>(endpoints_f_(sh));
  }

  [[nodiscard]] auto filtration_simplex_range() const {
    nanobind::object iterable = filtration_simplex_range_f_();
    return boost::make_iterator_range(nanobind::iter(iterable), nanobind::iterator::sentinel()) |
           boost::adaptors::transformed(
               [&](const nanobind::handle &item) -> Simplex_handle { return nanobind::cast<Simplex_handle>(item); });
  }

  [[nodiscard]] auto boundary_simplex_range(Simplex_handle sh) const {
    nanobind::object iterable = boundary_simplex_range_f_(sh);
    return boost::make_iterator_range(nanobind::iter(iterable), nanobind::iterator::sentinel()) |
           boost::adaptors::transformed(
               [&](const nanobind::handle &item) -> Simplex_handle { return nanobind::cast<Simplex_handle>(item); });
  }

 private:
  Complex cpx_;

  using F = decltype(cpx_.attr(""));

  Simplex_key null_key_;
  Simplex_handle null_simplex_;
  std::size_t num_simplices_;
  F filtration_f_;
  F dimension_f_;
  F endpoints_f_;
  F filtration_simplex_range_f_;
  F boundary_simplex_range_f_;

  std::vector<Simplex_key> handleToKey_;
  std::vector<Simplex_handle> keyToHandle_;
};

}  // namespace Gudhi

#endif  // PY_COMPLEX_PCOH_INTERFACE_H_INCLUDED
