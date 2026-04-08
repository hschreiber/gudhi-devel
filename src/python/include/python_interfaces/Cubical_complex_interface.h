/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_CUBICAL_INTERFACE_H_
#define INCLUDE_CUBICAL_INTERFACE_H_

#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>

namespace Gudhi {
namespace cubical_complex {

template <class CC>
inline std::vector<std::size_t> get_vertices(const CC& cpx, std::size_t sh) {
  if (cpx.dimension(sh) == 0) return {sh};
  return cpx.boundary_simplex_range(sh);
}

class Cubical_complex_interface : public Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>> {
 public:
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_base<double>>;
  using Base::Base;  // inheriting constructors

  explicit Cubical_complex_interface(const std::string& perseus_style_file) : Base(perseus_style_file.c_str()) {}

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // at another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<unsigned>& shape() { return this->sizes; };

  nanobind::ndarray<double, nanobind::numpy> get_numpy_array() {
    return nanobind::ndarray<double, nanobind::numpy>(Base::data.data(), {Base::data.size()});
  }

  [[nodiscard]] auto simplex_vertex_range(Simplex_handle sh) const {
    if (Base::dimension(sh) > 1)
      throw std::invalid_argument(
          "`simplex_vertex_range` is only valid for cells of dimension 0 and 1 in a cubical complex. Set argument "
          "max_dimension to 1 in the used python method.");

    return get_vertices(*this, sh);
  }
};

class Periodic_cubical_complex_interface
    : public Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>> {
 public:
  using Base = Bitmap_cubical_complex<Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>;
  using Base::Base;  // inheriting constructors

  explicit Periodic_cubical_complex_interface(const std::string& perseus_style_file)
      : Base(perseus_style_file.c_str()) {}

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // of another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<unsigned>& shape() { return this->sizes; };

  // TODO: nanobind is probably making a copy here (to verify), as it is only used privately we could think
  // of another strategy?
  // But as the vector is probably very small (number of dimensions), it is perhaps not worth it.
  const std::vector<bool>& periodicities() { return this->directions_in_which_periodic_b_cond_are_to_be_imposed; }

  nanobind::ndarray<double, nanobind::numpy> get_numpy_array() {
    return nanobind::ndarray<double, nanobind::numpy>(Base::data.data(), {Base::data.size()});
  }

  [[nodiscard]] auto simplex_vertex_range(Simplex_handle sh) const {
    if (Base::dimension(sh) > 1)
      throw std::invalid_argument(
          "`simplex_vertex_range` is only valid for cells of dimension 0 and 1 in a cubical complex. Set argument "
          "max_dimension to 1 in the used python method.");

    return get_vertices(*this, sh);
  }
};

}  // namespace cubical_complex
}  // namespace Gudhi

#endif  // INCLUDE_CUBICAL_INTERFACE_H_
