/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

// #include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <array>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include <gudhi/Debug_utils.h>
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Matrix.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

template <typename T>
static nanobind::object _to_nanobind_object(T&& val) {
  using U = std::decay_t<T>;
  // to avoid non trivial nanobind::cast when not necessary
  if constexpr (std::is_base_of_v<nanobind::object, U>) {
    return std::forward<T>(val);
  } else {
    return nanobind::cast(std::forward<T>(val));
  }
}

// more efficient than a list.append if n is known in advance
// (and the object does not have to be resized later obviously)
template <class F>
nanobind::tuple _build_tuple(std::size_t n, F&& construct_row) {
  PyObject* out = PyTuple_New(static_cast<Py_ssize_t>(n));
  if (!out) throw nanobind::python_error();

  for (std::size_t p = 0; p < n; ++p) {
    nanobind::object row = _to_nanobind_object(std::forward<F>(construct_row)(p));
    PyTuple_SET_ITEM(out, p, row.release().ptr());
  }

  return nanobind::steal<nanobind::tuple>(out);
}

template <typename D, typename I, bool is_RU>
struct Rep_matrix_options
    : Gudhi::persistence_matrix::Default_options<Gudhi::persistence_matrix::Column_types::NAIVE_VECTOR, true> {
  using Dimension = D;
  using Index = I;

  static const bool is_of_boundary_type = is_RU;
  static const Gudhi::persistence_matrix::Column_indexation_types column_indexation_type =
      Gudhi::persistence_matrix::Column_indexation_types::POSITION;
  static const bool has_column_pairings = true;
  static const bool has_vine_update = false;
  static const bool can_retrieve_representative_cycles = true;
};

template <typename Index, typename Dimension>
Index assign_keys_(Gudhi::Simplex_tree_interface& complex, Index start, Dimension maxDimension) {
  Index numberOfSimplices = start;
  // Vertex ID should correspond to position in the original point cloud if there was one.
  for (auto sh : complex.skeleton_simplex_range(0)) {
    complex.assign_key(sh, numberOfSimplices);
    ++numberOfSimplices;
  }
  for (auto sh : complex.filtration_simplex_range()) {
    if (complex.dimension(sh) != 0 && static_cast<Dimension>(complex.dimension(sh)) <= maxDimension) {
      complex.assign_key(sh, numberOfSimplices);
      ++numberOfSimplices;
    }
  }
  return numberOfSimplices;
}

template <typename Filtration_value, typename Index, typename Dimension>
void build_boundary_matrix_from_complex(Gudhi::Simplex_tree_interface& complex,
                                        std::vector<std::vector<Index>>& boundaries, std::vector<Dimension>& dimensions,
                                        std::vector<Filtration_value>& filtrationValues,
                                        Dimension maxDimension = std::numeric_limits<Dimension>::max()) {
  GUDHI_CHECK(boundaries.size() == dimensions.size() && boundaries.size() == filtrationValues.size(),
              std::invalid_argument("Output containers do not start with the same size"));

  if (maxDimension < 0) return;

  Index numberOfSimplices = assign_keys_(complex, boundaries.size(), maxDimension);

  boundaries.resize(numberOfSimplices);
  dimensions.resize(numberOfSimplices);
  filtrationValues.resize(numberOfSimplices);

  for (auto sh : complex.filtration_simplex_range()) {
    if (static_cast<Dimension>(complex.dimension(sh)) <= maxDimension) {
      auto index = complex.key(sh);
      dimensions[index] = complex.dimension(sh);
      filtrationValues[index] = complex.filtration(sh);
      std::vector<Index> boundary;
      for (auto b : complex.boundary_simplex_range(sh)) {
        boundary.push_back(complex.key(b));
      }
      std::sort(boundary.begin(), boundary.end());
      boundaries[index] = std::move(boundary);
    }
  }
}

template <class Index>
void compute_union(std::size_t i, const std::vector<std::vector<Index>>& boundaries,
                   std::vector<std::vector<Index>>& vertices, std::vector<Index>& tmp) {
  auto& vi = vertices[i];
  const auto& b = boundaries[i];
  vi = vertices[b[0]];
  for (std::size_t k = 1; k < b.size(); ++k) {
    const auto& other = vertices[b[k]];
    tmp.clear();
    std::set_union(vi.begin(), vi.end(), other.begin(), other.end(), std::back_inserter(tmp));
    vi.swap(tmp);
  }
}

template <class Index, class Dimension>
std::vector<std::vector<Index>> get_vertices_from_boundaries(const std::vector<std::vector<Index>>& boundaries,
                                                             const std::vector<Dimension>& dimensions,
                                                             int maxDim = -1) {
  GUDHI_CHECK(boundaries.size() == dimensions.size(),
              std::invalid_argument("There should be as many dimensions than boundaries."));

  if (dimensions.size() == 0) return {};

  const std::size_t numberOfCells = boundaries.size();
  std::vector<std::vector<Index>> vertices(numberOfCells);

  // initialize leafs/vertices
  Index c = 0;
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (boundaries[i].empty()) {
      vertices[i].assign(1, c++);
    }
  }

  // DFS in the tree of face/coface relation
  enum : unsigned char { NEW = 0, IN_PROGRESS = 1, DONE = 2 };
  std::vector<unsigned char> state(numberOfCells);
  for (std::size_t i = 0; i < numberOfCells; ++i)
    state[i] = (boundaries[i].empty() || (maxDim >= 0 && dimensions[i] > maxDim)) ? DONE : NEW;

  // TODO: replace this with the recursive version?
  // non-recursion is to avoid possible stack overflow if boundaries is malformed.
  // but not so sure that should really happen and the recursive version would be much shorter to write
  std::vector<Index> tmp;
  std::vector<std::pair<std::size_t, std::size_t>> stack;  // (cell, next child position)

  for (std::size_t root = 0; root < numberOfCells; ++root) {
    if (state[root] != NEW) continue;
    state[root] = IN_PROGRESS;
    stack.emplace_back(root, 0);

    while (!stack.empty()) {
      auto& [cell, pos] = stack.back();
      const auto& b = boundaries[cell];

      if (pos < b.size()) {
        auto child = static_cast<std::size_t>(b[pos++]);
        if (state[child] == IN_PROGRESS) throw std::invalid_argument("Cycle in boundaries.");
        if (state[child] == NEW) {
          state[child] = IN_PROGRESS;
          stack.emplace_back(child, 0);
        }
      } else {
        compute_union<Index>(cell, boundaries, vertices, tmp);
        state[cell] = DONE;
        stack.pop_back();
      }
    }
  }
  return vertices;
}

template <bool is_RU = true>
nb::tuple get_rep_cycles(Gudhi::Simplex_tree_interface& st) {
  using Index = Gudhi::Simplex_tree_interface::Simplex_key;
  using D = int;
  using F = Gudhi::Simplex_tree_interface::Filtration_value;
  using Matrix_options = Rep_matrix_options<D, Index, is_RU>;
  using Matrix = Gudhi::persistence_matrix::Matrix<Matrix_options>;
  using D_barcode = std::vector<F>;
  using Edge = std::array<Index, 2>;
  using Triangle = std::array<Index, 3>;
  using Cycle1 = std::vector<Edge>;
  using Cycle2 = std::vector<Triangle>;

  std::vector<std::vector<Index>> bc;
  std::vector<D> dc;
  std::vector<F> fc;

  if (st.dimension() != 1) throw std::invalid_argument("Expects a graph simplex tree.");

  st.expansion(3);
  st.initialize_filtration();

  build_boundary_matrix_from_complex(st, bc, dc, fc);
  std::vector<std::vector<Index>> vertices = get_vertices_from_boundaries(bc, dc, 2);

  Matrix matrix(bc);
  const auto& barcode = matrix.get_current_barcode();
  std::vector<D_barcode> outBC(2);
  std::vector<Cycle1> out1Cycles;
  std::vector<Cycle2> out2Cycles;
  for (const auto& bar : barcode) {
    // only essential 1 and 2-cycles
    if ((bar.dim == 1 || bar.dim == 2) && bar.death == Matrix::Bar::inf) {
      outBC[bar.dim - 1].push_back(fc[bar.birth]);
      matrix.update_representative_cycle(bar);
      const std::vector<Index>& cycle = matrix.get_representative_cycle(bar);
      if (bar.dim == 1) {
        out1Cycles.emplace_back();
        for (Index c : cycle) {
          out1Cycles.back().emplace_back(Edge{vertices[c][0], vertices[c][1]});
        }
      } else {
        out2Cycles.emplace_back();
        for (Index c : cycle) {
          out2Cycles.back().emplace_back(Triangle{vertices[c][0], vertices[c][1], vertices[c][2]});
        }
      }
    }
  }

  return nb::make_tuple(
      _wrap_as_numpy_array(std::move(outBC[0]), outBC[0].size()),
      _build_tuple(out1Cycles.size(),
                   [&out1Cycles](std::size_t i) { return _wrap_as_numpy_array(std::move(out1Cycles[i])); }),
      _wrap_as_numpy_array(std::move(outBC[1]), outBC[1].size()),
      _build_tuple(out2Cycles.size(),
                   [&out2Cycles](std::size_t i) { return _wrap_as_numpy_array(std::move(out2Cycles[i])); }));
}

NB_MODULE(_rep_cycles_ext, m) {
  m.attr("__license__") = "MIT";
  m.def("get_rep_cycles_v1", &get_rep_cycles<true>).def("get_rep_cycles_v2", &get_rep_cycles<false>);
}
