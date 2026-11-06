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
 * @file hilbert_signed_measure.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Contains the methods @ref Gudhi::multi_persistence::compute_hilbert_signed_measure and
 * @ref Gudhi::multi_persistence::compute_sparse_hilbert_signed_measure.
 */

#ifndef GUDHI_MP_HILBERT_SIGNED_H_
#define GUDHI_MP_HILBERT_SIGNED_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/tensor_helpers.h>
#include <gudhi/Multi_persistence/utils.h>
#include <gudhi/Multi_persistence/slice_sweep.h>

namespace Gudhi {
namespace multi_persistence {

namespace detail {

/////////// COMMON ///////////

/**
 * @private
 * @brief Family of paths of the Hilbert function: the lines of the grid which are parallel to the "barcode axis"
 * (the longest axis, which minimizes the number of lines). A line is identified by a "position", a container
 * (std::vector or std::array) with the same type and size as the grid shape, which holds the fixed value of every
 * axis except the barcode axis. position[barcodeAxis] is never read and stays 0.
 * Models the concept of family of @ref Slice_sweep.
 */
template <class Shape>
class Line_family {
 public:
  using Path = Shape;

  Line_family(const Shape &gridShape, bool zeroPad)
      : gridShape_(gridShape),
        activeShape_(_get_active_shape(gridShape, zeroPad)),
        barcodeAxis_(_longest_axis(activeShape_)),
        zeroPad_(zeroPad) {}

  /** @brief Shape of the output grid. */
  [[nodiscard]] const Shape &grid_shape() const { return gridShape_; }

  /** @brief Shape of the grid on which the lines are computed (one less on each axis with zero padding). */
  [[nodiscard]] const Shape &active_shape() const { return activeShape_; }

  /** @brief Axis along which the barcodes are computed. */
  [[nodiscard]] std::size_t barcode_axis() const { return barcodeAxis_; }

  [[nodiscard]] bool zero_pad() const { return zeroPad_; }

  /** @brief Number of grid points of a line. */
  [[nodiscard]] std::size_t width() const { return activeShape_[barcodeAxis_]; }

  /** @brief Number of lines parallel to the barcode axis in the active grid. */
  [[nodiscard]] std::size_t path_count() const {
    std::size_t count = 1;
    for (std::size_t axis = 0; axis < activeShape_.size(); ++axis) {
      if (axis != barcodeAxis_) count *= activeShape_[axis];
    }
    return count;
  }

  /**
   * @brief Returns the line number @p index. Lines are numbered like the digits of a mixed radix number: the lowest
   * axis other than the barcode axis is the fastest digit. Must stay consistent with @ref next_path.
   */
  [[nodiscard]] Path path_from_index(std::size_t index) const {
    Path position = activeShape_;
    position[barcodeAxis_] = 0;
    for (std::size_t axis = 0; axis < activeShape_.size(); ++axis) {
      if (axis == barcodeAxis_) continue;
      position[axis] = static_cast<typename Shape::value_type>(index % activeShape_[axis]);
      index /= activeShape_[axis];
    }
    return position;
  }

  /** @brief Advances @p position to the next line. Must stay consistent with @ref path_from_index. */
  void next_path(Path &position) const {
    for (std::size_t axis = 0; axis < position.size(); ++axis) {
      if (axis == barcodeAxis_) continue;
      if (++position[axis] < activeShape_[axis]) return;  // no carry
      position[axis] = 0;
    }
  }

  /**
   * @brief Returns the filtration value of a simplex restricted to the line described by @p position, that is the
   * smallest coordinate along the barcode axis of its generators which lie inside the slice.
   * Returns infinity if there is no such generator.
   */
  template <class Filtration_value>
  [[nodiscard]] typename Filtration_value::value_type slice_value(const Filtration_value &fil,
                                                                  const Path &position) const {
    typename Filtration_value::value_type birth = Filtration_value::T_inf;

    for (unsigned int g = 0; g < fil.num_generators(); ++g) {
      const auto candidate = fil(g, barcodeAxis_);

      if (candidate >= birth) continue;
      if (_is_generator_in_slice(fil, g, position)) birth = candidate;
    }

    return birth;
  }

 private:
  Shape gridShape_;
  Shape activeShape_;
  std::size_t barcodeAxis_;
  bool zeroPad_;

  /**
   * @private
   * @brief Returns @p shape, with every value decreased by one if @p zeroPad is true (the last grid point is then only
   * used as zero padding). In debug mode, checks that the shape is valid.
   */
  static Shape _get_active_shape(const Shape &shape, bool zeroPad) {
    GUDHI_CHECK_code(
      for (const auto extent : shape) {
        GUDHI_CHECK(extent > (zeroPad ? 1 : 0),
                    std::invalid_argument("Grid shape values must be positive, and greater than 1 when zero pad."));
      }
    )

        Shape activeShape = shape;
    if (zeroPad) {
      for (auto &extent : activeShape) --extent;
    }
    return activeShape;
  }

  /**
   * @private
   * @brief Returns the axis with the largest extent. Using it as barcode axis minimizes the number of lines.
   */
  static std::size_t _longest_axis(const Shape &shape) {
    std::size_t best = 0;
    for (std::size_t axis = 1; axis < shape.size(); ++axis) {
      if (shape[axis] > shape[best]) best = axis;
    }
    return best;
  }

  /**
   * @private
   * @brief Returns true if the generator @p g of @p fil appears within the slice described by @p position,
   * i.e. for all parameters other than the barcode axis.
   */
  template <class Filtration_value>
  bool _is_generator_in_slice(const Filtration_value &fil, unsigned int g, const Path &position) const {
    for (std::size_t p = 0; p < fil.num_parameters(); ++p) {
      if (p != barcodeAxis_ && fil(g, p) > position[p]) return false;
    }
    return true;
  }
};

/**
 * @private
 * @brief Computes the barcodes of the line at @p position with @p lineSlicer and calls
 * @p on_endpoint(degreeIndex, coordinate, sign) for every endpoint of every bar of the requested degrees: sign
 * is +1 at the birth of a bar and -1 at its death (or at the last coordinate if the bar does not die and there
 * is a zero padding). @p coordinate is the position of the endpoint along the barcode axis, and @p degreeIndex
 * the index of the degree in the degrees.
 */
template <class Slicer, class Shape, typename U, class LineSlicer, typename OnEndpoint>
inline void _for_each_endpoint(const Slice_sweep<Slicer, Line_family<Shape>, U> &sweep, LineSlicer &lineSlicer,
                               const Shape &position, OnEndpoint &&on_endpoint) {
  const typename LineSlicer::value_type width = sweep.family().width();
  const bool zeroPad = sweep.family().zero_pad();

  sweep.for_each_bar(lineSlicer, position, [&](std::size_t degreeIndex, auto birth, auto death) {
    // some births are the squeezed infinity sentinel
    if (birth >= width) return;

    on_endpoint(degreeIndex, static_cast<std::size_t>(birth), 1);
    if (death < width) {
      on_endpoint(degreeIndex, static_cast<std::size_t>(death), -1);
    } else if (zeroPad) {
      on_endpoint(degreeIndex, static_cast<std::size_t>(width) - 1, -1);
    }
  });
}

/////////// DENSE ONLY ///////////

/**
 * @private
 * @brief Computes the Hilbert surface into @p out (assumed to be a zero tensor of shape (degrees, grid shape)),
 * line by line. Along the barcode axis of @p sweep, the surface is already differentiated: the births and deaths
 * of the bars are written as +1 and -1.
 */
template <class Slicer, typename I, typename U, typename T>
inline void _get_hilbert_surface(const Slice_sweep<Slicer, Line_family<std::vector<I>>, U> &sweep,
                                 const Tensor<T> &out) {
  const std::size_t barcodeAxis = sweep.family().barcode_axis();

  sweep.for_each_path([&](auto &lineSlicer, const std::vector<I> &position) {
    // coordinates[0]: degree index, coordinates[c + 1]: parameter c. The barcode axis is set for each endpoint.
    std::vector<I> coordinates(position.size() + 1);
    std::copy(position.begin(), position.end(), coordinates.begin() + 1);

    _for_each_endpoint(sweep, lineSlicer, position, [&](std::size_t degreeIndex, std::size_t coordinate, int sign) {
      coordinates[0] = static_cast<I>(degreeIndex);
      coordinates[barcodeAxis + 1] = static_cast<I>(coordinate);
      out[coordinates] += static_cast<T>(sign);
    });
  });
}

/////////// SPARSE ONLY ///////////

template <std::size_t N, typename I>
inline std::array<I, 1 + N> _hilbert_sparse_key_shape(const std::array<I, N> &gridShape, I degreeCount) {
  std::array<I, 1 + N> keyShape{};
  keyShape[0] = degreeCount;
  for (std::size_t axis = 0; axis < N; ++axis) keyShape[1 + axis] = gridShape[axis];
  return keyShape;
}

/**
 * @private
 * @brief Adds to @p smPts the atoms of the signed measure which come from the endpoint at @p coordinate (along
 * @p barcodeAxis) of a bar of the line at @p position, with sign @p sign.
 *
 * The signed measure is the Möbius inversion of the Hilbert function. Along the barcode axis, it is given by the
 * endpoints of the bars. Along each other axis, the inversion turns a point into the point itself and its successor
 * (if it is inside the grid), with the opposite sign.
 */
template <std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_endpoint(Accumulator &smPts, std::size_t degreeIndex, const std::array<I, N> &position,
                           std::size_t barcodeAxis, std::size_t coordinate, int sign,
                           const std::array<I, N> &gridShape) {
  Atom_options<N, T, I> options;
  for (std::size_t axis = 0; axis < N; ++axis) {
    if (axis == barcodeAxis) {
      options.set_primary(axis, static_cast<I>(coordinate), static_cast<T>(sign));
      continue;
    }
    options.set_primary(axis, position[axis], static_cast<T>(1));
    if (position[axis] + 1 < gridShape[axis]) {
      options.set_secondary(axis, static_cast<I>(position[axis] + 1), static_cast<T>(-1));
    }
  }
  _emit_atom_product(smPts, degreeIndex, options, [](const auto &) { return true; });
}

template <std::size_t N, typename T, typename I, class Slicer, typename U>
inline Signed_measure<I, T> _compute_hilbert_signed_measure_sparse_by_slices(
    const Slice_sweep<Slicer, Line_family<std::array<I, N>>, U> &sweep) {
  using Accumulator = Packed_sparse_accumulator<1 + N, T, I>;

  const auto &family = sweep.family();
  const auto keyShape = _hilbert_sparse_key_shape<N>(family.grid_shape(), static_cast<I>(sweep.degree_count()));

  const Accumulator prototype(keyShape);
  Thread_local_storage<Accumulator> threadSmPts(prototype);

  sweep.for_each_path([&](auto &lineSlicer, const std::array<I, N> &position) {
    Accumulator &smPts = threadSmPts.local();

    _for_each_endpoint(sweep, lineSlicer, position, [&](std::size_t degreeIndex, std::size_t coordinate, int sign) {
      _emit_endpoint<N, T>(smPts, degreeIndex, position, family.barcode_axis(), coordinate, sign, family.grid_shape());
    });
  });

  return _merge_sparse_accumulators<1 + N, T, I>(threadSmPts, keyShape);
}

template <std::size_t N, typename T, typename I, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> _compute_sparse_hsm_for_n_param(Slicer &slicer, const std::vector<U1> &shape,
                                                            const std::vector<U2> &degrees, bool zeroPad,
                                                            bool ignoreInf, bool verbose, int n_jobs) {
  std::array<I, N> parameterShape;
  std::copy_n(shape.begin(), N, parameterShape.begin());  // to get the right type
  const Slice_sweep sweep(slicer, Line_family<std::array<I, N>>(parameterShape, zeroPad), degrees, ignoreInf, n_jobs);

  print_progress_wo_line_break(verbose, "Computing sparse Hilbert invariant ... ");

  auto out = _compute_hilbert_signed_measure_sparse_by_slices<N, T>(sweep);

  print_progress(verbose, "Done.");

  return out;
}

}  // namespace detail

// TODO: Once the complex is extracted from Slicer, we could pass here only the complex and pass its pointer
// to a Slicer we construct our-self? Except if this adds too many templates
/**
 * @brief Computes the signed measure of the Hilbert function of the persistence module of @p slicer.
 *
 * @param slicer Slicer of the filtered complex. The persistence related storage will be modified.
 * @param shape Size of the grid in each parameter. Its size has to be the number of parameters of the slicer.
 * @param degrees Homological degrees for which to compute the invariant.
 * @param zeroPad If true, the last grid point of each parameter is only used as zero padding. Default: false.
 * @param ignoreInf Ignore the infinite bars during the persistence computation. Default: true.
 * @param verbose Print progress messages. Default: false.
 * @param n_jobs Maximal number of threads (when compiled with TBB). Ignored by vine persistence, which is sequential.
 * Default: -1 (i.e. as much as possible).
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_hilbert_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                           const std::vector<U2> &degrees, bool zeroPad = false,
                                                           bool ignoreInf = true, bool verbose = false,
                                                           int n_jobs = -1) {
  GUDHI_CHECK(slicer.get_number_of_parameters() == shape.size(), std::invalid_argument("Invalid shape size."));
  GUDHI_CHECK(slicer.get_number_of_parameters() >= 2,
              std::invalid_argument(
                  "compute_hilbert_signed_measure is not implemented for number of parameters strictly lower than 2."));

  // The tensor has one more dimension than the grid: the degree.
  std::vector<U1> tensorShape;
  tensorShape.reserve(shape.size() + 1);
  tensorShape.push_back(static_cast<U1>(degrees.size()));
  tensorShape.insert(tensorShape.end(), shape.begin(), shape.end());

  std::size_t total = 1;
  for (U1 value : tensorShape) {
    GUDHI_CHECK(value >= 0, std::invalid_argument("The shape must be explicit and all its values positive."));
    total *= static_cast<std::size_t>(value);
  }
  if (total == 0) return {};
  std::vector<T> container(total, 0);
  Tensor<T> tensor(container, tensorShape);

  print_progress_wo_line_break(verbose, "Container shape : ");
  for (std::size_t axis = 0; axis < tensor.rank(); ++axis) {
    print_progress_wo_line_break(verbose, tensor.extent(axis), ", ");
  }
  print_progress_wo_line_break(verbose, "\nContainer size : ", tensor.size(), "\nComputing Hilbert invariant ...");

  const std::vector<I> gridShape(shape.begin(), shape.end());
  const detail::Slice_sweep sweep(slicer, detail::Line_family<std::vector<I>>(gridShape, zeroPad), degrees, ignoreInf,
                                  n_jobs);

  detail::_get_hilbert_surface(sweep, tensor);

  print_progress(verbose, "Done.");
  print_progress_wo_line_break(verbose, "Computing mobius inversion ...");

  // The inversion along the barcode axis was already done when writing the births and deaths.
  for (std::size_t axis = 0; axis < gridShape.size(); ++axis) {
    if (axis != sweep.family().barcode_axis()) differentiate(tensor, axis + 1);  // +1 for degree
  }

  print_progress(verbose, "Done.");
  print_progress_wo_line_break(verbose, "Sparsifying the measure ...");

  auto res = sparsify<I>(tensor);

  print_progress(verbose, "Done.");

  return res;
}

// TODO: Once the complex is extracted from Slicer, we could pass here only the complex and pass its pointer
// to a Slicer we construct our-self? Except if this adds too many templates
/**
 * @brief Computes the signed measure of the Hilbert function of the persistence module of @p slicer, without
 * building the dense tensor of the Hilbert function. Same result as @ref compute_hilbert_signed_measure, but only
 * implemented for 2, 3 and 4 parameters.
 *
 * @param slicer Slicer of the filtered complex. The persistence related storage will be modified.
 * @param shape Size of the grid in each parameter. Its size has to be the number of parameters of the slicer.
 * @param degrees Homological degrees for which to compute the invariant.
 * @param zeroPad If true, the last grid point of each parameter is only used as zero padding. Default: false.
 * @param ignoreInf Ignore the infinite bars during the persistence computation. Default: true.
 * @param verbose Print progress messages. Default: false.
 * @param n_jobs Maximal number of threads (when compiled with TBB). Ignored by vine persistence, which is sequential.
 * Default: -1 (i.e. as much as possible).
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_sparse_hilbert_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                                  const std::vector<U2> &degrees, bool zeroPad = false,
                                                                  bool ignoreInf = true, bool verbose = false,
                                                                  int n_jobs = -1) {
  GUDHI_CHECK(slicer.get_number_of_parameters() == shape.size(), std::invalid_argument("Invalid shape size."));

  if (degrees.size() == 0) return {};

  return detail::_dispatch_on_number_of_parameters(
      slicer.get_number_of_parameters(), "Sparse Hilbert signed measure is implemented for 2, 3, and 4 parameters.",
      [&](auto numberOfParameters) {
        return detail::_compute_sparse_hsm_for_n_param<decltype(numberOfParameters)::value, T, I>(
            slicer, shape, degrees, zeroPad, ignoreInf, verbose, n_jobs);
      });
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_HILBERT_SIGNED_H_