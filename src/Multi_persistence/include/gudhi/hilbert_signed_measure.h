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
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/tensor_helpers.h>
#include <gudhi/Multi_persistence/utils.h>

namespace Gudhi {
namespace multi_persistence {

namespace detail {

// The Hilbert function is computed line by line: a "line" is a set of grid points which only differ by their
// coordinate along the "barcode axis". A line is identified by a "position", a container (std::vector or std::array)
// with the same type and size as the grid shape, which holds the fixed value of every axis except the barcode axis.
// position[barcodeAxis] is never read and stays 0.

/////////// COMMON ///////////

/**
 * @private
 * @brief Holds what stays constant during one computation of a Hilbert signed measure (slicer, grid, degrees,
 * options) and enumerates the lines of the grid and the endpoints of their bars.
 */
template <class Slicer, class Shape, typename U>
class Line_sweep {
 public:
  // TODO: slicer could be constructed only here from complex
  Line_sweep(Slicer &slicer, const Shape &gridShape, const std::vector<U> &degrees, bool zeroPad, bool ignoreInf,
             int n_jobs)
      : slicer_(&slicer),
        gridShape_(&gridShape),
        degrees_(&degrees),
        activeShape_(_get_active_shape(gridShape, zeroPad)),
        barcodeAxis_(_longest_axis(activeShape_)),
        zeroPad_(zeroPad),
        ignoreInf_(ignoreInf),
        nJobs_(n_jobs) {}

  /** @brief Shape of the output grid. */
  [[nodiscard]] const Shape &grid_shape() const { return *gridShape_; }

  /** @brief Shape of the grid on which the lines are computed (one less on each axis with zero padding). */
  [[nodiscard]] const Shape &active_shape() const { return activeShape_; }

  /** @brief Axis along which the barcodes are computed. */
  [[nodiscard]] std::size_t barcode_axis() const { return barcodeAxis_; }

  [[nodiscard]] std::size_t degree_count() const { return degrees_->size(); }

  /**
   * @brief Calls @p on_line(slicer, position) for every line of the active grid which is parallel to the barcode
   * axis. The slicer given to the callback is the one to use to compute the barcodes of that line.
   *
   * With TBB, the lines are split into contiguous chunks which are processed in parallel, each thread working on
   * its own copy of the slicer. Otherwise, and when the persistence is a vine (which needs consecutive lines to be
   * close to each other), the lines are processed sequentially, in order, with the slicer itself.
   * Therefore, @p on_line has to be thread safe.
   */
  template <typename OnLine>
  void for_each_line(OnLine &&on_line) const {
    const std::size_t nbLines = _line_count();
    if (nbLines == 0) return;

    // Visits the lines [begin, end).
    auto visit = [&](auto &lineSlicer, std::size_t begin, std::size_t end) {
      Shape position = _position_from_index(begin);
      for (std::size_t line = begin; line != end; ++line) {
        on_line(lineSlicer, position);
        _increment_position(position);
      }
    };

#ifdef GUDHI_USE_TBB
    if (!Slicer::Persistence::is_vine && nJobs_ != 1) {
      using ThreadSafe = typename Slicer::Thread_safe;
      ThreadSafe prototype(*slicer_);
      tbb::enumerable_thread_specific<ThreadSafe> threadSlicers(prototype);

      tbb::task_arena arena(nJobs_);  // limits the number of threads
      arena.execute([&] {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nbLines),
                          [&](const tbb::blocked_range<std::size_t> &range) {
                            visit(threadSlicers.local(), range.begin(), range.end());
                          });
      });
      return;
    }
#endif

    visit(*slicer_, 0, nbLines);
  }

  /**
   * @brief Computes the barcodes of the line at @p position with @p lineSlicer and calls
   * @p on_endpoint(degreeIndex, coordinate, sign) for every endpoint of every bar of the requested degrees: sign
   * is +1 at the birth of a bar and -1 at its death (or at the last coordinate if the bar does not die and there
   * is a zero padding). @p coordinate is the position of the endpoint along the barcode axis, and @p degreeIndex
   * the index of the degree in the degrees.
   */
  template <class LineSlicer, typename OnEndpoint>
  void for_each_endpoint(LineSlicer &lineSlicer, const Shape &position, OnEndpoint &&on_endpoint) const {
    const std::size_t width = activeShape_[barcodeAxis_];

    auto &slice = lineSlicer.get_slice();
    const auto &filtrationValues = lineSlicer.get_filtration_values();

    for (std::size_t k = 0; k < filtrationValues.size(); ++k) {
      slice[k] = _get_line_filtration(filtrationValues[k], position);
    }

    if (!lineSlicer.persistence_computation_is_initialized()) {
      lineSlicer.initialize_persistence_computation(ignoreInf_);
    } else {
      lineSlicer.update_persistence_computation(ignoreInf_);
    }

    const auto barcodes = lineSlicer.template get_flat_barcode<true>();
    std::size_t degreeIndex = 0;
    for (U degree : *degrees_) {
      if (degree >= 0 && static_cast<std::size_t>(degree) < barcodes.size()) {
        for (const auto &bar : barcodes[degree]) {
          // some births are the squeezed infinity sentinel
          if (bar[0] >= static_cast<typename Slicer::value_type>(width)) continue;

          on_endpoint(degreeIndex, static_cast<std::size_t>(bar[0]), 1);
          if (bar[1] < static_cast<typename Slicer::value_type>(width)) {
            on_endpoint(degreeIndex, static_cast<std::size_t>(bar[1]), -1);
          } else if (zeroPad_) {
            on_endpoint(degreeIndex, width - 1, -1);
          }
        }
      }

      ++degreeIndex;
    }
  }

 private:
  Slicer *slicer_;
  const Shape *gridShape_;
  const std::vector<U> *degrees_;
  Shape activeShape_;
  std::size_t barcodeAxis_;
  bool zeroPad_;
  bool ignoreInf_;
  int nJobs_;

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
   * @brief Returns the number of lines parallel to @p barcodeAxis in a grid of shape @p shape.
   */
  std::size_t _line_count() const {
    std::size_t count = 1;
    for (std::size_t axis = 0; axis < activeShape_.size(); ++axis) {
      if (axis != barcodeAxis_) count *= activeShape_[axis];
    }
    return count;
  }

  /**
   * @private
   * @brief Returns the position of the line number @p index. Lines are numbered like the digits of a mixed radix
   * number: the lowest axis other than @p barcodeAxis is the fastest digit. Must stay consistent with
   * @ref _increment_position.
   */
  Shape _position_from_index(std::size_t index) const {
    Shape position = activeShape_;
    position[barcodeAxis_] = 0;
    for (std::size_t axis = 0; axis < activeShape_.size(); ++axis) {
      if (axis == barcodeAxis_) continue;
      position[axis] = static_cast<typename Shape::value_type>(index % activeShape_[axis]);
      index /= activeShape_[axis];
    }
    return position;
  }

  /**
   * @private
   * @brief Advances @p position to the next line, like an odometer. Must stay consistent with
   * @ref _position_from_index. After the last line, @p position wraps back to all zeros.
   */
  void _increment_position(Shape &position) const {
    for (std::size_t axis = 0; axis < position.size(); ++axis) {
      if (axis == barcodeAxis_) continue;
      if (++position[axis] < activeShape_[axis]) return;  // no carry
      position[axis] = 0;
    }
  }

  /**
   * @private
   * @brief Returns true if the generator @p g of @p fil appears within the slice described by @p position,
   * i.e. for all parameters other than the barcode axis.
   */
  template <class Filtration_value, class Position>
  bool _is_generator_in_slice(const Filtration_value &fil, unsigned int g, const Position &position) const {
    for (std::size_t p = 0; p < fil.num_parameters(); ++p) {
      if (p != barcodeAxis_ && fil(g, p) > position[p]) return false;
    }
    return true;
  }

  /**
   * @private
   * @brief Returns the filtration value of a simplex restricted to the line described by @p position, that is the
   * smallest coordinate along @p barcodeAxis of its generators which lie inside the slice.
   * Returns infinity if there is no such generator.
   */
  template <class Filtration_value, class Position>
  typename Filtration_value::value_type _get_line_filtration(const Filtration_value &fil,
                                                             const Position &position) const {
    typename Filtration_value::value_type birth = Filtration_value::T_inf;

    for (unsigned int g = 0; g < fil.num_generators(); ++g) {
      const auto candidate = fil(g, barcodeAxis_);

      if (candidate >= birth) continue;
      if (_is_generator_in_slice(fil, g, position)) birth = candidate;
    }

    return birth;
  }
};

/////////// DENSE ONLY ///////////

/**
 * @private
 * @brief Computes the Hilbert surface into @p out (assumed to be a zero tensor of shape (degrees, grid shape)),
 * line by line. Along the barcode axis of @p sweep, the surface is already differentiated: the births and deaths
 * of the bars are written as +1 and -1.
 */
template <typename I, class Slicer, class Shape, typename T, typename U>
inline void _get_hilbert_surface(const Line_sweep<Slicer, Shape, U> &sweep, const Tensor<T> &out) {
  sweep.for_each_line([&](auto &lineSlicer, const std::vector<I> &position) {
    // coordinates[0]: degree index, coordinates[c + 1]: parameter c. The barcode axis is set for each endpoint.
    std::vector<I> coordinates(position.size() + 1);
    std::copy(position.begin(), position.end(), coordinates.begin() + 1);

    sweep.for_each_endpoint(lineSlicer, position, [&](std::size_t degreeIndex, std::size_t coordinate, int sign) {
      coordinates[0] = static_cast<I>(degreeIndex);
      coordinates[sweep.barcode_axis() + 1] = static_cast<I>(coordinate);
      out[coordinates] += static_cast<T>(sign);
    });
  });
}

/////////// SPARSE ONLY ///////////

#ifdef GUDHI_USE_TBB
template <class T>
using Thread_local_storage = tbb::enumerable_thread_specific<T>;
#else
/**
 * @private
 * @brief Stand-in for tbb::enumerable_thread_specific when TBB is not used: a single value for the single thread.
 */
template <class T>
class Thread_local_storage {
 public:
  explicit Thread_local_storage(const T &exemplar) : value_(exemplar) {}

  T &local() { return value_; }

  T *begin() { return &value_; }

  T *end() { return &value_ + 1; }

  const T *begin() const { return &value_; }

  const T *end() const { return &value_ + 1; }

 private:
  T value_;
};
#endif

template <std::size_t N, typename I>
inline std::array<I, 1 + N> _hilbert_sparse_key_shape(const std::array<I, N> &gridShape, I degreeCount) {
  std::array<I, 1 + N> keyShape{};
  keyShape[0] = degreeCount;
  for (std::size_t axis = 0; axis < N; ++axis) keyShape[1 + axis] = gridShape[axis];
  return keyShape;
}

/**
 * @private
 * @brief Sums signed weights of the atoms of a measure on a (degree, grid) space. Atoms are stored in a hash map
 * under a 64-bit key, which packs the coordinates in mixed radix. Atoms whose weights cancel out are removed.
 */
template <std::size_t N, typename T, typename I>
class Packed_hilbert_sparse_accumulator {
 public:
  using Key = std::uint64_t;
  using Coordinate = std::array<I, 1 + N>;
  using Entry = std::pair<Key, T>;

  explicit Packed_hilbert_sparse_accumulator(const Coordinate &key_shape) : keyShape_(key_shape) {
    validate_key_shape();
  }

  [[nodiscard]] std::size_t size() const { return weights_.size(); }

  void add(const Coordinate &coordinates, const T weight) { add_key(pack(coordinates), weight); }

  void add_key(const Key key, const T weight) {
    if (weight == T{0}) return;
    auto [it, inserted] = weights_.emplace(key, weight);
    if (inserted) return;
    it->second += weight;
    if (it->second == T{0}) weights_.erase(it);
  }

  void append_entries_to(std::vector<Entry> &out) const {
    for (const auto &[key, weight] : weights_) out.emplace_back(key, weight);
  }

  [[nodiscard]] static Coordinate unpack(Key key, const Coordinate &keyShape) {
    Coordinate coordinates{};
    for (std::size_t axis = 1 + N; axis-- > 0;) {
      const auto extent = static_cast<Key>(keyShape[axis]);
      coordinates[axis] = static_cast<I>(key % extent);
      key /= extent;
    }
    return coordinates;
  }

 private:
  void validate_key_shape() const {
    Key capacity = 1;
    constexpr Key max_key = std::numeric_limits<Key>::max();
    for (const auto extent_value : keyShape_) {
      if (extent_value <= I{0}) {
        throw std::runtime_error("Internal error: invalid sparse Hilbert key shape.");
      }
      const auto extent = static_cast<Key>(extent_value);
      if (capacity > max_key / extent) {
        throw std::runtime_error("Sparse Hilbert key space exceeds 64-bit capacity.");
      }
      capacity *= extent;
    }
  }

  [[nodiscard]] Key pack(const Coordinate &coordinates) const {
    Key key = 0;
    for (std::size_t axis = 0; axis < 1 + N; ++axis) {
      const auto extent = static_cast<Key>(keyShape_[axis]);
      const auto coordinate = static_cast<Key>(coordinates[axis]);
      GUDHI_CHECK(coordinate < extent, std::invalid_argument("Internal error: invalid sparse Hilbert coordinate."));
      key = key * extent + coordinate;
    }
    return key;
  }

  Coordinate keyShape_{};
  std::unordered_map<Key, T> weights_;
};

/**
 * @private
 * @brief Adds to @p smPts the atoms of the signed measure which come from the endpoint at @p coordinate (along
 * @p barcodeAxis) of a bar of the line at @p position, with sign @p sign.
 *
 * The signed measure is the Möbius inversion of the Hilbert function. Along the barcode axis, it is given by the
 * endpoints of the bars. Along each other axis, the inversion turns a point into the point itself and its successor
 * (if it is inside the grid), with the opposite sign. So the endpoint yields one atom for each subset of the other
 * axes, shifted by one along the axes of the subset, with a sign flipped once per shifted axis.
 */
template <std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_endpoint(Accumulator &smPts, std::size_t degreeIndex, const std::array<I, N> &position,
                           std::size_t barcodeAxis, std::size_t coordinate, int sign,
                           const std::array<I, N> &outputShape) {
  for (std::size_t subset = 0; subset < (std::size_t{1} << N); ++subset) {
    if ((subset >> barcodeAxis) & 1) continue;  // no successor along the barcode axis

    std::array<I, 1 + N> atom{};
    atom[0] = static_cast<I>(degreeIndex);
    T weight = static_cast<T>(sign);
    bool inside = true;
    for (std::size_t axis = 0; axis < N && inside; ++axis) {
      const bool shifted = (subset >> axis) & 1;
      atom[1 + axis] =
          axis == barcodeAxis ? static_cast<I>(coordinate) : static_cast<I>(position[axis] + (shifted ? 1 : 0));
      inside = atom[1 + axis] < outputShape[axis];
      if (shifted) weight = static_cast<T>(-weight);
    }
    if (inside) smPts.add(atom, weight);
  }
}

template <std::size_t N, typename T, typename I>
inline Signed_measure<I, T> _merge_sparse_hilbert_accumulators(
    const Thread_local_storage<Packed_hilbert_sparse_accumulator<N, T, I>> &threadSmPts,
    const std::array<I, 1 + N> &keyShape) {
  using Accumulator = Packed_hilbert_sparse_accumulator<N, T, I>;

  std::size_t total = 0;
  for (const auto &local : threadSmPts) total += local.size();

  std::vector<typename Accumulator::Entry> entries;
  entries.reserve(total);
  for (const auto &local : threadSmPts) local.append_entries_to(entries);
  std::sort(entries.begin(), entries.end(),
            [](const auto &left, const auto &right) { return left.first < right.first; });

  Signed_measure<I, T> sm;
  sm.numberOfParameters = N + 1;
  // TODO: would be better to not over-reserve as this will survive for a long time
  sm.coordinates.reserve(entries.size() * sm.numberOfParameters);
  sm.values.reserve(entries.size());

  // The same key can come from several threads: combine them, and drop the atoms which cancel out.
  for (std::size_t i = 0; i < entries.size();) {
    const auto key = entries[i].first;
    T weight = 0;
    do {
      weight += entries[i].second;
      ++i;
    } while (i < entries.size() && entries[i].first == key);
    if (weight != 0) {
      const auto coordinate = Accumulator::unpack(key, keyShape);
      sm.coordinates.insert(sm.coordinates.end(), coordinate.begin(), coordinate.end());
      sm.values.push_back(weight);
    }
  }
  return sm;
}

template <std::size_t N, typename T, typename I, class Slicer, class Shape, typename U>
inline Signed_measure<I, T> _compute_hilbert_signed_measure_sparse_by_slices(
    const Line_sweep<Slicer, Shape, U> &sweep) {
  using Accumulator = Packed_hilbert_sparse_accumulator<N, T, I>;

  const auto keyShape = _hilbert_sparse_key_shape<N>(sweep.grid_shape(), static_cast<I>(sweep.degree_count()));

  const Accumulator prototype(keyShape);
  Thread_local_storage<Accumulator> threadSmPts(prototype);

  sweep.for_each_line([&](auto &lineSlicer, const std::array<I, N> &position) {
    Accumulator &smPts = threadSmPts.local();

    sweep.for_each_endpoint(lineSlicer, position, [&](std::size_t degreeIndex, std::size_t coordinate, int sign) {
      _emit_endpoint<N, T>(smPts, degreeIndex, position, sweep.barcode_axis(), coordinate, sign, sweep.grid_shape());
    });
  });

  return _merge_sparse_hilbert_accumulators<N, T, I>(threadSmPts, keyShape);
}

template <std::size_t N, typename T, typename I, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> _compute_sparse_hsm_for_n_param(Slicer &slicer, const std::vector<U1> &shape,
                                                            const std::vector<U2> &degrees, bool zeroPad,
                                                            bool ignoreInf, bool verbose, int n_jobs) {
  std::array<I, N> parameterShape;
  std::copy_n(shape.begin(), N, parameterShape.begin());
  detail::Line_sweep sweep(slicer, parameterShape, degrees, zeroPad, ignoreInf, n_jobs);

  print_progress_wo_line_break(verbose, "Computing sparse Hilbert invariant ... ");

  auto out = _compute_hilbert_signed_measure_sparse_by_slices<N, T, I>(sweep);

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
 * Default: 0.
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_hilbert_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                           const std::vector<U2> &degrees, bool zeroPad = false,
                                                           bool ignoreInf = true, bool verbose = false,
                                                           int n_jobs = 0) {
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
  for (auto r : tensor.extents()) print_progress_wo_line_break(verbose, r, ", ");
  print_progress_wo_line_break(verbose, "\nContainer size : ", tensor.size(), "\nComputing Hilbert invariant ...");

  const std::vector<I> gridShape(shape.begin(), shape.end()); // to have the right type.
  detail::Line_sweep sweep(slicer, gridShape, degrees, zeroPad, ignoreInf, n_jobs);

  detail::_get_hilbert_surface<I>(sweep, tensor);

  print_progress(verbose, "Done.");
  print_progress_wo_line_break(verbose, "Computing mobius inversion ...");

  // The inversion along the barcode axis was already done when writing the births and deaths.
  for (std::size_t axis = 0; axis < sweep.active_shape().size(); ++axis) {
    if (axis != sweep.barcode_axis()) differentiate(tensor, axis + 1, n_jobs);  // +1 for degree
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
 * @param slicer Slicer of the filtered complex. It is modified.
 * @param shape Size of the grid in each parameter. Its size has to be the number of parameters of the slicer.
 * @param degrees Homological degrees for which to compute the invariant.
 * @param zeroPad If true, the last grid point of each parameter is only used as zero padding. Default: false.
 * @param ignoreInf Ignore the infinite bars during the persistence computation. Default: true.
 * @param verbose Print progress messages. Default: false.
 * @param n_jobs Maximal number of threads (when compiled with TBB). Ignored by vine persistence, which is sequential.
 * Default: 0.
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_sparse_hilbert_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                                  const std::vector<U2> &degrees, bool zeroPad = false,
                                                                  bool ignoreInf = true, bool verbose = false,
                                                                  int n_jobs = 0) {
  GUDHI_CHECK(slicer.get_number_of_parameters() == shape.size(), std::invalid_argument("Invalid shape size."));

  if (degrees.size() == 0) return {};

  switch (slicer.get_number_of_parameters()) {
    case 2:
      return detail::_compute_sparse_hsm_for_n_param<2, T, I>(slicer, shape, degrees, zeroPad, ignoreInf, verbose,
                                                              n_jobs);
    case 3:
      return detail::_compute_sparse_hsm_for_n_param<3, T, I>(slicer, shape, degrees, zeroPad, ignoreInf, verbose,
                                                              n_jobs);
    case 4:
      return detail::_compute_sparse_hsm_for_n_param<4, T, I>(slicer, shape, degrees, zeroPad, ignoreInf, verbose,
                                                              n_jobs);
    default:
      throw std::runtime_error("Sparse Hilbert signed measure is implemented for 2, 3, and 4 parameters.");
  }
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_HILBERT_SIGNED_H_