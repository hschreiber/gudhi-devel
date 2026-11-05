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
 * @file rank_invariant.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Contains the methods @ref Gudhi::multi_persistence::compute_rank_invariant and
 * @ref Gudhi::multi_persistence::compute_sparse_rank_signed_measure.
 */

#ifndef GUDHI_MP_RANK_INVARIANT_H_
#define GUDHI_MP_RANK_INVARIANT_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/tensor_helpers.h>
#include <gudhi/Multi_persistence/slice_sweep.h>

namespace Gudhi {
namespace multi_persistence {

namespace detail {

/////////// COMMON ///////////

/**
 * @private
 * @brief Family of paths of the rank invariant. The rank of the map M(p) -> M(q) for p <= q is computed on a
 * monotone staircase path which goes through p and q: it starts at p, moves along the axis 0 until q_0, then along
 * the axis 1 until q_1, ..., and finally along the axis N - 1 until q_{N-1} (the axes can be taken in any other
 * order, see below). The path is described by the starting point of each of its segments (@c lower) and the end
 * point of each segment (@c upper). The first segment starts at lower[0] = 0, and the last one ends at
 * upper[N-1] = grid shape[N-1] - 1. There is exactly one path for each pair (p, q): p = (first point of the first
 * segment, lower[1], ..., lower[N-1]) and q = (upper[0], ..., upper[N-2], point of the last segment).
 *
 * The position of a point along the path is the sum of its coordinates. The pairs (p, q) of a path are the pairs
 * of a point of the first segment and a point of the last segment, they are "reached" by a bar of the path if
 * the bar contains the position of p and of q.
 *
 * The axes can be visited in any order, and the number of paths depends on it (each of the two end axes of the
 * order contributes a factor equal to its extent, each of the others a factor of e (e + 1) / 2). The order which
 * minimizes the number of paths is used. @c lower and @c upper are indexed by the position in this order, not by the
 * parameter.
 *
 * Models the concept of family of @ref Slice_sweep.
 */
template <std::size_t N, typename I>
class Slice_path_family {
  static_assert(N >= 2, "The rank invariant needs at least 2 parameters.");

 public:
  struct Path {
    std::array<I, N> lower{};
    std::array<I, N> upper{};
  };

  /** @brief Positions along a path of the segments which contain the birth points and the death points. */
  struct Bounds {
    std::size_t firstStart;  // first position of the first segment
    std::size_t firstStop;   // one after the last position of the first segment
    std::size_t lastStart;   // first position of the last segment
    std::size_t lastStop;    // one after the last position of the last segment
  };

  explicit Slice_path_family(const std::array<I, N> &gridShape) : gridShape_(gridShape), orderedShape_(gridShape) {
    std::iota(axisOrder_.begin(), axisOrder_.end(), std::size_t{0});
    if constexpr (N > 2) _choose_axis_order();
  }

  /** @brief Shape of the grid. */
  [[nodiscard]] const std::array<I, N> &grid_shape() const { return gridShape_; }

  /** @brief axis_order()[a] is the parameter which is visited at the step a of the paths. */
  [[nodiscard]] const std::array<std::size_t, N> &axis_order() const { return axisOrder_; }

  [[nodiscard]] std::size_t path_count() const { return _family_size(orderedShape_); }

  /**
   * @brief Returns the path number @p index. Paths are numbered like the digits of a mixed radix number: first
   * upper[0], then the pair (lower[a], upper[a]) of each middle step, and last lower[N - 1]. Must stay consistent
   * with @ref next_path.
   */
  [[nodiscard]] Path path_from_index(std::size_t index) const {
    Path path;
    for (std::size_t a = 0; a < N; ++a) path.upper[a] = static_cast<I>(orderedShape_[a] - 1);
    path.upper[0] = static_cast<I>(index % orderedShape_[0]);
    index /= orderedShape_[0];
    for (std::size_t a = 1; a + 1 < N; ++a) {
      const std::size_t extent = orderedShape_[a];
      const auto [lower, upper] = _pair_from_linear(index % (extent * (extent + 1) / 2), extent);
      path.lower[a] = static_cast<I>(lower);
      path.upper[a] = static_cast<I>(upper);
      index /= extent * (extent + 1) / 2;
    }
    path.lower[N - 1] = static_cast<I>(index % orderedShape_[N - 1]);
    return path;
  }

  /** @brief Advances @p path to the next path, like an odometer. Must stay consistent with @ref path_from_index. */
  void next_path(Path &path) const {
    if (++path.upper[0] < orderedShape_[0]) return;
    path.upper[0] = 0;
    for (std::size_t a = 1; a + 1 < N; ++a) {
      if (++path.upper[a] < orderedShape_[a]) return;
      if (++path.lower[a] < orderedShape_[a]) {
        path.upper[a] = path.lower[a];
        return;
      }
      path.lower[a] = 0;
      path.upper[a] = 0;
    }
    if (++path.lower[N - 1] < orderedShape_[N - 1]) return;
    path.lower[N - 1] = 0;
  }

  [[nodiscard]] Bounds bounds(const Path &path) const {
    Bounds bounds{0, 0, static_cast<std::size_t>(path.lower[N - 1]), static_cast<std::size_t>(orderedShape_[N - 1])};
    for (std::size_t a = 1; a < N; ++a) bounds.firstStart += path.lower[a];
    bounds.firstStop = bounds.firstStart + path.upper[0] + 1;
    for (std::size_t a = 0; a + 1 < N; ++a) {
      bounds.lastStart += path.upper[a];
      bounds.lastStop += path.upper[a];
    }
    return bounds;
  }

  /**
   * @brief Restricts the bar [@p birth, @p death) of the barcode of a path to the part which matters for the rank:
   * returns false if the bar is born after the first segment or dies before the last one, otherwise stores the
   * birth and the death as positions along the path, the death being clipped to the end of the last segment.
   * Born or dying at an infinite position is fine.
   */
  template <typename Value>
  static bool clip_bar(const Bounds &bounds, Value birth, Value death, std::size_t &clippedBirth,
                       std::size_t &clippedDeath) {
    // Compared before being converted: the values can be huge or infinite.
    if (!(birth < static_cast<Value>(bounds.firstStop))) return false;
    clippedBirth = static_cast<std::size_t>(birth);
    clippedDeath = static_cast<std::size_t>(std::min(death, static_cast<Value>(bounds.lastStop)));
    return clippedDeath > bounds.lastStart;
  }

  /**
   * @brief Returns the filtration value of a simplex restricted to @p path, that is the smallest position along the
   * path of its generators. Returns infinity if the simplex never appears on the path.
   */
  template <class Filtration_value>
  [[nodiscard]] typename Filtration_value::value_type slice_value(const Filtration_value &fil,
                                                                  const Path &path) const {
    typename Filtration_value::value_type best = Filtration_value::T_inf;
    for (unsigned int g = 0; g < fil.num_generators(); ++g) best = std::min(best, _generator_value(fil, g, path));
    return best;
  }

 private:
  std::array<I, N> gridShape_;
  std::array<I, N> orderedShape_;  // gridShape_ in the order of the paths
  std::array<std::size_t, N> axisOrder_;

  /** @brief Number of paths for a grid with the extents @p shape (in the order of the paths). */
  static std::size_t _family_size(const std::array<I, N> &shape) {
    std::size_t size = static_cast<std::size_t>(shape[0]) * shape[N - 1];
    for (std::size_t a = 1; a + 1 < N; ++a) {
      const std::size_t extent = shape[a];
      size *= extent * (extent + 1) / 2;
    }
    return size;
  }

  /** @brief Returns the @p linear -th pair (lower, upper) with lower <= upper < @p extent. */
  static std::pair<std::size_t, std::size_t> _pair_from_linear(std::size_t linear, std::size_t extent) {
    std::size_t lower = 0;
    while (linear >= extent - lower) {
      linear -= extent - lower;
      ++lower;
    }
    return {lower, lower + linear};
  }

  /** @brief Puts the two largest axes at the ends of the paths, and the others in the middle, if it is better. */
  void _choose_axis_order() {
    std::array<std::size_t, N> sortedAxes;
    std::iota(sortedAxes.begin(), sortedAxes.end(), std::size_t{0});
    std::sort(sortedAxes.begin(), sortedAxes.end(), [&](std::size_t left, std::size_t right) {
      return gridShape_[left] != gridShape_[right] ? gridShape_[left] > gridShape_[right] : left < right;
    });

    std::array<std::size_t, N> order;
    order[0] = sortedAxes[0];
    order[N - 1] = sortedAxes[1];
    for (std::size_t a = 1; a + 1 < N; ++a) order[a] = sortedAxes[a + 1];

    std::array<I, N> shape;
    for (std::size_t a = 0; a < N; ++a) shape[a] = gridShape_[order[a]];
    if (_family_size(shape) < _family_size(gridShape_)) {
      axisOrder_ = order;
      orderedShape_ = shape;
    }
  }

  /**
   * @brief Returns the position along @p path of the generator @p g of @p fil, infinity if it is not on the path.
   *
   * On the segment s, the path has finished the axes before s (they are at their upper value) and has not started
   * the axes after s (they are at their lower value). A generator is on the segment s if it appeared on all the other
   * axes by then, and if it appears on the axis s before the end of the segment; its position is then the position
   * of the path at the moment the generator appears on the axis s (or at the start of the segment if it appeared
   * before it).
   */
  template <class Filtration_value>
  typename Filtration_value::value_type _generator_value(const Filtration_value &fil, unsigned int g,
                                                         const Path &path) const {
    using Value = typename Filtration_value::value_type;

    for (std::size_t segment = 0; segment < N; ++segment) {
      bool onSegment = true;
      for (std::size_t a = 0; a < N && onSegment; ++a) {
        if (a == segment) continue;
        const Value bound = static_cast<Value>(a < segment ? path.upper[a] : path.lower[a]);
        onSegment = fil(g, axisOrder_[a]) <= bound;
      }
      if (!onSegment) continue;

      const Value coordinate = fil(g, axisOrder_[segment]);
      if (coordinate > static_cast<Value>(path.upper[segment])) continue;

      std::size_t position = std::max<std::size_t>(path.lower[segment], static_cast<std::size_t>(coordinate));
      for (std::size_t a = 0; a < segment; ++a) position += path.upper[a];
      for (std::size_t a = segment + 1; a < N; ++a) position += path.lower[a];
      return static_cast<Value>(position);
    }
    return Filtration_value::T_inf;
  }
};

/////////// INVARIANT ONLY ///////////

/**
 * @private
 * @brief Computes the rank invariant into @p out (assumed to be a zero tensor of shape
 * (degrees, grid shape, grid shape)), path by path: out[degree, p, q] is the rank of M(p) -> M(q) if p <= q.
 */
template <std::size_t N, typename T, typename I, class Slicer, typename U>
inline void _compute_rank_tensor(Slicer &slicer, const Tensor<T> &out, const std::array<I, N> &gridShape,
                                 const std::vector<U> &degrees, bool ignoreInf, int n_jobs) {
  const Slice_sweep sweep(slicer, Slice_path_family<N, I>(gridShape), degrees, ignoreInf, n_jobs);

  const auto &family = sweep.family();
  const auto &order = family.axis_order();

  sweep.for_each_path([&](auto &pathSlicer, const auto &path) {
    const auto bounds = family.bounds(path);

    // coordinates: degree index, then the point p, then the point q. Only the coordinates of p along the first
    // segment and of q along the last segment change.
    std::array<I, 1 + 2 * N> coordinates{};
    for (std::size_t a = 0; a < N; ++a) {
      coordinates[1 + order[a]] = path.lower[a];
      coordinates[1 + N + order[a]] = path.upper[a];
    }

    sweep.for_each_bar(pathSlicer, path, [&](std::size_t degreeIndex, auto birthValue, auto deathValue) {
      std::size_t birth, death;
      if (!family.clip_bar(bounds, birthValue, deathValue, birth, death)) return;

      coordinates[0] = static_cast<I>(degreeIndex);
      // The rank is increased for the pairs (p, q) such that birth <= position of p and position of q < death.
      for (std::size_t p = std::max(birth, bounds.firstStart); p < std::min(death, bounds.firstStop); ++p) {
        coordinates[1 + order[0]] = static_cast<I>(p - bounds.firstStart);
        for (std::size_t q = bounds.lastStart; q < death; ++q) {
          coordinates[1 + N + order[N - 1]] = static_cast<I>(q - bounds.lastStart + path.lower[N - 1]);
          out[coordinates]++;
        }
      }
    });
  });
}

/////////// SPARSE ONLY ///////////

template <std::size_t N, typename I>
inline std::array<I, 1 + 2 * N> _rank_sparse_key_shape(const std::array<I, N> &gridShape, I degreeCount,
                                                       bool zeroPad) {
  std::array<I, 1 + 2 * N> keyShape{};
  keyShape[0] = degreeCount;
  for (std::size_t axis = 0; axis < N; ++axis) {
    keyShape[1 + axis] = gridShape[axis];
    // the coordinates of the death points are shifted by one
    keyShape[1 + N + axis] = static_cast<I>(gridShape[axis] + (zeroPad ? 0 : 1));
  }
  return keyShape;
}

/**
 * @private
 * @brief Adds to @p acc the atoms of the signed measure which come from a bar of @p path.
 *
 * Restricted to the path, the bar adds 1 to the rank of the pairs (p, q) of a box of birth points times a box of
 * death points. The signed measure is the Möbius inversion of the rank invariant: along each axis, a box [low, high)
 * of birth points becomes +1 at low and -1 at high (if it is inside the grid), and a box [low, high) of death points
 * becomes +1 at high and -1 at low (if it is not 0). The death coordinates are therefore shifted by one: the atom
 * (p, q') stands for the pair (p, q' - 1), so only the atoms with p < q' are kept.
 */
template <std::size_t N, typename T, typename I, class Family, class Accumulator>
inline void _emit_rank_bar_atoms(Accumulator &acc, std::size_t degreeIndex, const Family &family,
                                 const typename Family::Path &path, const typename Family::Bounds &bounds,
                                 std::size_t birth, std::size_t death, bool zeroPad) {
  const auto &gridShape = family.grid_shape();
  const auto &order = family.axis_order();

  Atom_options<2 * N, T, I> options;
  for (std::size_t a = 0; a < N; ++a) {
    const std::size_t axis = order[a];

    // Boxes of birth and death points along this axis. They are single points, except for the birth points
    // along the first segment and the death points along the last one.
    std::size_t birthLow = path.lower[a];
    std::size_t birthHigh = path.lower[a] + 1;
    if (a == 0) {
      birthLow = std::max(birth, bounds.firstStart) - bounds.firstStart;
      birthHigh = std::min(death, bounds.firstStop) - bounds.firstStart;
    }
    std::size_t deathLow = path.upper[a];
    std::size_t deathHigh = path.upper[a] + 1;
    if (a == N - 1) {
      deathLow = path.lower[N - 1];
      deathHigh = path.lower[N - 1] + death - bounds.lastStart;
    }

    if (birthHigh <= birthLow) return;
    options.set_primary(axis, static_cast<I>(birthLow), static_cast<T>(1));
    if (birthHigh < static_cast<std::size_t>(gridShape[axis])) {
      options.set_secondary(axis, static_cast<I>(birthHigh), static_cast<T>(-1));
    }

    const std::size_t high = zeroPad ? std::min<std::size_t>(deathHigh, gridShape[axis] - 1) : deathHigh;
    if (high <= deathLow) return;
    options.set_primary(N + axis, static_cast<I>(high), static_cast<T>(1));
    if (deathLow > 0) options.set_secondary(N + axis, static_cast<I>(deathLow), static_cast<T>(-1));
  }

  _emit_atom_product(acc, degreeIndex, options, [](const auto &atom) {
    for (std::size_t axis = 0; axis < N; ++axis) {
      if (atom[1 + axis] >= atom[1 + N + axis]) return false;
    }
    return true;
  });
}

template <std::size_t N, typename T, typename I, class Slicer, typename U>
inline Signed_measure<I, T> _compute_sparse_rank_signed_measure(Slicer &slicer, const std::array<I, N> &gridShape,
                                                                const std::vector<U> &degrees, bool zeroPad,
                                                                bool ignoreInf, int n_jobs) {
  using Accumulator = Packed_sparse_accumulator<1 + 2 * N, T, I>;

  const Slice_sweep sweep(slicer, Slice_path_family<N, I>(gridShape), degrees, ignoreInf, n_jobs);
  const auto &family = sweep.family();

  const auto keyShape = _rank_sparse_key_shape<N>(gridShape, static_cast<I>(degrees.size()), zeroPad);
  const Accumulator prototype(keyShape);
  Thread_local_storage<Accumulator> threadAtoms(prototype);

  sweep.for_each_path([&](auto &pathSlicer, const auto &path) {
    Accumulator &acc = threadAtoms.local();
    const auto bounds = family.bounds(path);

    sweep.for_each_bar(pathSlicer, path, [&](std::size_t degreeIndex, auto birthValue, auto deathValue) {
      std::size_t birth, death;
      if (!family.clip_bar(bounds, birthValue, deathValue, birth, death)) return;
      _emit_rank_bar_atoms<N, T, I>(acc, degreeIndex, family, path, bounds, birth, death, zeroPad);
    });
  });

  return _merge_sparse_accumulators<1 + 2 * N, T, I>(threadAtoms, keyShape);
}

}  // namespace detail

/**
 * @brief Computes the signed measure of the rank invariant of the persistence module of @p slicer, that is the
 * Möbius inversion of the function which sends a pair of grid points p <= q to the rank of the map M(p) -> M(q).
 * Only implemented for 2, 3 and 4 parameters.
 *
 * The coordinates of an atom are (degree index, p_0, ..., p_{N-1}, q_0 + 1, ..., q_{N-1} + 1).
 *
 * @param slicer Slicer of the filtered complex. The persistence related storage will be modified.
 * @param shape Size of the grid in each parameter. Its size has to be the number of parameters of the slicer.
 * @param degrees Homological degrees for which to compute the invariant.
 * @param zeroPad If true, the coordinates of the death points are bounded by the last grid point (instead of one
 * after it).
 * @param ignoreInf Ignore the infinite bars during the persistence computation.
 * @param n_jobs Maximal number of threads (when compiled with TBB). Ignored by vine persistence, which is sequential.
 * Default: -1 (i.e. as much as possible).
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_sparse_rank_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                               const std::vector<U2> &degrees, bool zeroPad,
                                                               bool ignoreInf, int n_jobs = -1) {
  GUDHI_CHECK(shape.size() == slicer.get_number_of_parameters(), std::invalid_argument("Invalid rank grid shape."));

  if (degrees.size() == 0) return {};

  return detail::_dispatch_on_number_of_parameters(
      slicer.get_number_of_parameters(), "Sparse rank signed measure is implemented for 2, 3, and 4 parameters.",
      [&](auto numberOfParameters) {
        constexpr std::size_t N = decltype(numberOfParameters)::value;
        std::array<I, N> gridShape;
        std::copy_n(shape.begin(), N, gridShape.begin());
        return detail::_compute_sparse_rank_signed_measure<N, T>(slicer, gridShape, degrees, zeroPad, ignoreInf,
                                                                 n_jobs);
      });
}

/**
 * @brief Computes the rank invariant of the persistence module of @p slicer: the rank of the map M(p) -> M(q) for
 * every pair of grid points p <= q. Only implemented for 2, 3 and 4 parameters.
 *
 * @param slicer Slicer of the filtered complex. The persistence related storage will be modified.
 * @param gridShape Size of the grid in each parameter. Its size has to be the number of parameters of the slicer.
 * @param degrees Homological degrees for which to compute the invariant.
 * @param ignoreInf Ignore the infinite bars during the persistence computation.
 * @param n_jobs Maximal number of threads (when compiled with TBB). Ignored by vine persistence, which is sequential.
 * Default: -1 (i.e. as much as possible).
 * @return The values of the rank invariant, as a flat row-major tensor (0 for the pairs which are not p <= q), and
 * the shape of this tensor.
 */
template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline std::pair<std::vector<T>, std::vector<I>> compute_rank_invariant(Slicer &slicer,
                                                                        const std::vector<U1> &gridShape,
                                                                        const std::vector<U2> &degrees,
                                                                        bool ignoreInf, int n_jobs = -1) {
  GUDHI_CHECK(gridShape.size() == slicer.get_number_of_parameters(),
              std::invalid_argument("Invalid rank invariant tensor shape."));

  std::pair<std::vector<T>, std::vector<I>> out;
  out.second.reserve(1 + 2 * gridShape.size());
  out.second.push_back(static_cast<I>(degrees.size()));
  out.second.insert(out.second.end(), gridShape.begin(), gridShape.end());
  out.second.insert(out.second.end(), gridShape.begin(), gridShape.end());
  std::size_t total = 1;
  for (I value : out.second) total *= static_cast<std::size_t>(value);
  if (total == 0) return out;
  out.first.resize(total, 0);
  Tensor<T> tensor(out.first, out.second);

  detail::_dispatch_on_number_of_parameters(
      slicer.get_number_of_parameters(), "Rank invariant is implemented for 2, 3, and 4 parameters.",
      [&](auto numberOfParameters) {
        constexpr std::size_t N = decltype(numberOfParameters)::value;
        std::array<I, N> shape;
        std::copy_n(gridShape.begin(), N, shape.begin());
        detail::_compute_rank_tensor<N>(slicer, tensor, shape, degrees, ignoreInf, n_jobs);
      });
  return out;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_RANK_INVARIANT_H_