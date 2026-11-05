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
 * @file slice_sweep.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Tools shared by the invariants which are computed by restricting the module to a family of one parameter
 * slices of the grid (Hilbert function: lines, rank invariant: slice paths) and by summing the atoms of the
 * resulting signed measure (@ref Gudhi::multi_persistence::detail::Slice_sweep,
 * @ref Gudhi::multi_persistence::detail::Packed_sparse_accumulator, ...).
 */

#ifndef GUDHI_MP_SLICE_SWEEP_H_
#define GUDHI_MP_SLICE_SWEEP_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
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

namespace Gudhi {
namespace multi_persistence {

namespace detail {

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

/**
 * @private
 * @brief Calls @p fn with a std::integral_constant<std::size_t, N> for the value of N (2, 3 or 4) which is equal
 * to @p numberOfParameters, and returns what it returns. It lets the callers use the number of parameters as a
 * compile time constant. Throws std::runtime_error(@p errorMessage) for any other value.
 */
template <typename Fn>
inline decltype(auto) _dispatch_on_number_of_parameters(std::size_t numberOfParameters, const char *errorMessage,
                                                        Fn &&fn) {
  switch (numberOfParameters) {
    case 2:
      return fn(std::integral_constant<std::size_t, 2>{});
    case 3:
      return fn(std::integral_constant<std::size_t, 3>{});
    case 4:
      return fn(std::integral_constant<std::size_t, 4>{});
    default:
      throw std::runtime_error(errorMessage);
  }
}

/**
 * @private
 * @brief Holds what stays constant during the computation of an invariant which is obtained by restricting the
 * module to a family of one parameter slices ("paths") of the grid: the slicer, the family of paths, the degrees and
 * the options. It enumerates the paths, and the bars of the barcode of each path.
 *
 * @tparam Family Describes the paths, it has to provide:
 *  - `using Path = ...;` a value type which identifies one path.
 *  - `std::size_t path_count() const;` the number of paths.
 *  - `Path path_from_index(std::size_t index) const;` the path number @p index.
 *  - `void next_path(Path &path) const;` turns the path number i into the path number i + 1. Has to be consistent
 *    with @c path_from_index. Is never called on the last path.
 *  - `template <class Filtration_value> typename Filtration_value::value_type slice_value(const Filtration_value &,
 *    const Path &) const;` the filtration value of a simplex restricted to the path, `Filtration_value::T_inf` if it
 *    never appears on it.
 *
 * The object does not own the slicer or the degrees, they have to outlive it.
 * The persistence related storages in the slicer are modified.
 */
template <class Slicer, class Family, typename U>
class Slice_sweep {
 public:
  using Path = typename Family::Path;

  // TODO: slicer could be constructed only here from complex
  Slice_sweep(Slicer &slicer, Family family, const std::vector<U> &degrees, bool ignoreInf, int n_jobs)
      : slicer_(&slicer), family_(std::move(family)), degrees_(&degrees), ignoreInf_(ignoreInf), nJobs_(n_jobs) {}

  // The degrees are only referenced: refuse temporaries.
  Slice_sweep(Slicer &, Family, std::vector<U> &&, bool, int) = delete;

  [[nodiscard]] const Family &family() const { return family_; }

  [[nodiscard]] std::size_t degree_count() const { return degrees_->size(); }

  /**
   * @brief Calls @p on_path(slicer, path) for every path of the family. The slicer given to the callback is the
   * one to use to compute the barcode of that path (see @ref for_each_bar).
   *
   * With TBB, the paths are split into contiguous chunks which are processed in parallel, each thread working on
   * its own copy of the slicer. Otherwise, and when the persistence is a vine (which needs consecutive paths to be
   * close to each other), the paths are processed sequentially, in order, with the slicer itself.
   * Therefore, @p on_path has to be thread safe.
   */
  template <typename OnPath>
  void for_each_path(OnPath &&on_path) const {
    const std::size_t nbPaths = family_.path_count();
    if (nbPaths == 0) return;

    // Visits the paths [begin, end).
    auto visit = [&](auto &pathSlicer, std::size_t begin, std::size_t end) {
      Path path = family_.path_from_index(begin);
      for (std::size_t index = begin; index != end;) {
        on_path(pathSlicer, path);
        if (++index != end) family_.next_path(path);
      }
    };

#ifdef GUDHI_USE_TBB
    if (!Slicer::Persistence::is_vine && nJobs_ != 1) {
      using ThreadSafe = typename Slicer::Thread_safe;
      ThreadSafe prototype(*slicer_);
      tbb::enumerable_thread_specific<ThreadSafe> threadSlicers(prototype);

      // limits the number of threads (n_jobs <= 0: as many as possible)
      tbb::task_arena arena(nJobs_ > 0 ? nJobs_ : tbb::task_arena::automatic);
      arena.execute([&] {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nbPaths),
                          [&](const tbb::blocked_range<std::size_t> &range) {
                            visit(threadSlicers.local(), range.begin(), range.end());
                          });
      });
      return;
    }
#endif

    visit(*slicer_, 0, nbPaths);
  }

  /**
   * @brief Computes the barcode of @p path with @p pathSlicer and calls @p on_bar(degreeIndex, birth, death) for
   * every bar of the requested degrees. @p birth and @p death are the raw values of the slicer: they are positions
   * along the path, but can be larger than the path (or infinite for the death), it is up to the caller to clip them.
   * @p degreeIndex is the index of the degree in the degrees.
   */
  template <class PathSlicer, typename OnBar>
  void for_each_bar(PathSlicer &pathSlicer, const Path &path, OnBar &&on_bar) const {
    auto &slice = pathSlicer.get_slice();
    const auto &filtrationValues = pathSlicer.get_filtration_values();

    for (std::size_t k = 0; k < filtrationValues.size(); ++k) {
      slice[k] = family_.slice_value(filtrationValues[k], path);
    }

    if (!pathSlicer.persistence_computation_is_initialized()) {
      pathSlicer.initialize_persistence_computation(ignoreInf_);
    } else {
      pathSlicer.update_persistence_computation(ignoreInf_);
    }

    const auto barcodes = pathSlicer.template get_flat_barcode<true>();
    std::size_t degreeIndex = 0;
    for (U degree : *degrees_) {
      if (degree >= 0 && static_cast<std::size_t>(degree) < barcodes.size()) {
        for (const auto &bar : barcodes[degree]) on_bar(degreeIndex, bar[0], bar[1]);
      }
      ++degreeIndex;
    }
  }

 private:
  Slicer *slicer_;
  Family family_;
  const std::vector<U> *degrees_;
  bool ignoreInf_;
  int nJobs_;
};

/**
 * @private
 * @brief Sums signed weights of the atoms of a measure on a space with K coordinates (the first one being the
 * degree). Atoms are stored in a hash map under a 64-bit key, which packs the coordinates in mixed radix. Atoms
 * whose weights cancel out are removed.
 */
template <std::size_t K, typename T, typename I>
class Packed_sparse_accumulator {
 public:
  using Key = std::uint64_t;
  using Coordinate = std::array<I, K>;
  using Entry = std::pair<Key, T>;

  explicit Packed_sparse_accumulator(const Coordinate &key_shape) : keyShape_(key_shape) { validate_key_shape(); }

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
    for (std::size_t axis = K; axis-- > 0;) {
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
        throw std::runtime_error("Internal error: invalid sparse key shape.");
      }
      const auto extent = static_cast<Key>(extent_value);
      if (capacity > max_key / extent) {
        throw std::runtime_error("Sparse key space exceeds 64-bit capacity.");
      }
      capacity *= extent;
    }
  }

  [[nodiscard]] Key pack(const Coordinate &coordinates) const {
    Key key = 0;
    for (std::size_t axis = 0; axis < K; ++axis) {
      const auto extent = static_cast<Key>(keyShape_[axis]);
      const auto coordinate = static_cast<Key>(coordinates[axis]);
      GUDHI_CHECK(coordinate < extent, std::invalid_argument("Internal error: invalid sparse coordinate."));
      key = key * extent + coordinate;
    }
    return key;
  }

  Coordinate keyShape_{};
  std::unordered_map<Key, T> weights_;
};

/**
 * @private
 * @brief Merges the accumulators of all the threads into a signed measure, sorted by coordinates.
 */
template <std::size_t K, typename T, typename I>
inline Signed_measure<I, T> _merge_sparse_accumulators(
    const Thread_local_storage<Packed_sparse_accumulator<K, T, I>> &threadAccumulators,
    const std::array<I, K> &keyShape) {
  using Accumulator = Packed_sparse_accumulator<K, T, I>;

  std::size_t total = 0;
  for (const auto &local : threadAccumulators) total += local.size();

  std::vector<typename Accumulator::Entry> entries;
  entries.reserve(total);
  for (const auto &local : threadAccumulators) local.append_entries_to(entries);
  std::sort(entries.begin(), entries.end(),
            [](const auto &left, const auto &right) { return left.first < right.first; });

  Signed_measure<I, T> sm;
  sm.numberOfParameters = K;
  // TODO: would be better to not over-reserve as this will survive for a long time
  sm.coordinates.reserve(entries.size() * K);
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

/**
 * @private
 * @brief For each of K axes, one or two alternatives (a coordinate and a weight), the second one being optional.
 * Every axis has to get a primary alternative before it is used.
 */
template <std::size_t K, typename T, typename I>
struct Atom_options {
  std::array<std::array<I, 2>, K> coordinate{};
  std::array<std::array<T, 2>, K> weight{};
  std::array<std::size_t, K> count{};

  void set_primary(std::size_t axis, I coord, T w) {
    coordinate[axis][0] = coord;
    weight[axis][0] = w;
    count[axis] = 1;
  }

  void set_secondary(std::size_t axis, I coord, T w) {
    coordinate[axis][1] = coord;
    weight[axis][1] = w;
    count[axis] = 2;
  }
};

/**
 * @private
 * @brief Adds to @p acc the atoms (@p degreeIndex, c_0, ..., c_{K-1}) for every combination of one alternative per
 * axis of @p options for which @p accept(atom) is true. The weight of an atom is the product of the weights of the
 * chosen alternatives.
 *
 * A signed measure obtained by a Möbius inversion is a product over the axes of "a point and its neighbour with
 * the opposite sign": this is the sum of such a product.
 */
template <std::size_t K, typename T, typename I, class Accumulator, typename Accept>
inline void _emit_atom_product(Accumulator &acc, std::size_t degreeIndex, const Atom_options<K, T, I> &options,
                               Accept &&accept) {
  std::array<I, 1 + K> atom{};
  atom[0] = static_cast<I>(degreeIndex);
  std::array<std::size_t, K> choice{};  // mixed radix counter over the alternatives

  std::size_t axis = 0;
  while (axis < K) {
    T weight = 1;
    for (std::size_t ax = 0; ax < K; ++ax) {
      atom[1 + ax] = options.coordinate[ax][choice[ax]];
      weight = static_cast<T>(weight * options.weight[ax][choice[ax]]);
    }
    if (accept(atom)) acc.add(atom, weight);

    axis = 0;
    while (axis < K && ++choice[axis] == options.count[axis]) choice[axis++] = 0;
  }
}

}  // namespace detail

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_SLICE_SWEEP_H_
