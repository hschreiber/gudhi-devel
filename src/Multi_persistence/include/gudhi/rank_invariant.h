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
 * @brief Contains the methods @ref Gudhi::multi_persistence::compute_rank_invariant and
 * @ref Gudhi::multi_persistence::compute_sparse_rank_signed_measure.
 */

#ifndef GUDHI_MP_RANK_INVARIANT_H_
#define GUDHI_MP_RANK_INVARIANT_H_

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

template <std::size_t N, typename I>
inline std::array<I, 1 + 2 * N> _rank_sparse_key_shape(const std::array<I, N> &gridShape, const I degreeCount,
                                                       const bool zeroPad) {
  std::array<I, 1 + 2 * N> keyShape{};
  keyShape[0] = degreeCount;
  for (std::size_t axis = 0; axis < N; ++axis) {
    keyShape[1 + axis] = gridShape[axis];
    keyShape[1 + N + axis] = gridShape[axis] + (zeroPad ? I{0} : I{1});
  }
  return keyShape;
}

template <std::size_t N, typename T, typename I>
class Packed_rank_sparse_accumulator {
 public:
  using key_type = std::uint64_t;
  using coordinate_type = std::array<I, 1 + 2 * N>;

  explicit Packed_rank_sparse_accumulator(const coordinate_type &key_shape) : keyShape_(key_shape) {
    key_type capacity = 1;
    constexpr key_type max_key = std::numeric_limits<key_type>::max();
    for (const auto extent_value : keyShape_) {
      if (extent_value <= I{0}) [[unlikely]] {
        throw std::runtime_error("Internal error: invalid sparse rank key shape.");
      }
      const auto extent = static_cast<key_type>(extent_value);
      if (capacity > max_key / extent) [[unlikely]] {
        throw std::runtime_error("Sparse rank key space exceeds 64-bit capacity.");
      }
      capacity *= extent;
    }
    weights_.reserve(1024);
  }

  void reserve(const std::size_t size) { weights_.reserve(size); }

  [[nodiscard]] std::size_t size() const { return weights_.size(); }

  void add(const coordinate_type &coordinates, const T weight) { add_key(pack(coordinates), weight); }

  void merge_from(const Packed_rank_sparse_accumulator &other) {
    for (const auto &[key, weight] : other.weights_) add_key(key, weight);
  }

  [[nodiscard]] std::vector<key_type> sorted_keys() const {
    std::vector<key_type> keys;
    keys.reserve(weights_.size());
    for (const auto &[key, weight] : weights_) keys.push_back(key);
    std::sort(keys.begin(), keys.end());
    return keys;
  }

  [[nodiscard]] T weight(const key_type key) const { return weights_.at(key); }

  [[nodiscard]] coordinate_type unpack(key_type key) const {
    coordinate_type coordinates{};
    for (std::size_t axis = 1 + 2 * N; axis-- > 0;) {
      const auto extent = static_cast<key_type>(keyShape_[axis]);
      coordinates[axis] = static_cast<I>(key % extent);
      key /= extent;
    }
    return coordinates;
  }

 private:
  [[nodiscard]] key_type pack(const coordinate_type &coordinates) const {
    key_type key = 0;
    for (std::size_t axis = 0; axis < 1 + 2 * N; ++axis) {
      const auto extent = static_cast<key_type>(keyShape_[axis]);
      const auto coordinate = static_cast<key_type>(coordinates[axis]);
      if (coordinate >= extent) [[unlikely]] {
        throw std::runtime_error("Internal error: invalid sparse rank coordinate.");
      }
      key = key * extent + coordinate;
    }
    return key;
  }

  void add_key(const key_type key, const T weight) {
    if (weight == T{0}) return;
    auto [it, inserted] = weights_.emplace(key, weight);
    if (inserted) return;
    it->second += weight;
    if (it->second == T{0}) weights_.erase(it);
  }

  coordinate_type keyShape_{};
  std::unordered_map<key_type, T> weights_;
};

template <typename Idx, typename value_type>
inline value_type _get_slice_rank_filtration(const value_type x, const value_type y, const Idx I, const Idx J) {
  // Use F::T_inf instead
  if (x > static_cast<value_type>(I))
    return std::numeric_limits<value_type>::has_infinity ? std::numeric_limits<value_type>::infinity()
                                                         : std::numeric_limits<value_type>::max();
  if (y > static_cast<value_type>(J)) return I + static_cast<Idx>(y);
  return J + static_cast<Idx>(x);
}

template <typename Idx>
inline std::pair<Idx, Idx> _get_coordinates(Idx inSliceValue, Idx I, Idx J) {
  if (inSliceValue <= J) return {0, J};
  if (inSliceValue <= I + J) return {inSliceValue - J, J};
  return {I, inSliceValue - I};
}

template <typename Output, typename... Indices>
inline void _increment_output(const Output &out, Indices... coordinates) {
  using index_type = std::common_type_t<Indices...>;
  const std::array<index_type, sizeof...(Indices)> coordinates_array{coordinates...};
  out[coordinates_array]++;
}

template <typename T, class Extents, class LayoutPolicy, typename... Indices>
inline void _increment_output(const Gudhi::Simple_mdspan<T, Extents, LayoutPolicy> &out, Indices... coordinates) {
  using index_type = typename Gudhi::Simple_mdspan<T, Extents, LayoutPolicy>::index_type;
  out(static_cast<index_type>(coordinates)...)++;
}

template <std::size_t N, typename I>
inline std::array<I, 1 + 2 * N> _vector_to_full_rank_shape(const std::vector<I> &values) {
  if (values.size() != 1 + 2 * N) [[unlikely]] {
    throw std::runtime_error("Internal error: invalid fixed rank tensor shape dispatch.");
  }
  std::array<I, 1 + 2 * N> out{};
  std::copy_n(values.begin(), 1 + 2 * N, out.begin());
  return out;
}

template <std::size_t N, typename Output, typename I>
inline void _increment_rank_output(const Output &out, I degree_index, const std::array<I, N> &birth,
                                   const std::array<I, N> &death) {
  std::array<I, 1 + 2 * N> coordinates{};
  coordinates[0] = degree_index;
  for (std::size_t axis = 0; axis < N; ++axis) {
    coordinates[1 + axis] = birth[axis];
    coordinates[1 + N + axis] = death[axis];
  }
  out[coordinates]++;
}

template <std::size_t N, typename I>
inline I _monotone_pair_count(I extent) {
  return extent * (extent + 1) / 2;
}

template <typename I>
inline std::pair<I, I> _monotone_pair_from_linear(I linear, I extent) {
  for (I lower = 0; lower < extent; ++lower) {
    const I count = extent - lower;
    if (linear < count) return {lower, lower + linear};
    linear -= count;
  }
  return {extent - 1, extent - 1};
}

template <std::size_t N, typename I>
inline I _slice_path_family_size(const std::array<I, N> &gridShape) {
  I size = gridShape[0] * gridShape[N - 1];
  for (std::size_t axis = 1; axis + 1 < N; ++axis) size *= _monotone_pair_count<N>(gridShape[axis]);
  return size;
}

template <std::size_t N, typename I>
inline std::array<std::size_t, N> _sparse_slice_path_axis_order(const std::array<I, N> &gridShape) {
  std::array<std::size_t, N> identity{};
  for (std::size_t axis = 0; axis < N; ++axis) identity[axis] = axis;
  if constexpr (N > 2) {
    std::array<std::size_t, N> sortedAxes = identity;
    std::sort(sortedAxes.begin(), sortedAxes.end(), [&](std::size_t left, std::size_t right) {
      if (gridShape[left] != gridShape[right]) return gridShape[left] > gridShape[right];
      return left < right;
    });
    std::array<std::size_t, N> order{};
    order[0] = sortedAxes[0];
    order[N - 1] = sortedAxes[1];
    for (std::size_t axis = 1; axis + 1 < N; ++axis) order[axis] = sortedAxes[axis + 1];
    std::array<I, N> orderedShape{};
    for (std::size_t axis = 0; axis < N; ++axis) orderedShape[axis] = gridShape[order[axis]];
    if (_slice_path_family_size<N>(orderedShape) >= _slice_path_family_size<N>(gridShape)) return identity;
    return order;
  }
  return identity;
}

template <std::size_t N>
inline bool _axis_order_is_identity(const std::array<std::size_t, N> &axisOrder) {
  for (std::size_t axis = 0; axis < N; ++axis)
    if (axisOrder[axis] != axis) return false;
  return true;
}

template <std::size_t N, typename I>
inline std::array<I, N> _ordered_grid_shape(const std::array<I, N> &gridShape,
                                            const std::array<std::size_t, N> &axisOrder) {
  std::array<I, N> out{};
  for (std::size_t axis = 0; axis < N; ++axis) out[axis] = gridShape[axisOrder[axis]];
  return out;
}

template <std::size_t N, typename I>
inline void _slice_path_bounds_from_linear(I linear, const std::array<I, N> &gridShape, std::array<I, N> &lower,
                                           std::array<I, N> &upper) {
  lower.fill(I{0});
  upper = gridShape;
  for (auto &extent : upper) --extent;

  upper[0] = linear % gridShape[0];
  linear /= gridShape[0];

  for (std::size_t axis = 1; axis + 1 < N; ++axis) {
    const I count = _monotone_pair_count<N>(gridShape[axis]);
    const auto [axis_lower, axis_upper] = _monotone_pair_from_linear(linear % count, gridShape[axis]);
    linear /= count;
    lower[axis] = axis_lower;
    upper[axis] = axis_upper;
  }

  lower[N - 1] = linear % gridShape[N - 1];
}

template <std::size_t N, typename I>
inline I _path_scalar_min(const std::array<I, N> &lower) {
  I value = 0;
  for (std::size_t axis = 1; axis < N; ++axis) value += lower[axis];
  return value;
}

template <std::size_t N, typename I>
inline I _path_last_segment_start(const std::array<I, N> &upper, const std::array<I, N> &lower) {
  I value = lower[N - 1];
  for (std::size_t axis = 0; axis + 1 < N; ++axis) value += upper[axis];
  return value;
}

template <std::size_t N, typename I>
inline I _path_last_segment_stop(const std::array<I, N> &upper, const std::array<I, N> &gridShape) {
  I value = gridShape[N - 1];
  for (std::size_t axis = 0; axis + 1 < N; ++axis) value += upper[axis];
  return value;
}

template <std::size_t N, class Filtration_value, typename I>
inline typename Filtration_value::value_type _get_slice_rank_filtration_on_slice_path(
    const Filtration_value &filtration, unsigned int generator, const std::array<I, N> &lower,
    const std::array<I, N> &upper) {
  using value_type = typename Filtration_value::value_type;
  if constexpr (N == 2) {
    return _get_slice_rank_filtration<I, value_type>(filtration(generator, 0), filtration(generator, 1), upper[0],
                                                     lower[1]);
  }
  for (std::size_t segment = 0; segment < N; ++segment) {
    bool appearsOnSegment = true;
    for (std::size_t axis = 0; axis < segment; ++axis) {
      if (filtration(generator, static_cast<unsigned int>(axis)) > static_cast<value_type>(upper[axis])) {
        appearsOnSegment = false;
        break;
      }
    }
    if (!appearsOnSegment) continue;
    for (std::size_t axis = segment + 1; axis < N; ++axis) {
      if (filtration(generator, static_cast<unsigned int>(axis)) > static_cast<value_type>(lower[axis])) {
        appearsOnSegment = false;
        break;
      }
    }
    if (!appearsOnSegment) continue;
    const auto coordinate = static_cast<I>(filtration(generator, static_cast<unsigned int>(segment)));
    if (coordinate > upper[segment]) continue;
    I scalar = std::max(lower[segment], coordinate);
    for (std::size_t axis = 0; axis < segment; ++axis) scalar += upper[axis];
    for (std::size_t axis = segment + 1; axis < N; ++axis) scalar += lower[axis];
    return static_cast<value_type>(scalar);
  }
  return Filtration_value::T_inf;
}

template <std::size_t N, class Filtration_value, typename I>
inline typename Filtration_value::value_type _get_slice_rank_filtration_on_ordered_slice_path(
    const Filtration_value &filtration, unsigned int generator, const std::array<I, N> &lower,
    const std::array<I, N> &upper, const std::array<std::size_t, N> &axisOrder) {
  using value_type = typename Filtration_value::value_type;
  for (std::size_t segment = 0; segment < N; ++segment) {
    bool appearsOnSegment = true;
    for (std::size_t axis = 0; axis < segment; ++axis) {
      if (filtration(generator, static_cast<unsigned int>(axisOrder[axis])) > static_cast<value_type>(upper[axis])) {
        appearsOnSegment = false;
        break;
      }
    }
    if (!appearsOnSegment) continue;
    for (std::size_t axis = segment + 1; axis < N; ++axis) {
      if (filtration(generator, static_cast<unsigned int>(axisOrder[axis])) > static_cast<value_type>(lower[axis])) {
        appearsOnSegment = false;
        break;
      }
    }
    if (!appearsOnSegment) continue;
    const auto coordinate = static_cast<I>(filtration(generator, static_cast<unsigned int>(axisOrder[segment])));
    if (coordinate > upper[segment]) continue;
    I scalar = std::max(lower[segment], coordinate);
    for (std::size_t axis = 0; axis < segment; ++axis) scalar += upper[axis];
    for (std::size_t axis = segment + 1; axis < N; ++axis) scalar += lower[axis];
    return static_cast<value_type>(scalar);
  }
  return Filtration_value::T_inf;
}

template <std::size_t Axis, std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_endpoint_atom_combinations(Accumulator &atoms,
                                             const std::array<std::array<I, 2>, 2 * N> &endpointCoordinates,
                                             const std::array<std::array<T, 2>, 2 * N> &endpointWeights,
                                             const std::array<std::size_t, 2 * N> &endpointCounts,
                                             std::array<I, 1 + 2 * N> &coordinates, T weight) {
  if constexpr (Axis == 2 * N) {
    bool validRankCoordinate = true;
    for (std::size_t axis = 0; axis < N; ++axis) {
      if (coordinates[1 + axis] >= coordinates[1 + N + axis]) {
        validRankCoordinate = false;
        break;
      }
    }
    if (validRankCoordinate) atoms.add(coordinates, weight);
  } else {
    for (std::size_t choice = 0; choice < endpointCounts[Axis]; ++choice) {
      coordinates[1 + Axis] = endpointCoordinates[Axis][choice];
      _emit_endpoint_atom_combinations<Axis + 1, N>(atoms, endpointCoordinates, endpointWeights, endpointCounts,
                                                    coordinates,
                                                    static_cast<T>(weight * endpointWeights[Axis][choice]));
    }
  }
}

template <std::size_t N, typename I, typename Output>
inline void _add_slice_path_bar_contribution(const Output &out, I degreeIndex, I birth, I death,
                                             const std::array<I, N> &lower, const std::array<I, N> &upper,
                                             const std::array<I, N> &gridShape) {
  const I firstStart = _path_scalar_min<N>(lower);
  const I firstStop = firstStart + upper[0] + 1;
  const I lastStart = _path_last_segment_start<N>(upper, lower);
  const I lastStop = _path_last_segment_stop<N>(upper, gridShape);
  if (birth >= firstStop || death <= lastStart) return;

  const I birthStop = std::min(death, firstStop);
  const I deathStop = std::min(death, lastStop);
  std::array<I, N> birthCoordinate = lower;
  std::array<I, N> deathCoordinate = upper;
  for (I birthScalar = std::max(birth, firstStart); birthScalar < birthStop; ++birthScalar) {
    birthCoordinate[0] = birthScalar - firstStart;
    for (I deathScalar = lastStart; deathScalar < deathStop; ++deathScalar) {
      deathCoordinate[N - 1] = deathScalar - (lastStart - lower[N - 1]);
      _increment_rank_output<N>(out, degreeIndex, birthCoordinate, deathCoordinate);
    }
  }
}

template <std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_interval_endpoint_atoms(Accumulator &atoms, I degree_index, const std::array<I, N> &birthLow,
                                          const std::array<I, N> &birthHigh, const std::array<I, N> &deathLow,
                                          const std::array<I, N> &deathHigh, const std::array<I, N> &gridShape,
                                          bool zeroPad) {
  std::array<std::array<I, 2>, 2 * N> endpointCoordinates{};
  std::array<std::array<T, 2>, 2 * N> endpointWeights{};
  std::array<std::size_t, 2 * N> endpointCounts{};

  for (std::size_t axis = 0; axis < N; ++axis) {
    if (birthHigh[axis] <= birthLow[axis]) return;
    endpointCoordinates[axis][0] = birthLow[axis];
    endpointWeights[axis][0] = T{1};
    endpointCounts[axis] = 1;
    if (birthHigh[axis] < gridShape[axis]) {
      endpointCoordinates[axis][1] = birthHigh[axis];
      endpointWeights[axis][1] = T{-1};
      endpointCounts[axis] = 2;
    }
  }

  for (std::size_t axis = 0; axis < N; ++axis) {
    const auto high = zeroPad ? std::min(deathHigh[axis], gridShape[axis] - 1) : deathHigh[axis];
    if (high <= deathLow[axis]) return;
    const auto endpointAxis = N + axis;
    endpointCoordinates[endpointAxis][0] = high;
    endpointWeights[endpointAxis][0] = T{1};
    endpointCounts[endpointAxis] = 1;
    if (deathLow[axis] > 0) {
      endpointCoordinates[endpointAxis][1] = deathLow[axis];
      endpointWeights[endpointAxis][1] = T{-1};
      endpointCounts[endpointAxis] = 2;
    }
  }

  std::array<I, 1 + 2 * N> coordinates{};
  coordinates[0] = degree_index;
  _emit_endpoint_atom_combinations<0, N>(atoms, endpointCoordinates, endpointWeights, endpointCounts, coordinates,
                                         T{1});
}

template <std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_slice_path_bar_signed_measure_atoms(Accumulator &atoms, I degreeIndex, I birth, I death,
                                                      const std::array<I, N> &lower, const std::array<I, N> &upper,
                                                      const std::array<I, N> &gridShape, bool zeroPad) {
  if constexpr (N == 2) {
    const I firstStart = lower[1];
    const I firstStop = firstStart + upper[0] + 1;
    const I lastStart = upper[0] + lower[1];
    const I lastStop = upper[0] + gridShape[1];
    if (birth >= firstStop || death <= lastStart) return;

    const I birthStop = std::min(death, firstStop);
    const I deathStop = std::min(death, lastStop);
    const I birth_x_low = std::max(birth, firstStart) - firstStart;
    const I birth_x_high = birthStop - firstStart;
    const I birth_y_low = lower[1];
    const I birth_y_high = lower[1] + 1;
    const I death_x_low = upper[0];
    const I death_x_high = upper[0] + 1;
    const I death_y_low = lower[1];
    const I death_y_high = lower[1] + deathStop - lastStart;

    if (birth_x_high <= birth_x_low || death_y_high <= death_y_low) return;

    std::array<std::array<I, 2>, 3> endpointCoordinates{};
    std::array<std::array<T, 2>, 3> endpointWeights{};
    std::array<std::size_t, 3> endpointCounts{};

    endpointCoordinates[0][0] = birth_x_low;
    endpointWeights[0][0] = T{1};
    endpointCounts[0] = 1;
    if (birth_x_high < gridShape[0]) {
      endpointCoordinates[0][1] = birth_x_high;
      endpointWeights[0][1] = T{-1};
      endpointCounts[0] = 2;
    }

    endpointCoordinates[1][0] = birth_y_low;
    endpointWeights[1][0] = T{1};
    endpointCounts[1] = 1;
    if (birth_y_high < gridShape[1]) {
      endpointCoordinates[1][1] = birth_y_high;
      endpointWeights[1][1] = T{-1};
      endpointCounts[1] = 2;
    }

    const I death_x_high_endpoint = zeroPad ? std::min(death_x_high, gridShape[0] - 1) : death_x_high;
    if (death_x_high_endpoint <= death_x_low) return;
    endpointCoordinates[2][0] = death_x_high_endpoint;
    endpointWeights[2][0] = T{1};
    endpointCounts[2] = 1;
    if (death_x_low > 0) {
      endpointCoordinates[2][1] = death_x_low;
      endpointWeights[2][1] = T{-1};
      endpointCounts[2] = 2;
    }

    const I death_y_high_endpoint = zeroPad ? std::min(death_y_high, gridShape[1] - 1) : death_y_high;
    if (death_y_high_endpoint <= death_y_low) return;
    std::array<I, 5> coordinates{};
    coordinates[0] = degreeIndex;
    coordinates[4] = death_y_high_endpoint;
    for (std::size_t bx = 0; bx < endpointCounts[0]; ++bx) {
      coordinates[1] = endpointCoordinates[0][bx];
      const T bx_weight = endpointWeights[0][bx];
      for (std::size_t by = 0; by < endpointCounts[1]; ++by) {
        coordinates[2] = endpointCoordinates[1][by];
        const T birth_weight = static_cast<T>(bx_weight * endpointWeights[1][by]);
        for (std::size_t dx = 0; dx < endpointCounts[2]; ++dx) {
          coordinates[3] = endpointCoordinates[2][dx];
          const T death_x_weight = static_cast<T>(birth_weight * endpointWeights[2][dx]);
          if (coordinates[1] < coordinates[3] && coordinates[2] < coordinates[4]) {
            atoms.add(coordinates, death_x_weight);
          }
        }
      }
    }
    return;
  }

  const I first_start = _path_scalar_min<N>(lower);
  const I first_stop = first_start + upper[0] + 1;
  const I last_start = _path_last_segment_start<N>(upper, lower);
  const I last_stop = _path_last_segment_stop<N>(upper, gridShape);
  if (birth >= first_stop || death <= last_start) return;

  const I birth_stop = std::min(death, first_stop);
  const I death_stop = std::min(death, last_stop);

  std::array<I, N> birth_low = lower;
  std::array<I, N> birth_high = lower;
  std::array<I, N> death_low = upper;
  std::array<I, N> death_high = upper;

  birth_low[0] = std::max(birth, first_start) - first_start;
  birth_high[0] = birth_stop - first_start;
  for (std::size_t axis = 1; axis < N; ++axis) birth_high[axis] = lower[axis] + 1;

  for (std::size_t axis = 0; axis + 1 < N; ++axis) death_high[axis] = upper[axis] + 1;
  death_low[N - 1] = lower[N - 1];
  death_high[N - 1] = lower[N - 1] + death_stop - last_start;

  _emit_interval_endpoint_atoms<N, T, I, Accumulator>(atoms, degreeIndex, birth_low, birth_high, death_low, death_high,
                                                      gridShape, zeroPad);
}

template <std::size_t N, typename T, typename I, typename Accumulator>
inline void _emit_ordered_slice_path_bar_signed_measure_atoms(
    Accumulator &atoms, I degreeIndex, I birth, I death, const std::array<I, N> &lower, const std::array<I, N> &upper,
    const std::array<I, N> &gridShape, const std::array<I, N> &orderedShape,
    const std::array<std::size_t, N> &axisOrder, bool zeroPad) {
  const I firstStart = _path_scalar_min<N>(lower);
  const I firstStop = firstStart + upper[0] + 1;
  const I lastStart = _path_last_segment_start<N>(upper, lower);
  const I lastStop = _path_last_segment_stop<N>(upper, orderedShape);
  if (birth >= firstStop || death <= lastStart) return;

  const I birthStop = std::min(death, firstStop);
  const I deathStop = std::min(death, lastStop);

  std::array<I, N> birth_low{};
  std::array<I, N> birth_high{};
  std::array<I, N> death_low{};
  std::array<I, N> death_high{};
  for (std::size_t axis = 0; axis < N; ++axis) {
    const auto original_axis = axisOrder[axis];
    birth_low[original_axis] = lower[axis];
    birth_high[original_axis] = lower[axis] + 1;
    death_low[original_axis] = upper[axis];
    death_high[original_axis] = upper[axis] + 1;
  }

  birth_low[axisOrder[0]] = std::max(birth, firstStart) - firstStart;
  birth_high[axisOrder[0]] = birthStop - firstStart;

  death_low[axisOrder[N - 1]] = lower[N - 1];
  death_high[axisOrder[N - 1]] = lower[N - 1] + deathStop - lastStart;

  _emit_interval_endpoint_atoms<N, T, I, Accumulator>(atoms, degreeIndex, birth_low, birth_high, death_low, death_high,
                                                      gridShape, zeroPad);
}

template <typename Idx, typename Output>
inline void _add_bar_contribution(const Output &out, Idx degreeIndex, Idx birth, Idx death, Idx I, Idx J,
                                  bool flipDeath) {
  const Idx corner = I + J;
  if (birth > corner || death <= corner) return;

  const Idx lastBirth = std::min(death, corner + 1);
  if (flipDeath) {
    for (Idx intermediateBirth = birth; intermediateBirth < lastBirth; ++intermediateBirth) {
      const auto [i, j] = _get_coordinates(intermediateBirth, I, J);
      for (Idx l = J; l < death - I; ++l) {
        _increment_output(out, degreeIndex, i, j, I - 1 - I, J - 1 - l);
      }
    }
  } else {
    for (Idx intermediateBirth = birth; intermediateBirth < lastBirth; ++intermediateBirth) {
      const auto [i, j] = _get_coordinates(intermediateBirth, I, J);
      for (Idx l = J; l < death - I; ++l) {
        _increment_output(out, degreeIndex, i, j, I, l);
      }
    }
  }
}

template <class Slicer, typename Idx, typename Output>
inline void _compute_2d_rank_invariant_of_elbow(typename Slicer::Thread_safe &slicer, const Output &out, const Idx I,
                                                const Idx J, const std::array<Idx, 5> &gridShape,
                                                const std::vector<Idx> &degrees, const bool flipDeath = false,
                                                const bool ignoreInf = true) {
  using value_type = typename Slicer::value_type;
  using F = typename Slicer::Filtration_value;
  const auto &filtrationsValues = slicer.get_filtration_values();
  auto numGenerators = filtrationsValues.size();
  const auto Y = gridShape[2];

  for (auto i = 0u; i < numGenerators; ++i) {
    const auto &f = filtrationsValues[i];
    value_type filtrationInSlice = F::T_inf;
    for (unsigned int g = 0; g < f.num_generators(); ++g) {
      value_type x = f(g, 0);
      value_type y = f(g, 1);

      filtrationInSlice = std::min(filtrationInSlice, _get_slice_rank_filtration(x, y, I, J));
    }
    slicer.get_slice()[i] = filtrationInSlice;
  }

  Idx degreeIndex = 0;

  if (!slicer.persistence_computation_is_initialized()) [[unlikely]] {
    slicer.initialize_persistence_computation(ignoreInf);
  } else {
    slicer.update_persistence_computation(ignoreInf);
  }
  auto barcodes = slicer.template get_flat_barcode<true>();

  // note one_pers not necessary when vine, but does the same computation

  for (auto degree : degrees) {
    // this assumes barcodes degrees starts from 0
    if (degree >= Idx{0} && static_cast<std::size_t>(degree) < barcodes.size()) {
      const auto &barcode = barcodes[degree];
      for (const auto &bar : barcode) {
        if (bar[0] > Y + I) continue;
        auto birth = static_cast<Idx>(bar[0]);
        auto death = static_cast<Idx>(std::min(bar[1],
                                               static_cast<value_type>(Y + I)));  // I,J reached, not X nor Y
        _add_bar_contribution(out, degreeIndex, birth, death, I, J, flipDeath);
      }
    }
    degreeIndex++;
  }
};

template <class Slicer, typename Idx, typename Output>
inline void _compute_2d_rank_invariant(Slicer &slicer, const Output &out, const std::array<Idx, 5> &gridShape,
                                       const std::vector<Idx> &degrees, const bool flipDeath, const bool ignoreInf) {
  if (degrees.size() == 0) return;
  Idx X = gridShape[1];
  Idx Y = gridShape[2];  // First axis is degree

  using ThreadSafe = typename Slicer::Thread_safe;
  ThreadSafe slicer_thread(slicer);
  tbb::enumerable_thread_specific<ThreadSafe> thread_locals(slicer_thread);
  tbb::parallel_for(0, X, [&](Idx I) {
    tbb::parallel_for(0, Y, [&](Idx J) {
      ThreadSafe &slicer = thread_locals.local();
      _compute_2d_rank_invariant_of_elbow<Slicer>(slicer, out, I, J, gridShape, degrees, flipDeath, ignoreInf);
    });
  });
}

template <std::size_t N, class Slicer, typename I, typename Output>
inline void _compute_nd_rank_invariant_by_slice_paths(Slicer &slicer, const Output &out,
                                                      const std::array<I, N> &gridShape, const std::vector<I> &degrees,
                                                      const bool ignoreInf) {
  using value_type = typename Slicer::value_type;
  using ThreadSafe = typename Slicer::Thread_safe;
  using F = typename Slicer::Filtration_value;

  if (degrees.size() == 0) return;
  const I totalSlicePaths = _slice_path_family_size<N>(gridShape);
  ThreadSafe slicerThread(slicer);
  tbb::enumerable_thread_specific<ThreadSafe> threadLocals(slicerThread);

  tbb::parallel_for(I{0}, totalSlicePaths, [&](I slice_path_index) {
    std::array<I, N> lower{};
    std::array<I, N> upper{};
    _slice_path_bounds_from_linear<N>(slice_path_index, gridShape, lower, upper);

    ThreadSafe &localSlicer = threadLocals.local();
    auto &sliceFiltration = localSlicer.get_slice();
    const auto &filtrationValues = localSlicer.get_filtration_values();
    for (std::size_t i = 0; i < filtrationValues.size(); ++i) {
      const auto &filtration = filtrationValues[i];
      value_type filtrationInSlice = F::T_inf;
      for (unsigned int generator = 0; generator < filtration.num_generators(); ++generator) {
        filtrationInSlice = std::min(filtrationInSlice,
                                     _get_slice_rank_filtration_on_slice_path<N>(filtration, generator, lower, upper));
      }
      sliceFiltration[i] = filtrationInSlice;
    }

    if (!localSlicer.persistence_computation_is_initialized()) [[unlikely]] {
      localSlicer.initialize_persistence_computation(ignoreInf);
    } else {
      localSlicer.update_persistence_computation(ignoreInf);
    }
    auto barcodes = localSlicer.template get_flat_barcode<true>();

    I degree_index = 0;
    for (auto degree : degrees) {
      if (degree >= I{0} && static_cast<std::size_t>(degree) < barcodes.size()) {
        const auto &barcode = barcodes[degree];
        for (const auto &bar : barcode) {
          const auto birth = static_cast<I>(bar[0]);
          const auto death =
              static_cast<I>(std::min(bar[1], static_cast<value_type>(_path_last_segment_stop<N>(upper, gridShape))));
          _add_slice_path_bar_contribution<N>(out, degree_index, birth, death, lower, upper, gridShape);
        }
      }
      ++degree_index;
    }
  });
}

template <std::size_t N, typename T, typename I>
inline Signed_measure<I, T> _merge_sparse_rank_accumulators(
    tbb::enumerable_thread_specific<Packed_rank_sparse_accumulator<N, T, I>> &threadAtoms,
    const std::array<I, 1 + 2 * N> &keyShape) {
  Packed_rank_sparse_accumulator<N, T, I> merged(keyShape);
  std::size_t total = 0;
  for (const auto &localAtoms : threadAtoms) total += localAtoms.size();
  merged.reserve(total);
  for (const auto &localAtoms : threadAtoms) merged.merge_from(localAtoms);

  const auto keys = merged.sorted_keys();

  Signed_measure<I, T> sm;
  sm.numberOfParameters = 1 + 2 * N;
  sm.coordinates.reserve(keys.size() * sm.numberOfParameters);
  sm.values.reserve(keys.size());

  for (const auto key : keys) {
    const auto coordinate = merged.unpack(key);
    const auto weight = merged.weight(key);
    sm.coordinates.insert(sm.coordinates.end(), coordinate.begin(), coordinate.end());
    sm.values.push_back(weight);
  }
  return sm;
}

template <std::size_t N, typename T, typename Idx, class Slicer>
inline Signed_measure<Idx, T> _compute_rank_signed_measure_by_slice_paths(Slicer &slicer,
                                                                          const std::array<Idx, N> &gridShape,
                                                                          const std::vector<Idx> &degrees,
                                                                          const bool zeroPad, const bool ignoreInf) {
  using value_type = typename Slicer::value_type;
  using ThreadSafe = typename Slicer::Thread_safe;
  using F = typename Slicer::Filtration_value;

  if (degrees.size() == 0) return {};

  if constexpr (N == 2) {
    const auto keyShape = _rank_sparse_key_shape<N>(gridShape, static_cast<Idx>(degrees.size()), zeroPad);
    ThreadSafe slicerThread(slicer);
    tbb::enumerable_thread_specific<ThreadSafe> threadLocals(slicerThread);
    Packed_rank_sparse_accumulator<N, T, Idx> atomAccumulator(keyShape);
    tbb::enumerable_thread_specific<Packed_rank_sparse_accumulator<N, T, Idx>> thread_atoms(atomAccumulator);

    const Idx X = gridShape[0];
    const Idx Y = gridShape[1];
    tbb::parallel_for(Idx{0}, X * Y, [&](Idx slicePathIndex) {
      const Idx I = slicePathIndex % X;
      const Idx J = slicePathIndex / X;
      const std::array<Idx, N> lower{Idx{0}, J};
      const std::array<Idx, N> upper{I, static_cast<Idx>(Y - 1)};

      ThreadSafe &localSlicer = threadLocals.local();
      auto &slice = localSlicer.get_slice();
      const auto &filtrationValues = localSlicer.get_filtration_values();
      for (std::size_t i = 0; i < filtrationValues.size(); ++i) {
        const auto &filtration = filtrationValues[i];
        value_type filtrationInSlice = F::T_inf;
        for (unsigned int generator = 0; generator < filtration.num_generators(); ++generator) {
          filtrationInSlice = std::min(
              filtrationInSlice,
              _get_slice_rank_filtration<Idx, value_type>(filtration(generator, 0), filtration(generator, 1), I, J));
        }
        slice[i] = filtrationInSlice;
      }

      if (!localSlicer.persistence_computation_is_initialized()) [[unlikely]] {
        localSlicer.initialize_persistence_computation(ignoreInf);
      } else {
        localSlicer.update_persistence_computation(ignoreInf);
      }
      auto barcodes = localSlicer.template get_flat_barcode<true>();

      auto &smPts = thread_atoms.local();
      Idx degreeIndex = 0;
      for (auto degree : degrees) {
        if (degree >= Idx{0} && static_cast<std::size_t>(degree) < barcodes.size()) {
          const auto &barcode = barcodes[degree];
          for (const auto &bar : barcode) {
            const auto birth = static_cast<Idx>(bar[0]);
            const auto death = static_cast<Idx>(std::min(bar[1], static_cast<value_type>(I + Y)));
            _emit_slice_path_bar_signed_measure_atoms<N, T>(smPts, degreeIndex, birth, death, lower, upper, gridShape,
                                                            zeroPad);
          }
        }
        ++degreeIndex;
      }
    });

    return _merge_sparse_rank_accumulators<N, T, Idx>(thread_atoms, keyShape);
  }

  const auto axisOrder = _sparse_slice_path_axis_order<N>(gridShape);
  const bool reorderAxes = !_axis_order_is_identity<N>(axisOrder);
  const auto orderedShape = _ordered_grid_shape<N>(gridShape, axisOrder);
  const auto slicePathGridShape = reorderAxes ? orderedShape : gridShape;
  const Idx totalSlicePaths = _slice_path_family_size<N>(slicePathGridShape);
  const auto keyShape = _rank_sparse_key_shape<N>(gridShape, static_cast<Idx>(degrees.size()), zeroPad);

  ThreadSafe slicerThread(slicer);
  tbb::enumerable_thread_specific<ThreadSafe> threadLocals(slicerThread);
  Packed_rank_sparse_accumulator<N, T, Idx> atomAccumulator(keyShape);
  tbb::enumerable_thread_specific<Packed_rank_sparse_accumulator<N, T, Idx>> threadAtoms(atomAccumulator);

  tbb::parallel_for(Idx{0}, totalSlicePaths, [&](Idx slice_path_index) {
    std::array<Idx, N> lower{};
    std::array<Idx, N> upper{};
    _slice_path_bounds_from_linear<N>(slice_path_index, slicePathGridShape, lower, upper);

    ThreadSafe &localSlicer = threadLocals.local();
    auto &slice = localSlicer.get_slice();
    const auto &filtrationValues = localSlicer.get_filtration_values();
    for (std::size_t i = 0; i < filtrationValues.size(); ++i) {
      const auto &filtration = filtrationValues[i];
      value_type filtrationInSlice = F::T_inf;
      for (unsigned int generator = 0; generator < filtration.num_generators(); ++generator) {
        filtrationInSlice = std::min(
            filtrationInSlice,
            reorderAxes
                ? _get_slice_rank_filtration_on_ordered_slice_path<N>(filtration, generator, lower, upper, axisOrder)
                : _get_slice_rank_filtration_on_slice_path<N>(filtration, generator, lower, upper));
      }
      slice[i] = filtrationInSlice;
    }

    if (!localSlicer.persistence_computation_is_initialized()) [[unlikely]] {
      localSlicer.initialize_persistence_computation(ignoreInf);
    } else {
      localSlicer.update_persistence_computation(ignoreInf);
    }
    auto barcodes = localSlicer.template get_flat_barcode<true>();

    auto &smPts = threadAtoms.local();
    Idx degreeIndex = 0;
    for (auto degree : degrees) {
      if (degree >= Idx{0} && static_cast<std::size_t>(degree) < barcodes.size()) {
        const auto &barcode = barcodes[degree];
        for (const auto &bar : barcode) {
          const auto birth = static_cast<Idx>(bar[0]);
          const auto death = static_cast<Idx>(
              std::min(bar[1], static_cast<value_type>(_path_last_segment_stop<N>(upper, slicePathGridShape))));
          if (reorderAxes) {
            _emit_ordered_slice_path_bar_signed_measure_atoms<N, T>(smPts, degreeIndex, birth, death, lower, upper,
                                                                    gridShape, slicePathGridShape, axisOrder, zeroPad);
          } else {
            _emit_slice_path_bar_signed_measure_atoms<N, T>(smPts, degreeIndex, birth, death, lower, upper, gridShape,
                                                            zeroPad);
          }
        }
      }
      ++degreeIndex;
    }
  });

  return _merge_sparse_rank_accumulators<N, T, Idx>(threadAtoms, keyShape);
}

}  // namespace detail

template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline Signed_measure<I, T> compute_sparse_rank_signed_measure(Slicer &slicer, const std::vector<U1> &shape,
                                                               const std::vector<U2> &degrees, const bool zeroPad,
                                                               const bool ignoreInf, int n_jobs = -1) {
  GUDHI_CHECK(shape.size() == slicer.get_number_of_parameters(), std::invalid_argument("Invalid rank grid shape."));
  if (degrees.size() == 0) return {};
  tbb::task_arena arena(Slicer::Persistence::is_vine ? 1 : n_jobs);
  switch (slicer.get_number_of_parameters()) {
    case 2: {
      const std::array<I, 2> parameter_shape = {shape[0], shape[1]};
      Signed_measure<I, T> out;
      arena.execute([&] {
        out = detail::_compute_rank_signed_measure_by_slice_paths<2, T, I>(slicer, parameter_shape, degrees, zeroPad,
                                                                           ignoreInf);
      });
      return out;
    }
    case 3: {
      const std::array<I, 3> parameter_shape = {shape[0], shape[1], shape[2]};
      Signed_measure<I, T> out;
      arena.execute([&] {
        out = detail::_compute_rank_signed_measure_by_slice_paths<3, T, I>(slicer, parameter_shape, degrees, zeroPad,
                                                                           ignoreInf);
      });
      return out;
    }
    case 4: {
      const std::array<I, 4> parameter_shape = {shape[0], shape[1], shape[2], shape[3]};
      Signed_measure<I, T> out;
      arena.execute([&] {
        out = detail::_compute_rank_signed_measure_by_slice_paths<4, T, I>(slicer, parameter_shape, degrees, zeroPad,
                                                                           ignoreInf);
      });
      return out;
    }
    default:
      throw std::runtime_error("Sparse rank signed measure is implemented for 2, 3, and 4 parameters.");
  }
}

template <typename T = int, typename I = std::uint16_t, class Slicer, typename U1, typename U2>
inline std::pair<std::vector<T>, std::vector<I>> compute_rank_invariant(Slicer &slicer,
                                                                        const std::vector<U1> &gridShape,
                                                                        const std::vector<U2> &degrees, bool ignoreInf,
                                                                        I n_jobs = -1) {
  const auto num_parameters = slicer.get_number_of_parameters();
  GUDHI_CHECK(gridShape.size() == num_parameters, std::invalid_argument("Invalid rank invariant tensor shape."));

  std::pair<std::vector<T>, std::vector<I>> out;
  out.second.reserve(1 + 2 * gridShape.size());
  out.second.push_back(degrees.size());
  out.second.insert(out.second.end(), gridShape.begin(), gridShape.end());
  out.second.insert(out.second.end(), gridShape.begin(), gridShape.end());
  std::size_t total = 1;
  for (I value : out.second) {
    total *= static_cast<std::size_t>(value);
  }
  if (total == 0) return out;
  out.first.resize(total, 0);
  Tensor<T> tensor(out.first, out.second);

  oneapi::tbb::task_arena arena(Slicer::Persistence::is_vine ? 1 : n_jobs);  // limits the number of threads
  switch (num_parameters) {
    case 2: {
      const auto shape = detail::_vector_to_full_rank_shape<2>(gridShape);
      arena.execute([&] { detail::_compute_2d_rank_invariant(slicer, tensor, shape, degrees, false, ignoreInf); });
      return out;
    }
    case 3: {
      const auto shape = detail::_vector_to_full_rank_shape<3>(gridShape);
      const std::array<I, 3> parameter_shape = {shape[1], shape[2], shape[3]};
      arena.execute([&] {
        detail::_compute_nd_rank_invariant_by_slice_paths<3>(slicer, tensor, parameter_shape, degrees, ignoreInf);
      });
      return out;
    }
    case 4: {
      const auto shape = detail::_vector_to_full_rank_shape<4>(gridShape);
      const std::array<I, 4> parameter_shape = {shape[1], shape[2], shape[3], shape[4]};
      arena.execute([&] {
        detail::_compute_nd_rank_invariant_by_slice_paths<4>(slicer, tensor, parameter_shape, degrees, ignoreInf);
      });
      return out;
    }
    default:
      throw std::runtime_error("Rank invariant is implemented for 2, 3, and 4 parameters.");
  }
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_RANK_INVARIANT_H_