/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @private
 * @file tensor_helpers.h
 * @author Hannah Schreiber, David Loiseaux
 */

#ifndef GUDHI_MP_TENSOR_HELPERS_H_INCLUDED
#define GUDHI_MP_TENSOR_HELPERS_H_INCLUDED

#include <cstddef>      //std::size_t
#include <stdexcept>    //std::out_of_range
#include <type_traits>  //std::is_const_v, std::is_pointer_v
#include <algorithm>    //std::min, std::max
#include <initializer_list>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>

namespace Gudhi {
namespace multi_persistence {

// Dynamic version of mdspan where ndim does not have to be known at compile time and data has to be continuous
template <typename T>
class Tensor {
 public:
  using element_type = T;
  using Container = std::vector<element_type>;
  using size_type = typename Container::size_type;
  using index_type = typename Container::size_type;
  using rank_type = std::size_t;
  using Shape = std::vector<size_type>;
  using Strides = std::vector<size_type>;
  using data_handle_type = element_type*;
  using reference = element_type&;

  // TODO: remove?
  template <class Range = std::initializer_list<element_type>>
  Tensor(element_type* data, const Range& shape)
      : ptr_(data), shape_(shape.begin(), shape.end()), strides_(shape.size()) {
    _initialize_strides();
  }

  template <class Range = std::initializer_list<element_type>>
  Tensor(const Container& data, const Range& shape)
      : ptr_(data.data()), shape_(shape.begin(), shape.end()), strides_(shape.size()) {
    _initialize_strides();
    GUDHI_CHECK(data.size() == size(), std::invalid_argument("Given shape does not corresponds to container size."));
  }

  Tensor(const Container& data, Shape&& shape) : ptr_(data.data()), shape_(std::move(shape)), strides_(shape_.size()) {
    _initialize_strides();
    GUDHI_CHECK(data.size() == size(), std::invalid_argument("Given shape does not corresponds to container size."));
  }

  template <class... IndexTypes>
  reference operator()(IndexTypes... indices) const {
    return operator[]({static_cast<index_type>(indices)...});
  }

  template <class IndexRange = std::initializer_list<index_type>>
  reference operator[](const IndexRange& indices) const {
    return *(ptr_ + _map(indices));
  }

  rank_type rank() const noexcept { return shape_.size(); }

  rank_type rank_dynamic() const noexcept { return shape_.size(); }

  static constexpr std::size_t static_extent([[maybe_unused]] rank_type r) noexcept {
    return std::numeric_limits<std::size_t>::max();
  }

  index_type extent(rank_type r) const {
    GUDHI_CHECK(r < shape_.size(), std::out_of_range("Out of bound index."));
    return shape_[r];
  }

  size_type size() const noexcept {
    if (shape_.size() == 0)
      return 0;
    else
      return strides_[0] * shape_[0];
  }

  bool empty() const noexcept { return size() == 0; }

  index_type stride(rank_type r) const {
    GUDHI_CHECK(r < strides_.size(), std::out_of_range("Out of bound index."));
    return strides_[r];
  }

  // not extents_type using here as shape_ does not has the same interface at all
  const Shape& extents() const noexcept { return shape_; }

  constexpr const data_handle_type& data_handle() const noexcept { return ptr_; }

  // if is_unique() is true for all possible instantiations of this class
  static constexpr bool is_always_unique() noexcept { return true; }

  // if is_exhaustive() is true for all possible instantiations of this class
  static constexpr bool is_always_exhaustive() noexcept { return true; }

  // if is_strided() is true for all possible instantiations of this class
  static constexpr bool is_always_strided() noexcept { return true; }

  // unicity of the mapping (i,j,k,...) -> real index
  static constexpr bool is_unique() noexcept { return true; }

  // if all real indices have a preimage in form (i,j,k,...)
  static constexpr bool is_exhaustive() noexcept { return true; }

  // if distance in memory is constant between two values in same rank
  static constexpr bool is_strided() noexcept { return true; }

  friend constexpr void swap(Tensor& x, Tensor& y) noexcept {
    std::swap(x.ptr_, y.ptr_);
    x.shape_.swap(y.shape_);
    x.strides_.swap(y.strides_);
  }

 private:
  Container::pointer ptr_;
  Shape shape_;
  Strides strides_;

  void _initialize_strides() {
    strides_[rank() - 1] = 1;
    for (auto i = rank() - 1; i > 0; --i) {
      strides_[i - 1] = strides_[i] * extent(i);
    }
  }

  template <class IndexRange = std::initializer_list<index_type>>
  index_type _map(const IndexRange& indices) const {
    GUDHI_CHECK(indices.size() == rank(), std::invalid_argument("Wrong number of parameters."));

    index_type newIndex = 0;
    auto it = indices.begin();
    GUDHI_CHECK_code(unsigned int i = 0);
    for (auto stride : strides_) {
      GUDHI_CHECK_code(GUDHI_CHECK(*it < extent(i), std::out_of_range("Out of bound index.")));
      newIndex += (stride * (*it));
      ++it;
      GUDHI_CHECK_code(++i);
    }

    return newIndex;
  }
};

namespace detail {

// Backward difference along dimension 1 of `v`, restricted to
//   o in [o_begin, o_end)  and  k in [k_begin, k_end).
// Distinct (o, k) pairs touch disjoint elements, so ranges can run concurrently
template <typename T, class Extents>
inline void _differentiate_range(const Gudhi::Simple_mdspan<T, Extents>& data, std::size_t o_begin, std::size_t o_end,
                                 std::size_t k_begin, std::size_t k_end) {
  for (std::size_t o = o_begin; o < o_end; ++o)
    for (std::size_t h = data.extent(1) - 1; h >= 1; --h)
      for (std::size_t k = k_begin; k < k_end; ++k) data(o, h, k) -= data(o, h - 1, k);
}

template <typename Index, class Tensor_t>
inline void _fill_data_index_coordinates(const Tensor_t& tensor, typename Tensor_t::index_type dataIndex,
                                         const std::vector<bool>& flipAxes, std::vector<Index>& coordinates) {
  const std::size_t startIndex = coordinates.size();
  coordinates.resize(startIndex + tensor.rank());
  for (int axis = static_cast<int>(tensor.rank()) - 1; axis >= 0; --axis) {
    const auto q = dataIndex / tensor.extent(axis);
    const auto r = dataIndex - (q * tensor.extent(axis));
    if (static_cast<int>(flipAxes.size()) > axis && flipAxes[axis])
      coordinates[startIndex + axis] = static_cast<Index>(tensor.extent(axis) - r);
    else
      coordinates[startIndex + axis] = static_cast<Index>(r);
    dataIndex = q;
  }
}

}  // namespace detail

template <class Tensor_t>
inline void differentiate(const Tensor_t& tensor, std::size_t axis) {
  static_assert(!std::is_const_v<typename Tensor_t::element_type>, "differentiate: the tensor must be writable.");
  static_assert(std::is_pointer_v<typename Tensor_t::data_handle_type>,
                "differentiate: requires a raw-pointer data handle.");
  static_assert(Tensor_t::is_always_exhaustive() && Tensor_t::is_always_strided(),
                "differentiate: requires a contiguous layout.");

  if (axis >= tensor.rank()) throw std::out_of_range("differentiate: invalid axis");

  const std::size_t axis_size = static_cast<std::size_t>(tensor.extent(axis));
  if (axis_size <= 1 || tensor.empty()) return;

  // The data is a flat block seen as [outer][axis_size][stride], where stride is the distance
  // between consecutive indices along `axis` (valid because contiguous layout).
  const std::size_t stride = static_cast<std::size_t>(tensor.stride(axis));
  const std::size_t outer = tensor.size() / (axis_size * stride);
  const auto view = Gudhi::Simple_mdspan(tensor.data_handle(), outer, axis_size, stride);

#ifdef GUDHI_USE_TBB
  const std::size_t block = axis_size * stride;
  // Those are average general values, we could do benchmarks to affine them to this special case
  constexpr std::size_t min_parallel_size = 1u << 15;  // below this, threading costs more than it saves
  constexpr std::size_t target_chunk = 1u << 12;       // ~elements of work per TBB task

  if (tensor.size() >= min_parallel_size) {
    // 2D split over (outer, column): stride == 1 -> only outer is split,
    // outer == 1 -> only columns are split, otherwise both.
    const std::size_t grain_k = std::min<std::size_t>(stride, target_chunk);
    const std::size_t grain_o = std::max<std::size_t>(1, target_chunk / block);
    tbb::parallel_for(tbb::blocked_range2d<std::size_t>(0, outer, grain_o, 0, stride, grain_k),
                      [=](const tbb::blocked_range2d<std::size_t>& r) {
                        detail::_differentiate_range(view, r.rows().begin(), r.rows().end(), r.cols().begin(),
                                                     r.cols().end());
                      });
    return;
  }
#endif

  detail::_differentiate_range(view, 0, outer, 0, stride);
}

template <typename Index, class Tensor_t>
inline std::pair<std::vector<Index>, std::vector<typename Tensor_t::element_type>> sparsify(
    const Tensor_t& tensor, const std::vector<bool>& flipAxes = {}) {
  static_assert(std::is_pointer_v<typename Tensor_t::data_handle_type>,
                "differentiate: requires a raw-pointer data handle.");
  static_assert(Tensor_t::is_always_exhaustive() && Tensor_t::is_always_strided(),
                "sparsify: requires a contiguous layout.");

  std::vector<Index> coordinates;
  std::vector<typename Tensor_t::element_type> values;
  for (typename Tensor_t::index_type i = 0; i < tensor.size(); ++i) {
    const auto value = *(tensor.data_handle() + i);
    if (value != 0) {
      detail::_fill_data_index_coordinates(tensor, i, flipAxes, coordinates);
      values.push_back(value);
    }
  }
#ifdef DEBUG_TRACES
  for (std::size_t i = 0; i < (coordinates.size() / tensor.rank()); ++i) {
    for (std::size_t j = 0; j < tensor.rank(); ++j) std::cout << coordinates[j + tensor.rank() * i] << " ";
    std::cout << "| " << values[i] << std::endl;
  }
#endif  // DEBUG_TRACES
  // coordinates has shape number_of_non_null_values * tensor.rank()
  return {coordinates, values};
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // GUDHI_MP_TENSOR_HELPERS_H_INCLUDED
