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
 * @private
 * @file utils.h
 * @author Hannah Schreiber
 */

#ifndef MP_UTILS_H_
#define MP_UTILS_H_

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace multi_persistence {

template <class Index>
constexpr Index null_vertex() {
  return static_cast<Index>(-1);
}

namespace detail {

template <class T, typename = void>
struct is_forward_iterator : std::false_type {};

template <class T>
struct is_forward_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category>>
    : std::bool_constant<
          std::is_base_of_v<std::forward_iterator_tag, typename std::iterator_traits<T>::iterator_category>> {};

template <typename T>
struct type_identity {
  using type = T;
};

template <class Index, class Boundaries>
void _compute_union(std::size_t i, const Boundaries& boundaries, std::vector<std::vector<Index>>& vertices,
                    std::vector<Index>& tmp) {
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

}  // namespace detail

template <class T>
constexpr bool is_forward_iterator_v = detail::is_forward_iterator<T>::value;

// std::make_signed_t does not compile for T signed and std::conditional evaluates both possibilities
// so this trick is necessary if we want to avoid using `if constexpr` everywhere
template <typename T>
using maybe_make_signed_t =
    typename std::conditional_t<std::is_unsigned_v<T>, std::make_signed<T>, detail::type_identity<T>>::type;

// boundaries ordered by dimension and empty boundaries are of lowest dimension
template <class Index, class Boundaries, class Dimensions>
std::vector<std::vector<Index>> get_vertices_from_ordered_boundaries(const Boundaries& boundaries,
                                                                     const Dimensions& dimensions, int maxDim = -1) {
  const std::size_t numberOfCells =
      (dimensions.back() <= maxDim || maxDim < 0)
          ? boundaries.size()
          : static_cast<std::size_t>(
                std::distance(dimensions.begin(), std::upper_bound(dimensions.begin(), dimensions.end(), maxDim)));
  if (numberOfCells == 0) return {};

  GUDHI_CHECK(boundaries[0].empty(), std::invalid_argument("At least one cell with empty boundary is required."));

  std::vector<std::vector<Index>> vertices(numberOfCells);

  // No vertices of dimension 0: every leaf is null, so every cell's set is just {-1}.
  if (dimensions[0] != 0) {
    for (auto& v : vertices) v.assign(1, null_vertex<Index>());
    return vertices;
  }

  std::size_t curr = 0;
  while (curr < numberOfCells && boundaries[curr].empty()) {
    vertices[curr].assign(1, static_cast<Index>(curr));
    ++curr;
  }

  std::vector<Index> tmp;
  GUDHI_CHECK_code(int prevDim = 0;);
  while (curr < numberOfCells) {
    GUDHI_CHECK(dimensions[curr] > prevDim, std::invalid_argument("Boundaries are not ordered by dimension!"));
    GUDHI_CHECK(!boundaries[curr].empty(),
                std::invalid_argument("There should be no empty boundary which is not of dimension 0"));
    // facets have smaller indices
    detail::_compute_union<Index>(curr, boundaries, vertices, tmp);
    ++curr;
    GUDHI_CHECK_code(if (curr < numberOfCells && dimensions[curr] != dimensions[curr - 1])
                         prevDim = dimensions[curr - 1]);
  }
  return vertices;
}

// general case
template <class Index, class Boundaries, class Dimensions>
std::vector<std::vector<Index>> get_vertices_from_boundaries(const Boundaries& boundaries, const Dimensions& dimensions,
                                                             int maxDim = -1) {
  const std::size_t numberOfCells = boundaries.size();
  std::vector<std::vector<Index>> vertices(numberOfCells);

  // initialize leafs/vertices
  Index c = 0;
  GUDHI_CHECK_code(std::size_t leaves = 0;);
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (boundaries[i].empty()) {
      vertices[i].assign(1, dimensions[i] == 0 ? c++ : null_vertex<Index>());
      GUDHI_CHECK_code(++leaves;);
    }
  }
  GUDHI_CHECK(leaves > 0, std::invalid_argument("At least one boundaries has to be empty."));

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
        // TODO: Put the throw in a GUDHI_CHECK?
        // This will trigger an endless loop, so it is probably better to throw in any compile mode?
        if (state[child] == IN_PROGRESS) throw std::invalid_argument("Cycle in boundaries.");
        if (state[child] == NEW) {
          state[child] = IN_PROGRESS;
          stack.emplace_back(child, 0);
        }
      } else {
        detail::_compute_union<Index>(cell, boundaries, vertices, tmp);
        state[cell] = DONE;
        stack.pop_back();
      }
    }
  }
  return vertices;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_UTILS_H_
