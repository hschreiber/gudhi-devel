/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
#define INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_

#include <array>
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <utility>    // for std::pair
#include <algorithm>  // for sort
#include <unordered_map>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>  // for Extended_simplex_type
#include <python_interfaces/numpy_utils.h>

namespace Gudhi {

/**
 * @private
 */
template <typename T>
class ComplexTraits {
 private:
  static auto check_simplex_vertex_range(...) -> std::false_type;
  template <typename U>
  static auto check_simplex_vertex_range(const U& x)
      -> decltype(x.simplex_vertex_range(typename U::Simplex_handle()), std::true_type{});

  static auto check_find_vertex(...) -> std::false_type;
  template <typename U>
  static auto check_find_vertex(const U& x) -> decltype(x.find_vertex(0), std::true_type{});

 public:
  static constexpr bool has_simplex_vertex_range = decltype(check_simplex_vertex_range(std::declval<T>()))::value;
  static constexpr bool has_find_vertex = decltype(check_find_vertex(std::declval<T>()))::value;
};

template <class FilteredComplex>
class Persistent_cohomology_interface
    : public persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp>
{
 private:
  typedef persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> Base;

  /*
   * Compare two intervals by dimension, then by length.
   */
  struct cmp_intervals_by_dim_then_length {
    template <typename Persistent_interval>
    bool operator()(const Persistent_interval& p1, const Persistent_interval& p2)
    {
      if (std::get<0>(p1) == std::get<0>(p2)) {
        auto& i1 = std::get<1>(p1);
        auto& i2 = std::get<1>(p2);
        return std::get<1>(i1) - std::get<0>(i1) > std::get<1>(i2) - std::get<0>(i2);
      } else
        return (std::get<0>(p1) > std::get<0>(p2));
      // Why does this sort by decreasing dimension?
    }
  };

 public:
  Persistent_cohomology_interface(FilteredComplex& stptr, bool persistence_dim_max = false)
      : Base(stptr, persistence_dim_max), stptr_(&stptr)
  {}

  // TODO: move to the constructors?
  void compute_persistence(int homology_coeff_field, double min_persistence)
  {
    Base::init_coefficients(homology_coeff_field);
    Base::compute_persistent_cohomology(min_persistence);
  }

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence()
  {
    std::vector<std::pair<int, std::pair<double, double>>> persistence;
    auto const& persistent_pairs = Base::get_persistent_pairs();
    persistence.reserve(persistent_pairs.size());
    for (auto pair : persistent_pairs) {
      persistence.emplace_back(stptr_->dimension(get<0>(pair)),
                               std::make_pair(stptr_->filtration(get<0>(pair)), stptr_->filtration(get<1>(pair))));
    }
    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp;
    std::sort(std::begin(persistence), std::end(persistence), cmp);
    return persistence;
  }

  nanobind::tuple get_intervals()
  {
    std::pair<std::vector<std::vector<std::array<double, 2> > >, std::vector<std::vector<double> > > intervals;
    intervals.first.resize(Base::dim_max_);
    intervals.second.resize(Base::dim_max_);

    for (const auto& pair : Base::get_persistent_pairs()) {
      if (get<1>(pair) == stptr_->null_simplex()) {
        intervals.second[stptr_->dimension(get<0>(pair))].push_back(stptr_->filtration(get<0>(pair)));
      } else {
        intervals.first[stptr_->dimension(get<0>(pair))].push_back({stptr_->filtration(get<0>(pair)), stptr_->filtration(get<1>(pair))});
      }
    }

    nanobind::list finites;
    nanobind::list infinites;

    for (int d = 0; d < Base::dim_max_; ++d) {
      // TODO: sort like for the others?
      finites.append(_wrap_as_numpy_array(std::move(intervals.first[d])));
      infinites.append(_wrap_as_numpy_array(std::move(intervals.second[d]), intervals.second[d].size()));
    }

    return nanobind::make_tuple(finites, infinites);
  }

  // This function computes the top-dimensional cofaces associated to the positive and negative
  // simplices of a cubical complex. The output format is a vector of vectors of three integers,
  // which are [homological dimension, index of top-dimensional coface of positive simplex,
  // index of top-dimensional coface of negative simplex]. If the topological feature is essential,
  // then the index of top-dimensional coface of negative simplex is arbitrarily set to -1.
  std::vector<std::vector<int>> cofaces_of_cubical_persistence_pairs()
  {
    // Warning: this function is meant to be used with CubicalComplex only!!

    auto&& pairs = Base::get_persistent_pairs();

    // Compute the ordering function of the top-dimensional cells simplex handles
    // This function allows to go directly from the simplex handle to the position of the corresponding top-dimensional
    // cell in the input data
    std::unordered_map<std::size_t, int> order;
    unsigned idx = 0;
    for (auto splx : stptr_->top_dimensional_cells_range()) {
      order.emplace(splx, idx);
      idx++;
    }

    std::vector<std::vector<int>> persistence_pairs;

    for (auto pair : pairs) {
      int h = static_cast<int>(stptr_->dimension(get<0>(pair)));
      // Recursively get the top-dimensional cell / coface associated to the persistence generator
      std::size_t face0 = stptr_->get_top_dimensional_coface_of_a_cell(get<0>(pair));
      // Retrieve the index of the corresponding top-dimensional cell in the input data
      int splx0 = order[face0];

      int splx1 = -1;
      if (get<1>(pair) != stptr_->null_simplex()) {
        // Recursively get the top-dimensional cell / coface associated to the persistence generator
        std::size_t face1 = stptr_->get_top_dimensional_coface_of_a_cell(get<1>(pair));
        // Retrieve the index of the corresponding top-dimensional cell in the input data
        splx1 = order[face1];
      }
      persistence_pairs.push_back({h, splx0, splx1});
    }

    return persistence_pairs;
  }

  // This function computes the vertices associated to the positive and negative
  // simplices of a cubical complex. The output format is a vector of vectors of three integers,
  // which are [homological dimension, index of vertex of positive simplex,
  // index of vertex of negative simplex]. If the topological feature is essential,
  // then the index of vertex of negative simplex is arbitrarily set to -1.
  std::vector<std::vector<int>> vertices_of_cubical_persistence_pairs()
  {
    // Warning: this function is meant to be used with CubicalComplex only!!
    auto&& pairs = Base::get_persistent_pairs();

    // Compute the ordering function of the vertices simplex handles
    // This function allows to go directly from the simplex handle to the position of the corresponding vertex in the
    // input data
    std::unordered_map<std::size_t, int> order;
    unsigned idx = 0;
    for (auto splx : stptr_->vertices_range()) {
      order.emplace(splx, idx);
      idx++;
    }

    std::vector<std::vector<int>> persistence_pairs;

    for (auto pair : pairs) {
      int h = static_cast<int>(stptr_->dimension(get<0>(pair)));
      // Recursively get the vertex associated to the persistence generator
      std::size_t face0 = stptr_->get_vertex_of_a_cell(get<0>(pair));
      // Retrieve the index of the corresponding vertex in the input data
      int splx0 = order[face0];

      int splx1 = -1;
      if (get<1>(pair) != stptr_->null_simplex()) {
        // Recursively get the vertex associated to the persistence generator
        std::size_t face1 = stptr_->get_vertex_of_a_cell(get<1>(pair));
        // Retrieve the index of the corresponding vertex in the input data
        splx1 = order[face1];
      }
      persistence_pairs.push_back({h, splx0, splx1});
    }

    return persistence_pairs;
  }

  std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs()
  {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs;
    auto const& pairs = Base::get_persistent_pairs();
    persistence_pairs.reserve(pairs.size());
    std::vector<int> birth;
    std::vector<int> death;
    for (auto pair : pairs) {
      birth.clear();
      if (get<0>(pair) != stptr_->null_simplex()) {
        for (auto vertex : stptr_->simplex_vertex_range(get<0>(pair))) {
          birth.push_back(vertex);
        }
      }

      death.clear();
      if (get<1>(pair) != stptr_->null_simplex()) {
        death.reserve(birth.size() + 1);
        for (auto vertex : stptr_->simplex_vertex_range(get<1>(pair))) {
          death.push_back(vertex);
        }
      }

      persistence_pairs.emplace_back(birth, death);
    }
    return persistence_pairs;
  }

  // only for simplices as it is not really usefull for general cells
  nanobind::tuple get_simplicial_intervals_as_vertices(double min_persistence, int maxDim = -1,
                                                       bool only_finite = false, bool only_infinite = false) {
    // won't compile even without the static_assert, it is just to have a clearer message as
    // `get_simplicial_intervals_as_vertices` is not a necessary function
    static_assert(ComplexTraits<FilteredComplex>::has_simplex_vertex_range,
                  "Contrary to all other methods, `get_simplicial_intervals_as_vertices` needs the complex class to "
                  "have a `simplex_vertex_range` method.");
    
    using Simplex_handle = typename FilteredComplex::Simplex_handle;

    auto add_simplex = [&](Simplex_handle sh, std::size_t dim, std::vector<std::vector<int>>& cont) {
      for (auto vertex : stptr_->simplex_vertex_range(sh)) {
        cont[dim].push_back(vertex);
      }
    };
    auto add_dummy_simplex = [&](std::size_t dim, std::vector<std::vector<int>>& cont) {
      cont[dim].resize(cont[dim].size() + dim + 2, -1);
    };

    maxDim = maxDim < 0 ? Base::dim_max_ : maxDim + 1;
    std::pair<std::vector<std::vector<int> >, std::vector<std::vector<int> > > intervals;
    intervals.first.resize(only_infinite ? 0 : maxDim);
    intervals.second.resize(only_finite ? 0 : maxDim);

    for (const auto& pair : Base::get_persistent_pairs()) {
      if (get<0>(pair) != stptr_->null_simplex()) {
        int dim = stptr_->dimension(get<0>(pair));
        if (dim < maxDim) {
          if (!only_finite && get<1>(pair) == stptr_->null_simplex()) {
            add_simplex(get<0>(pair), dim, intervals.second);
          } else if (!only_infinite && get<1>(pair) != stptr_->null_simplex()) {
            auto b = get<0>(pair);
            auto d = get<1>(pair);
            if (stptr_->filtration(d) - stptr_->filtration(b) >= min_persistence) {
              add_simplex(b, dim, intervals.first);
              intervals.first[dim].push_back(-1);
              if (dim == maxDim - 1) add_dummy_simplex(dim, intervals.first);
              else add_simplex(d, dim, intervals.first);
            }
          }
        }
      }
    }

    nanobind::list finites;
    nanobind::list infinites;

    for (int d = 0; d < maxDim; ++d) {
      int size = d + 2;
      if (!only_infinite)
        finites.append(
            _wrap_as_numpy_array(std::move(intervals.first[d]), intervals.first[d].size() / (size * 2), 2, size));
      if (!only_finite)
        infinites.append(_wrap_as_numpy_array(std::move(intervals.second[d]), intervals.second[d].size() / (size - 1), size - 1));
    }

    return nanobind::make_tuple(finites, infinites);
  }

  // TODO: (possibly at the python level)
  // - an option to return only some of those vectors?
  typedef std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> Generators;

  Generators lower_star_generators()
  {
    Generators out;
    // diags[i] should be interpreted as vector<array<int,2>>
    auto& diags = out.first;
    // diagsinf[i] should be interpreted as vector<int>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      int dim = stptr_->dimension(s);
      auto v = stptr_->vertex_with_same_filtration(s);
      if (t == stptr_->null_simplex()) {
        while (static_cast<int>(diagsinf.size()) < dim + 1) diagsinf.emplace_back();
        diagsinf[dim].push_back(v);
      } else {
        while (static_cast<int>(diags.size()) < dim + 1) diags.emplace_back();
        auto w = stptr_->vertex_with_same_filtration(t);
        auto& d = diags[dim];
        d.insert(d.end(), {v, w});
      }
    }
    return out;
  }

  // An alternative, to avoid those different sizes, would be to "pad" vertex generator v as (v, v) or (v, -1). When
  // using it as index, this corresponds to adding the vertex filtration values either on the diagonal of the distance
  // matrix, or as an extra row or column. We could also merge the vectors for different dimensions into a single one,
  // with an extra column for the dimension (converted to type double).
  Generators flag_generators()
  {
    Generators out;
    // diags[0] should be interpreted as vector<array<int,3>> and other diags[i] as vector<array<int,4>>
    auto& diags = out.first;
    // diagsinf[0] should be interpreted as vector<int> and other diagsinf[i] as vector<array<int,2>>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      int dim = stptr_->dimension(s);
      bool infinite = t == stptr_->null_simplex();
      if (infinite) {
        if (dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if (diagsinf.size() == 0) diagsinf.emplace_back();
          diagsinf[0].push_back(v);
        } else {
          auto e = stptr_->edge_with_same_filtration(s);
          auto&& e_vertices = stptr_->simplex_vertex_range(e);
          auto i = std::begin(e_vertices);
          auto v1 = *i;
          auto v2 = *++i;
          GUDHI_CHECK(++i == std::end(e_vertices), "must be an edge");
          while (static_cast<int>(diagsinf.size()) < dim + 1) diagsinf.emplace_back();
          auto& d = diagsinf[dim];
          d.insert(d.end(), {v1, v2});
        }
      } else {
        auto et = stptr_->edge_with_same_filtration(t);
        auto&& et_vertices = stptr_->simplex_vertex_range(et);
        auto it = std::begin(et_vertices);
        auto w1 = *it;
        auto w2 = *++it;
        GUDHI_CHECK(++it == std::end(et_vertices), "must be an edge");
        if (dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if (diags.size() == 0) diags.emplace_back();
          auto& d = diags[0];
          d.insert(d.end(), {v, w1, w2});
        } else {
          auto es = stptr_->edge_with_same_filtration(s);
          auto&& es_vertices = stptr_->simplex_vertex_range(es);
          auto is = std::begin(es_vertices);
          auto v1 = *is;
          auto v2 = *++is;
          GUDHI_CHECK(++is == std::end(es_vertices), "must be an edge");
          while (static_cast<int>(diags.size()) < dim + 1) diags.emplace_back();
          auto& d = diags[dim];
          d.insert(d.end(), {v1, v2, w1, w2});
        }
      }
    }
    return out;
  }

  using Filtration_value = typename FilteredComplex::Filtration_value;
  using Birth_death = std::pair<Filtration_value, Filtration_value>;
  using Persistence_subdiagrams = std::vector<std::vector<std::pair<int, Birth_death>>>;

  Persistence_subdiagrams compute_extended_persistence_subdiagrams(Filtration_value min_persistence)
  {
    Persistence_subdiagrams pers_subs(4);
    auto const& persistent_pairs = Base::get_persistent_pairs();
    for (auto pair : persistent_pairs) {
      std::pair<Filtration_value, Extended_simplex_type> px =
          stptr_->decode_extended_filtration(stptr_->filtration(get<0>(pair)), stptr_->efd);
      std::pair<Filtration_value, Extended_simplex_type> py =
          stptr_->decode_extended_filtration(stptr_->filtration(get<1>(pair)), stptr_->efd);
      std::pair<int, Birth_death> pd_point =
          std::make_pair(stptr_->dimension(get<0>(pair)), std::make_pair(px.first, py.first));
      if (std::abs(px.first - py.first) > min_persistence) {
        // Ordinary
        if (px.second == Extended_simplex_type::UP && py.second == Extended_simplex_type::UP) {
          pers_subs[0].push_back(pd_point);
        }
        // Relative
        else if (px.second == Extended_simplex_type::DOWN && py.second == Extended_simplex_type::DOWN) {
          pers_subs[1].push_back(pd_point);
        } else {
          // Extended+
          if (px.first < py.first) {
            pers_subs[2].push_back(pd_point);
          }
          // Extended-
          else {
            pers_subs[3].push_back(pd_point);
          }
        }
      }
    }
    return pers_subs;
  }

 private:
  // A copy
  FilteredComplex* stptr_;
};

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
