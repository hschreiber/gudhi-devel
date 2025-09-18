/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2025/09 Hannah Schreiber: Split of the Persistent_cohomology class into smaller modules
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_H_
#define PERSISTENT_COHOMOLOGY_H_

#include <stdexcept>  // std::out_of_range, std::invalid_argument
#include <iostream>   // std::ostream, std::cout, std::endl
#include <fstream>    // std::ofstream
#include <algorithm>  // std::sort, std::max
#include <limits>     // std::numeric_limits
#include <utility>    // std::pair
#include <string>     // std::string
#include <tuple>      // std::tuple
#include <vector>     //std::vector, std::begin, std::end

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistent_cohomology/Persistent_cohomology_base.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>  //only here for external use....

namespace Gudhi {

namespace persistent_cohomology {

/** \brief Computes the persistent cohomology of a filtered complex.
 *
 * \ingroup persistent_cohomology
 *
 * The computation is implemented with a Compressed Annotation Compressed_annotation_matrix
 * (CAM)\cite DBLP:conf/esa/BoissonnatDM13,
 * and is adapted to the computation of Multi-Field Persistent Homology (MF)
 * \cite boissonnat:hal-00922572 .
 *
 * \implements PersistentHomology
 *
 */
// TODO(CM): Memory allocation policy: classic, use a mempool, etc.
template <class FilteredComplex, class CoefficientField>
class Persistent_cohomology
{
 public:
  // Data attached to each simplex to interface with a Property Map.

  /** \brief Data stored for each simplex. */
  using Simplex_key = typename FilteredComplex::Simplex_key;
  /** \brief Handle to specify a simplex. */
  using Simplex_handle = typename FilteredComplex::Simplex_handle;
  /** \brief Type for the value of the filtration function. */
  using Filtration_value = typename FilteredComplex::Filtration_value;
  /** \brief Type of element of the field. */
  using Arith_element = typename CoefficientField::Element;
  /** \brief Type for birth and death FilteredComplex::Simplex_handle.
   * The Arith_element field is used for the multi-field framework. */
  using Persistent_interval = std::tuple<Simplex_handle, Simplex_handle, typename CoefficientField::Element>;

  /** \brief Initializes the Persistent_cohomology class.
   *
   * @param[in] cpx Complex for which the persistent homology is computed.
   * cpx is a model of FilteredComplex
   *
   * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the
   *                                complex is computed. If false, it is ignored. Default is false.
   *
   * @exception std::out_of_range In case the number of simplices is more than Simplex_key type numeric limit.
   */
  explicit Persistent_cohomology(FilteredComplex& cpx, bool persistence_dim_max = false)
      : cpx_(&cpx),
        dim_max_(cpx.dimension()),             // upper bound on the dimension of the simplices
        coeff_field_(),                        // initialize the field coefficient structure.
        num_simplices_(cpx_->num_simplices())  // num_simplices save to avoid to call thrice the function
  {
    if (num_simplices_ > std::numeric_limits<Simplex_key>::max()) {
      // num_simplices must be strictly lower than the limit, because a value is reserved for null_key.
      throw std::out_of_range("The number of simplices is more than Simplex_key type numeric limit.");
    }
    if (persistence_dim_max) {
      ++dim_max_;
    }
  }

  /** \brief Initializes the coefficient field.*/
  void init_coefficients(int charac)
  {
    if (charac <= 0) throw std::invalid_argument("Characteristic has to be a positive prime number.");
    coeff_field_.set_characteristic(charac);
  }

  /** \brief Initializes the coefficient field for multi-field persistent homology.*/
  void init_coefficients(int charac_min, int charac_max) { coeff_field_.set_characteristic(charac_min, charac_max); }

  /** \brief Compute the persistent homology of the filtered simplicial
   * complex.
   *
   * @param[in] min_interval_length the computation discards all intervals of length
   *                                less or equal than min_interval_length
   *
   * Assumes that the filtration provided by the simplicial complex is
   * valid. Undefined behavior otherwise. */
  void compute_persistent_cohomology(Filtration_value min_interval_length = 0)
  {
    if (coeff_field_.get_characteristic() == CoefficientField::nullCharacteristic)
      throw std::logic_error("Coefficient field was not initialized! Please call `init_coefficients` first.");

    if (dim_max_ <= 0) return;

    Persistent_cohomology_base<FilteredComplex, CoefficientField, length_interval, true> base(
        cpx_, coeff_field_, length_interval(cpx_, min_interval_length), num_simplices_, cpx_->num_vertices());

    Simplex_key idx_fil = -1;
    Simplex_key idx_fil_v = -1;
    Simplex_handle u, v;
    // Compute all finite intervals
    for (auto sh : cpx_->filtration_simplex_range()) {
      int dim_simplex = cpx_->dimension(sh);
      switch (dim_simplex) {
        case 0:
          cpx_->assign_key(sh, ++idx_fil_v);
          ++idx_fil;
          base.update_cohomology_groups_vertex(sh);
          break;
        case 1:
          cpx_->assign_key(sh, ++idx_fil);
          boost::tie(u, v) = cpx_->endpoints(sh);
          base.update_cohomology_groups_edge(sh, u, v, dim_max_ == 1);
          break;
        default:
          cpx_->assign_key(sh, ++idx_fil);
          base.update_cohomology_groups(sh, cpx_->boundary_simplex_range(sh), dim_simplex == dim_max_);
          break;
      }
    }
    base.finalize_computation();
    persistent_pairs_ = base.get_persistent_pairs();
  }

  /**
   * @private
   * Temporary patch to make the computation of more general FilteredComplex possible
   * Will be removed when a better solution was agreed on.
   */
  void compute_persistent_cohomology_without_optimizations(Filtration_value min_interval_length = 0)
  {
    if (coeff_field_.get_characteristic() == CoefficientField::nullCharacteristic)
      throw std::logic_error("Coefficient field was not initialized! Please call `init_coefficients` first.");

    Persistent_cohomology_base<FilteredComplex, CoefficientField, length_interval, false> base(
        cpx_, coeff_field_, length_interval(cpx_, min_interval_length), num_simplices_, 0);

    Simplex_key idx_fil = -1;
    for (auto sh : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(sh, ++idx_fil);
      base.update_cohomology_groups(sh, cpx_->boundary_simplex_range(sh), cpx_->dimension(sh) == dim_max_);
    }
    base.finalize_computation();
    persistent_pairs_ = base.get_persistent_pairs();
  }

  /** \brief Output the persistence diagram in ostream.
   *
   * The file format is the following:
   *    p1*...*pr   dim b d
   *
   * where "dim" is the dimension of the homological feature,
   * b and d are respectively the birth and death of the feature and
   * p1*...*pr is the product of prime numbers pi such that the homology
   * feature exists in homology with Z/piZ coefficients.
   */
  void output_diagram(std::ostream& ostream = std::cout)
  {
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(persistent_pairs_), std::end(persistent_pairs_), cmp);
    for (auto pair : persistent_pairs_) {
      ostream << get<2>(pair) << "  " << cpx_->dimension(get<0>(pair)) << " " << cpx_->filtration(get<0>(pair)) << " "
              << cpx_->filtration(get<1>(pair)) << " " << std::endl;
    }
  }

  void write_output_diagram(const std::string& diagram_name)
  {
    std::ofstream diagram_out(diagram_name.c_str());
    diagram_out.exceptions(std::ofstream::failbit);
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(persistent_pairs_), std::end(persistent_pairs_), cmp);
    for (auto pair : persistent_pairs_) {
      diagram_out << cpx_->dimension(get<0>(pair)) << " " << cpx_->filtration(get<0>(pair)) << " "
                  << cpx_->filtration(get<1>(pair)) << std::endl;
    }
  }

  /** @brief Returns Betti numbers.
   * @return A vector of Betti numbers.
   */
  [[nodiscard]] std::vector<int> betti_numbers() const
  {
    // Init Betti numbers vector with zeros until Simplicial complex dimension and don't allocate a vector of negative
    // size for an empty complex
    std::vector<int> betti_numbers(std::max(dim_max_, 0));

    for (auto pair : persistent_pairs_) {
      // Count never ended persistence intervals
      if (cpx_->null_simplex() == get<1>(pair)) {
        // Increment corresponding betti number
        betti_numbers[cpx_->dimension(get<0>(pair))] += 1;
      }
    }
    return betti_numbers;
  }

  /** @brief Returns the Betti number of the dimension passed by parameter.
   * @param[in] dimension The Betti number dimension to get.
   * @return Betti number of the given dimension
   *
   */
  [[nodiscard]] int betti_number(int dimension) const
  {
    int betti_number = 0;

    for (auto pair : persistent_pairs_) {
      // Count never ended persistence intervals
      if (cpx_->null_simplex() == get<1>(pair)) {
        if (cpx_->dimension(get<0>(pair)) == dimension) {
          // Increment betti number found
          ++betti_number;
        }
      }
    }
    return betti_number;
  }

  /** @brief Returns the persistent Betti numbers.
   * @param[in] from The persistence birth limit to be added in the number \f$(persistent birth \leq from)\f$.
   * @param[in] to The persistence death limit to be added in the number  \f$(persistent death > to)\f$.
   * @return A vector of persistent Betti numbers.
   */
  [[nodiscard]] std::vector<int> persistent_betti_numbers(Filtration_value from, Filtration_value to) const
  {
    // Init Betti numbers vector with zeros until Simplicial complex dimension and don't allocate a vector of negative
    // size for an empty complex
    std::vector<int> betti_numbers(std::max(dim_max_, 0));
    for (auto pair : persistent_pairs_) {
      // Count persistence intervals that covers the given interval
      // null_simplex test : if the function is called with to=+infinity, we still get something useful. And it will
      // still work if we change the complex filtration function to reject null simplices.
      if (cpx_->filtration(get<0>(pair)) <= from &&
          (get<1>(pair) == cpx_->null_simplex() || cpx_->filtration(get<1>(pair)) > to)) {
        // Increment corresponding betti number
        betti_numbers[cpx_->dimension(get<0>(pair))] += 1;
      }
    }
    return betti_numbers;
  }

  /** @brief Returns the persistent Betti number of the dimension passed by parameter.
   * @param[in] dimension The Betti number dimension to get.
   * @param[in] from The persistence birth limit to be added in the number \f$(persistent birth \leq from)\f$.
   * @param[in] to The persistence death limit to be added in the number  \f$(persistent death > to)\f$.
   * @return Persistent Betti number of the given dimension
   */
  [[nodiscard]] int persistent_betti_number(int dimension, Filtration_value from, Filtration_value to) const
  {
    int betti_number = 0;

    for (auto pair : persistent_pairs_) {
      // Count persistence intervals that covers the given interval
      // null_simplex test : if the function is called with to=+infinity, we still get something useful. And it will
      // still work if we change the complex filtration function to reject null simplices.
      if (cpx_->filtration(get<0>(pair)) <= from &&
          (get<1>(pair) == cpx_->null_simplex() || cpx_->filtration(get<1>(pair)) > to)) {
        if (cpx_->dimension(get<0>(pair)) == dimension) {
          // Increment betti number found
          ++betti_number;
        }
      }
    }
    return betti_number;
  }

  /** @brief Returns a list of persistence birth and death FilteredComplex::Simplex_handle pairs.
   * @return A list of Persistent_cohomology::Persistent_interval
   */
  const std::vector<Persistent_interval>& get_persistent_pairs() const { return persistent_pairs_; }

  /** @brief Returns persistence intervals for a given dimension.
   * @param[in] dimension Dimension to get the birth and death pairs from.
   * @return A vector of persistence intervals (birth and death) on a fixed dimension.
   */
  [[nodiscard]] std::vector<std::pair<Filtration_value, Filtration_value>> intervals_in_dimension(int dimension) const
  {
    std::vector<std::pair<Filtration_value, Filtration_value>> result;
    // auto && pair, to avoid unnecessary copying
    for (auto&& pair : persistent_pairs_) {
      if (cpx_->dimension(get<0>(pair)) == dimension) {
        result.emplace_back(cpx_->filtration(get<0>(pair)), cpx_->filtration(get<1>(pair)));
      }
    }
    return result;
  }

 private:
  struct length_interval {
    length_interval(FilteredComplex* cpx, Filtration_value min_length) : cpx_(cpx), min_length_(min_length) {}

    bool operator()(Simplex_handle sh1, Simplex_handle sh2)
    {
      return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_;
    }

    FilteredComplex* cpx_;
    Filtration_value min_length_;
  };

  /*
   * Compare two intervals by length.
   */
  struct cmp_intervals_by_length {
    explicit cmp_intervals_by_length(FilteredComplex* sc) : sc_(sc) {}

    bool operator()(const Persistent_interval& p1, const Persistent_interval& p2)
    {
      return (sc_->filtration(get<1>(p1)) - sc_->filtration(get<0>(p1)) >
              sc_->filtration(get<1>(p2)) - sc_->filtration(get<0>(p2)));
    }

    FilteredComplex* sc_;
  };

  FilteredComplex* cpx_;
  int dim_max_;
  CoefficientField coeff_field_;
  size_t num_simplices_;

  /* Persistent intervals. */
  std::vector<Persistent_interval> persistent_pairs_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_H_
