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

#ifndef PCH_BASE_H_
#define PCH_BASE_H_

#include <initializer_list>
#include <stdexcept>    // std::out_of_range, std::invalid_argument
#include <type_traits>  // std::conditional_t
#include <algorithm>    // std::sort
#include <limits>       // std::numeric_limits
#include <utility>      // std::pair, std::make_pair
#include <tuple>        // std::tuple, std::tie
#include <map>
#include <unordered_map>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistent_cohomology/Compressed_annotation_matrix.h>
#include <gudhi/Persistent_cohomology/Simple_union_find.h>

namespace Gudhi {

namespace persistent_cohomology {

class Persistent_cohomology_no_optimizations
{
 public:
  template <class KeyHandleMap, class Cell_key>
  Persistent_cohomology_no_optimizations([[maybe_unused]] KeyHandleMap* key_handle_map,
                                         [[maybe_unused]] Cell_key num_vertices)
  {}
};

template <class KeyHandleMap>
class Persistent_cohomology_optimizations
{
 public:
  /** \brief Data stored for each simplex. */
  using Cell_key = typename KeyHandleMap::Simplex_key;
  /** \brief Handle to specify a simplex. */
  using Cell_handle = typename KeyHandleMap::Simplex_handle;

  /** \brief Initializes the Persistent_cohomology_base class.
   *
   * @param[in] cpx Complex for which the persistent homology is computed.
   * cpx is a model of KeyHandleMap
   *
   * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the
   *                                complex is computed. If false, it is ignored. Default is false.
   *
   * @exception std::out_of_range In case the number of simplices is more than Cell_key type numeric limit.
   */
  explicit Persistent_cohomology_optimizations(KeyHandleMap* key_handle_map, Cell_key num_vertices)
      : key_handle_map_(key_handle_map),
        current_vertex_index_(0),
        vertex_sets_(num_vertices),
        vertices_(num_vertices),
        zero_cocycles_()
  {}

  void update_cohomology_groups_vertex(Cell_handle vertex)
  {
    vertex_sets_.make_set(current_vertex_index_);
    vertices_[current_vertex_index_] = vertex;
    ++current_vertex_index_;
  }

  template <class F1, class F2>
  void update_cohomology_groups_edge(Cell_handle sigma, Cell_key u, Cell_key v, F1&& add_pair, F2&& create_cocycle)
  {
    Cell_key ku = vertex_sets_.find_set(u);
    Cell_key kv = vertex_sets_.find_set(v);

    if (ku != kv) {  // Destroy a connected component
      vertex_sets_.merge_sets(ku, kv);
      // Keys of the simplices which created the connected components containing
      // respectively u and v.
      Cell_key idx_coc_u, idx_coc_v;
      auto map_it_u = zero_cocycles_.find(ku);
      // If the index of the cocycle representing the class is already ku.
      if (map_it_u == zero_cocycles_.end()) {
        idx_coc_u = ku;
      } else {
        idx_coc_u = map_it_u->second;
      }

      auto map_it_v = zero_cocycles_.find(kv);
      // If the index of the cocycle representing the class is already kv.
      if (map_it_v == zero_cocycles_.end()) {
        idx_coc_v = kv;
      } else {
        idx_coc_v = map_it_v->second;
      }

      if (idx_coc_u < idx_coc_v) {  // Kill cocycle [idx_coc_v], which is younger.
        std::forward<F1>(add_pair)(vertices_[idx_coc_v], sigma);
        // Maintain the index of the 0-cocycle alive.
        if (kv != idx_coc_v) {
          zero_cocycles_.erase(map_it_v);
        }
        if (kv == vertex_sets_.find_set(kv)) {
          if (ku != idx_coc_u) {
            zero_cocycles_.erase(map_it_u);
          }
          zero_cocycles_[kv] = idx_coc_u;
        }
      } else {  // Kill cocycle [idx_coc_u], which is younger.
        std::forward<F1>(add_pair)(vertices_[idx_coc_u], sigma);
        // Maintain the index of the 0-cocycle alive.
        if (ku != idx_coc_u) {
          zero_cocycles_.erase(map_it_u);
        }
        if (ku == vertex_sets_.find_set(ku)) {
          if (kv != idx_coc_v) {
            zero_cocycles_.erase(map_it_v);
          }
          zero_cocycles_[ku] = idx_coc_v;
        }
      }
      key_handle_map_->assign_key(sigma, key_handle_map_->null_key());
    } else {  // If ku == kv, same connected component: create a 1-cocycle class.
      std::forward<F2>(create_cocycle)(sigma);
    }
  }

  template <class F>
  void finalize_computation(F&& add_pair)
  {
    Cell_key i = 0;
    for (const Cell_handle& key : vertices_) {
      if (vertex_sets_.get_parent(i) == i  // root of its tree
          && zero_cocycles_.find(i) == zero_cocycles_.end()) {
        std::forward<F>(add_pair)(key);
      }
      ++i;
    }
    for (const auto& zero_idx : zero_cocycles_) {
      std::forward<F>(add_pair)(vertices_[zero_idx.second]);
    }
  }

 private:
  KeyHandleMap* key_handle_map_;

  Cell_key current_vertex_index_;
  /*  Disjoint sets data structure to link the model of KeyHandleMap
   * with the compressed annotation matrix.
   * ds_rank_ is a property map Cell_key -> int, ds_parent_ is a property map
   * Cell_key -> simplex_key_t */
  Simple_union_find<Cell_key> vertex_sets_;
  // maps vertex position wrt other vertices to position in the whole filtration.
  // remains empty of no optimization is used.
  std::vector<Cell_handle> vertices_;
  /*  Dictionary establishing the correspondence between the Cell_key of
   * the root vertex in the union-find ds and the Cell_key of the vertex which
   * created the connected component as a 0-dimension homology feature.*/
  std::unordered_map<Cell_key, Cell_key> zero_cocycles_;
};

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
template <class KeyHandleMap, class CoefficientField, class LengthInterval, bool with_optimizations = true>
class Persistent_cohomology_base : std::conditional_t<with_optimizations,
                                                      Persistent_cohomology_optimizations<KeyHandleMap>,
                                                      Persistent_cohomology_no_optimizations>
{
 public:
  /** \brief Data stored for each simplex. */
  using Cell_key = typename KeyHandleMap::Simplex_key;
  /** \brief Handle to specify a simplex. */
  using Cell_handle = typename KeyHandleMap::Simplex_handle;
  /** \brief Type of element of the field. */
  using Arith_element = typename CoefficientField::Element;
  /** \brief Type for birth and death KeyHandleMap::Cell_handle.
   * The Arith_element field is used for the multi-field framework. */
  using Persistent_interval = std::tuple<Cell_handle, Cell_handle, Arith_element>;

 private:
  using Base = std::conditional_t<with_optimizations,
                                  Persistent_cohomology_optimizations<KeyHandleMap>,
                                  Persistent_cohomology_no_optimizations>;
  using Cam = Compressed_annotation_matrix<Cell_key, CoefficientField>;
  // Column type
  using Column = typename Cam::Column;
  // Sparse column type for the annotation of the boundary of an element.
  using A_ds_type = std::vector<std::pair<Cell_key, Arith_element>>;

 public:
  /** \brief Initializes the Persistent_cohomology_base class.
   *
   * @param[in] cpx Complex for which the persistent homology is computed.
   * cpx is a model of KeyHandleMap
   *
   * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the
   *                                complex is computed. If false, it is ignored. Default is false.
   *
   * @exception std::out_of_range In case the number of simplices is more than Cell_key type numeric limit.
   */
  explicit Persistent_cohomology_base(KeyHandleMap* key_handle_map,
                                      const CoefficientField& coeff_field,
                                      const LengthInterval& interval_length_policy,
                                      Cell_key num_cells,
                                      Cell_key num_vertices)
      : Base(key_handle_map, num_vertices),
        coeff_field_(coeff_field),
        key_handle_map_(key_handle_map),
        current_index_(0),
        cam_(num_cells, coeff_field),    // collection of annotation vectors
        transverse_idx_(),  // key -> row
        interval_length_policy_(interval_length_policy),
        persistent_pairs_()
  {
    if (num_cells > std::numeric_limits<Cell_key>::max()) {
      // num_simplices must be strictly lower than the limit, because a value is reserved for null_key.
      throw std::out_of_range("The number of simplices is more than Cell_key type numeric limit.");
    }
  }

  void update_cohomology_groups_vertex(Cell_handle vertex)
  {
    if constexpr (with_optimizations) {
      Base::update_cohomology_groups_vertex(vertex);
    } else {
      update_cohomology_groups(vertex, {}, false);
    }
  }

  /** \brief Update the cohomology groups under the insertion of an edge.
   *
   * The 0-homology is maintained with a simple Union-Find data structure, which
   * explains the existence of a specific function of edge insertions. */
  void update_cohomology_groups_edge(Cell_handle edge, Cell_handle u, Cell_handle v, const bool kill_only)
  {
    if constexpr (with_optimizations) {
      Base::update_cohomology_groups_edge(
          edge,
          key_handle_map_->key(u),
          key_handle_map_->key(v),
          [&](Cell_handle birth, Cell_handle death) {
            if (interval_length_policy_(birth, death)) {
              persistent_pairs_.emplace_back(birth, death, coeff_field_.get_characteristic());
            }
          },
          [&](Cell_handle birth) {
            if (!kill_only) {  // If ku == kv, same connected component: create a 1-cocycle class.
              create_cocycle(birth, coeff_field_.get_multiplicative_identity(), coeff_field_.get_characteristic());
            }
          });
    } else {
      update_cohomology_groups(edge, {u, v}, kill_only);
    }
  }

  /*
   * Update the cohomology groups under the insertion of a simplex.
   */
  template <class Simplex_boundary_range = std::initializer_list<Cell_handle>>
  void update_cohomology_groups(Cell_handle sigma, const Simplex_boundary_range& boundary, const bool kill_only)
  {
    // Compute the annotation of the boundary of sigma:
    std::map<Cell_key, Arith_element> map_a_ds;
    annotation_of_the_boundary(map_a_ds, boundary);
    // Update the cohomology groups:
    if (map_a_ds.empty()) {  // sigma is a creator in all fields represented in coeff_field_
      if (!kill_only) {
        create_cocycle(sigma, coeff_field_.get_multiplicative_identity(), coeff_field_.get_characteristic());
      }
    } else {  // sigma is a destructor in at least a field in coeff_field_
      // Convert map_a_ds to a vector
      A_ds_type a_ds;  // admits reverse iterators
      for (const auto& map_a_ds_ref : map_a_ds) {
        a_ds.push_back(std::pair<Cell_key, Arith_element>(map_a_ds_ref.first, map_a_ds_ref.second));
      }

      Arith_element inv_x, charac;
      Arith_element prod = coeff_field_.get_characteristic();  // Product of characteristic of the fields
      for (auto a_ds_rit = a_ds.rbegin();
           (a_ds_rit != a_ds.rend()) && (prod != coeff_field_.get_multiplicative_identity());
           ++a_ds_rit) {
        std::tie(inv_x, charac) = coeff_field_.get_partial_inverse(a_ds_rit->second, prod);

        if (inv_x != coeff_field_.get_additive_identity()) {
          destroy_cocycle(sigma, a_ds, a_ds_rit->first, inv_x, charac);
          prod /= charac;
        }
      }
      if (prod != coeff_field_.get_multiplicative_identity() && !kill_only) {
        create_cocycle(sigma, coeff_field_.get_partial_multiplicative_identity(prod), prod);
      }
    }
  }

  void finalize_computation()
  {
    if constexpr (with_optimizations) {
      Base::finalize_computation([&](Cell_handle birth) {
        persistent_pairs_.emplace_back(birth, key_handle_map_->null_simplex(), coeff_field_.get_characteristic());
      });
    }
    for (const auto& cocycle : transverse_idx_) {
      persistent_pairs_.emplace_back(
          key_handle_map_->simplex(cocycle.first), key_handle_map_->null_simplex(), cocycle.second);
    }
  }

  /** @brief Returns a list of persistence birth and death KeyHandleMap::Cell_handle pairs.
   * @return A list of Persistent_cohomology_base::Persistent_interval
   */
  const std::vector<Persistent_interval>& get_persistent_pairs() const { return persistent_pairs_; }

 private:
  CoefficientField coeff_field_;
  KeyHandleMap* key_handle_map_;

  Cell_key current_index_;
  /* The compressed annotation matrix fields.*/
  Cam cam_;
  /*  Key -> row. */
  std::map<Cell_key, Arith_element> transverse_idx_;

  LengthInterval interval_length_policy_;
  /* Persistent intervals. */
  std::vector<Persistent_interval> persistent_pairs_;

  /*
   * Compute the annotation of the boundary of a simplex.
   */
  template <class Simplex_boundary_range>
  void annotation_of_the_boundary(std::map<Cell_key, Arith_element>& map_a_ds, const Simplex_boundary_range& boundary)
  {
    // traverses the boundary of sigma, keeps track of the annotation vectors,
    // with multiplicity. We used to sum the coefficients directly in
    // annotations_in_boundary by using a map, we now do it later.
    using annotation_t = std::pair<Column const*, int>;
    thread_local std::vector<annotation_t> annotations_in_boundary;
    annotations_in_boundary.clear();
    int sign = -1;  // \in {-1,1} provides the sign in the alternate sum in the boundary.
    Cell_key key;

    for (const auto& sh : boundary) {
      key = key_handle_map_->key(sh);
      if (key != key_handle_map_->null_key()) {  // A simplex with null_key is a killer, and have null annotation
        // Find its annotation vector
        const auto& col = cam_.get_column(key);
        if (!col.is_empty()) {
          // and insert it in annotations_in_boundary with multiplicative factor "sign".
          annotations_in_boundary.emplace_back(&col, sign);
        }
      }
      sign = -sign;
    }
    // Place identical annotations consecutively so we can easily sum their multiplicities.
    std::sort(annotations_in_boundary.begin(),
              annotations_in_boundary.end(),
              [](annotation_t const& a, annotation_t const& b) { return a.first < b.first; });

    // Sum the annotations with multiplicity, using a map<key,coeff>
    // to represent a sparse vector.
    std::pair<typename std::map<Cell_key, Arith_element>::iterator, bool> result_insert_a_ds;

    for (auto ann_it = annotations_in_boundary.begin(); ann_it != annotations_in_boundary.end(); /**/) {
      Column const* col = ann_it->first;
      int mult = ann_it->second;
      ++ann_it;
      while (ann_it != annotations_in_boundary.end() && ann_it->first == col) {
        mult += ann_it->second;
        ++ann_it;
      }
      Arith_element mult_coeff = coeff_field_.get_value(mult);
      // The following test is just a heuristic, it is not required, and it is fine that is misses p == 0.
      if (mult_coeff != coeff_field_.get_additive_identity()) {  // For all columns in the boundary,
        for (const auto& entry_ref : col->col_) {                // insert every entry in map_a_ds with multiplicity
          Arith_element w_y = coeff_field_.multiply(entry_ref.coefficient_, mult_coeff);  // coefficient * multiplicity

          if (w_y != coeff_field_.get_additive_identity()) {  // if != 0
            result_insert_a_ds = map_a_ds.insert(std::pair<Cell_key, Arith_element>(entry_ref.key_, w_y));
            if (!(result_insert_a_ds.second)) {  // if entry_ref.key_ already a Key in map_a_ds
              result_insert_a_ds.first->second = coeff_field_.add(result_insert_a_ds.first->second, w_y);
              if (result_insert_a_ds.first->second == coeff_field_.get_additive_identity()) {
                map_a_ds.erase(result_insert_a_ds.first);
              }
            }
          }
        }
      }
    }
  }

  /*  \brief Create a new cocycle class.
   *
   * The class is created by the insertion of the simplex sigma.
   * The methods adds a cocycle, representing the new cocycle class,
   * to the matrix representing the cohomology groups.
   * The new cocycle has value 0 on every simplex except on sigma
   * where it worths 1.*/
  void create_cocycle(Cell_handle sigma, Arith_element x, Arith_element charac)
  {
    Cell_key key = key_handle_map_->key(sigma);
    cam_.insert_column(key, x);
    transverse_idx_[key] = charac;  // insert the new row
  }

  /*  \brief Destroy a cocycle class.
   *
   * The cocycle class is destroyed by the insertion of sigma.
   * The methods proceeds to a reduction of the matrix representing
   * the cohomology groups using Gauss pivoting. The reduction zeros-out
   * the row containing the entry with highest key in
   * a_ds, the annotation of the boundary of simplex sigma. This key
   * is "death_key".*/
  void destroy_cocycle(Cell_handle sigma,
                       A_ds_type const& a_ds,
                       Cell_key death_key,
                       Arith_element inv_x,
                       Arith_element charac)
  {
    // Create a finite persistent interval for which the interval exists
    if (interval_length_policy_(key_handle_map_->simplex(death_key), sigma)) {
      persistent_pairs_.emplace_back(key_handle_map_->simplex(death_key)  // creator
                                     ,
                                     sigma  // destructor
                                     ,
                                     charac);  // fields
    }

    const auto& death_key_row = cam_.get_row(death_key);  // Find the beginning of the row.
    auto death_cocycle_coeff = transverse_idx_.find(death_key);

    auto row_entry_it = death_key_row.begin();

    while (row_entry_it != death_key_row.end()) {  // Traverse all entries in the row at index death_key.
      const auto& entry = *row_entry_it;
      Arith_element w = coeff_field_.multiply(inv_x, -entry.get_element());
      // Arith_element w = coeff_field_.times_minus(inv_x, entry.get_element());

      if (w != coeff_field_.get_additive_identity()) {
        ++row_entry_it;  // has to be done before column addition to avoid invalidating the iterator
        cam_.multiply_source_and_add_to(w, a_ds, entry.get_column_index());
      } else {
        ++row_entry_it;
      }  // If w == 0, pass.
    }

    // Because it is a killer simplex, set the data of sigma to null_key().
    if (charac == coeff_field_.get_characteristic()) {
      key_handle_map_->assign_key(sigma, key_handle_map_->null_key());
    }
    if (death_cocycle_coeff->second == charac) {
      cam_.erase_empty_row(death_key);
      transverse_idx_.erase(death_cocycle_coeff);
    } else {
      death_cocycle_coeff->second /= charac;
    }
  }
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PCH_BASE_H_
