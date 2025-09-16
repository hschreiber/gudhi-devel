/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_H_
#define PERSISTENT_COHOMOLOGY_H_

#include <initializer_list>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>
#include <list>
#include <vector>
#include <set>
#include <fstream>  // std::ofstream
#include <limits>   // for numeric_limits<>
#include <tuple>
#include <algorithm>
#include <string>
#include <stdexcept>  // for std::out_of_range

#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/intrusive/list.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistent_cohomology/Persistent_cohomology_column.h>
#include <gudhi/Simple_object_pool.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h> //only here for external use....

namespace Gudhi {

namespace persistent_cohomology {

template <typename Simplex_key>
struct Cell_union_find {
 public:
  Cell_union_find(int num_simplices)
      : ds_rank_(num_simplices), ds_parent_(num_simplices), dsets_(ds_rank_.data(), ds_parent_.data())
  {}

  void make_set(Simplex_key key) { dsets_.make_set(key); }

  Simplex_key find_set(Simplex_key key) { return dsets_.find_set(key); }

  void merge_sets(Simplex_key k1, Simplex_key k2) { dsets_.link(k1, k2); }

  Simplex_key get_parent(Simplex_key key) const { return ds_parent_[key]; }

 private:
  std::vector<int> ds_rank_;
  std::vector<Simplex_key> ds_parent_;
  boost::disjoint_sets<int*, Simplex_key*> dsets_;
};

template <typename Simplex_key, class CoefficientField>
class Compressed_annotation_matrix

{
 public:
  using Characteristic = typename CoefficientField::Element;
  using Index = Simplex_key;
  using ID_index = unsigned int;
  using Column = Persistent_cohomology_column<Index, Characteristic>;
  // Entry type
  using Entry = typename Column::Entry;  // contains 2 list_hooks
  // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
  using Row = boost::intrusive::
      list<Entry, boost::intrusive::constant_time_size<false>, boost::intrusive::base_hook<base_hook_cam_h>>;

  using Cam = boost::intrusive::set<Column, boost::intrusive::constant_time_size<false>>;

  /**
   * @brief Constructs a new empty matrix and reserves space for the given number of columns.
   *
   * @param numberOfColumns Number of columns to reserve space for.
   * @param characteristic Characteristic of the coefficient field. If not specified and
   * @ref PersistenceMatrixOptions::is_z2 is false, the characteristic has to be set later with the use of
   * @ref set_characteristic before calling for the first time a method needing it. Ignored if
   * @ref PersistenceMatrixOptions::is_z2 is true.
   */
  Compressed_annotation_matrix(unsigned int numberOfColumns)
      : coeff_field_(),
        ds_repr_(numberOfColumns, nullptr),
        transverse_idx_(),
        cam_(),
        cell_sets_(numberOfColumns),
        column_pool_(),
        entry_pool_()
  {}

  /**
   * @brief Copy constructor.
   *
   * @param matrixToCopy %Compressed_annotation_matrix to copy.
   */
  Compressed_annotation_matrix(const Compressed_annotation_matrix& matrixToCopy) = delete;

  /**
   * @brief Move constructor.
   * After the move, the given matrix will be empty.
   *
   * @param other %Compressed_annotation_matrix to move.
   */
  Compressed_annotation_matrix(Compressed_annotation_matrix&& other) noexcept
      : coeff_field_(std::move(other.coeff_field_)),
        ds_repr_(std::move(other.ds_repr_)),
        transverse_idx_(std::move(other.transverse_idx_)),
        cam_(std::move(other.cam_)),
        cell_sets_(std::move(other.cell_sets_)),
        column_pool_(std::move(other.column_pool_)),
        entry_pool_(std::move(other.entry_pool_))
  {}

  ~Compressed_annotation_matrix()
  {
    // Clean the transversal lists
    for (auto& row : transverse_idx_) {
      // Destruct all the entries
      row.second.clear_and_dispose([&](Entry* p) { p->~Entry(); });
    }
  }

  // TODO: compatibility with multi fields:
  //   - set_characteristic(Characteristic min, Characteristic max)
  //   - readapt reduction?
  /**
   * @brief Sets the characteristic of the coefficient field if @ref PersistenceMatrixOptions::is_z2 is false,
   * does nothing otherwise.
   * Should be used if no characteristic could be specified at the creation of the empty matrix.
   * Do not change the value of the characteristic once used.
   *
   * @warning The coefficient values stored in the matrix are stored after computing the corresponding modulo.
   * Therefore, changing the characteristic after is very likely to invalidate all entry values.
   *
   * @param characteristic The characteristic to set.
   */
  void set_characteristic(Characteristic characteristic) { coeff_field_.init(characteristic); }

  void set_characteristic(int charac_min, int charac_max) { coeff_field_.init(charac_min, charac_max); }

  // (TODO: if there is no row access and the column type corresponds to the internal column type of the matrix,
  // moving the column instead of copying it should be possible. Is it worth implementing it?)
  /**
   * @brief Inserts a new ordered column at the end of the matrix by copying the given range of
   * @ref Entry_representative. The content of the range is assumed to be sorted by increasing ID value.
   *
   * Only available for @ref basematrix "base matrices".
   * Otherwise use @ref insert_boundary which will deduce a new column from the boundary given.
   *
   * @tparam Container Range of @ref Entry_representative. Assumed to have a begin(), end() and size() method.
   * @param column Column to be inserted.
   */
  template <class Container = std::initializer_list<std::pair<Index, Characteristic>>>
  void insert_column(const Container& column)
  {
    // Should only be called on a single entry column
    GUDHI_CHECK(column.size() == 1, "Internal problem: column should contain only one element.");
    auto element = *column.begin();
    Index key = element.first;
    Characteristic x = element.second;
    // Create a column containing only one entry,
    Column* new_col = column_pool_.construct(key);
    Entry* new_entry = entry_pool_.construct(key, x, new_col);
    new_col->col_.push_back(*new_entry);
    // and insert it in the matrix, in constant time thanks to the hint cam_.end().
    // Indeed *new_col has the biggest lexicographic value because key is the
    // biggest key used so far.
    cam_.insert(cam_.end(), *new_col);
    // Update the disjoint sets data structure.
    transverse_idx_[key] = Row();  // insert the new row
    transverse_idx_[key].push_back(*new_entry);
    cell_sets_.make_set(key);
    ds_repr_[key] = new_col;
  }

  /**
   * @brief Only available for @ref chainmatrix "chain matrices". Returns the column at the given @ref MatIdx index.
   * The type of the column depends on the chosen options, see @ref PersistenceMatrixOptions::column_type.
   *
   * @param columnIndex @ref MatIdx index of the column to return.
   * @return Const reference to the column.
   */
  const Column& get_column(Index columnIndex) { return *ds_repr_[cell_sets_.find_set(columnIndex)]; }

  const Column& get_column(Column* col) { return *col; }

  /**
   * @brief Only available for @ref chainmatrix "chain matrices" and matrices with column compression.
   * Returns the row at the given @ref rowindex "row index".
   * The type of the row depends on the chosen options, see @ref PersistenceMatrixOptions::has_intrusive_rows.
   *
   * @param rowIndex @ref rowindex "Row index" of the row to return: @ref IDIdx for @ref chainmatrix "chain matrices"
   * or updated @ref IDIdx for @ref boundarymatrix "boundary matrices" if swaps occurred.
   * @return Const reference to the row.
   */
  const Row& get_row(ID_index rowIndex) const { return transverse_idx_.find(rowIndex)->second; }

  // TODO: rename method to be less confusing.
  /**
   * @brief The effect varies depending on the matrices and the options:
   * - @ref basematrix "base matrix" and @ref boundarymatrix "boundary matrix":
   *    - @ref PersistenceMatrixOptions::has_map_column_container and @ref
   *      PersistenceMatrixOptions::has_column_and_row_swaps are true: cleans up maps used for the lazy row swaps.
   *    - @ref PersistenceMatrixOptions::has_row_access and @ref PersistenceMatrixOptions::has_removable_rows are true:
   *      assumes that the row is empty and removes it.
   *    - Otherwise, does nothing.
   *
   * @warning The removed rows are always assumed to be empty. If it is not the case, the deleted row entries are not
   * removed from their columns. And in the case of intrusive rows, this will generate a segmentation fault when
   * the column entries are destroyed later. The row access is just meant as a "read only" access to the rows and the
   * @ref erase_empty_row method just as a way to specify that a row is empty and can therefore be removed from
   * dictionaries. This allows to avoid testing the emptiness of a row at each column entry removal, what can be quite
   * frequent.
   *
   * @param rowIndex @ref rowindex "Row index" of the empty row to remove.
   */
  void erase_empty_row(ID_index rowIndex)
  {
    // row should be empty
    GUDHI_CHECK(transverse_idx_.find(rowIndex)->second.size() == 0, "Internal problem: row should be empty.");
    transverse_idx_.erase(rowIndex);
  }

  /**
   * @brief Multiplies the source column with the coefficient before adding it to the target column.
   * That is: `targetColumn += (coefficient * sourceColumn)`. The source column will **not** be modified.
   *
   * The representatives of redundant columns are summed together, which means that
   * all column compressed together with the target column are affected by the change, not only the target.
   *
   * @tparam Entry_range_or_column_index Either a range of @ref Entry with a begin() and end() method,
   * or any integer type.
   * @param coefficient Value to multiply.
   * @param sourceColumn Either a @ref Entry range or the @ref MatIdx index of the column to add.
   * @param targetColumnIndex @ref MatIdx index of the target column.
   */
  template <class Entry_range>
  void multiply_source_and_add_to(Characteristic coefficient, const Entry_range& sourceColumn, Column const* constTargetColumn)
  {
    Column * targetColumn = const_cast<Column *>(constTargetColumn);

    // Disconnect the column from the rows in the CAM.
    for (auto& col_entry : targetColumn->col_) {
      col_entry.base_hook_cam_h::unlink();
    }

    // Remove the column from the CAM before modifying its value
    cam_.erase(cam_.iterator_to(*targetColumn));
    // Proceed to the reduction of the column
    plus_equal_column(*targetColumn, sourceColumn, coefficient);

    if (targetColumn->col_.empty()) {  // If the column is null
      ds_repr_[targetColumn->class_key_] = nullptr;
      column_pool_.destroy(targetColumn);  // delete targetColumn;
    } else {
      std::pair<typename Cam::iterator, bool> result_insert_cam;

      // Find whether the column obtained is already in the CAM
      result_insert_cam = cam_.insert(*targetColumn);
      if (result_insert_cam.second) {  // If it was not in the CAM before: insertion has succeeded
        for (auto& col_entry : targetColumn->col_) {
          // re-establish the row links
          transverse_idx_[col_entry.key_].push_front(col_entry);
        }
      } else {  // There is already an identical column in the CAM:
        // merge two disjoint sets.
        cell_sets_.merge_sets(targetColumn->class_key_, result_insert_cam.first->class_key_);

        Simplex_key key_tmp = cell_sets_.find_set(targetColumn->class_key_);
        ds_repr_[key_tmp] = &(*(result_insert_cam.first));
        result_insert_cam.first->class_key_ = key_tmp;
        // intrusive containers don't own their elements, we have to release them manually
        targetColumn->col_.clear_and_dispose([&](Entry* p) { entry_pool_.destroy(p); });
        column_pool_.destroy(targetColumn);  // delete targetColumn;
      }
    }
  }

  /**
   * @brief Indicates if the column at given index has value zero.
   *
   * For @ref boundarymatrix "RU matrices", equivalent to
   * @ref is_zero_column(Index columnIndex, bool inR) "is_zero_column(columnIndex, true)".
   *
   * Note that for @ref chainmatrix "chain matrices", this method should always return false, as a valid
   * @ref chainmatrix "chain matrix" never has empty columns.
   *
   * @param columnIndex @ref MatIdx index of the column.
   * @return true If the column has value zero.
   * @return false Otherwise.
   */
  bool is_zero_column(Index columnIndex)
  {
    auto col = ds_repr_[cell_sets_.find_set(columnIndex)];
    if (col == nullptr) return true;
    return col->is_empty();
  }

  /**
   * @brief Assign operator.
   *
   * @param other %Compressed_annotation_matrix to copy
   * @return Reference to this object.
   */
  Compressed_annotation_matrix& operator=(Compressed_annotation_matrix other) = delete;

  Compressed_annotation_matrix& operator=(Compressed_annotation_matrix&& other) noexcept
  {
    if (&ds_repr_ == &(other.ds_repr_)) return *this;

    for (auto& row : transverse_idx_) {
      // Destruct all potential entries
      row.second.clear_and_dispose([&](Entry* p) { p->~Entry(); });
    }

    coeff_field_ = std::move(other.coeff_field_);
    ds_repr_ = std::move(other.ds_repr_);
    transverse_idx_ = std::move(other.transverse_idx_);
    cam_ = std::move(other.cam_);
    cell_sets_ = std::move(other.cell_sets_);
    column_pool_ = std::move(other.column_pool_);
    entry_pool_ = std::move(other.entry_pool_);
  }

  void print();  // for debug

 private:
  CoefficientField coeff_field_;
  std::vector<Column*> ds_repr_;
  std::map<Simplex_key, Row> transverse_idx_;
  Cam cam_;
  Cell_union_find<Simplex_key> cell_sets_;

  Simple_object_pool<Column> column_pool_;
  Simple_object_pool<Entry> entry_pool_;

  /*
   * Assign:    target <- target + w * other.
   */
  template <class Entry_range>
  void plus_equal_column(Column& target,
                         Entry_range const& other  // value_type is pair<Simplex_key,Arith_element>
                         ,
                         Characteristic w)
  {
    auto target_it = target.col_.begin();
    auto other_it = other.begin();
    while (target_it != target.col_.end() && other_it != other.end()) {
      if (target_it->key_ < other_it->first) {
        ++target_it;
      } else {
        if (target_it->key_ > other_it->first) {
          Entry* entry_tmp = entry_pool_.construct(Entry(other_it->first  // key
                                                     ,
                                                     coeff_field_.additive_identity(),
                                                     &target));

          entry_tmp->coefficient_ = coeff_field_.plus_times_equal(entry_tmp->coefficient_, other_it->second, w);

          target.col_.insert(target_it, *entry_tmp);

          ++other_it;
        } else {  // it1->key == it2->key
          // target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
          target_it->coefficient_ = coeff_field_.plus_times_equal(target_it->coefficient_, other_it->second, w);
          if (target_it->coefficient_ == coeff_field_.additive_identity()) {
            auto tmp_it = target_it;
            ++target_it;
            ++other_it;  // iterators remain valid
            Entry* tmp_entry_ptr = &(*tmp_it);
            target.col_.erase(tmp_it);  // removed from column

            entry_pool_.destroy(tmp_entry_ptr);  // delete from memory
          } else {
            ++target_it;
            ++other_it;
          }
        }
      }
    }
    while (other_it != other.end()) {
      Entry* entry_tmp = entry_pool_.construct(Entry(other_it->first, coeff_field_.additive_identity(), &target));
      entry_tmp->coefficient_ = coeff_field_.plus_times_equal(entry_tmp->coefficient_, other_it->second, w);
      target.col_.insert(target.col_.end(), *entry_tmp);

      ++other_it;
    }
  }
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
  using Persistent_interval = std::tuple<Simplex_handle, Simplex_handle, Arith_element>;

 private:
  using Cam = Compressed_annotation_matrix<Simplex_key, CoefficientField>;
  // Column type
  using Column = typename Cam::Column;
  // Sparse column type for the annotation of the boundary of an element.
  using A_ds_type = std::vector<std::pair<Simplex_key, Arith_element>>;

 public:
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
        dim_max_(cpx.dimension()),              // upper bound on the dimension of the simplices
        coeff_field_(),                         // initialize the field coefficient structure.
        num_simplices_(cpx_->num_simplices()),  // num_simplices save to avoid to call thrice the function
        vertex_sets_(num_simplices_), // TODO: super wasteful
        cam_(num_simplices_),  // collection of annotation vectors
        zero_cocycles_(),      // union-find -> Simplex_key of creator for 0-homology
        transverse_idx_(),     // key -> row
        persistent_pairs_(),
        interval_length_policy(&cpx, 0)
  {
    if (num_simplices_ > std::numeric_limits<Simplex_key>::max()) {
      // num_simplices must be strictly lower than the limit, because a value is reserved for null_key.
      throw std::out_of_range("The number of simplices is more than Simplex_key type numeric limit.");
    }
    if (persistence_dim_max) {
      ++dim_max_;
    }
  }

 private:
  struct length_interval {
    length_interval(FilteredComplex* cpx, Filtration_value min_length) : cpx_(cpx), min_length_(min_length) {}

    bool operator()(Simplex_handle sh1, Simplex_handle sh2)
    {
      return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_;
    }

    void set_length(Filtration_value new_length) { min_length_ = new_length; }

    FilteredComplex* cpx_;
    Filtration_value min_length_;
  };

 public:
  /** \brief Initializes the coefficient field.*/
  void init_coefficients(int charac)
  {
    coeff_field_.init(charac);
    cam_.set_characteristic(charac);
  }

  /** \brief Initializes the coefficient field for multi-field persistent homology.*/
  void init_coefficients(int charac_min, int charac_max)
  {
    coeff_field_.init(charac_min, charac_max);
    cam_.set_characteristic(charac_min, charac_max);
  }

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
    if (coeff_field_.characteristic() == CoefficientField::nullCharacteristic)
      throw std::logic_error("Coefficient field was not initialized! Please call `init_coefficients` first.");

    if (dim_max_ <= 0) return;  // --------->>

    interval_length_policy.set_length(min_interval_length);
    Simplex_key idx_fil = -1;
    std::vector<Simplex_key> vertices;  // so we can check the connected components at the end
    // Compute all finite intervals
    for (auto sh : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(sh, ++idx_fil);
      int dim_simplex = cpx_->dimension(sh);
      switch (dim_simplex) {
        case 0:
          vertex_sets_.make_set(idx_fil);
          vertices.push_back(idx_fil);
          break;
        case 1:
          update_cohomology_groups_edge(sh);
          break;
        default:
          update_cohomology_groups(sh, dim_simplex);
          break;
      }
    }
    // Compute infinite intervals of dimension 0
    for (Simplex_key key : vertices) {         // for all 0-dimensional simplices
      if (vertex_sets_.get_parent(key) == key  // root of its tree
          && zero_cocycles_.find(key) == zero_cocycles_.end()) {
        persistent_pairs_.emplace_back(cpx_->simplex(key), cpx_->null_simplex(), coeff_field_.characteristic());
      }
    }
    for (auto zero_idx : zero_cocycles_) {
      persistent_pairs_.emplace_back(
          cpx_->simplex(zero_idx.second), cpx_->null_simplex(), coeff_field_.characteristic());
    }
    // Compute infinite interval of dimension > 0
    for (const auto& cocycle : transverse_idx_) {
      persistent_pairs_.emplace_back(cpx_->simplex(cocycle.first), cpx_->null_simplex(), cocycle.second);
    }
  }

 private:
  /** \brief Update the cohomology groups under the insertion of an edge.
   *
   * The 0-homology is maintained with a simple Union-Find data structure, which
   * explains the existence of a specific function of edge insertions. */
  void update_cohomology_groups_edge(Simplex_handle sigma)
  {
    Simplex_handle u, v;
    boost::tie(u, v) = cpx_->endpoints(sigma);

    Simplex_key ku = vertex_sets_.find_set(cpx_->key(u));
    Simplex_key kv = vertex_sets_.find_set(cpx_->key(v));

    if (ku != kv) {  // Destroy a connected component
      vertex_sets_.merge_sets(ku, kv);
      // Keys of the simplices which created the connected components containing
      // respectively u and v.
      Simplex_key idx_coc_u, idx_coc_v;
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

      if (cpx_->filtration(cpx_->simplex(idx_coc_u)) <
          cpx_->filtration(cpx_->simplex(idx_coc_v))) {  // Kill cocycle [idx_coc_v], which is younger.
        if (interval_length_policy(cpx_->simplex(idx_coc_v), sigma)) {
          persistent_pairs_.emplace_back(cpx_->simplex(idx_coc_v), sigma, coeff_field_.characteristic());
        }
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
        if (interval_length_policy(cpx_->simplex(idx_coc_u), sigma)) {
          persistent_pairs_.emplace_back(cpx_->simplex(idx_coc_u), sigma, coeff_field_.characteristic());
        }
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
      cpx_->assign_key(sigma, cpx_->null_key());
    } else if (dim_max_ > 1) {  // If ku == kv, same connected component: create a 1-cocycle class.
      create_cocycle(sigma, coeff_field_.multiplicative_identity(), coeff_field_.characteristic());
    }
  }

  /*
   * Compute the annotation of the boundary of a simplex.
   */
  void annotation_of_the_boundary(std::map<Simplex_key, Arith_element>& map_a_ds, Simplex_handle sigma, int dim_sigma)
  {
    // traverses the boundary of sigma, keeps track of the annotation vectors,
    // with multiplicity. We used to sum the coefficients directly in
    // annotations_in_boundary by using a map, we now do it later.
    using annotation_t = std::pair<Column const*, int>;
    thread_local std::vector<annotation_t> annotations_in_boundary;
    annotations_in_boundary.clear();
    int sign = 1 - (2 * (dim_sigma % 2));  // \in {-1,1} provides the sign in the
                                           // alternate sum in the boundary.
    Simplex_key key;

    for (auto sh : cpx_->boundary_simplex_range(sigma)) {
      key = cpx_->key(sh);
      if (key != cpx_->null_key()) {  // A simplex with null_key is a killer, and have null annotation
        // Find its annotation vector
        if (!cam_.is_zero_column(key)) {
          // and insert it in annotations_in_boundary with multiplicative factor "sign".
          annotations_in_boundary.emplace_back(&cam_.get_column(key), sign);
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
    std::pair<typename std::map<Simplex_key, Arith_element>::iterator, bool> result_insert_a_ds;

    for (auto ann_it = annotations_in_boundary.begin(); ann_it != annotations_in_boundary.end(); /**/) {
      Column const* col = ann_it->first;
      int mult = ann_it->second;
      ++ann_it;
      while (ann_it != annotations_in_boundary.end() && ann_it->first == col) {
        mult += ann_it->second;
        ++ann_it;
      }
      // The following test is just a heuristic, it is not required, and it is fine that is misses p == 0.
      if (mult != coeff_field_.additive_identity()) {  // For all columns in the boundary,
        for (const auto& entry_ref : col->col_) {              // insert every entry in map_a_ds with multiplicity
          Arith_element w_y = coeff_field_.times(entry_ref.coefficient_, mult);  // coefficient * multiplicity

          if (w_y != coeff_field_.additive_identity()) {  // if != 0
            result_insert_a_ds = map_a_ds.insert(std::pair<Simplex_key, Arith_element>(entry_ref.key_, w_y));
            if (!(result_insert_a_ds.second)) {  // if entry_ref.key_ already a Key in map_a_ds
              result_insert_a_ds.first->second = coeff_field_.plus_equal(result_insert_a_ds.first->second, w_y);
              if (result_insert_a_ds.first->second == coeff_field_.additive_identity()) {
                map_a_ds.erase(result_insert_a_ds.first);
              }
            }
          }
        }
      }
    }
  }

  /*
   * Update the cohomology groups under the insertion of a simplex.
   */
  void update_cohomology_groups(Simplex_handle sigma, int dim_sigma)
  {
    // Compute the annotation of the boundary of sigma:
    std::map<Simplex_key, Arith_element> map_a_ds;
    annotation_of_the_boundary(map_a_ds, sigma, dim_sigma);
    // Update the cohomology groups:
    if (map_a_ds.empty()) {  // sigma is a creator in all fields represented in coeff_field_
      if (dim_sigma < dim_max_) {
        create_cocycle(sigma, coeff_field_.multiplicative_identity(), coeff_field_.characteristic());
      }
    } else {  // sigma is a destructor in at least a field in coeff_field_
      // Convert map_a_ds to a vector
      A_ds_type a_ds;  // admits reverse iterators
      for (auto map_a_ds_ref : map_a_ds) {
        a_ds.push_back(std::pair<Simplex_key, Arith_element>(map_a_ds_ref.first, map_a_ds_ref.second));
      }

      Arith_element inv_x, charac;
      Arith_element prod = coeff_field_.characteristic();  // Product of characteristic of the fields
      for (auto a_ds_rit = a_ds.rbegin(); (a_ds_rit != a_ds.rend()) && (prod != coeff_field_.multiplicative_identity());
           ++a_ds_rit) {
        std::tie(inv_x, charac) = coeff_field_.inverse(a_ds_rit->second, prod);

        if (inv_x != coeff_field_.additive_identity()) {
          destroy_cocycle(sigma, a_ds, a_ds_rit->first, inv_x, charac);
          prod /= charac;
        }
      }
      if (prod != coeff_field_.multiplicative_identity() && dim_sigma < dim_max_) {
        create_cocycle(sigma, coeff_field_.multiplicative_identity(prod), prod);
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
  void create_cocycle(Simplex_handle sigma, Arith_element x, Arith_element charac)
  {
    Simplex_key key = cpx_->key(sigma);
    cam_.insert_column({std::make_pair(key, x)});
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
  void destroy_cocycle(Simplex_handle sigma,
                       A_ds_type const& a_ds,
                       Simplex_key death_key,
                       Arith_element inv_x,
                       Arith_element charac)
  {
    // Create a finite persistent interval for which the interval exists
    if (interval_length_policy(cpx_->simplex(death_key), sigma)) {
      persistent_pairs_.emplace_back(cpx_->simplex(death_key)  // creator
                                     ,
                                     sigma  // destructor
                                     ,
                                     charac);  // fields
    }

    const auto& death_key_row = cam_.get_row(death_key);  // Find the beginning of the row.
    auto& death_cocycle_coeff = transverse_idx_[death_key];

    auto row_entry_it = death_key_row.begin();

    while (row_entry_it != death_key_row.end()) {  // Traverse all entries in the row at index death_key.
      const auto& entry = *row_entry_it;
      Arith_element w = coeff_field_.times_minus(inv_x, entry.get_element());

      if (w != coeff_field_.additive_identity()) {
        ++row_entry_it;  // has to be done before column addition to avoid invalidating the iterator

        cam_.multiply_source_and_add_to(w, a_ds, entry.get_column_index());
      } else {
        ++row_entry_it;
      }  // If w == 0, pass.
    }

    // Because it is a killer simplex, set the data of sigma to null_key().
    if (charac == coeff_field_.characteristic()) {
      cpx_->assign_key(sigma, cpx_->null_key());
    }
    if (death_cocycle_coeff == charac) {
      cam_.erase_empty_row(death_key);
      transverse_idx_.erase(death_key);
    } else {
      death_cocycle_coeff /= charac;
    }
  }

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

 public:
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
  std::vector<int> betti_numbers() const
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
  int betti_number(int dimension) const
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
  std::vector<int> persistent_betti_numbers(Filtration_value from, Filtration_value to) const
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
  int persistent_betti_number(int dimension, Filtration_value from, Filtration_value to) const
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
  std::vector<std::pair<Filtration_value, Filtration_value>> intervals_in_dimension(int dimension)
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
  FilteredComplex* cpx_;
  int dim_max_;
  CoefficientField coeff_field_;
  size_t num_simplices_;

  /*  Disjoint sets data structure to link the model of FilteredComplex
   * with the compressed annotation matrix.
   * ds_rank_ is a property map Simplex_key -> int, ds_parent_ is a property map
   * Simplex_key -> simplex_key_t */
  Cell_union_find<Simplex_key> vertex_sets_;

  /* The compressed annotation matrix fields.*/
  Cam cam_;
  /*  Dictionary establishing the correspondence between the Simplex_key of
   * the root vertex in the union-find ds and the Simplex_key of the vertex which
   * created the connected component as a 0-dimension homology feature.*/
  std::unordered_map<Simplex_key, Simplex_key> zero_cocycles_;
  /*  Key -> row. */
  std::map<Simplex_key, Arith_element> transverse_idx_;
  /* Persistent intervals. */
  std::vector<Persistent_interval> persistent_pairs_;
  length_interval interval_length_policy;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_H_
