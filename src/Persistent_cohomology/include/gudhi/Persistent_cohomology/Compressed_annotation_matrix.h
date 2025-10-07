/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria, Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PCH_CAM_H_
#define PCH_CAM_H_

#include <stdexcept>  // for std::logic_error
#include <utility>    // std::move, std::pair
#include <map>
#include <vector>
#include <initializer_list>

#include <boost/intrusive/set.hpp>
#include <boost/intrusive/list.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistent_cohomology/Persistent_cohomology_column.h>
#include <gudhi/Simple_object_pool.h>
#include <gudhi/Persistent_cohomology/Simple_union_find.h>

namespace Gudhi {

namespace persistent_cohomology {

template <typename Cell_key, class CoefficientField>
class Compressed_annotation_matrix
{
 public:
  using Characteristic = typename CoefficientField::Element;
  using Index = Cell_key;
  using ID_index = Cell_key;
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
  Compressed_annotation_matrix(unsigned int numberOfColumns, const CoefficientField& coeff_field)
      : coeff_field_(coeff_field),
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
  void set_characteristic(Characteristic characteristic) { coeff_field_.set_characteristic(characteristic); }

  void set_characteristic(int charac_min, int charac_max) { coeff_field_.set_characteristic(charac_min, charac_max); }

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
  // template <class Container = std::initializer_list<std::pair<Index, Characteristic>>>
  void insert_column(Index key, Characteristic x)
  {
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
  const Column& get_column(Index columnIndex)
  {
    auto* col = ds_repr_[cell_sets_.find_set(columnIndex)];
    if (col == nullptr) return empty_column_;
    return *col;
  }

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
    GUDHI_CHECK(transverse_idx_.find(rowIndex)->second.size() == 0,
                std::logic_error("Internal problem: row should be empty."));
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
  void multiply_source_and_add_to(Characteristic coefficient,
                                  const Entry_range& sourceColumn,
                                  Column const* constTargetColumn)
  {
    Column* targetColumn = const_cast<Column*>(constTargetColumn);

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

        Cell_key key_tmp = cell_sets_.find_set(targetColumn->class_key_);
        ds_repr_[key_tmp] = &(*(result_insert_cam.first));
        result_insert_cam.first->class_key_ = key_tmp;
        // intrusive containers don't own their elements, we have to release them manually
        targetColumn->col_.clear_and_dispose([&](Entry* p) { entry_pool_.destroy(p); });
        column_pool_.destroy(targetColumn);  // delete targetColumn;
      }
    }
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
  std::map<Cell_key, Row> transverse_idx_;
  Cam cam_;
  Simple_union_find<Cell_key> cell_sets_;

  Simple_object_pool<Column> column_pool_;
  Simple_object_pool<Entry> entry_pool_;

  inline static const Column empty_column_; /**< Representative for empty columns. */

  /*
   * Assign:    target <- target + w * other.
   */
  template <class Entry_range>
  void plus_equal_column(Column& target,
                         Entry_range const& other  // value_type is pair<Cell_key,Arith_element>
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
                                                         coeff_field_.get_additive_identity(),
                                                         &target));

          entry_tmp->coefficient_ = coeff_field_.multiply_and_add(other_it->second, w, entry_tmp->coefficient_);

          target.col_.insert(target_it, *entry_tmp);

          ++other_it;
        } else {  // it1->key == it2->key
          // target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
          target_it->coefficient_ = coeff_field_.multiply_and_add(other_it->second, w, target_it->coefficient_);
          if (target_it->coefficient_ == coeff_field_.get_additive_identity()) {
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
      Entry* entry_tmp = entry_pool_.construct(Entry(other_it->first, coeff_field_.get_additive_identity(), &target));
      entry_tmp->coefficient_ = coeff_field_.multiply_and_add(other_it->second, w, entry_tmp->coefficient_);
      target.col_.insert(target.col_.end(), *entry_tmp);

      ++other_it;
    }
  }
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PCH_CAM_H_
