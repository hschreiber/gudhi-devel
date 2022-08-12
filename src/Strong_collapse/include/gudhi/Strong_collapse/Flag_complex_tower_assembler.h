/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STRONG_COLLAPSE_FLAG_COMPLEX_TOWER_ASSEMBLER_H_
#define STRONG_COLLAPSE_FLAG_COMPLEX_TOWER_ASSEMBLER_H_

#include <gudhi/Strong_collapse/Flag_complex_sparse_matrix.h>

#include <set>
#include <fstream>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>

namespace Gudhi {

namespace strong_collapse {

// assumptions : (1) K1 and K2 have the same vertex set
//               (2) The set of simplices of K1 is a subset of set of simplices of K2
// K1  ->  K2    [Original Simplicial Complexes]
// |       |
// |       |
// K1c ->  K2c   [Strongly Collapsed Flag Complexes]

template<class SimplicialComplex>
class Flag_complex_tower_assembler {
 private:
  using Filtration_value = typename SimplicialComplex::Filtration_value;
  using Reduction_map = typename Flag_complex_sparse_matrix<SimplicialComplex>::Reduction_map;
  using Edge = typename Flag_complex_sparse_matrix<SimplicialComplex>::Edge;
  using Sparse_row_matrix = typename Flag_complex_sparse_matrix<SimplicialComplex>::Sparse_row_matrix;
  using Sparse_row_iterator = typename Flag_complex_sparse_matrix<SimplicialComplex>::Sparse_row_iterator;

  Reduction_map renamed_vertices_;
  std::size_t current_rename_counter_;
  Flag_complex_sparse_matrix<SimplicialComplex> flag_filtration_;

 public:
  Flag_complex_tower_assembler(const std::size_t num_vertices)
      : flag_filtration_((num_vertices * log2(num_vertices)) +1) {
    for (std::size_t i = 0; i <= num_vertices; ++i) {
      renamed_vertices_[i] = i;
    }
    current_rename_counter_ = num_vertices + 1;
  }

  // mat_1 and mat_2 are simplex_trees of K1c and K2c (the
  // collapsed ones), redmap_2 is the map of K2 -> K2c
  void build_tower_for_two_complexes(Flag_complex_sparse_matrix<SimplicialComplex>& mat_1,
                                     const Flag_complex_sparse_matrix<SimplicialComplex>& mat_2,
                                     const Filtration_value filtration_value, const std::string& outFile = "")
  {
    std::ofstream myfile(outFile, std::ios::app);
    if (myfile.is_open() || outFile.empty()) {
      Reduction_map reduction_map_2 = mat_2.reduction_map();
      for (auto& v : mat_1.vertex_set()) {
        // If v collapsed to something?
        auto collapsed_to = reduction_map_2.find(v);
        // Collapse happened, because there is a vertex in the map
        if (collapsed_to != reduction_map_2.end()) {
          // Collapsed to an existing vertex in mat_1.
          if (mat_1.membership(collapsed_to->second)) {
            flag_filtration_.active_strong_expansion(renamed_vertices_.at(v), renamed_vertices_.at(collapsed_to->second),
                                                     filtration_value);
            renamed_vertices_.at(v) = current_rename_counter_;
            current_rename_counter_++;
          } else {
            if (!outFile.empty()) {
              myfile << filtration_value << " i " << renamed_vertices_.at(collapsed_to->second) << std::endl;
              myfile << filtration_value << " c " << renamed_vertices_.at(v) << " "
                     << renamed_vertices_.at(collapsed_to->second) << std::endl;
            }
            flag_filtration_.active_strong_expansion(renamed_vertices_.at(v), renamed_vertices_.at(collapsed_to->second),
                                                     filtration_value);
            renamed_vertices_.at(v) = current_rename_counter_;
            current_rename_counter_++;
          }
          // If the vertex "collapsed_to->second" is not a member of mat_1, the contraction function will simply add and
          // then collapse
          mat_1.contraction(v, collapsed_to->second);
        }
      }

      // The core K1c (mat_1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of
      // K2c, the remaining simplices need to be included
      // Writing the inclusion of all remaining simplices...
      for (const Edge& e : mat_2.all_edges()) {
        auto u = std::get<0>(e);
        auto v = std::get<1>(e);
        if (!mat_1.membership(u)) {
          flag_filtration_.insert_vertex(renamed_vertices_.at(u), filtration_value);
          if (!outFile.empty()) {
            myfile << filtration_value << " i";
            myfile << " " << renamed_vertices_.at(u);
            myfile << std::endl;
          }
          mat_1.insert_vertex(u, 1);
        }

        if (!mat_1.membership(v)) {
          flag_filtration_.insert_vertex(renamed_vertices_.at(v), filtration_value);

          if (!outFile.empty()) {
            myfile << filtration_value << " i";
            myfile << " " << renamed_vertices_.at(v);
            myfile << std::endl;
          }
          mat_1.insert_vertex(v, 1);
        }
        if (!mat_1.membership(e)) {
          flag_filtration_.insert_new_edges(renamed_vertices_.at(u), renamed_vertices_.at(v), filtration_value);

          if (!outFile.empty()) {
            myfile << filtration_value << " i";
            myfile << " " << renamed_vertices_.at(u) << ", " << renamed_vertices_.at(v);
            myfile << std::endl;
          }
          mat_1.insert_new_edges(u, v, 1);
        }
      }
      if (!outFile.empty()) {
        myfile << "# Tower updated for the additional subcomplex.\n";
        myfile.close();
      }
    } else {
      std::cerr << "Unable to open file " << outFile;
      exit(-1);
    }
  }

  typename Flag_complex_sparse_matrix<SimplicialComplex>::Distance_matrix distance_matrix() {
    std::size_t non_zero_rw = flag_filtration_.num_vertices();
    Filtration_value inf = std::numeric_limits<Filtration_value>::max();
    Sparse_row_matrix mat = flag_filtration_.uncollapsed_matrix();

    typename Flag_complex_sparse_matrix<SimplicialComplex>::Distance_matrix distance_mat;
    for (std::size_t indx = 0; indx < non_zero_rw; indx++) {
      std::vector<Filtration_value> distances;
      Sparse_row_iterator it(mat, indx);
      // Iterate over the non-zero columns
      for (std::size_t j = 0; j <= indx; j++) {
        if (it.index() == j && j != indx) {
          // inner index, here it is equal to it.columns()
          distances.push_back(it.value());
          ++it;
        } else if (j == indx) {
          distances.push_back(0);
        } else {
          distances.push_back(inf);
        }
      }
      distance_mat.push_back(distances);
    }
    return distance_mat;
  }

  void print_sparse_matrix() {
    flag_filtration_.print_sparse_skeleton();
  }
};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // STRONG_COLLAPSE_FLAG_COMPLEX_TOWER_ASSEMBLER_H_
