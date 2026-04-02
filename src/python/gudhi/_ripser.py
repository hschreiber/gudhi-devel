# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2022 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


import math
import numpy as np
from typing import Literal
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform

from gudhi._ripser_ext import _lower, _full, _sparse, _lower_to_coo, _lower_cone_radius
from gudhi.flag_filtration.edge_collapse import reduce_graph
from gudhi import SimplexTree


def _compute_ripser_input(
    input,
    max_dimension: int,
    use_simplex_tree: bool,
    num_collapses: int,
    input_type: Literal[
        "point cloud", "full distance matrix", "lower distance matrix", "distance coo_matrix"
    ],
    threshold: float,
    n: int
):
    # TODO: give the user more control over the strategy
    # Should we use threshold in the sparse case?

    # Points -> distance matrix
    if input_type == "point cloud":
        if threshold < float("inf"):
            # Hope that the user gave a useful threshold
            tree = cKDTree(input)

            ## V1: Returns self-loops and every edge twice (symmetry)
            # inp = tree.sparse_distance_matrix(tree, max_distance=threshold, output_type="coo_matrix")
            # mask = inp.row < inp.col
            # inp = coo_matrix((inp.data[mask], (inp.row[mask], inp.col[mask])), shape=inp.shape)

            # V2: Gets the right edges, but forgets the distances
            pairs = tree.query_pairs(r=threshold, output_type="ndarray")
            data = np.ravel(np.linalg.norm(np.diff(input[pairs], axis=1), axis=-1))
            input = coo_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(n,) * 2)

            input_type = "distance coo_matrix"
        else:
            input = squareform(pdist(input))
            input_type = "full distance matrix"

    # Dense -> sparse
    if input_type in ("full distance matrix", "lower distance matrix"):
        # After this filtration value, all complexes are cones, nothing happens
        if input_type == "full distance matrix":
            input = np.asarray(input)
            cone_radius = input.max(-1).min()
        else:
            cone_radius = _lower_cone_radius(input)
        sparsify = (
            use_simplex_tree or num_collapses > 0 or threshold < cone_radius
        )  # heuristic
        threshold = min(threshold, cone_radius)
        if sparsify:
            # For 'full' we could use i, j = np.triu_indices_from(inp, k=1), etc
            i, j, f = _lower_to_coo(input, threshold)
            input = coo_matrix((f, (i, j)), shape=(n,) * 2)
            input_type = "distance coo_matrix"

    if num_collapses > 0:
        if input_type != "distance coo_matrix":
            raise ValueError(
                "Input type has to be 'distance coo_matrix' when 'num_collapses' is not 0"
            )
        input = reduce_graph(input, num_collapses)

    return input, input_type


def _compute_persistence_with_ripser(
    input,
    max_dimension: int,
    homology_coeff_field: int,
    input_type: Literal[
        "full distance matrix", "lower distance matrix", "distance coo_matrix"
    ],
    threshold: float,
):
    if input_type == "full distance matrix":
        ## Possibly transpose for performance?
        # if input.strides[0] > input.strides[1]: # or the reverse?
        #     input = input.T
        dgm = _full(
            input,
            max_dimension=max_dimension,
            max_edge_length=threshold,
            homology_coeff_field=homology_coeff_field,
        )
    elif input_type == "lower distance matrix":
        dgm = _lower(
            input,
            max_dimension=max_dimension,
            max_edge_length=threshold,
            homology_coeff_field=homology_coeff_field,
        )
    elif input_type == "distance coo_matrix":
        # Switch to coo_array (danger: row/col seem deprecated)?
        dgm = _sparse(
            input.row,
            input.col,
            input.data,
            input.shape[0],
            max_dimension=max_dimension,
            max_edge_length=threshold,
            homology_coeff_field=homology_coeff_field,
        )
    else:
        raise ValueError(
            "Only 'point cloud', 'lower distance matrix', 'full distance matrix' and 'distance coo_matrix' are valid input_type"
        )  # move to __init__?

    return dgm


def _build_simplex_tree_from_ripser_matrix(input, num_vertices: int, max_dimension: int):
    st = SimplexTree()
    # Use create_from_array in case of full matrix?
    # (not important since this fallback mostly matters in high dimension, where we use edge-collapse anyway)
    st.insert_batch(np.arange(num_vertices).reshape(1, -1), np.zeros(num_vertices))
    st.insert_edges_from_coo_matrix(input)
    st.expansion(max_dimension + 1)
    return st


def _compute_ripser_arguments(
    input,
    homology_dim_list: np.ndarray[int] | None = None,
    max_dim: int | None = None,
    homology_coeff_field: int = 11,
    num_collapses: int | Literal["auto"] = "auto",
    input_type: Literal[
        "point cloud", "full distance matrix", "lower distance matrix", "distance coo_matrix"
    ] = "point cloud",
    threshold: float = float("inf"),
):
    # TODO: give the user more control over the strategy
    # Should we use threshold in the sparse case?
    num_vertices = input.shape[0] if input_type == "distance coo_matrix" else len(input)
    if max_dim is None:
        # assumes homology_dim_list not None if max_dim is None
        # as this is private, i don't think that an assert is necessary here?
        max_dim = max(homology_dim_list)
    max_dimension = min(max_dim, max(0, num_vertices - 3))
    # Ripser needs to encode simplices and coefficients in 128 bits, which may not always fit
    # Instead of a 256 bit version which may not always suffice either, fall back to SimplexTree
    use_simplex_tree = math.comb(num_vertices, min(num_vertices // 2, max_dimension + 2)) >= (
        1 << (128 - (homology_coeff_field - 2).bit_length())
    )
    if num_collapses == "auto":
        num_collapses = 1 if max_dimension > (not use_simplex_tree) else 0
        # or num_collapses=max_dimension-1 maybe?
    elif max_dimension == 0:
        num_collapses = 0

    input, input_type = _compute_ripser_input(
        input=input,
        max_dimension=max_dimension,
        use_simplex_tree=use_simplex_tree,
        num_collapses=num_collapses,
        input_type=input_type,
        threshold=threshold,
        n=num_vertices,
    )

    if use_simplex_tree:
        input = _build_simplex_tree_from_ripser_matrix(
            input=input, num_vertices=num_vertices, max_dimension=max_dimension
        )
        input_type = use_simplex_tree

    return (input_type, input, max_dimension)

