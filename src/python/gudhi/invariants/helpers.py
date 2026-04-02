# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from typing import Literal, Protocol, runtime_checkable
from collections.abc import Iterable
from numpy.typing import NDArray

import numpy as np

from gudhi.invariants import PersistenceObject
from gudhi.invariants import _pers_cohomology_ext as t
from gudhi import SimplexTree, CubicalComplex, PeriodicCubicalComplex
from gudhi._ripser import _compute_ripser_arguments, _compute_persistence_with_ripser


@runtime_checkable  # enables isinstance()
class ComplexForCohomology(Protocol):
    def num_cells(self) -> int: ...

    def filtration(self, sh: np.uint32) -> np.double: ...

    def dimension(self, sh: int | None = None) -> int: ...

    def null_handle(self) -> np.uint32: ...

    def endpoints(self, sh: np.uint32) -> tuple[np.uint32, np.uint32]: ...

    def filtration_cell_range(self) -> Iterable[np.uint32]: ...

    def boundary_cell_range(self, sh: np.uint32) -> Iterable[np.uint32]: ...


def _compute_persistence_annotation_cohomology(
    cpx: SimplexTree | CubicalComplex | PeriodicCubicalComplex | ComplexForCohomology | None,
    coef: int,
    max_dim: int | None,
    persistence_dim_max: bool,
):
    if cpx is None:
        raise ValueError("'cpx' argument has to be given for 'annotation_cohomology' backend")
    elif isinstance(cpx, SimplexTree):
        if max_dim is not None:
            cpx = t._St_dim_overlay(cpx, max_dim)
            pers = t._Simplex_tree_max_dim_persistence_interface(cpx, persistence_dim_max)
        else:
            pers = t._Simplex_tree_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, CubicalComplex):
        if max_dim is not None:
            cpx = t._Cc_dim_overlay(cpx, max_dim)
            pers = t._Cubical_complex_max_dim_persistence_interface(cpx, persistence_dim_max)
        else:
            pers = t._Cubical_complex_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, PeriodicCubicalComplex):
        if max_dim is not None:
            cpx = t._Pcc_complex_dim_overlay(cpx, max_dim)
            pers = t._Periodic_cubical_complex_max_dim_persistence_interface(
                cpx, persistence_dim_max
            )
        else:
            pers = t._Periodic_cubical_complex_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, ComplexForCohomology):
        cpx = t._complex_coho_overlay(cpx)
        if max_dim is not None:
            # careful to not release cpx as it is passed by pointer
            d_cpx = t._Complex_dim_overlay(cpx, max_dim)
            pers = t._Complex_max_dim_persistence_interface(d_cpx, persistence_dim_max)
        else:
            pers = t._Complex_persistence_interface(cpx, persistence_dim_max)
    else:
        raise ValueError("Complex type not supported yet.")

    pers._compute_persistence(coef, 0)
    return pers._get_intervals()


def _compute_persistence_ripser_(
    data,
    rips_max_dimension: int,
    homology_coeff_field: int,
    input_type: Literal[
        "full distance matrix", "lower distance matrix", "distance coo_matrix"
    ],
    threshold: float,
):
    dgm = _compute_persistence_with_ripser(
        input=data,
        max_dimension=rips_max_dimension,
        homology_coeff_field=homology_coeff_field,
        input_type=input_type,
        threshold=threshold,
    )

    finites = [np.asarray([p for p in d_dgm if p[1] != float("inf")]) for d_dgm in dgm]
    infinites = [np.asarray([p[0] for p in d_dgm if p[1] == float("inf")]) for d_dgm in dgm]
    return (finites, infinites)


def _progress(verbose: bool, msg: str):
    if verbose:
        print(msg)


def compute_persistence(
    cpx: (
        SimplexTree | CubicalComplex | PeriodicCubicalComplex | ComplexForCohomology | None
    ) = None,
    data: list[list[float]] | list[np.ndarray] | None = None,
    data_input_type: (
        Literal[
            "point cloud",
            "full distance matrix",
            "lower distance matrix",
            "distance coo_matrix",
        ]
        | None
    ) = None,
    coef: int = 2,
    max_dim: int | None = None,
    persistence_dim_max: bool = False,
    threshold: float = float("inf"),
    num_collapses: int | None = None,
    force_backend: Literal["annotation_cohomology", "ripser"] | None = None,
    verbose: bool = False,
) -> PersistenceObject:
    if force_backend is not None:
        backend = force_backend
    else:
        if (cpx is None and data is None) or (cpx is not None and data is not None):
            raise ValueError(
                "Exactly one argument has to be specified among `cpx` and `data`."
            )
        if cpx is None:
            if data_input_type is None or max_dim is None:
                raise ValueError(
                    "If `data` is specified instead of `cpx`, the arguments `data_input_type`"
                    + " and `max_dim` have to be given."
                )
            backend = "ripser"
        else:
            backend = "annotation_cohomology"

    _progress(verbose, "Chosen backend: `" + backend + "`.")

    if backend == "ripser":
        if data is None or data_input_type is None or max_dim is None:
            raise ValueError(
                "Arguments `data`, `data_input_type` and `max_dim` have to be given for 'ripser' backend."
            )
        if num_collapses is None:
            num_collapses = "auto"
        input_data, input, rips_max_dimension = _compute_ripser_arguments(
            input=data,
            max_dim=max_dim,
            homology_coeff_field=coef,
            num_collapses=num_collapses,
            input_type=data_input_type,
            threshold=threshold,
        )
        if input_data is True:
            _progress(verbose, "Fallback to backend `annotation_cohomology` because of size.")
            backend = "annotation_cohomology"
            cpx = input
            persistence_dim_max = rips_max_dimension >= cpx.dimension()

    _progress(verbose, "Computing persistence pairs.")

    match backend:
        case "annotation_cohomology":
            pers = _compute_persistence_annotation_cohomology(
                cpx, coef, max_dim, persistence_dim_max
            )
        case "ripser":
            pers = _compute_persistence_ripser_(
                input, rips_max_dimension, coef, input_data, threshold
            )
        case _:
            raise ValueError(
                "argument `force_backend` does not contain a valid literal: "
                + force_backend
                + ". Possibilities are: `annotation_cohomology`, `ripser`"
            )

    _progress(verbose, "Done.")

    return PersistenceObject(pers)
