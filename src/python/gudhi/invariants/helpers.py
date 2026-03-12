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
    coef: int = 2,
    persistence_dim_max: bool = False,
):
    if cpx is None:
        raise ValueError("'cpx' argument has to be given for 'annotation_cohomology' backend")
    elif isinstance(cpx, SimplexTree):
        pers = t._Simplex_tree_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, CubicalComplex):
        pers = t._Cubical_complex_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, PeriodicCubicalComplex):
        pers = t._Periodic_cubical_complex_persistence_interface(cpx, persistence_dim_max)
    elif isinstance(cpx, ComplexForCohomology):
        cpx = t._complex_coho_overlay(cpx)
        pers = t._Complex_persistence_interface(cpx, persistence_dim_max)
    else:
        raise ValueError("Complex type not supported yet.")

    pers._compute_persistence(coef, 0)
    return pers._get_intervals()


def compute_persistence(
    cpx: (
        SimplexTree | CubicalComplex | PeriodicCubicalComplex | ComplexForCohomology | None
    ) = None,
    coef: int = 2,
    persistence_dim_max: bool = False,
    force_backend: Literal["annotation_cohomology"] | None = None,
) -> PersistenceObject:
    if force_backend is not None:
        match force_backend:
            case "annotation_cohomology":
                pers = _compute_persistence_annotation_cohomology(
                    cpx, coef, persistence_dim_max
                )
            case _:
                raise ValueError(
                    "argument `force_backend` does not contain a valid literal: "
                    + force_backend
                    + ". Possibilities are: `annotation_cohomology`, ..."
                )
    else:
        if cpx is None:
            return
        else:
            pers = _compute_persistence_annotation_cohomology(cpx, coef, persistence_dim_max)

    return PersistenceObject(pers)


def compute_persistence2(
    cpx: (
        SimplexTree | CubicalComplex | PeriodicCubicalComplex | ComplexForCohomology | None
    ) = None,
    coef: int = 2,
    persistence_dim_max: bool = False,
    force_backend: Literal["annotation_cohomology"] | None = None,
) -> PersistenceObject:
    if force_backend is not None:
        backend = force_backend
    else:
        if cpx is None:
            return
        else:
            backend = "annotation_cohomology"

    match backend:
        case "annotation_cohomology":
            pers = _compute_persistence_annotation_cohomology(cpx, coef, persistence_dim_max)
        case _:
            raise ValueError(
                "argument `force_backend` does not contain a valid literal: "
                + force_backend
                + ". Possibilities are: `annotation_cohomology`, ..."
            )

    return PersistenceObject(pers)
