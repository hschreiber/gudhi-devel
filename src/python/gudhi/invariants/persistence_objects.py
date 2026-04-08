# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from typing import Optional, Union
from numpy.typing import NDArray

import copy
import numpy as np


class IntervalObject:
    def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]]):
        # pair(finite, infinite) of dim x bar number x (birth, death)
        # same numbers of dim
        # no death for infinite bars
        self._intervals = intervals

    def _get_union_io(self, f: NDArray[np.double], i: NDArray[np.double], min_persistence: float):
        return np.concatenate(
            (
                f[f[:, 1] - f[:, 0] >= min_persistence],
                np.column_stack([i, np.full(len(i), np.inf)]),
            ),
            axis=0,
        )

    def get_intervals(
        self, min_persistence: float = 0, dimension: Optional[int] = None
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [
                self._get_union_io(f, i, min_persistence)
                for f, i in zip(self._intervals[0], self._intervals[1])
            ]
        if dimension >= len(self._intervals[0]):
            return np.empty(shape=[0, 2])
        return self._get_union_io(
            self._intervals[0][dimension], self._intervals[1][dimension], min_persistence
        )

    def get_finite_intervals(
        self, min_persistence: float = 0, dimension: Optional[int] = None
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [f[f[:, 1] - f[:, 0] >= min_persistence] for f in self._intervals[0]]
        if dimension >= len(self._intervals[0]):
            return np.empty(shape=[0, 2])
        f = self._intervals[0][dimension]
        return f[f[:, 1] - f[:, 0] >= min_persistence]

    def get_infinite_intervals(
        self, dimension: Optional[int] = None, no_inf_birth: bool = False
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [i[i != np.inf] if no_inf_birth else i for i in self._intervals[1]]
        if dimension >= len(self._intervals[1]):
            return np.empty(shape=[0, 1])
        i = self._intervals[1][dimension]
        return i[i != np.inf] if no_inf_birth else np.array(i)


class SimplexIntervalObject:
    def __init__(self, pers):
        self._pers = pers

    def _get_union_sio(self, f: NDArray, i: NDArray):
        return np.concatenate(
            (
                f,
                np.pad(i[:, np.newaxis, :], ((0, 0), (0, 1), (0, 0)), constant_values=-1),
            ),
            axis=0,
        )

    def build_simplicial_intervals_as_vertices(
        self,
        min_persistence: float = 0,
        max_dimension: Optional[int] = None,
        dimension: Optional[int] = None,
    ) -> Union[list[NDArray], NDArray]:
        if not self._pers:
            raise RuntimeError(
                "`build_simplicial_intervals_as_vertices` cannot be executed as it is " +
                "not compatible with the chosen persistence computation method."
            )

        if (dimension is not None) and (max_dimension is None or max_dimension > dimension):
            max_dimension = dimension
        if dimension is not None and max_dimension is not None and dimension > max_dimension:
            return np.empty(shape=[0, 2, dimension + 1])

        intervals = self._pers._get_simplicial_intervals_as_vertices(
            min_persistence, max_dimension if max_dimension is not None else -1, False, False
        )
        if dimension is None:
            return [self._get_union_sio(f, i) for f, i in zip(intervals[0], intervals[1])]
        if dimension >= len(intervals[0]):
            return np.empty(shape=[0, 2, dimension + 1])
        return self._get_union_sio(intervals[0][dimension], intervals[1][dimension])

    def build_finite_simplicial_intervals_as_vertices(
        self,
        min_persistence: float = 0,
        max_dimension: Optional[int] = None,
        dimension: Optional[int] = None,
    ) -> Union[list[NDArray], NDArray]:
        if not self._pers:
            raise RuntimeError(
                "`build_finite_simplicial_intervals_as_vertices` cannot be executed as it is " +
                "not compatible with the chosen persistence computation method."
            )

        if (dimension is not None) and (max_dimension is None or max_dimension > dimension):
            max_dimension = dimension
        if dimension is not None and max_dimension is not None and dimension > max_dimension:
            return np.empty(shape=[0, 2, dimension + 1])

        intervals = self._pers._get_simplicial_intervals_as_vertices(
            min_persistence, max_dimension if max_dimension is not None else -1, True, False
        )
        if dimension is None:
            return intervals[0]
        if dimension >= len(intervals[0]):
            return np.empty(shape=[0, 2, dimension + 1])
        return intervals[0][dimension]

    def build_infinite_simplicial_intervals_as_vertices(
        self, max_dimension: Optional[int] = None, dimension: Optional[int] = None
    ) -> Union[list[NDArray], NDArray]:
        if not self._pers:
            raise RuntimeError(
                "`build_infinite_simplicial_intervals_as_vertices` cannot be executed as it is " +
                "not compatible with the chosen persistence computation method."
            )

        if (dimension is not None) and (max_dimension is None or max_dimension > dimension):
            max_dimension = dimension
        if dimension is not None and max_dimension is not None and dimension > max_dimension:
            return np.empty(shape=[0, dimension + 1])

        intervals = self._pers._get_simplicial_intervals_as_vertices(
            min_persistence, max_dimension if max_dimension is not None else -1, False, True
        )
        if dimension is None:
            return intervals[1]
        if dimension >= len(intervals[1]):
            return np.empty(shape=[0, dimension + 1])
        return intervals[1][dimension]


class BettiObject:
    def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]]):
        # pair(finite, infinite) of dim x bar number x (birth, death)
        # same numbers of dim
        # no death for infinite bars
        self._intervals = intervals

    def get_betti_numbers(
        self, min: float = np.inf, max: float = np.inf, dimension: Optional[int] = None
    ) -> Union[NDArray[int], int]:
        if dimension is None:
            if min == np.inf:
                if max != np.inf:
                    mask_fd = [(f[:, 1] > max) for f in self._intervals[0]]
                mask_id = [np.ones_like(f, dtype=bool) for f in self._intervals[1]]
            else:
                if max != np.inf:
                    mask_fd = [(f[:, 0] <= min & f[:, 1] > max) for f in self._intervals[0]]
                mask_id = [(f <= min) for f in self._intervals[1]]
            i_betti = [np.count_nonzero(mask_i) for mask_i in mask_id]
            if max != np.inf:
                f_betti = [np.count_nonzero(mask_f) for mask_f in mask_fd]
                return np.sum([f_betti, i_betti], axis=0)
            return np.array(i_betti)
        if dimension >= len(self._intervals[0]):
            return 0
        f = self._intervals[0][dimension]
        i = self._intervals[1][dimension]
        if min == np.inf:
            if max != np.inf:
                mask_fd = f[:, 1] > max
            mask_id = np.ones_like(i, dtype=bool)
        else:
            if max != np.inf:
                mask_fd = f[:, 0] <= min & f[:, 1] > max
            mask_id = i <= min
        i_betti = np.count_nonzero(mask_id)
        if max != np.inf:
            f_betti = np.count_nonzero(mask_fd)
            return f_betti + i_betti
        return i_betti


class PersistenceObject(IntervalObject, BettiObject, SimplexIntervalObject):
    def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]], pers = None):
        # pair(finite, infinite) of dim x bar number x (birth, death)
        # same numbers of dim
        # no death for infinite bars
        self._intervals = copy.deepcopy(intervals)
        IntervalObject.__init__(self, self._intervals)
        BettiObject.__init__(self, self._intervals)
        # if pers is None, calling the methods of SimplexIntervalObject will raise a RuntimeError
        if pers is not None:
            self._pers_pair = pers
            SimplexIntervalObject.__init__(self, self._pers_pair[0])
