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

# from gudhi.invariants import _persistence_objects_ext as t


class IntervalObject:
    def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]]):
        # pair(finite, infinite) of dim x bar number x (birth, death)
        # same numbers of dim
        # no death for infinite bars
        self._intervals = intervals

    def get_intervals(
        self, min_persistence: float = 0, dimension: Optional[int] = None
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [
                np.concatenate(
                    (
                        f[f[:, 1] - f[:, 0] >= min_persistence],
                        np.column_stack([i, np.full(len(i), np.inf)]),
                    ),
                    axis=0,
                )
                for f, i in zip(self._intervals[0], self._intervals[1])
            ]
        f = self._intervals[0][dimension]
        i = self._intervals[1][dimension]
        return np.concatenate(
            (
                f[f[:, 1] - f[:, 0] >= min_persistence],
                np.column_stack([i, np.full(len(i), np.inf)]),
            ),
            axis=0,
        )

    def get_finite_intervals(
        self, min_persistence: float = 0, dimension: Optional[int] = None
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [f[f[:, 1] - f[:, 0] >= min_persistence] for f in self._intervals[0]]
        f = self._intervals[0][dimension]
        return f[f[:, 1] - f[:, 0] >= min_persistence]

    def get_infinite_intervals(
        self, dimension: Optional[int] = None, no_inf_birth: bool = False
    ) -> Union[list[NDArray[np.double]], NDArray[np.double]]:
        if dimension is None:
            return [i[i != np.inf] if no_inf_birth else i for i in self._intervals[1]]
        i = self._intervals[1][dimension]
        return i[i != np.inf] if no_inf_birth else np.array(i)


# class BettiObject2:
#     def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]]):
#         # pair(finite, infinite) of dim x bar number x (birth, death)
#         # same numbers of dim
#         # no death for infinite bars
#         self._intervals = intervals

#     def get_betti_numbers(
#         self, min: float = np.inf, max: float = np.inf, dimension: Optional[int] = None
#     ) -> Union[NDArray[int], int]:
#         betti_numbers = t._get_betti_numbers(self._intervals[0], self._intervals[1], min, max)
#         if dimension is None:
#             return betti_numbers
#         return 0 if dimension >= betti_numbers.shape[0] else betti_numbers[dimension]


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


class PersistenceObject(IntervalObject, BettiObject):
    def __init__(self, intervals: tuple[list[NDArray[np.double]], list[NDArray[np.double]]]):
        # pair(finite, infinite) of dim x bar number x (birth, death)
        # same numbers of dim
        # no death for infinite bars
        self._intervals = copy.deepcopy(intervals)
        IntervalObject.__init__(self, self._intervals)
        BettiObject.__init__(self, self._intervals)
