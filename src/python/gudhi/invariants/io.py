# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from typing import Union, Iterable, Literal, Optional, Any
import numpy as np

from gudhi.invariants import PersistenceObject


def write_diagram_in_file(
    diagram: PersistenceObject,
    persistence_file: str,
    min_persistence: float = 0,
    delimiter: str = " ",
):
    values = diagram.get_intervals(min_persistence=min_persistence)

    total_rows = sum(arr.shape[0] for arr in values)
    out = np.empty((total_rows, 3), dtype=np.double)

    row = 0
    for d, arr in enumerate(values):
        n = arr.shape[0]
        out[row : row + n, 0] = d
        out[row : row + n, 1:] = arr[:, :]
        row += n

    np.savetxt(persistence_file, out, fmt=["%d", "%.18f", "%.18f"], delimiter=delimiter)


def read_diagram_from_file(
    persistence_file: str,
    delimiter: str = " ",
) -> PersistenceObject:
    raw = np.loadtxt(persistence_file, dtype=np.double, delimiter=delimiter)
    # to order by dimension and have inf deaths last
    raw = raw[np.lexsort((raw[:, 1], raw[:, 2], raw[:, 0]))]
    max_dim = int(raw[-1][0])
    raw_finite = raw[raw[:, 2] != np.inf]
    raw_infinite = raw[raw[:, 2] == np.inf]
    raw_infinite = raw_infinite[:, :2]

    split_points = np.where(np.diff(raw_finite[:, 0]))[0] + 1
    dimensions = raw_finite[np.r_[0, split_points], 0]
    finite = np.split(raw_finite[:, 1:], split_points)
    finite = dict(zip(dimensions, finite))
    finite = [finite.get(np.double(d), np.empty((0, 2))) for d in range(max_dim + 1)]

    split_points = np.where(np.diff(raw_infinite[:, 0]))[0] + 1
    dimensions = raw_infinite[np.r_[0, split_points], 0]
    infinite = np.split(raw_infinite[:, 1], split_points)
    infinite = dict(zip(dimensions, infinite))
    infinite = [infinite.get(np.double(d), np.empty((0))) for d in range(max_dim + 1)]

    return PersistenceObject(intervals=(finite, infinite))
