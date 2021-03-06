from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2017 Inria

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2017 Inria"
__license__ = "GPL v3"

cdef extern from "Reader_utils_interface.h" namespace "Gudhi":
    vector[vector[double]] read_matrix_from_csv_file(string off_file, char separator)
    map[int, vector[pair[double, double]]] read_pers_intervals_grouped_by_dimension(string filename)
    vector[pair[double, double]] read_pers_intervals_in_dimension(string filename, int only_this_dim)

def read_lower_triangular_matrix_from_csv_file(csv_file='', separator=';'):
    """Read lower triangular matrix from a CSV style file.

    :param csv_file: A CSV file style name.
    :type csv_file: string
    :param separator: The value separator in the CSV file. Default value is ';'
    :type separator: char

    :returns:  The lower triangular matrix.
    :rtype: vector[vector[double]]
    """
    if csv_file is not '':
        if os.path.isfile(csv_file):
            return read_matrix_from_csv_file(str.encode(csv_file), ord(separator[0]))
    print("file " + csv_file + " not set or not found.")
    return []

def read_persistence_intervals_grouped_by_dimension(persistence_file=''):
    """Reads a file containing persistence intervals.
    Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
    The return value is an `map[dim, vector[pair[birth, death]]]`
    where `dim` is an `int`, `birth` a `double`, and `death` a `double`.
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string

    :returns:  The persistence pairs grouped by dimension.
    :rtype: map[int, vector[pair[double, double]]]
    """
    if persistence_file is not '':
        if os.path.isfile(persistence_file):
            return read_pers_intervals_grouped_by_dimension(str.encode(persistence_file))
    print("file " + persistence_file + " not set or not found.")
    return []

def read_persistence_intervals_in_dimension(persistence_file='', only_this_dim=-1):
    """Reads a file containing persistence intervals.
    Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
    If `only_this_dim` = -1, dimension is ignored and all lines are returned.
    If `only_this_dim` is >= 0, only the lines where dimension = `only_this_dim`
    (or where dimension is not specified) are returned.
    The return value is an `vector[pair[birth, death]]`
    where `birth` a `double`, and `death` a `double`.
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string

    :returns:  The persistence pairs grouped by dimension.
    :rtype: map[int, vector[pair[double, double]]]
    """
    if persistence_file is not '':
        if os.path.isfile(persistence_file):
            return read_pers_intervals_in_dimension(str.encode(persistence_file), only_this_dim)
    print("file " + persistence_file + " not set or not found.")
    return []
