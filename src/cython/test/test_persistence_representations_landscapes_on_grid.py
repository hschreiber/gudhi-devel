import gudhi
import sys
import math

"""
    This file is part of the Gudhi Library. The Gudhi library
    (Geometric Understanding in Higher Dimensions) is a generic C++
    library for computational topology.

    Author(s):       Pawel Dlotko

    Copyright (C) 2017  Swansea University

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

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"

epsilon = 0.0000005;

def test_check_construction_of_landscape():
    l = gudhi.PersistenceLandscapesOnGrid(file_with_intervals='data/file_with_diagram_1')
    l.print_to_file("landscape_from_file_with_diagram_1")
    g = gudhi.PersistenceLandscapesOnGrid()
    g.load_from_file("landscape_from_file_with_diagram_1")
    assert l.check_if_the_same(g)
    
def test_check_computations_of_integrals():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", number_of_points_=100 )
    integral = p.compute_integral_of_landscape()
    assert math.fabs(integral - 27.343) <= 0.00005
    
def test_check_construction_of_landscape_using_only_ten_levels():	
	l = gudhi.PersistenceLandscapesOnGrid(file_with_intervals='data/file_with_diagram_1', number_of_points_=100 )
	g = gudhi.PersistenceLandscapesOnGrid(file_with_intervals='data/file_with_diagram_1', number_of_points_=100 )
	for level in range(1,9):
		v1 = l.vectorize(level)
		v2 = g.vectorize(level)
		assert v1 == v2        

def test_check_computations_of_integrals_for_each_level_separatelly():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", number_of_points_=100)
    integrals_fir_different_levels = [0.241168,0.239276,0.237882,0.235193,0.230115,0.227626,0.226132,0.223643,0.221651,0.220556,0.21727,0.215976,0.213685,0.211993,0.2102,0.208707,0.207014,0.205122,0.204226,0.202633]
    for level in range(0,len(integrals_fir_different_levels)):
        integral = p.compute_integral_of_a_level_of_a_landscape(level);    
        assert math.fabs(integral - integrals_fir_different_levels[level]) <= 0.00005
  
def test_check_computations_of_values_on_different_points():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", number_of_points_=100)
    results_level_0 = [0.00997867,0.0521921,0.104312,0.156432,0.208552,0.260672,0.312792,0.364912,0.417032,0.429237]
    results_level_10 = [7.21433e-05,0.0422135,0.0943335,0.146453,0.198573,0.240715,0.272877,0.324997,0.359232,0.379344]
    x = 0.0012321;
    dx = 0.05212;
    for i in range(0,10):
        assert math.fabs(p.compute_value_at_a_given_point(0, x) - results_level_0[i]) < 0.00005;
        assert math.fabs(p.compute_value_at_a_given_point(10, x) - results_level_10[i]) < 0.00005;
        x += dx;
  
def test_check_default_parameters_of_distances():
    diag = gudhi.read_persistence_intervals_in_dimension(persistence_file="data/file_with_diagram")
    p = gudhi.PersistenceLandscapesOnGrid(vector_of_intervals=diag, grid_min_=0., grid_max_=1., number_of_points_=100)
    diag1 = gudhi.read_persistence_intervals_in_dimension(persistence_file="data/file_with_diagram_1")
    q = gudhi.PersistenceLandscapesOnGrid(vector_of_intervals=diag1, grid_min_=0., grid_max_=1., number_of_points_=100)
    dist_numeric_limit_max = p.distance(q, sys.maxsize);
    dist_infinity = p.distance(q, sys.maxsize);
    assert dist_numeric_limit_max == dist_infinity

def test_check_computations_of_distances():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram", grid_min_=0., grid_max_=1., number_of_points_=10000)
    q = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", grid_min_=0., grid_max_=1., number_of_points_=10000)
    assert math.fabs(p.distance(q,1) - 25.5779) <= 0.00005
    assert math.fabs(p.distance(q, 2) - 2.048980326807754) <= 0.00001
    assert math.fabs(p.distance(q, sys.float_info.max) - 0.359) <= 0.00001

def test_check_computations_of_scalar_product():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram", grid_min_=0., grid_max_=1., number_of_points_=10000)
    q = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", grid_min_=0., grid_max_=1., number_of_points_=10000)
    assert math.fabs(p.compute_scalar_product(q)- 0.754367) < 0.000001
  
def test_check_computations_of_maxima_and_norms():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", grid_min_=0., grid_max_=1., number_of_points_=100)
    assert math.fabs(p.compute_minimum_maximum()[1] - 0.46) <= 0.00001
    assert math.fabs(p.compute_norm_of_landscape(1) - 27.3373) <= 0.00001
    assert math.fabs(p.compute_norm_of_landscape(2) - 1.84143) <= 0.00001
    assert math.fabs(p.compute_norm_of_landscape(3) - 0.927067) <= 0.00001   

def test_check_computations_of_averages():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram", grid_min_=0., grid_max_=1., number_of_points_=100)
    q = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", grid_min_=0., grid_max_=1., number_of_points_=100)
    av = gudhi.PersistenceLandscapesOnGrid()
    av.compute_average([p,q])

    template_average = gudhi.PersistenceLandscapesOnGrid()
    template_average.load_from_file("data/average_on_a_grid")
    assert template_average.check_if_the_same(av)
 
def test_check_computations_of_integrals_of_powers_of_landscape():
    p = gudhi.PersistenceLandscapesOnGrid(file_with_intervals="data/file_with_diagram_1", number_of_points_=100)
    integrals_fir_different_powers = [146.567,17.7334,3.43279,0.812015,0.21863]
    for power in range(0,5):
        integral = p.compute_integral_of_power_of_landscape(power)
        assert math.fabs(integral - integrals_fir_different_powers[power]) <= 0.00001