/**
 * This file is part of the point set analysis tool psa
 * 
 * Copyright 2011, Thomas Schl\"{o}mer, thomas.schloemer@uni-konstanz.de
 * 
 * psa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <list>
#include "util.h"


/**
 *  Computes mean, standard deviation, minimum, and maximum for the
 *  values of a given array.
 *
 *  Algorithm based on the recurrence formulas by Welford, B. P. from
 *  [Welford 1962], Technometrics 4, pp. 419--420.
 *
 *  Also cf. [Knuth 68], TAOCP, Vol. 2, Seminumerical Algorithms,
 *  pp. 231--232.
 */
void BasicStats(const std::list<float> &values, struct stats &s);


/**
 *  L_{2}-norm based star discrepancy for two-dimensional point sets.
 *
 *  Algorithm from Warnock, T. T., "Computational Investigations of
 *  Low-Discrepancy Point Sets", Applications of Number Theory to
 *  Numerical Analysis, pp. 319--343, 1972.
 *
 *  For a faster computation in O(n log(n)^2), cf. Heinrich S.,
 *  "Efficient algorithms for computing the L_{2} discrepancy", 
 *  Math. Comput. 65, pp. 1621--1633, 1996.
 *
 *  For an introduction related to computer graphics, cf. Shirley P.,
 *  "Discrepancy as a Quality Measure for Sample Distributions",
 *  Proceedings of Eurographics, pp. 183--194, 1991.
 */
float L2normStarDiscrepancy(const float * const points, int npoints);


/**
 *  Distance measures based on the local minimum distance
 *    d_x = \min_{y\in X\backslash\{x\}} d_T(x,y)
 *  which denotes the mininum distance from the point x to any other point
 *  in the point set under consideration of the unit torus.
 *
 *  The computed measures are
 *    mindist    = \min_{x\in X} d_x         (global mindist)
 *    avgmindist = 1/|X| \sum{x\in X} d_x    (average mindist)
 *
 *  This is a naive implementation without any acceleration structure.
 */
void DistanceMeasures(const float * const points, int npoints,
                      float &mindist, float &avgmindist);


#endif    // MEASUREMENTS_H

