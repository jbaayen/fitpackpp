// Copyright (C) 2015 Deltares
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License version 2.1 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

/**
 * @file BSplineSurface.h
 * @author Jorn Baayen
 * @version 1.0
 * @date 2015
 */

#ifndef BSPLINE_SURFACE_H
#define BSPLINE_SURFACE_H

#include <vector>

namespace fitpackpp
{

/**
 * @brief B-Spline surface
 * @details A wrapper for the surfit(), bispev(), and parder() routines from FITPACK by  P. Dierckx:
 * http://www.netlib.org/fitpack/index.html
 * More specifically, we use the double-precision FITPACK version included with scipy:
 * https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 */
class BSplineSurface
{
public:
	BSplineSurface(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, int preferredDegree=3, double smoothing=0.0);

	~BSplineSurface();

	double eval(double x, double y);
	double der(double x, double y, int xOrder, int yOrder);

private:
	int     k;    // Spline degree
	int     nx;   // Number of knots in X direction
	int     ny;   // Number of knots in Y direction
	double *tx;   // Knot X coordinates
	double *ty;   // Knot Y coordinates
	double *c;    // Spline coefficients
	int     lwrk; // Working memory for spline (and derivative) evaluation: length
};

};

#endif