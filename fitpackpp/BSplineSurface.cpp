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
 * @file BSplineSurface.cpp
 * @author Jorn Baayen
 * @version 1.0
 * @date 2015
 */

#include <assert.h>
#include <limits>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "BSplineSurface.h"

using namespace fitpackpp;

#ifdef _WIN32
#define surfit SURFIT
#define bispev BISPEV
#define parder PARDER
#else
#define surfit surfit_
#define bispev bispev_
#define parder parder_
#endif

extern "C" {
	void surfit(int *iopt, int *m, double *x, double *y, double *z, double *w, double *xb, double *xe, double *yb, double *ye, int *kx, int *ky,
		        double *s, int *nxest, int *nyest, int *nmax, double *eps, int *nx, double *tx, int *ny, double *ty, double *c, double *fp,
		        double *wrk1, int *lwrk1, double *wrk2, int *lwrk2, int *iwrk, int *kwrk, int *ier);
	void bispev(double *tx, int *nx, double *ty, int *ny, double *c, int *kx, int *ky, double *x, int *mx, double *y, int *my, double *z,
		        double *wrk, int *lwrk, int *iwrk, int *kwrk, int *ier);
	void parder(double *tx, int *nx, double *ty, int *ny, double *c, int *kx, int *ky, int *nux, int *nuy, double *x, int *mx, double *y, int *my,
		        double *z, double *wrk, int *lwrk, int *iwrk, int *kwrk, int *ier);
}

/**
 * @brief Constructor
 * @details Construct a B-Spline surface interpolation for the points specified by the list of data points.
 * 
 * @param x Data points:  X coordinates
 * @param y Data points:  Y coordinates
 * @param z Data points:  Z coordinates
 * @param preferredDegree Preferred degree of the interpolating spline. 
 * The actual degree is chosen such as to be one less than the number of data points, but no higher than preferredDegree.
 * @param smoothing Smoothing factor.  Must be non-negative. Set to 0.0, i.e., no smoothing, by default.
 */
BSplineSurface::BSplineSurface(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, int preferredDegree, double smoothing)
{
	// Number of data points
	int m = (int) x.size();

	// Compute feasible spline degree
	k = preferredDegree;
	if ((k + 1) * (k + 1) > m)
		k = (int) floor(sqrt((double) m) - 1);

	// Configure surfit() parameters
	int iopt = 0;                                     // Compute a smoothing spline
	int nest = m + k + 1;                             // Over-estimate the number of knots

	// Allocate weighting vector
	double *w = new double[m];
	for (int i = 0; i < m; i++)
		w[i] = 1.0;

	// Allocate memory for knots and coefficients
	tx = new double[nest];                            // X knots
	ty = new double[nest];                            // Y knots
	c  = new double[(nest - k - 1) * (nest - k - 1)]; // Coefficients

	double fp; // Weighted sum of squared residuals

	// Allocate working memory required by surfit
	int  u = nest - k - 1;
	int km = k + 1;
	int b1 = k * u + k + 1;
	int b2 = b1 + u - k;

	int     lwrk1 = u * u * (2 + b1 + b2) + 2 * (u + u + km * (m + nest) + nest - k - k) + b2 + 1;
	double *wrk1  = new double[lwrk1];

	int     lwrk2 = u * u * (b2 + 1) + b2;
	double *wrk2  = new double[lwrk2];

	int     kwrk  = m + (nest - 2 * k - 1) * (nest - 2 * k - 1);
	int    *iwrk  = new int   [kwrk];

	double eps = std::numeric_limits<double>::epsilon();

	int ier;
	surfit(&iopt, &m, (double*) &x[0], (double*) &y[0], (double*) &z[0], w, &x[0], &x[m - 1], &y[0], &y[m - 1], &k, &k, &smoothing, &nest, &nest, &nest, &eps, &nx, tx, &ny, ty, c, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
	if (ier >= 10) {
		std::stringstream s;
		s << "Error fitting B-Spline surface using surfit(): " << ier;
		throw std::runtime_error(s.str());
	}

	// De-allocate temporary memory
	delete[] w;
	delete[] wrk1;
	delete[] wrk2;
	delete[] iwrk;

	// Allocate work vectors for spline (and derivative) evaluation
	lwrk = 2 * (k + 1) + (nx - k - 1) * (ny - k - 1);
	wrk  = new double[lwrk];
}

BSplineSurface::~BSplineSurface(void)
{
	// Free memory
	delete[] tx;
	delete[] ty;
	delete[] c;
	delete[] wrk;
}

/**
 * @brief Evaluate the surface at point x
 * 
 * @param x Evaluation point:  X coordinate
 * @param y Evaluation point:  Y coordinate
 * @return Surface Z coordinate at point x
 */
double BSplineSurface::eval(double x, double y)
{
	double z;
	int m = 1; // Evaluate a single point
	int kwrk = 2;
	int iwrk[2];
	int ier;
	bispev(tx, &nx, ty, &ny, c, &k, &k, &x, &m, &y, &m, &z, wrk, &lwrk, iwrk, &kwrk, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating B-Spline surface using bispev() at point " << x << ": " << ier;
		throw std::runtime_error(s.str());
	}

	return z;
}

/**
 * @brief Evaluate the partial derivative of the surface at point x,
 * with X direction order xOrder and Y direction order yOrder.
 * 
 * @param x      Evaluation point:  X coordinate
 * @param y      Evaluation point:  Y coordinate
 * @param xOrder Derivative order in the X direction
 * @param yOrder Derivative order in the Y direction
 * @return Partial derivative of the specified order, evaluated at point x
 */
double BSplineSurface::der(double x, double y, int xOrder, int yOrder)
{
	double z;
	int m = 1; // Evaluate a single point
	int kwrk = 2;
	int iwrk[2];
	int ier;
	parder(tx, &nx, ty, &ny, c, &k, &k, &xOrder, &yOrder, &x, &m, &y, &m, &z, wrk, &lwrk, iwrk, &kwrk, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating (" << xOrder << ", " << yOrder << ")th partial B-Spline surface derivative using parder() at point " << x << ": " << ier;
		throw std::runtime_error(s.str());
	}

	return z;
}