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
 * @file BSplineCurve.cpp
 * @author Jorn Baayen
 * @version 1.0
 * @date 2015
 */

#include <assert.h>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "BSplineCurve.h"

#include "FCMangle.h"

using namespace fitpackpp;

extern "C" {
	void curfit(int *iopt, int *m, double *x, double *y, double *w, double *xb, double *xe, int *k, double *s, int *nest, int *n, double *t, double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier);
	void splev (double *t, int *n, double *c, int *k, double *x, double *y, int *m, int *e, int *ier);
	void splder(double *t, int *n, double *c, int *k, int *nu, double *x, double *y, int *m, int *e, double *wrk, int *ier);
}

/**
 * @brief Constructor
 * @details Construct a B-Spline curve interpolation for the points specified by the list of abscissae x and ordinates y.
 * 
 * @param x Abscissae
 * @param y Ordinates
 * @param preferredDegree Preferred degree of the interpolating spline. 
 * The actual degree is chosen such as to be one less than the number of data points, but no higher than preferredDegree.
 * @param smoothing Smoothing factor.  Must be non-negative. Set to 0.0, i.e., no smoothing, by default.
 */
BSplineCurve::BSplineCurve(std::vector<double> &x, std::vector<double> &y, int preferredDegree, double smoothing)
{
	// Number of data points
	int m = (int) x.size();

	// The actual degree of the spline must be less than m
	k = preferredDegree;
	if (k >= m)
		k = m - 1;

	// Configure curfit() parameters
	int iopt = 0;                       // Compute a smoothing spline
	int nest = m + k + 1;               // Over-estimate the number of knots

	// Allocate weighting vector
	double *w = new double[m];
	for (int i = 0; i < m; i++)
		w[i] = 1.0;

	// Allocate memory for knots and coefficients
	t = new double[nest];               // Knots
	c = new double[nest];               // Coefficients

	double fp; // Weighted sum of squared residuals

	// Allocate working memory required by curfit
	int     lwrk = (m * (k + 1) + nest * (7 + 3 * k));
	double *wrk  = new double[lwrk];
	int    *iwrk = new int   [nest];

	int ier;
	curfit(&iopt, &m, (double*) &x[0], (double*) &y[0], w, &x[0], &x[m - 1], &k, &smoothing, &nest, &n, t, c, &fp, wrk, &lwrk, iwrk, &ier);
	if (ier >= 10) {
		std::stringstream s;
		s << "Error fitting B-Spline curve using curfit(): " << ier;
		throw std::runtime_error(s.str());
	}

	// De-allocate temporary memory
	delete[] w;
	delete[] wrk;
	delete[] iwrk;

	// Allocate work vector for derivative computation
	wder = new double[n]; 
}

BSplineCurve::~BSplineCurve(void)
{
	// Free memory
	delete[] t;
	delete[] c;
	delete[] wder;
}

/**
 * @brief Evaluate the curve at point x
 * 
 * @param x Evaluation point
 * @return Curve ordinate at point x
 */
double BSplineCurve::eval(double x)
{
	double y;
	int m = 1; // Evaluate a single point
	int e = 0; // Don't clip argument to range
	int ier;
	splev(t, &n, c, &k, &x, &y, &m, &e, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating B-Spline curve using splev() at point " << x << ": " << ier;
		throw std::runtime_error(s.str());
	}

	return y;
}

/**
 * @brief Evaluate the order'th derivative of the curve at point x
 * 
 * @param x     Evaluation point
 * @param order Derivative order.  Defaults to 1
 * @return Derivative of the specified order, evaluated at point x
 */
double BSplineCurve::der(double x, int order)
{
	double y;
	int m = 1; // Evaluate a single point
	int e = 0; // Don't clip argument to range
	int ier;
	splder(t, &n, c, &k, &order, &x, &y, &m, &e, wder, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating " << order << "nth B-Spline curve derivative using splder() at point " << x << ": " << ier;
		throw std::runtime_error(s.str());
	}

	return y;
}