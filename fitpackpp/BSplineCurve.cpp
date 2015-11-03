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
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ios>

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
	if (k >= m) {
		k = m - 1;

		std::cerr << "WARNING:  Too few data points (" << m << ") to create B-Spline curve of order " << preferredDegree << ". Reducing order to " << k << "." << std::endl;
	}

	// Configure curfit() parameters
	int iopt = 0;                       // Compute a smoothing spline
	int nest = m + k + 1;               // Over-estimate the number of knots

	// Allocate weighting vector
	double *w = new double[m];
	std::fill(w, w + m, 1.0);

	// Allocate memory for knots and coefficients
	t = new double[nest];               // Knots
	std::fill(t, t + nest, 0.0);

	c = new double[nest];               // Coefficients
	std::fill(c, c + nest, 0.0);

	double fp = 0.0; // Weighted sum of squared residuals

	// Allocate working memory required by curfit
	int     lwrk = (m * (k + 1) + nest * (7 + 3 * k));
	double *wrk  = new double[lwrk];
	std::fill(wrk, wrk + lwrk, 0.0);

	int    *iwrk = new int   [nest];
	std::fill(iwrk, iwrk + nest, 0);

	int ier = 0;
	curfit(&iopt, &m, (double*) &x[0], (double*) &y[0], w, &x[0], &x[m - 1], &k, &smoothing, &nest, &n, t, c, &fp, wrk, &lwrk, iwrk, &ier);
	if (ier > 0) {
		if (ier >= 10) {
			std::stringstream s;
			s << "Error fitting B-Spline curve using curfit(): " << ier;
			throw std::runtime_error(s.str());
		} else {
			std::cerr << "WARNING:  Non-fatal error while fitting B-Spline curve using curfit(): " << ier << std::endl;
		}
	}

	// De-allocate temporary memory
	delete[] w;
	delete[] wrk;
	delete[] iwrk;
}

/**
 * @brief Constructor
 * @details Construct a B-Spline curve interpolation for the given knots and coefficients.
 * 
 * @param knotX Knot X coordinates
 * @param coefs B-Spline coefficients
 * @param degree Preferred degree of the interpolating spline. 
 */
BSplineCurve::BSplineCurve(std::vector<double> &knotX, std::vector<double> &coefs, int degree)
{
	// Store parameters
	k = degree;

	n = (int) knotX.size();

	t = new double[knotX.size()];
	std::copy(knotX.begin(), knotX.end(), t);

	c = new double[coefs.size()];
	std::copy(coefs.begin(), coefs.end(), c);
}

/**
 * @brief Constructor
 * @details Construct a B-Spline curve interpolation from a previously serialized BSplineCurve object
 * 
 * @param filename File to load 
 */
BSplineCurve::BSplineCurve(const std::string &filename)
{
	std::ifstream f;
	f.open(filename, std::ios::binary);
	if (!f.good()) {
		std::stringstream s;
		s << "Failed to open B-Spline curve cache file: " << filename;
		throw std::runtime_error(s.str());
	}

	// Read spline from cache
	f.read((char*) &n, sizeof(int));
	t = new double[n];
	for (int i = 0; i < n; i++) {
		f.read((char*) &t[i], sizeof(double));
	}

	int nc;
	f.read((char*) &nc, sizeof(int));
	c = new double[nc];
	for (int i = 0; i < nc; i++) {
		f.read((char*) &c[i], sizeof(double));
	}

	f.read((char*) &k, sizeof(int));

	f.close();
}

BSplineCurve::~BSplineCurve(void)
{
	// Free memory
	delete[] t;
	delete[] c;
}

/**
 * @brief Knot X coordinates
 * @return Knot X coordinates
 */
std::vector<double> BSplineCurve::knotX() 
{
	std::vector<double> knotX;
	knotX.assign(t, t + n);
	return knotX;
}

/**
 * @brief Coefficients
 * @return Coefficients
 */
std::vector<double> BSplineCurve::coefs() 
{
	std::vector<double> coefs;
	coefs.assign(c, c + n);
	return coefs;
}	

/**
 * @brief Degree
 * @return Degree
 */
int BSplineCurve::degree()
{
	return k;
}

/**
 * @brief Serialize the BSplineCurve object
 * 
 * @param filename Destination file
 */
void BSplineCurve::serialize(const std::string &filename)
{
	std::ofstream f;
	f.open(filename, std::ios::binary);

	f.write((char*) &n, sizeof(int));
	for (int i = 0; i < n; i++) {
		f.write((char*) &t[i], sizeof(double));
	}

	f.write((char*) &n, sizeof(int));
	for (int i = 0; i < n; i++) {
		f.write((char*) &c[i], sizeof(double));
	}

	f.write((char*) &k, sizeof(int));

	f.close();
}

/**
 * @brief Evaluate the curve at point x
 * 
 * @param x Evaluation point
 * @return Curve ordinate at point x
 */
double BSplineCurve::eval(double x)
{
	double y = 0.0;
	int m = 1; // Evaluate a single point
	int e = 0; // Don't clip argument to range
	int ier = 0;

	// splev also evaluates points in the exterior
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
	if (order < 0 || order >= k) {
		std::stringstream s;
		s << "Cannot evaluate order " << order << " derivative of B-Spline curve of order " << k;
		throw std::runtime_error(s.str());
	}

	// Allocate working memory on the stack to keep this function thread-safe
	double *wder = (double*) alloca(sizeof(double) * n);
	std::fill(wder, wder + n, 0.0);

	double y = 0.0;
	int m = 1; // Evaluate a single point
	int e = 0; // Don't clip argument to range
	int ier = 0;

	// splder also evaluates points in the exterior
	splder(t, &n, c, &k, &order, &x, &y, &m, &e, wder, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating order " << order << " B-Spline curve derivative using splder() at point " << x << ": " << ier;
		throw std::runtime_error(s.str());
	}

	return y;
}