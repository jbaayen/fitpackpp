fitpackpp: C++ bindings for FITPACK
===

Fitpackpp wraps the package of Fortran subroutines for smoothing splines by P. Dierckx, [FITPACK](http://www.netlib.org/dierckx). This wrapper uses the double precision version of FITPACK distributed with [scipy](http://www.scipy.org).

The BSplineCurve class wraps the 1D routines "curfit", "splev", and "splder".

The BSplineSurface class wraps the 2D routines "surfit", "bispev", and "parder".

References
===

Paul Dierckx, *Curve and Surface Fitting with Splines*, Oxford University Press, 1993.
