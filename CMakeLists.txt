enable_language(Fortran)

set(fitpack_SRC
    fitpack/bispev.f
    fitpack/curfit.f
    fitpack/fpback.f
    fitpack/fpchec.f
    fitpack/fpcurf.f
    fitpack/fpbspl.f
    fitpack/fpdisc.f
    fitpack/fpgivs.f
    fitpack/fpknot.f
    fitpack/fprati.f
    fitpack/fprota.f
    fitpack/parder.f
    fitpack/splder.f
    fitpack/splev.f
    fitpack/surfit.f
)

add_library(fitpack STATIC 
            ${fitpack_SRC}
)

enable_language(CXX)

set(fitpackpp_SRC
	fitpackpp/BSplineCurve.cpp)

add_library(fitpackpp SHARED
	        ${fitpackpp_SRC})

target_link_libraries(fitpackpp fitpack)