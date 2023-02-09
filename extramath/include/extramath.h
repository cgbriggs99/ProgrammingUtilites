/*
 * extramath.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#define EXTRAMATH_H_

#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <stdint.h>

/*
 * If we want to precompte values and interpolate, here are the macros for that. The EXTRAMATH_DECL_NOPRE
 * will allow for a function to be defined without a precomputed version. If precomputed functions are turned off,
 * then EXTRAMATH_DECL and EXTRAMATH_DECL_NOPRE should give the same result.
 */
#define __FORTRAN__(fname) __CONCAT(fname, _fortran)
#define EXTRAMATH_DECL(fname, args) extern __TYPE_NAME__ __FNAME__(fname) args __attribute__((nothrow))
#define EXTRAMATH_DECL_FORT(fname, args) extern __TYPE_NAME__ __FNAME__(fname) args __attribute__((nothrow)); \
	extern __TYPE_NAME__ __FORTRAN__(__FNAME__(fname)) args
#define EXTRAMATH_ARRDECL(fname,args) extern int __FNAME__(fname) args __attribute__((nothrow))
#define EXTRAMATH_DECL_NOPRE(fname, args) extern __TYPE_NAME__ __FNAME__(fname) args __attribute__((nothrow)); \
	extern __TYPE_NAME__ __FORTRAN__(__FNAME__(fname)) args

#define __TYPE_NAME__ double
#define __SCALAR_TYPE__ double
#define __FNAME__(fname) fname
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __SCALAR_TYPE__
#undef __TYPE_NAME__
#undef __FNAME__

#define __TYPE_NAME__ long double
#define __SCALAR_TYPE__ long double
#define __FNAME__(fname) fname##l
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ float
#define __SCALAR_TYPE__ float
#define __FNAME__(fname) fname##f
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ double _Complex
#define __SCALAR_TYPE__ double
#define __FNAME__(fname) c##fname
#define __IS_COMPLEX_HEADER__
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __IS_COMPLEX_HEADER__
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ long double _Complex
#define __SCALAR_TYPE__ long double
#define __FNAME__(fname) c##fname##l
#define __IS_COMPLEX_HEADER__
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __IS_COMPLEX_HEADER__
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ float _Complex
#define __SCALAR_TYPE__ float
#define __FNAME__(fname) c##fname##f
#define __IS_COMPLEX_HEADER__
#include "roots.h"
#include "bessels.h"
#include "polynomials.h"
#include "integrated.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "elliptic.h"
#include "zeta.h"
#include "spheroidal.h"
#include "numbertheory_frac.h"
#include "statistics.h"
#include "solver.h"
#include "utils.h"
#include "longlinalg.h"
#undef __IS_COMPLEX_HEADER__
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ int32_t
#define __SCALAR_TYPE__ int32_t
#define __FNAME__(fname) fname
#include "integer.h"
#include "numbertheory_int.h"
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ int64_t
#define __SCALAR_TYPE__ int64_t
#define __FNAME__(fname) fname##l
#include "integer.h"
#include "numbertheory_int.h"
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define __TYPE_NAME__ long long
#define __SCALAR_TYPE__ long long
#define __FNAME__(fname) fname##L
#include "integer.h"
#include "numbertheory_int.h"
#undef __TYPE_NAME__
#undef __SCALAR_TYPE__
#undef __FNAME__

#define MASCHERONI 0.57721566490153286060651209008240243
#define APERY 1.2020569031595942853997381615114499907649862923404988817922715553

#endif /* EXTRAMATH_H_ */
