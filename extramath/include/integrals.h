/*
 * integrals.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

#define __KERNEL__ __FNAME__(kernel_)

typedef __TYPE_NAME__ (*__KERNEL__)(__SCALAR_TYPE__, const void *);

#define __INTEGRAL_NONE__ 0
#define __INTEGRAL_UPPER__ 1
#define __INTEGRAL_LOWER__ 2
#define __INTEGRAL_BOTH__ 3

/*
 * Riemann sums.
 */
EXTRAMATH_DECL(lriemannint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));
EXTRAMATH_DECL(rriemannint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));
EXTRAMATH_DECL(mriemannint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_DECL(trapezoidint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Simpson's rule.
 */
EXTRAMATH_DECL(simpsonint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_DECL(gausslegendreint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussjacobiint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, __SCALAR_TYPE__ __alpha, __SCALAR_TYPE__ __beta,
		unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint1, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint2, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_DECL(gausslaguerreint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussgenlaguerreint, (__KERNEL__ __func,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __alpha, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausshermiteint, (__KERNEL__ __func,
		unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
EXTRAMATH_DECL(tanhsinhint, (__KERNEL__ __func, __SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points,
		const void *__extra));

// Multidimensional integrals.
#define __KERNEL_ND__ __FNAME__(kernel_nd_)

typedef double (*__KERNEL_ND__)(const __SCALAR_TYPE__ *in, int indim, const void *pass_data);

/*
 * Riemann sums.
 */
EXTRAMATH_DECL(lriemannintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const int *include_ends,
		const void *__extra));
EXTRAMATH_DECL(rriemannintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const int *include_ends,
		const void *__extra));
EXTRAMATH_DECL(mriemannintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_DECL(trapezoidintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const int *include_ends,
		const void *__extra));

/*
 * Simpson's rule.
 */
//EXTRAMATH_DECL(simpsonintnd, (__KERNEL_ND__ __func, int indim,
//		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const int *include_ends,
//		const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_DECL(gausslegendreintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussjacobiintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, const __SCALAR_TYPE__ *__alphas,
		const __SCALAR_TYPE__ *__betas, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint1nd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint2nd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_DECL(gausslaguerreintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussgenlaguerreintnd, (__KERNEL_ND__ __func, int indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__alphas, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausshermiteintnd, (__KERNEL_ND__ __func, int indim,
		unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
//EXTRAMATH_DECL(tanhsinhintnd, (__KERNEL_ND__ __func, int indim,
//		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

/*
 * Combination. If one of the requested integrals is improper, it will not consider the infinite endpoint.
 */
#ifndef __INTEGRAL_SELECTOR__
#define __INTEGRAL_SELECTOR__
typedef enum {
	RIEMANN_LEFT,
	RIEMANN_RIGHT,
	RIEMANN_MID,
	TRAPEZOID,
	SIMPSON,
	GAUSSLEGENDRE,
	GAUSSJACOBI,
	GAUSSCHEBYCHEV1,
	GAUSSCHEBYCHEV2,
	GAUSSLAGUERRE,
	GAUSSGENLAGUERRE,
	GAUSSHERMITE,
	TANSINH
} integral_selector_t;
#endif

EXTRAMATH_DECL(genintegralnd, (__KERNEL_ND__ __func, int __indim,
		const __SCALAR_TYPE__ *__starts, const __SCALAR_TYPE__ *__ends,
		const integral_selector_t *__kinds, const __SCALAR_TYPE__ *__alphas, const __SCALAR_TYPE__ *__betas,
		unsigned int __points, const unsigned int *__includes, const void *__extra));

// Multidimensional integrals with varying boundaries.
#define __BOUNDARY_FUNCTION__ __FNAME__(bounds_)
typedef void (*__BOUNDARY_FUNCTION__)(const __SCALAR_TYPE__ *__curr_pos, int __ncurrdim, __SCALAR_TYPE__ *low,
		__SCALAR_TYPE__ *high, const void *data);

/*
 * Riemann sums.
 */
EXTRAMATH_DECL(lriemannintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const int *include_ends,
		const void *__extra));
EXTRAMATH_DECL(rriemannintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const int *include_ends,
		const void *__extra));
EXTRAMATH_DECL(mriemannintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_DECL(trapezoidintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const int *include_ends,
		const void *__extra));

/*
 * Simpson's rule.
 */
EXTRAMATH_DECL(simpsonintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const int *include_ends,
		const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_DECL(gausslegendreintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussjacobiintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, const __SCALAR_TYPE__ *__alphas,
		const __SCALAR_TYPE__ *__betas, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint1nd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gausschebychevint2nd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_DECL(gausslaguerreintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, unsigned int __points, const void *__extra));
EXTRAMATH_DECL(gaussgenlaguerreintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start0, const __SCALAR_TYPE__ *__alphas, unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
EXTRAMATH_DECL(tanhsinhintnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int indim,
		__SCALAR_TYPE__ __start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

/*
 * Combination. If one of the requested integrals is improper, it will not consider the infinite endpoint.
 */

EXTRAMATH_DECL(genintegralnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds, int __indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0,
		const integral_selector_t *__kinds, const __SCALAR_TYPE__ *__alphas, const __SCALAR_TYPE__ *__betas,
		unsigned int __points, const unsigned int *__includes, const void *__extra));

/*
 * Vector integrals.
 */
#define __KERNEL_VEC__ __FNAME__(kernel_vec_)
typedef void (*__KERNEL_VEC__)(__SCALAR_TYPE__ __x, __TYPE_NAME__ *out, int nouts, const void *);

/*
 * Riemann sums.
 */
EXTRAMATH_ARRDECL(lriemannintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));
EXTRAMATH_ARRDECL(rriemannintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));
EXTRAMATH_ARRDECL(mriemannintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_ARRDECL(trapezoidintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Simpson's rule.
 */
EXTRAMATH_ARRDECL(simpsonintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_ARRDECL(gausslegendreintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussjacobiintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, __SCALAR_TYPE__ __alpha, __SCALAR_TYPE__ __beta,
		unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint1vec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint2vec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_ARRDECL(gausslaguerreintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussgenlaguerreintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __alpha, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausshermiteintvec, (__KERNEL_VEC__ __func, __TYPE_NAME__ *__out, int __outdim,
		unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
EXTRAMATH_ARRDECL(tanhsinhintvec, (__KERNEL_VEC__ __func,__TYPE_NAME__ *__out, int __outdim,
		__SCALAR_TYPE__ __start, __SCALAR_TYPE__ __end, unsigned int __points, const void *__extra));


/*
 * Multi-dimensional vector integrals.
 */
#define __KERNEL_VEC_ND__ __FNAME__(kernel_vec_nd_)
typedef void (*__KERNEL_VEC_ND__)(const __SCALAR_TYPE__ *__x, int indim, __TYPE_NAME__ *out, int nouts, const void *);

/*
 * Riemann sums.
 */
EXTRAMATH_ARRDECL(lriemannintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends,
		const void *__extra));
EXTRAMATH_ARRDECL(rriemannintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends,
		const void *__extra));
EXTRAMATH_ARRDECL(mriemannintvecnd, (__KERNEL_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_ARRDECL(trapezoidintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Simpson's rule.
 */
EXTRAMATH_ARRDECL(simpsonintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends,
		const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_ARRDECL(gausslegendreintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussjacobiintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, const __SCALAR_TYPE__ *__alphas,
		const __SCALAR_TYPE__ *__betas, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint1vecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint2vecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_ARRDECL(gausslaguerreintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussgenlaguerreintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__alphas, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausshermiteintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
EXTRAMATH_ARRDECL(tanhsinhintvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

EXTRAMATH_ARRDECL(genintegralvecnd, (__KERNEL_VEC_ND__ __func, __TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__starts, const __SCALAR_TYPE__ *__ends,
		const integral_selector_t *__kinds, const __SCALAR_TYPE__ *__alphas, const __SCALAR_TYPE__ *__betas,
		unsigned int __points, const unsigned int *__includes, const void *__extra));

/*
 * Multidimensional vector integrals with generalized bounds.
 */

/*
 * Riemann sums.
 */
EXTRAMATH_ARRDECL(lriemannintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, int include_ends,
		const void *__extra));
EXTRAMATH_ARRDECL(rriemannintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, int include_ends,
		const void *__extra));
EXTRAMATH_ARRDECL(mriemannintvecnd_genbounds, (__KERNEL_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		__SCALAR_TYPE__ __start0, __SCALAR_TYPE__ __end0, unsigned int __points, const void *__extra));

/*
 * Trapezoidal rule.
 */
EXTRAMATH_ARRDECL(trapezoidintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends, const void *__extra));

/*
 * Simpson's rule.
 */
EXTRAMATH_ARRDECL(simpsonintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, int include_ends,
		const void *__extra));

/*
 * Gaussian quadrature.
 */
//Propers
EXTRAMATH_ARRDECL(gausslegendreintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussjacobiintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, const __SCALAR_TYPE__ *__alphas,
		const __SCALAR_TYPE__ *__betas, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint1vecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausschebychevint2vecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));
//Impropers
EXTRAMATH_ARRDECL(gausslaguerreintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gaussgenlaguerreintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__alphas, unsigned int __points, const void *__extra));
EXTRAMATH_ARRDECL(gausshermiteintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		unsigned int __points, const void *__extra));

/*
 * tanh-sinh quadrature.
 */
EXTRAMATH_ARRDECL(tanhsinhintvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__start, const __SCALAR_TYPE__ *__end, unsigned int __points, const void *__extra));

EXTRAMATH_ARRDECL(genintegralvecnd_genbounds, (__KERNEL_VEC_ND__ __func, __BOUNDARY_FUNCTION__ __bounds,
		__TYPE_NAME__ *__out, int __outdim, int __indim,
		const __SCALAR_TYPE__ *__starts, const __SCALAR_TYPE__ *__ends,
		const integral_selector_t *__kinds, const __SCALAR_TYPE__ *__alphas, const __SCALAR_TYPE__ *__betas,
		unsigned int __points, const unsigned int *__includes, const void *__extra));


#undef __KERNEL__
#undef __KERNEL_ND__
#undef __KERNEL_VEC__
#undef __BOUNDARY_FUNCTION__
