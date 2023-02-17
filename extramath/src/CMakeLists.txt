
cmake_minimum_required(VERSION 3.13)

project(ProgrammingUtilities)

configure_file(polyroots.pre.c polyroots.c)

add_library(float OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(float PUBLIC __TYPENAME__=float
   __SCALARTYPE__=float
   __TYPESUFF__=f
   __TYPEVAL__=__TYPEVAL_FLOAT__)

add_library(double OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(double PUBLIC __TYPENAME__=double
   __SCALARTYPE__=double
   __TYPESUFF__=
   __TYPEVAL__=__TYPEVAL_DOUBLE__)

add_library(longdouble OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(longdouble PUBLIC __TYPENAME__=long\ double
   __SCALARTYPE__=long\ double
   __TYPESUFF__=l
   __TYPEVAL__=__TYPEVAL_LONGDOUBLE__)

add_library(complex_float OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(complex_float PUBLIC __TYPENAME__=float\ _Complex
   __SCALARTYPE__=float
   __TYPESUFF__=f
   __TYPEVAL__=__TYPEVAL_COMPLEXFLOAT__
   __IS_COMPLEX__)

add_library(complex_double OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(complex_double PUBLIC __TYPENAME__=double\ _Complex
   __SCALARTYPE__=double
   __TYPESUFF__=
   __TYPEVAL__=__TYPEVAL_COMPLEXDOUBLE__
   __IS_COMPLEX__)

add_library(complex_longdouble OBJECT bessels.c beta.c decompose.c
   eigenvals.c erf.c gamma.c gauss_quad.c hypergeom.c integer.c
   integrals_nd.c integrated.c matarith.c matprops.c polynomials.c
   polyroots.c riemann_int.c tanhsinh.c utils.c zeta.c)
target_compile_definitions(complex_longdouble PUBLIC __TYPENAME__=long\ double\ _Complex
   __SCALARTYPE__=long\ double
   __TYPESUFF__=l
   __TYPEVAL__=__TYPEVAL_COMPLEXLONGDOUBLE__
   __IS_COMPLEX__)

 
target_sources(extramath PUBLIC $<TARGET_OBJECTS:float>
  $<TARGET_OBJECTS:double>
  $<TARGET_OBJECTS:longdouble>
  $<TARGET_OBJECTS:complex_float>
  $<TARGET_OBJECTS:complex_double>
  $<TARGET_OBJECTS:complex_longdouble>)
