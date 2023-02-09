cmake_minimum_required(VERSION 3.13)

project(ProgrammingUtilities)

include(CheckLibraryExists)
include(FindBLAS)
include(FindLAPACK)
include(CTest)

check_library_exists(m sin "" MATH)

find_library(LAPACKE_LIBRARIES lapacke)

find_path(LAPACKE_INC lapacke.h)

if(${MATH})
  set(MATH m)
endif()

add_library(extramath SHARED)
add_dependencies(extramath float double longdouble complex_float
  complex_double complex_longdouble)
#  $<TARGET_OBJECTS:float>
#  $<TARGET_OBJECTS:double>
#  $<TARGET_OBJECTS:longdouble>
#  $<TARGET_OBJECTS:complex_float>
#  $<TARGET_OBJECTS:complex_double>
#  $<TARGET_OBJECTS:complex_longdouble>)

target_link_libraries(extramath INTERFACE ${MATH} ${LAPACK_LIBRARIES}
  ${BLAS_LIBRAIRES} ${LAPACKE_LIBRARIES})

install(TARGETS extramath LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(extramath)
