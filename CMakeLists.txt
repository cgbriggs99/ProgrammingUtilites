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

add_library(connormath SHARED)
target_link_libraries(connormath INTERFACE ${MATH} ${LAPACK_LIBRARIES}
  ${BLAS_LIBRAIRES} ${LAPACKE_LIBRARIES})

install(TARGETS connormath LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(connormath)

