 
#include "connormath.h"
#include <stdlib.h>
#include <string.h>

// Divide polynomials.
static int divide(double *num, int ndim, double *den, int ddim, double **out,
		  int *outdim) {
  double *rem = calloc(ndim, sizeof(double));
  double q;
  *out = malloc((ndim - ddim) * sizeof(double));
  *outdim = ndim - ddim;
  memcpy(rem, num, ndim * sizeof(double));
  
  for(int i = 0; i < *outdim; i++) {
    q = rem[ndim - i - 1] / den[ddim - 1];
    *out[*outdim - i - 1] = q;
    for(int j = 0; j < ddim; j++) {
      rem[ndim - i - 1 - j] -= q * den[ddim - j - 1];
    }
  }
  free(rem);
  return (0);
}


int poly_roots(double *coefs, int len, double *out) {
  // First, if this polynomial is solvable, do that.
  if(len < 2) {
    // Can't do these.
    return (-1);
  } else if(len == 2) {
    out[0] = -coefs[0] / coefs[1];
    return (0);
  } else if(len == 3) {
    double a = coefs[2], b = coefs[1], c = coefs[0];
    if(b * b - 4 * a * c < 0) {
      return (-1);
    }
    out[0] = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    out[1] = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    return (0);
  }
  
  // First, find the Gershgorin circles.
  // Centers at 0 and -c_(n-1). Radii are 1 + max(coef) up to n-2, and 1.

  // If one circle is disjoint, then use it to find a root.
  double high_sep = coefs[dim - 2] / coefs[dim - 1] + 1,
    low_sep = coefs[dim - 2] / coesf[dim - 1] - 1;
  double rad_origin = 0;
  for(int i = 0; i < dim - 2; i++) {
    if(1 + abs(coefs[i] / coefs[dim - 1]) > rad_origin) {
      rad_origin = 1 + abs(coefs[i] / coefs[dim - 1]);
    }
  }
  
			   
  
