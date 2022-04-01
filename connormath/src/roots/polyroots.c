 
#include "connormath.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
  if(len < 2) {
    // Can't do these.
    return (0);
  }

  // If the leading coefficient is zero, move to a different length.
  if(coefs[len - 1] == 0) {
    return (poly_roots(coefs, len - 1, out));
  }
  // If the trailing coefficient is zero, the polynomial can be reduced.
  if(coefs[0] == 0) {
    out[0] = 0;
    return (poly_roots(coefs + 1, len - 1, out + 1));
  }

  // First, if this polynomial is solvable, do that.
  if(len == 2) { // Linear
    out[0] = -coefs[0] / coefs[1];
    return (1);
  } else if(len == 3) { // Quadratic
    double a = coefs[2], b = coefs[1], c = coefs[0];
    if(b * b - 4 * a * c < 0) {
      return (0);
    }
    out[0] = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    out[1] = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    return (2);
  } else if(len == 4) { // Cubic
    double a = coefs[3], b = coefs[2], c = coefs[1], d = coefs[0],
      del0 = b * b - 3 * a * c,
      del1 = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d;
    double Cr, Ci;
    if(del0 == 0 && del1 == 0) {
      out[0] = out[1] = out[2] = -b / (3 * a);
      return (3);
    } else if(del0 == 0) {
      Cr = cbrt(del1);
      Ci = 0;
    } else if(del1 * del1 - 4 * del0 * del0 * del0 < 0) {
      // Handle complex stuff.
      double R0 = del1 * del1 - 4 * del0 * del0 * del0;
      
      Cr = cbrt((del1 / sqrt(del1 * del1 + abs(R0)) - 3) / 4);
      Ci = cbrt((3 - abs(R0) / (sqrt(del1 * del1 + abs(R0)))) / 4);
    } else {
      Cr = cbrt((del1 + sqrt(del1 * del1 - 4 * del0 * del0 * del0)) / 2);
      Ci = 0;
    }
    double C = hypot(Cr, Ci);
    double r0 = -(b + Cr + del0 * Cr / C) / (3 * a);
    double i0 = -(Ci - del0 * Ci / C) / (3 * a);
    double r1 = -(b + Cr/2 - sqrt(3) * Ci/2 + del0 * Cr / (2 * C) +
		  del0 * Ci * sqrt(3) / (2 * C)) / (3 * a);
    double i1 = -(-Ci / 2 + sqrt(3) * Cr / 2 + del0 * Ci / (2 * C) -
		  del0 * Cr * sqrt(3) / (2 * C)) / (3 * a);
    double r2 = -(b + Cr/2 + sqrt(3) * Ci/2 - del0 * Cr / (2 * C) -
		  del0 * Ci * sqrt(3) / (2 * C)) / (3 * a);
    double i2 = -(-Ci / 2 - sqrt(3) * Cr / 2 - del0 * Ci / (2 * C) +
		  del0 * Cr * sqrt(3) / (2 * C)) / (3 * a);
    int pos = 0;
    if(i0 == 0) {
      out[pos] = r0;
      pos++;
    }
    if(i1 == 0) {
      out[pos] = r1;
      pos++;
    }
    if(i2 == 0) {
      out[pos] = r2;
      pos++;
    }
    return (pos - 1);
  } else if(len == 5) { // Quartic
    // First, check if the polynomial can be written as a quadratic.
    if(coefs[1] == 0 and coefs[3] == 0) {
      if(b * b - 4 * a * c < 0) {
	return (0);
      }
      double a = coefs[4], b = coefs[2], c = coefs[0];
      out[0] = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
      out[1] = -(-b + sqrt(b * b - 4 * a * c)) / (2 * a);
      out[2] = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
      out[3] = -(-b - sqrt(b * b - 4 * a * c)) / (2 * a);
      return (4);
    } else {
      _Complex
    
      
  
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
  
			   
  
