 
#include "connormath.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <errno.h>

// Divide polynomials.
static int divide(const double *num, int ndim, const double *den, int ddim,
		  double *out, int *outdim) {
  double *rem = calloc(ndim, sizeof(double));
  double q;
  *outdim = ndim - ddim;
  memcpy(rem, num, ndim * sizeof(double));
  
  for(int i = 0; i < *outdim; i++) {
    q = rem[ndim - i - 1] / den[ddim - 1];
    out[*outdim - i - 1] = q;
    for(int j = 0; j < ddim; j++) {
      rem[ndim - i - 1 - j] -= q * den[ddim - j - 1];
    }
  }
  int remdeg = 0;
  for(int i = 0; i < ndim; i++) {
    if(fabs(rem[0]) > 1e-6) {
      remdeg++;
    }
  }
  free(rem);
  return (remdeg);
}

static _Complex double eval_poly(const double *coefs, int len, _Complex double pos) {
  _Complex double sum = 0, xp = 1;
  for(int i = 0; i < len; i++) {
    sum += coefs[i] * xp;
    xp *= pos;
  }
  return (sum);
}

static _Complex double deriv_poly(const double *coefs, int len, _Complex double pos, int order) {
  /*
   * 1 x x^2 x^3 ...
   * 0 1 2x 3x^2 ...
   * 0 0 2 2*3x ...
   */
  _Complex double sum = 0, xp = 1;
  if(order == 1) {
    for(int i = 1; i < len; i++) {
      sum += i * coefs[i] * xp;
      xp *= pos;
    }
    return (sum);
  }
  
  for(int i = 0; i < order; i++) {
    xp *= order - i;
  }
  for(int i = order; i < len; i++) {
    sum += coefs[i] * xp;
    xp *= pos * (i + 1) / (i + 1 - order);
  }
  return (sum);
}


int poly_roots(double *coefs, int len, double *out, double conv) {
  if(len < 2) {
    // Can't do these.
    errno = EDOM;
    return (0);
  }

  // If the leading coefficient is zero, move to a different length.
  if(coefs[len - 1] == 0) {
    return (poly_roots(coefs, len - 1, out, conv));
  }
  // If the trailing coefficient is zero, the polynomial can be reduced.
  if(coefs[0] == 0) {
    out[0] = 0;
    return (1 + poly_roots(coefs + 1, len - 1, out + 1, conv));
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
    _Complex double C;
    if(del0 == 0 && del1 == 0) {
      out[0] = out[1] = out[2] = -b / (3 * a);
      return (3);
    } else if(del0 == 0) {
      C = cbrt(del1);
    } else {
      _Complex double R = del1 + csqrt(del1 * del1 - 4 * del0 * del0 * del0);
      C = cbrt(cabs(R / 2)) * cexp(I * carg(R/2) / 3);
    }
    _Complex double r0, r1, r2;
    _Complex double x0 = CMPLX(-0.5, sqrt(3) / 2),
      x1 = CMPLX(-0.5, -sqrt(3) / 2);
    r0 = -(b + C + del0 / C) / (3 * a);
    r1 = -(b + C * x0 + del0 / (C * x0)) / (3 * a);
    r2 = -(b + C * x1 + del0 / (C * x1)) / (3 * a);
    int pos = 0;
    if(fabs(cimag(r0)) < conv) {
      out[pos] = creal(r0);
      pos++;
    }
    if(fabs(cimag(r1)) < conv) {
      out[pos] = creal(r1);
      pos++;
    }
    if(fabs(cimag(r2)) < conv) {
      out[pos] = creal(r2);
      pos++;
    }
    
    return (pos);
  } else if(len == 5) { // Quartic
    // First, check if the polynomial can be written as a quadratic.
    if(coefs[1] == 0 && coefs[3] == 0) {
      double a = coefs[4], b = coefs[2], c = coefs[0];
      if(b * b - 4 * a * c < 0) {
	return (0);
      }
      _Complex double r0, r1, r2, r3;
      r0 = (-b + csqrt(b * b - 4 * a * c)) / (2 * a);
      r1 = -(-b + csqrt(b * b - 4 * a * c)) / (2 * a);
      r2 = (-b - csqrt(b * b - 4 * a * c)) / (2 * a);
      r3 = -(-b - csqrt(b * b - 4 * a * c)) / (2 * a);
      int pos = 0;

      if(fabs(cimag(r0)) < conv) {
	out[pos] = creal(r0);
	pos++;
      }
      if(fabs(cimag(r1)) < conv) {
	out[pos] = creal(r1);
	pos++;
      }
      if(fabs(cimag(r2)) < conv) {
	out[pos] = creal(r2);
	pos++;
      }
      if(fabs(cimag(r3)) < conv) {
	out[pos] = creal(r3);
	pos++;
      }
      return (pos);
    } else {
      double a = coefs[4], b = coefs[3], c = coefs[2],
	d = coefs[1], e = coefs[0];
      double offset = -b / 4 / a;
      double alp, bet, gam;
      alp = -3 * b * b / (8 * a * a) + c / a;
      bet = b * b * b / (8 * a * a * a) - b * c / (2 * a * a) + d / a;
      gam = -3 * b * b * b * b / (256 * a * a * a * a) +
	c * b * b / (16 * a * a * a) - b * d / (4 * a * a) + e / a;
      if(bet == 0) {
	double r0, r1;
	if(alp * alp - 4 * gam < 0) {
	  return (0);
	}
	r0 = (-alp + sqrt(alp * alp - 4 * gam)) / 2;
	r1 = (-alp - sqrt(alp * alp - 4 * gam)) / 2;
	int pos = 0;
	if(r0 >= 0) {
	  out[pos] = offset + sqrt(r0);
	  out[pos + 1] = offset - sqrt(r0);
	  pos += 2;
	}
	if(r1 >= 0) {
	  out[pos] = offset + sqrt(r1);
	  out[pos + 1] = offset - sqrt(r1);
	  pos += 2;
	}
	return (pos);
      }
      if(gam == 0) {
	double *cub = calloc(4, sizeof(double));
	cub[0] = bet;
	cub[1] = alp;
	cub[2] = 0;
	cub[3] = 1;
	out[0] = 0;
	int retval = 1 + poly_roots(cub, len - 1, out + 1, conv);
	free(cub);
	return (retval);
      }
      double p, q;
      _Complex double r, u, y, w;
      _Complex double r0, r1, r2, r3;
      p = -alp * alp / 12 - gam;
      q = -alp * alp * alp / 108 + alp * gam / 3 - bet * bet / 8;
      r = -q / 2 + csqrt(q * q / 4 + p * p * p / 27);
      u = cbrt(cabs(r)) * CMPLX(cos(carg(r) / 3), sin(carg(r) / 3));
      if(u == 0) {
	y = -5 * alp / 6 - cbrt(q);
      } else {
	y = -5 * alp / 6 + u - p / (3 * u);
      }
      w = csqrt(alp + 2 * y);
      r0 = -b / (4 * a) + (w + csqrt(-(3 * alp + 2 * y + 2 * bet / w))) / 2;
      r1 = -b / (4 * a) + (w - csqrt(-(3 * alp + 2 * y + 2 * bet / w))) / 2;
      r2 = -b / (4 * a) + (-w + csqrt(-(3 * alp + 2 * y - 2 * bet / w))) / 2;
      r3 = -b / (4 * a) + (-w - csqrt(-(3 * alp + 2 * y - 2 * bet / w))) / 2;
      int pos = 0;
      if(fabs(cimag(r0)) < conv) {
	out[pos] = creal(r0);
	pos++;
      }
      if(fabs(cimag(r1)) < conv) {
	out[pos] = creal(r1);
	pos++;
      }
      if(fabs(cimag(r2)) < conv) {
	out[pos] = creal(r2);
	pos++;
      }
      if(fabs(cimag(r3)) < conv) {
	out[pos] = creal(r3);
	pos++;
      }
      return (pos);
    }
  }
      
  
  // First, try to get roots from the companion matrix.
  double *companion = calloc((len - 1) * (len - 1), sizeof(double));
  double *eye = calloc((len - 1) * (len - 1), sizeof(double));
  for(int i = 0; i < len - 1; i++) {
    companion[i * (len - 1) + len - 2] = -coefs[i] / coefs[len - 1];
    eye[i * (len - 1) + i] = 1;
  }
  for(int i = 1; i < len - 1; i++) {
    companion[i * (len - 1) + i - 1] = 1;
  }
  double *revals = calloc(len - 1, sizeof(double)),
    *ievals = calloc(len - 1, sizeof(double)),
    *dens = calloc(len - 1, sizeof(double));
  // Use the results of this calculation as seeds for Newton's method.
  LAPACKE_dggev(LAPACK_COL_MAJOR, 'N', 'N', len - 1, companion,
		len - 1, eye, len - 1, revals, ievals, dens, NULL,
		1, NULL, 1);
  free(companion);
  free(eye);
  int pos = 0;
  for(int i = 0; i < len - 1; i++) {
    _Complex double x1 = CMPLX(revals[i] / dens[i], ievals[i] / dens[i]), x0;

    do {
      x0 = x1;
      x1 = x0 - eval_poly(coefs, len, x0) / deriv_poly(coefs, len, x0, 1);
    } while(cabs(x0 - x1) > conv);

    if(fabs(cimag(x1)) < sqrt(conv)) {
      out[pos] = creal(x1);
      pos++;
    }
  }
  free(revals);
  free(ievals);
  free(dens);
  return (pos);
}