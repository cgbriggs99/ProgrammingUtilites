/*
 * matprops.c
 *
 *  Created on: Feb 6, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

EXTRAMATH_FUNDEF(det,(const __TYPENAME__ *in, int dim)) {
	__TYPENAME__ determinant = 1;
	__TYPENAME__ *hold = calloc(dim * dim, sizeof(__TYPENAME__));

	memcpy(hold, in, dim * dim * sizeof(__TYPENAME__));

	for(int i = 0; i < dim; i++) {
		// Find the largest element in the column.
		__TYPENAME__ max = 0;
		int maxrow;
		for(int k = i; k < dim; k++) {
			if(__FNAMESRC_PREF__(abs)(hold[k * dim + i]) > __FNAMESRC_PREF__(abs)(max)) {
				max = hold[k * dim + i];
				maxrow = k;
			}
		}
		if(max == 0) {	// Matrix is not invertible, so zero.
			free(hold);
			return 0;
		}
		if(maxrow != i) {
			// Swap rows.
			for(int k = 0; k < dim; k++) {
				__TYPENAME__ temp = hold[i * dim + k];
				hold[i * dim + k] = hold[maxrow * dim + k];
				hold[maxrow * dim + k] = temp;
			}
		}
		// Update the determinant.
		determinant *= ((maxrow - i) % 2)? -1: 1;
		determinant /= max;
		// Eliminate. These steps do not affect the determinant.
		for(int j = 0; j < dim; j++) {
			if(j == i) {
				continue;
			}
			__TYPENAME__ scale = hold[j * dim + i] / max;
			hold[j * dim + i] = 0;
			for(int k = i + 1; k < dim; k++) {
				hold[j * dim + k] -= hold[i * dim + k] * scale;
			}
		}
		// Reduce the row. Already accounted for in the determinant.
		for(int j = 0; j < dim; j++) {
			hold[i * dim + j] /= max;
		}
	}
	free(hold);
	return determinant;
}

EXTRAMATH_FUNDEF(trace,(const __TYPENAME__ *in, int dim)) {
	__TYPENAME__ sum = 0;
	for(int i = 0; i < dim; i++) {
		sum += in[i * dim + i];
	}
	return sum;
}

EXTRAMATH_ARRFUNDEF(linearsolve,(__TYPENAME__ *coef_inout, __TYPENAME__ *right_inout, int eqs, int vars, int rightcols)) {
	int curr_row = 0, unsolved = 0;
	for(int i = 0; i < vars; i++) {
		if(curr_row >= eqs) {
			return 0;
		}
		// Find the largest element in the column.
		__TYPENAME__ max = 0;
		int maxrow;
		for(int k = curr_row; k < eqs; k++) {
			if(__FNAMESRC_PREF__(abs)(coef_inout[k * vars + i]) > __FNAMESRC_PREF__(abs)(max)) {
				max = coef_inout[k * vars + i];
				maxrow = k;
			}
		}
		if(max == 0) {	// System is not solvable, but continue anyway.
			unsolved++;
			continue;
		}
		if(maxrow != i) {
			// Swap rows.
			for(int k = 0; k < vars; k++) {
				__TYPENAME__ temp = coef_inout[curr_row * vars + k];
				coef_inout[curr_row * vars + k] = coef_inout[maxrow * vars + k];
				coef_inout[maxrow * vars + k] = temp;
			}
			for(int k = 0; k < rightcols; k++) {
				__TYPENAME__ temp = right_inout[curr_row * rightcols + k];
				right_inout[curr_row * rightcols + k] = right_inout[maxrow * rightcols + k];
				right_inout[maxrow * rightcols + k] = temp;
			}
		}
		// Eliminate.
		for(int j = 0; j < eqs; j++) {
			if(j == curr_row) {
				continue;
			}
			__TYPENAME__ scale = coef_inout[j * vars + i] / max;
			coef_inout[j * vars + i] = 0;
			for(int k = i + 1; k < vars; k++) {
				coef_inout[j * vars + k] -= coef_inout[curr_row * vars + k] * scale;
			}
			for(int k = 0; k < rightcols; k++) {
				right_inout[j * rightcols + k] -= right_inout[curr_row * rightcols + k] * scale;
			}
		}
		// Reduce the row.
		for(int j = 0; j < vars; j++) {
			coef_inout[curr_row * vars + j] /= max;
		}
		for(int j = 0; j < rightcols; j++) {
			right_inout[curr_row * rightcols + j] /= max;
		}
		// Go to the next row.
		curr_row++;
	}
	return unsolved;
}
EXTRAMATH_ARRFUNDEF(matinv,(__TYPENAME__ *inout, int dim)) {
	__TYPENAME__ *right = calloc(dim * dim, sizeof(__TYPENAME__));

	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			right[i * dim + j] = (i == j)? 1: 0;
		}
	}

	for(int i = 0; i < dim; i++) {
		// Find the largest element in the column.
		__TYPENAME__ max = 0;
		int maxrow;
		for(int k = i; k < dim; k++) {
			if(__FNAMESRC_PREF__(abs)(inout[k * dim + i]) > __FNAMESRC_PREF__(abs)(max)) {
				max = inout[k * dim + i];
				maxrow = k;
			}
		}
		if(max == 0) {	// Matrix is not invertible.
			free(right);
			return -1;
		}
		if(maxrow != i) {
			// Swap rows.
			for(int k = 0; k < dim; k++) {
				__TYPENAME__ temp = inout[i * dim + k];
				inout[i * dim + k] = inout[maxrow * dim + k];
				inout[maxrow * dim + k] = temp;
			}
			for(int k = 0; k < dim; k++) {
				__TYPENAME__ temp = right[i * dim + k];
				right[i * dim + k] = right[maxrow * dim + k];
				right[maxrow * dim + k] = temp;
			}
		}
		// Eliminate.
		for(int j = 0; j < dim; j++) {
			if(j == i) {
				continue;
			}
			__TYPENAME__ scale = inout[j * dim + i] / max;
			inout[j * dim + i] = 0;
			for(int k = i + 1; k < dim; k++) {
				inout[j * dim + k] -= inout[i * dim + k] * scale;
			}
			for(int k = 0; k < dim; k++) {
				right[j * dim + k] -= right[i * dim + k] * scale;
			}
		}
		// Reduce the row.
		for(int j = 0; j < dim; j++) {
			inout[i * dim + j] /= max;
		}
		for(int j = 0; j < dim; j++) {
			right[i * dim + j] /= max;
		}
	}
	memcpy(inout, right, dim * dim * sizeof(__TYPENAME__));
	free(right);
	return 0;
}


