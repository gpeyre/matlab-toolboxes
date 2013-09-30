#ifndef INVERT
#define INVERT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix_vector.h"
#include "defs_and_types.h"

#define EPSILON 0.000001

double PYTHAG(double a, double b);
int dsvd(double **a, int m, int n, double *w, double **v);
matrix_t pseudo_invert(matrix_t a);
matrix_t gram_invert(matrix_t a);

#endif
