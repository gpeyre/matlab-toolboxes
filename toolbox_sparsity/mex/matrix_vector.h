#ifndef MATRIX_VECTOR
#define MATRIX_VECTOR

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


typedef struct {
  double *vector;
  int h;
} vector_t;

typedef struct {
  double *matrix;
  int w, h;
} matrix_t;


vector_t random_sparse_vector(int h, int k);
matrix_t random_matrix(int h, int w, double mean, double var);
matrix_t mult_matrix(matrix_t a, matrix_t b);
vector_t mult_matrix_vector(matrix_t a, vector_t x);
matrix_t transp(matrix_t a);
double norm_zero(vector_t v);
double norm_l1(vector_t v);
double norm_l2(vector_t v);
double scal(vector_t v1, vector_t v2);
void del_matrix(matrix_t a);
void del_vector(vector_t x);

#endif
