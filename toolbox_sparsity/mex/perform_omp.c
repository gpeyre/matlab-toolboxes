#include <stdio.h>
#include <stdlib.h>

#include "invert.h"
#include "matrix_vector.h"

#define MAX_TOL 0.001

vector_t perform_omp(matrix_t a, vector_t b, int k) {
  int i, i_max, j, ii, j_max;
  double e0, s, m;
  vector_t x_sol, indices, residual, pinv_times_b;
  matrix_t pinv_a_extr, a_extr;

  indices.h=0;
  indices.vector=(double *)malloc(a.w*sizeof(double));
  if (indices.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<indices.h;i++) indices.vector[i]=0;
  
  residual.h=b.h;
  residual.vector=(double *)malloc(residual.h*sizeof(double));
  if (residual.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  x_sol.h=a.w;
  x_sol.vector=(double *)malloc(x_sol.h*sizeof(double));
  if (x_sol.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }

  a_extr.w=a.w;
  a_extr.h=a.h;
  a_extr.matrix=(double *)malloc(a_extr.w*a_extr.h*sizeof(double));
  if (a_extr.matrix==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }

  pinv_times_b.h=a.w;
  pinv_times_b.vector=(double *)malloc(pinv_times_b.h*sizeof(double));
  if (pinv_times_b.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }

  for (i=0;i<residual.h;i++) {
    residual.vector[i]=b.vector[i];
  }

  e0=norm_l2(residual);

  i_max=a.h;
  if ((k<i_max)&&(k>0)) {
    i_max=k;
  }

  for (i=0;i<i_max;i++) {
    
    /* on cherche l'indice du max de A'*r */
    m=0;
    j_max=0;
    for (ii=0;ii<a.w;ii++) {
      s=0;
      for (j=0;j<a.h;j++) {
	s+=a.matrix[j+ii*a.h]*residual.vector[j];
      }

      if (s>m) {
	m=s;
	j_max=ii;
      }
      if (-s>m) {
	m=-s;
	j_max=ii;
      }
    }

    indices.vector[i]=j_max;
    indices.h++;

    /* calcul de r =  x - a(:,indices(1:j)) *( pinv(a(:,indices(1:j))) * b ) */
    for (ii=0;ii<a.h;ii++) {
      a_extr.matrix[ii+i*a.h]=a.matrix[ii+j_max*a.h];
    }
    a_extr.w=i+1;

    pinv_a_extr=pseudo_invert(a_extr);

    for (ii=0;ii<pinv_a_extr.h;ii++) {
      pinv_times_b.vector[ii]=0;
      for (j=0;j<pinv_a_extr.w;j++) {
	pinv_times_b.vector[ii]+=pinv_a_extr.matrix[ii+pinv_a_extr.h*j]*b.vector[j];
      }
    }
    pinv_times_b.h=pinv_a_extr.h;

    for (ii=0;ii<residual.h;ii++) {
      residual.vector[ii]=b.vector[ii];
      for (j=0;j<a_extr.w;j++) {
	residual.vector[ii]-=a_extr.matrix[ii+a_extr.h*j]*pinv_times_b.vector[j];
      }
    }

    if ((norm_l2(residual)/e0)<MAX_TOL) {
      /* sortie de la boucle */
      i=a.w;
    }

    del_matrix(pinv_a_extr);
  }

  for (j=0;j<x_sol.h;j++) {
    x_sol.vector[j]=0;
  }
  for (j=0;j<indices.h;j++) {
    x_sol.vector[(int)(indices.vector[j])]=pinv_times_b.vector[j];
  }

  del_vector(indices);
  del_vector(residual);
  del_vector(pinv_times_b);
  del_matrix(a_extr);
  
  return x_sol;
}
