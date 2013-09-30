#include "matrix_vector.h"

#define TWOPI (6.2831853071795864769252867665590057683943387987502)
#define RAND (rand())/((double) RAND_MAX)
#define RANDN (sqrt(-2.0*log(RAND))*cos(TWOPI*RAND))

double rand_normal(double mean, double sigma) {
  return (mean+sigma*RANDN);
}

vector_t random_sparse_vector(int h, int k) {
  int n, r;
  vector_t v;

  v.h=h;
  v.vector=(double *)malloc(h*sizeof(double));
  if (v.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (n=0;n<h;n++)
    v.vector[n]=0;
  n=0;
  while (n<k) {
    r=rand()%h;
    if (!v.vector[r]) {
      v.vector[r]=1/(double)k;
      n++;
    }   
  }
  return v;
}

matrix_t random_matrix(int h, int w, double mean, double sigma) {
  matrix_t m;
  int i, j;
  double s;

  m.w=w;
  m.h=h;
  m.matrix=(double *)malloc(h*w*sizeof(double));
  if (m.matrix==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<h;i++)
    for (j=0;j<w;j++)
      m.matrix[i+h*j]=rand_normal(mean,sigma);
  for (j=0;j<w;j++) {
    s=0;
    for (i=0;i<h;i++)
      s+=(m.matrix[i+h*j])*(m.matrix[i+h*j]);
    for (i=0;i<h;i++)
      m.matrix[i+h*j]=m.matrix[i+h*j]/sqrt(s);
  }
  return m;
}

matrix_t mult_matrix(matrix_t a, matrix_t b) {
  int i, j, k;
  matrix_t m;

  m.h=a.h;
  m.w=b.w;
  m.matrix=(double *)malloc(m.h*m.w*sizeof(double));
  if (m.matrix==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<m.h;i++) {
    for (j=0;j<m.w;j++) {
      m.matrix[i+m.h*j]=0;
      for (k=0;k<a.w;k++) {
	m.matrix[i+m.h*j]+=a.matrix[i+a.h*k]*b.matrix[k+b.h*j];
      }
    }
  }
  return m;
}

vector_t mult_matrix_vector(matrix_t a, vector_t x) {
  int i, j;
  vector_t ax;
  
  ax.h=a.h;
  ax.vector=(double *)malloc(ax.h*sizeof(double));
  if (ax.vector==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<ax.h;i++) {
    ax.vector[i]=0;
    for (j=0;j<a.w;j++) {
      ax.vector[i]+=a.matrix[i+a.h*j]*x.vector[j];
    }
  }
  return ax;
}

matrix_t transp(matrix_t a) {
  int i, j;
  matrix_t ta;

  ta.h=a.w;
  ta.w=a.h;
  ta.matrix=(double *)malloc(ta.h*ta.w*sizeof(double));
  if (ta.matrix==NULL) {
    fprintf(stderr,"malloc() error.\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<ta.h;i++) {
    for (j=0;j<ta.w;j++) {
      ta.matrix[i+ta.h*j]=a.matrix[j+a.h*i];
    }
  }
  return ta;
}

double norm_zero(vector_t v) {
  int i;
  double s;
  
  s=0;
  for (i=0;i<v.h;i++) {
    if (v.vector[i]!=0.0) {
      s+=1;
    }
  }
  return s;
}

double norm_l1(vector_t v) {
  int i;
  double s;

  s=0;
  for (i=0;i<v.h;i++) {
    if (v.vector[i]<0)
      s+=-v.vector[i];
    else
      s+=v.vector[i];
  }
  return s;
}

double norm_l2(vector_t v) {
  int i;
  double s;

  s=0;
  for (i=0;i<v.h;i++) {
      s+=v.vector[i]*v.vector[i];
  }
  return sqrt(s);
}

double scal(vector_t v1, vector_t v2) {
  int i;
  double s;
  s=0;
  for (i=0;i<v1.h;i++) {
    s+=v1.vector[i]*v2.vector[i];
  }
  return s;
}

void del_matrix(matrix_t a) {
  free(a.matrix);
}

void del_vector(vector_t x) {
  free(x.vector);
}
