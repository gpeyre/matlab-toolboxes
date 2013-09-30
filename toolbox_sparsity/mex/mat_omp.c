#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "config.h"
#include "matrix_vector.h"
#include "mex.h"
// #include "perform_omp.c"
extern vector_t perform_omp(matrix_t a, vector_t b, int k);

/* Global variables */
double *x=NULL;
int ah=-1;
int aw=-1;
int bh=-1;
int bw=-1;
int kh=-1;
int kw=-1;
double *k=NULL;
double *a=NULL;
double *b=NULL;

void test_omp() {
  int i;
  matrix_t amatrix;
  vector_t bvector, xvector;

  srand(time(0));

  if (ah>aw) {
    mexErrMsgTxt("a has more rows than columns.");
  }
  if (ah!=bh) {
    mexErrMsgTxt("a and b dimensions do not match.");
  }
  if (bw!=1) {
    mexErrMsgTxt("b must be a vector");
  }

  if ((kh!=1) || (kw!=1)) {
    mexErrMsgTxt("k must be an integer");
  }
  
  amatrix.h=ah;
  amatrix.w=aw;
  amatrix.matrix=a;
  
  bvector.h=bh;
  bvector.vector=b;
  
  xvector=perform_omp(amatrix, bvector, (int)k[0]);

  for (i=0;i<xvector.h;i++) {
    x[i]=xvector.vector[i];
  }
  
  free(xvector.vector);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if(nrhs<2)
    mexErrMsgTxt("2 input arguments required."); 
  if(nlhs!=1)
    mexErrMsgTxt("1 output arguments required.");
		
  /* -- input 1 : a -- */
  ah=mxGetDimensions(prhs[0])[0];
  aw=mxGetDimensions(prhs[0])[1];
  a=mxGetPr(prhs[0]);
	
  /* -- input 2 : b -- */
  bh=mxGetDimensions(prhs[1])[0];
  bw=mxGetDimensions(prhs[1])[1];
  b=mxGetPr(prhs[1]);

  /* -- input 3 : k -- */
  kh=mxGetDimensions(prhs[2])[0];
  kw=mxGetDimensions(prhs[2])[1];
  k=mxGetPr(prhs[2]);

  /* -- outpout 1 : x --  */
  plhs[0]=mxCreateDoubleMatrix(aw, 1, mxREAL);	
  x=mxGetPr(plhs[0]);

  test_omp();
}
