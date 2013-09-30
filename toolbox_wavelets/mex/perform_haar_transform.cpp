/*=================================================================
% perform_haar_transform - compute a wavelet haar transform
%
% y = perform_haar_transform(x,Jmin,dir);
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/
#include "mex.h"
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "perform_haar_transform.h"

/* Global variables */
int n, p;			// size
int dir = 1;
int Jmin = 0;
double* y = NULL;
double* x = NULL;

void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrieve arguments */
	if( nrhs<1 || nrhs>3 ) 
		mexErrMsgTxt("1-3 input arguments are required."); 
	if( nlhs!=1 ) 
		mexErrMsgTxt("1 output arguments are required."); 

	// first argument
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[0]);
	if( n==1 )
	{
		n = p;
		p = 1;
	}
	else if( p!=1 && n!=p )
		mexErrMsgTxt("2D Haar transform works only for square matrix."); 
	
	x = mxGetPr(prhs[0]);

    // second input : Jmin
    if( nrhs>=2 )
        Jmin = (int) *mxGetPr(prhs[1]);
    else
        Jmin = 0;

	// third input : dir
	if( nrhs>=3 )
		dir = (int) *mxGetPr(prhs[2]);
	else
		dir = 1;

    if( dir!=1 && dir!=-1 )
        mexErrMsgTxt("dir should be either +1 or -1."); 

	// first ouput : y
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL); 
	y = mxGetPr(plhs[0]);

	// perform the transform
	if( p==1 )
		perform_haar_transform_1d(x, y, n, Jmin, dir);
	else
		perform_haar_transform_2d(x, y, n, Jmin, dir);
}