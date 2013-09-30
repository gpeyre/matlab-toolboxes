/*=================================================================
% perform_79_transform - compute a wavelet biorthogonal 79 transform
%
% y = perform_79_transform(x,Jmin,dir);
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/
#include "mex.h"
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "perform_79_transform.h"

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
    
    if( p>1 )
        mexErrMsgTxt("7-9 transform not yet implemented for 2D matrix."); 

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

    memcpy( y, x, n*p*sizeof(double) );
    // perform the transform
    if( p==1 )
        perform_79_transform_1d( y, n, Jmin, dir );
    else
        perform_79_transform_2d( y, n, Jmin, dir );
}