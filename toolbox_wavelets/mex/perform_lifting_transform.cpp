/*=================================================================
% perform_lifting_transform - compute a wavelet biorthogonal 79 transform
%
% y = perform_lifting_transform( x, step_type, step_param, Jmin, dir );
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
    if( nrhs<3 || nrhs>5 ) 
        mexErrMsgTxt("3-5 input arguments are required."); 
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
        mexErrMsgTxt("lifting transform not yet implemented for 2D matrix."); 

    x = mxGetPr(prhs[0]);

    // input 2 : step_type
    int nbr_step = mxGetM(prhs[1]);
    if( nbr_step==1 )
        nbr_step = mxGetN(prhs[1]);
    double* step_type = (double*) mxGetPr(prhs[1]);

    // input 3 : step_param
    int tmp = mxGetM(prhs[2]);
    if( tmp==1 )
        tmp = mxGetN(prhs[2]);
    if( tmp!=nbr_step )
        mexErrMsgTxt("step_type and step_param must be of same size."); 
    double* step_param = (double*) mxGetPr(prhs[2]);

    // input 4 : Jmin
    if( nrhs>=4 )
        Jmin = (int) *mxGetPr(prhs[3]);
    else
        Jmin = 0;

    // input 5 : dir
    if( nrhs>=4 )
        dir = (int) *mxGetPr(prhs[4]);
    else
        dir = 1;

    if( dir!=1 && dir!=-1 )
        mexErrMsgTxt("dir should be either +1 or -1."); 

    // first ouput : y
    plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL); 
    y = mxGetPr(plhs[0]);

    memcpy( y, x, n*p*sizeof(double) );

    int* step_typei = new int[nbr_step];
    for( int i=0; i<nbr_step; ++i )
        step_typei[i] = (int) step_type[i];

    // perform the transform
    if( p==1 )
        perform_lifting_transform_1d( y, n, nbr_step, step_typei, step_param, Jmin, dir );
    else
        perform_lifting_transform_2d( y, n, p, nbr_step, step_typei, step_param, Jmin, dir );

    delete [] step_typei;
}