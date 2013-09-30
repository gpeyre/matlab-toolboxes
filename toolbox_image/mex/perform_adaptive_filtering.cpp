/*=================================================================
% perform_adaptive_filtering - perform adaptive filtering
%
%   B = perform_adaptive_filtering(A,H,I);
%
%   A is an n1 x n2 input image.
%   H is a p1 x p2 x m matrix, each H(:,:,k) being a filter.
%       Note that p1 and p2 should be odd integers.
%       Note that during filter, H(:,:,k) is automatically normalized
%       to sum to 1.
%   I is a p1 x p2 matrix of integer in {1,...,m}. I(i,j) is the filter to 
%       use in pixel (i,j)
%   B is an n x p output image.
%
%       B(i,j) = sum_x A(i+x1,j+x1) H(x1,x2,I(i,j))
%   
%   Copyright (c) 2006 Gabriel Peyré
*=================================================================*/

#include "mex.h"


#define		GW_ABS(a)       ((a) > 0 ? (a) : -(a))			//!<	Returns the absolute value a


#define A_(i,j) A[(i)+n1*(j)]
#define B_(i,j) B[(i)+n1*(j)]
#define I_(i,j) I[(i)+n1*(j)]
#define H_(i,j,k) H[(i)+p1*(j)+p1*p2*(k)]


void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<3 ) 
		mexErrMsgTxt("3 input arguments are required."); 
    if( nlhs!=1 ) 
        mexErrMsgTxt("1 output arguments are required."); 

	// first argument : input array
	int n1 = mxGetM(prhs[0]); 
	int n2 = mxGetN(prhs[0]);
    double* A = mxGetPr(prhs[0]);

    // secong argument : input filters
    int p1 = mxGetDimensions(prhs[1])[0];
    int p2 = mxGetDimensions(prhs[1])[1];
    int m = mxGetDimensions(prhs[1])[2];
    double* H = mxGetPr(prhs[1]);
    if( (p1%2)!=1 || (p2%2)!=1 )
        mexErrMsgTxt("Filters should be of odd size."); 

    // third argument : input index
    int a1 = mxGetM(prhs[2]); 
    int a2 = mxGetN(prhs[2]);
    double* I = mxGetPr(prhs[2]);
    if( a1!=n1 || a2!=n2 )
        mexErrMsgTxt("Array A and I should be of the same size."); 

    // output results
	plhs[0] = mxCreateDoubleMatrix(n1, n2, mxREAL);
	double* B = mxGetPr(plhs[0]);

    int q1 = (p1-1)/2;
    int q2 = (p2-1)/2;

    // perform filtering
    for( int j=0; j<n2; ++j )
    {
        for( int i=0; i<n1; ++i )
        {
            // number of the filter
            int h = (int) (I_(i,j)-1);
            B_(i,j) = 0;
            double v = 0;   // normalization value
            for( int s1=-q1; s1<=q1; ++s1 )
            for( int s2=-q2; s2<=q2; ++s2 )
            {
                if( (i+s1>=0) && (i+s1<n1) && (j+s2>=0) && (j+s2<n2)  )
                {
                    B_(i,j) += A_( i+s1,j+s2 ) * H_( s1+q1,s2+q2,h );
                    v += GW_ABS(H_( s1+q1,s2+q2,h ));
                }
            }
            B_(i,j) /= v;
        }
    }
}
