// mex driver function
#include "mex.h"
#include "perform_moment_transform.h"

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	int n,k,k2,P;			// size
	int dir;

	/* retrieve arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 or 5 input arguments are required."); 
	if( nlhs!=1 ) 
		mexErrMsgTxt("1 output arguments are required."); 

	// first argument : v
	n = mxGetM(prhs[0]);
	if( n==1 )
		n = mxGetN(prhs[0]);
	v = mxGetPr(prhs[0]);
	// second argument : pos
	int a  = mxGetM(prhs[1]);
	int b  = mxGetN(prhs[1]);
	if( a!=2 || b!=n )
		mexErrMsgTxt("pos must be of size 2xn."); 
	pos = mxGetPr(prhs[1]);
	// third argument : monomials
	a = mxGetM(prhs[2]);
	if( a!=2 )
		mexErrMsgTxt("monomials must be of size 2x(2*k2)."); 
	k = mxGetN(prhs[2]);
	k2 = k >> 1;
	if( k!=2*k2 )
		mexErrMsgTxt("monomials must be of even size."); 
	monomials = mxGetPr(prhs[2]);
	// fourth argument : dir
	dir = (int) *mxGetPr(prhs[3]);
	// fifth argument : part
	if( nrhs>=5 )
	{
		const mxArray* mx_part = prhs[4];
		if( !mxIsCell(mx_part) )
			mexErrMsgTxt("part must be a cell array."); 
		P = mxGetNumberOfElements(mx_part);
		part = new double*[P];
		// construct G, part, and si
		G = new int[P];
		si = new int[P+1]; si[0] = 0;
		for( int i=0; i<P; ++i )
		{
			mxArray* parti = mxGetCell(mx_part, i);
			if( parti==NULL )
				mexErrMsgTxt("error in part.");  
			G[i] = mxGetNumberOfElements(parti);
			si[i+1] = si[i]+G[i];
			part[i] = mxGetPr(parti);
		}
	}
	else
	{
		P = -1;
		part = NULL;
	}

	// first ouput : y
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); 
	w = mxGetPr(plhs[0]);

	// perform the transform
	perform_moment_transform(n, P, k, dir);

	if( nrhs>=5 )
	{
		// we are responsible for G and si
		GW_DELETE( G );
		GW_DELETE( si );
	}
}