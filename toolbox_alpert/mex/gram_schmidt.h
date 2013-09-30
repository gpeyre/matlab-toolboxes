#include "config.h"
#include <math.h>
#ifdef GW_DEBUG
	#include <crtdbg.h>
#endif

#define A_(i,j) A[i+(j)*n]
#define B_(i,j) B[i+(j)*n]

#define uk(i) B_(i,k)
#define uj(i) B_(i,j)
#define xk(i) A_(i,k)

INLINE
double dotp(const double* uu, double* vv, int n)
{
	int i;
	double d = 0;
	for( i=0; i<n; ++i )
		d += uu[i]*vv[i];
	return d;
}

INLINE
void gram_schmidt(const double* A, double *B, int n, int p)
{
	int i,j,k;
	double norm;
	if( p>n )
		p = n;
	// init to zeros
	for( k=0; k<p; ++k )
	{
		// search for u_k of the form u_k = x_k + sum_{j<i}{lambda_j*u_j}
		// init u_k to x_k
		for( i=0; i<n; ++i )
			uk(i) = xk(i);
		for( j=0; j<k; ++j )
		{
			// we have lambda_j=-<x_k,u_j>
			double lambda = -dotp( A + k*n, B + j*n, n );
			// add lambda*u_j to uk
			for( i=0; i<n; ++i )
				uk(i) += lambda*uj(i);
		}
		// renormalize
		norm = dotp( &B[k*n], &B[k*n], n );
		norm = sqrt(norm);
		/*
#if GW_DEBUG
		if( norm<1e-19 )
			GW_ASSERT(GW_False);		// the moment matrix is singular
#endif
		*/
		if( k>n )
			norm = 1.0;
		for( i=0; i<n; ++i )
			uk(i) /= norm;
	}
}
