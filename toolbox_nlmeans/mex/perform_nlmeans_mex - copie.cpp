/*=================================================================
% denoise - denoise an image.
%
%   [M1,Wx,Wy] = perform_nlmeans_vectorized(Ma,H,Ha,Vx,Vy,T,max_dist,do_median);
%
%	Ma is the image used to perform denoising.
%	H is a high dimensional representation (generaly patch wise) of the image to denoise.
%	Ha is a h.d. representation for Ma.
%	Vx is the x position of center of the search in Ma (same size as H).
%	Vy is the y position of center of the search in Ma (same size as H).
%	T is the variance of the gaussian used to compute the weights.
%	max_dist restricts the search size around position given by Vx,Vy.
%	
%	M1 is the denoised image
%	Wx is the new center x position for next seach (center of best fit)
%	Wy is the new center y position for next seach (center of best fit)
%   
%   Copyright (c) 2006 Gabriel Peyré
*=================================================================*/


#include <math.h>
#include "mex.h"
#include "config.h"
#include <stdio.h>
#include <string.h>

#define access(M,a,b,c) M[(a)+m*(b)+m*n*(c)]
#define accessa(Ma,a,b,c) Ma[(a)+ma*(b)+ma*na*(c)]

#define Ma_(a,b,c) accessa(Ma,a,b,c)
#define Ha_(a,b,c) accessa(Ha,a,b,c)
#define w_(a,b) accessa(w,a,b,0)

#define H_(a,b,c) access(H,a,b,c)
#define M1_(a,b,c) access(M1,a,b,c)
#define Vx_(a,b) access(Vx,a,b,0)
#define Vy_(a,b) access(Vy,a,b,0)
#define Wx_(a,b) access(Wx,a,b,0)
#define Wy_(a,b) access(Wy,a,b,0)
#define Cac_(a,b,c) access(Cac,a,b,c)

/* Global variables */
int n = -1;	// width of M
int m = -1;	// height of M
int na = -1;	// width of Ma
int ma = -1;	// height of Ma
double* M1 = NULL;	// output
double* Ma = NULL;	// exemplar image to transfer 
double* H = NULL;	// vectorized patches to denoise
double* Ha = NULL;	// vectorized exemplar patches to transfer
double* Vx = NULL;
double* Vy = NULL;
double* Wx = NULL;
double* Wy = NULL;
double* w = NULL;
int k = -1;			// dimensionality of the vectorized patches H,Ha
int s = -1;			// number of color chanels of M,Ma
double T = 0.05f;	// width of the gaussian
int max_dist = 10;	// max distance for searching
bool do_median= false; // use L1 fit
bool do_patchwise = false;
bool use_lun= false; // use L1 or L2 for distances

inline void display_message(const char* mess, int v)
{
	char str[128];
	sprintf(str, mess, v);
	mexWarnMsgTxt(str);
}
inline void display_messagef(const char* mess, float v)
{
	char str[128];
	sprintf(str, mess, v);
	mexWarnMsgTxt(str);
}


inline void get_dimensions( const mxArray* arr, int& a, int& b, int& c )
{
	int nd = mxGetNumberOfDimensions(arr);
	if( nd==3 )
	{
		a = mxGetDimensions(arr)[0];
		b = mxGetDimensions(arr)[1];
		c = mxGetDimensions(arr)[2];
	}
	else if( nd==2 )
	{
		a = mxGetM(arr);
		b = mxGetN(arr);
		c = 1;
	}
	else
		mexErrMsgTxt("only 2D and 3D arrays supported."); 
}



/* the distance function used */
inline double weight_gaussian(double x)
{
	// x is the distance
	return exp(-(x*x)/(2*T*T))/T;
}
inline double weight_laplacian(double x)
{
	// x is the distance
	return exp(-(x)/(2*T))/T;
}
// #define weight_func weight_laplacian
#define weight_func weight_gaussian

/* 
	compute the SQUARE distance between two windows.
	(i,j) is the center for the pixel to denoise (in H)
	(i1,j1) is the center for the potential pixel (in Ha), 
	k is the dimensionality of the vectors
*/
inline
double dist_windows( int i, int j, int i1, int j1 )
{
	double dist = 0;
	for( int a=0; a<k; ++a )
	{
		double d = H_(i,j,a) - Ha_(i1,j1,a);
		if( use_lun ) 
			dist += GW_ABS(d);
		else
			dist += d*d;
	}
	if( use_lun )
		return dist/k;
	else
		return sqrt( dist/k );
}

// structure to hold result of qsort
struct pixel
{
	double v;
	int i;
	int j;
};

int pixel_cmp(const void *a1, const void *a2)
{
	const pixel* c1 = (const pixel*) a1;
	const pixel* c2 = (const pixel*) a2;
	if(c1->v < c2->v)
		return -1;
	else if(c1->v > c2->v)
		return 1;
  	return 0;
}

double compute_weights(int i, int j, int i_min,int i_max,int j_min,int j_max)
{		
	double dmin = 1e9; // value of minimum distance
	double w_sum = 0;	// sum of all weights
	for( int i1=i_min; i1<=i_max; ++i1 )	// pixels of Ma
	for( int j1=j_min; j1<=j_max; ++j1 )
	{
		double ww = dist_windows( i,j, i1,j1 );
		if( ww<dmin )
		{
			// update best fit
			Wx_(i,j) = i1; Wy_(i,j) = j1;
			dmin = ww;
		}
		ww = weight_func(ww);
		w_sum += ww;
		w_(i1,j1) = ww;
	}
	if( w_sum<1e-9 )
	{
		// too low weights : using best fit
		// display_messagef("w_sum=%.8f", w_sum);
		w_sum = 1;
		w_((int) Wx_(i,j),(int) Wy_(i,j)) = 1;
	}
	return w_sum;
}

void denoise()
{		
	int i,j;
	// allocate some space to do the sorting operation
	pixel* vals = NULL;
	if( do_median )
		vals = (pixel*) malloc( (2*max_dist+1)*(2*max_dist+1)*sizeof(pixel) );
	for( i=0; i<m; ++i ) // pixels of M
	for( j=0; j<n; ++j )
	{
		// center for the search
		int ic = Vx_(i,j);
		int jc = Vy_(i,j);
		/* compute the weight for each windows */
		int i_min = GW_MAX(0,ic-max_dist);
		int i_max = GW_MIN(ma-1,ic+max_dist);
		int j_min = GW_MAX(0,jc-max_dist);
		int j_max = GW_MIN(na-1,jc+max_dist);
//		display_message("max_dist=%d", max_dist);
		double w_sum = compute_weights(i,j, i_min,i_max,j_min,j_max);
		/* perform reconstruction */
		for( int a=0; a<s; ++a )
		{
			M1_(i,j,a) = 0;
			if( !do_median )
			{
				// traditional mean
				for( int i1=i_min; i1<=i_max; ++i1 )
				for( int j1=j_min; j1<=j_max; ++j1 )
				{
					M1_(i,j,a) = M1_(i,j,a) + w_(i1,j1)/w_sum*Ma_(i1,j1,a);
				}
			}
			else	// median
			{
				// create the list for ranking
				int count = -1;
				for( int i1=i_min; i1<=i_max; ++i1 )
				for( int j1=j_min; j1<=j_max; ++j1 )
				{
					// if( i1!=i || j1!=j )
					{
						count++;
						vals[count].v =  Ma_(i1,j1,a);
						vals[count].i = i1;
						vals[count].j = j1;	
					}
				}
				w_sum -= w_(i,j);
				// do the sorting
				qsort(vals, count+1, sizeof(pixel), pixel_cmp);
				// extract medial rank
				double wcum = 0; int count1 = -1;
				while( wcum<=w_sum/2 && count1<count )
				{
					count1++;
					wcum = wcum + w_(vals[count1].i,vals[count1].j);
					//display_messagef("v=%f",vals[count1].v);
				}
				//display_messagef("opt=%f",(double) count1);
				// the medial value is count1
				M1_(i,j,a) = vals[count1].v;
			}
		}
	}
	if( do_median )
		free(vals);
}

int wdist = 3; // width of the patches
double lambda = 0.5;
int niter = 1; 

void denoise_patchwise()
{		
	double* Cac = (double*) malloc( m*n*s*sizeof(double) );
	// clear the accumulation buffers
	memset(Cac,0,m*n*s*sizeof(double));
	for( int i=0; i<m; ++i ) // pixels of M
	for( int j=0; j<n; ++j )
	{
		// center for the search
		int ic = Vx_(i,j);
		int jc = Vy_(i,j);
		// compute explorating region around (ic,jc)
		int i_min = GW_MAX(0,ic-max_dist);
		int i_max = GW_MIN(ma-1,ic+max_dist);
		int j_min = GW_MAX(0,jc-max_dist);
		int j_max = GW_MIN(na-1,jc+max_dist);
		// compute weight
		double w_sum = compute_weights(i,j, i_min,i_max,j_min,j_max);
		for( int i1=i_min; i1<=i_max; ++i1 )
		for( int j1=j_min; j1<=j_max; ++j1 )
		{
			// all the correct points at distance wdist
			int ti_min = -wdist;
			ti_min = GW_MAX( GW_MAX( ti_min, -i ), -i1 );
			ti_min = GW_MIN( GW_MIN( ti_min, m-1-i ), ma-1-i1 );
			int ti_max = wdist;
			ti_max = GW_MAX( GW_MAX( ti_max, -i ), -i1 );
			ti_max = GW_MIN( GW_MIN( ti_max, m-1-i ), ma-1-i1 );
			int tj_min = -wdist;
			tj_min = GW_MAX( GW_MAX( tj_min, -j ), -j1 );
			tj_min = GW_MIN( GW_MIN( tj_min, n-1-j ), na-1-j1 );
			int tj_max = wdist;
			tj_max = GW_MAX( GW_MAX( tj_max, -j ), -j1 );
			tj_max = GW_MIN( GW_MIN( tj_max, n-1-j ), na-1-j1 );
			for( int a=0; a<s; ++a )
			for( int ti = ti_min; ti<=ti_max; ++ti )
			for( int tj = tj_min; tj<=tj_max; ++tj )
			{
				M1_(i+ti,j+tj,a) += w_(i1,j1)*Ma_(i1+ti,j1+tj,a);
				Cac_(i+ti,j+tj,a) += w_(i1,j1);
			}
		}
	}	
	// normalize the result
	for( int a=0; a<s; ++a )
	for( int i=0; i<m; ++i ) // pixels of M
	for( int j=0; j<n; ++j )
	{
		M1_(i,j,a) /= Cac_(i,j,a);
	}
	free(Cac);
}

void mexFunction(	int nlhs, mxArray *plhs[], 
					int nrhs, const mxArray*prhs[] ) 
{ 
	if( nrhs<4 ) 
		mexErrMsgTxt("4 input arguments required."); 
	if( nlhs!=3 )
		mexErrMsgTxt("3 output arguments required.");
		
	// -- input 1 : Ma --
	get_dimensions( prhs[0], ma,na,s );
	Ma = mxGetPr(prhs[0]);
	
	// -- input 2 : H --
	get_dimensions( prhs[1], m,n,k );
	H = mxGetPr(prhs[1]);
		
	// -- input 3 : Ha --
	int m1,n1,k1;
	get_dimensions( prhs[2], m1,n1,k1 );
	if( na!=n1 || ma!=m1 || k!=k1  )
		mexErrMsgTxt("Ha should be of same size as Ma."); 
	Ha = mxGetPr(prhs[2]);
		
	// -- input 4 : Vx (same size as M) --
	get_dimensions( prhs[3], m1,n1,k1 );
	if( n!=n1 || m!=m1 || k1!=1  )
		mexErrMsgTxt("Vx should be of same size as H."); 
	Vx = mxGetPr(prhs[3]);
	
	// -- input 5 : Vy (same size as M) --
	get_dimensions( prhs[4], m1,n1,k1 );
	if( n!=n1 || m!=m1 || k1!=1  )
		mexErrMsgTxt("Vx should be of same size as H."); 
	Vy = mxGetPr(prhs[4]);
	
	// -- input 6 : T --
	if( nrhs>=6 )  
		T = *mxGetPr(prhs[5]);
	// -- input 7 : max_dist --
	if( nrhs>=7 )  
		max_dist = (int) *mxGetPr(prhs[6]);
	// -- input 8 : do_median
	do_median = false;
	if( nrhs>=8 )  
		do_median = *mxGetPr(prhs[7])>0.5;
	// -- input 9 : do_patchwise
	do_patchwise = false;
	if( nrhs>=8 )  
		do_patchwise = *mxGetPr(prhs[8])>0.5;
	if( nlhs>8 ) 
		mexErrMsgTxt("Too many input arguments.");
		
	// use robust distance with median
//	use_lun = do_median;
	use_lun = false;

	
	// -- outpout 1 : M1 -- 
	int dims[3] = {m,n,s};
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL );	
	M1 = mxGetPr(plhs[0]);
	
	// -- outpout 2 : Vx -- 
	plhs[1] = mxCreateDoubleMatrix(m,n, mxREAL);	
	plhs[2] = mxCreateDoubleMatrix(m,n, mxREAL);	
	Wx = mxGetPr(plhs[1]);
	Wy = mxGetPr(plhs[2]);


	// matrix of the weights for each pixel (same size as Ma)
	w = (double*) malloc( ma*na*sizeof(double) );
	memset(w,0,ma*na*sizeof(double));
    /* Do the actual computations in a subroutine */
	if( !do_patchwise )
		denoise();
	else
		denoise_patchwise();
	free( w );
}


