/*=================================================================
% perform_arithmetic_coding_fixed - perform adaptive arithmetic coding
%
%   code = perform_arithmetic_coding_fixed(x, dir, counts, n);
%
%   This is a wrapper to the code of Numerical Recipes.
%		www.nr.com
%
%	This coder works only for token with values >0.
%
%	dir=1 for coding, dir=-1 for decoding.
%	counts(i) is the number of token that have value i.
%	n is the size of the signal (you MUST provide it).
%
%	You should use the matlab driver for this Mex function
%	perform_arithmetid_coding together with options.coder_type=8.
%   
%   Copyright (c) 2006 Gabriel Peyré
*=================================================================*/

#include "mex.h"
#include "nr/nrutil.h"

typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned int uint;

#define MC 2024 // Maximum anticipated value of nchh in arcmak.
#define NWK 20 // Keep this value the same as in arcode, below.
typedef struct {
    ulong *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;

extern void arcode(ulong *ich, uchar **codep, ulong *lcode,
            ulong *lcd, int isign, arithcode *acode);
extern void arcmak(ulong nfreq[], ulong nchh, ulong nradd,
                   arithcode *acode);
#define MAXLINE 2048

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrieve arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 input arguments are required."); 
    if( nlhs!=1 ) 
        mexErrMsgTxt("1 output arguments are required."); 

	// first argument : input array
	int n1 = mxGetM(prhs[0]); 
	int p = mxGetN(prhs[0]);
    if( p>1 )
        mexErrMsgTxt("Works only for vector arrays."); 
    double* x1 = mxGetPr(prhs[0]);


	if( n1<=0 ) 
		mexErrMsgTxt("Works only for non empty arrays.."); 

    // second: direction
    int dir = (int) *mxGetPr(prhs[1]);

    // third: frequency count, add special token EOF
    ulong nchh = ((ulong) mxGetM(prhs[2])) + 1; 
    p = mxGetN(prhs[2]);
    if( p>1 )
        mexErrMsgTxt("Works only for vector arrays.");
    double* nfreq1 = mxGetPr(prhs[2]);
    ulong *nfreq = new ulong[nchh+10];
    for( uint i=0; i<nchh+10; ++i )
        nfreq[i] = 0;
    for( uint i=1; i<nchh; ++i )
        nfreq[i] = (ulong) nfreq1[i-1];
    nfreq[nchh] = 1; // special token occurence : 1

    ulong n;      // size of the signal
    ulong nc;     // size of the code
    ulong *x = NULL;     // signal 
    uchar *code = NULL; // code
    if( dir==1 )
    {
        n = n1; 
        nc = MAXLINE;
        x = new ulong[n];
        code=cvector(0,nc);
        // copy signal in int format
        for( uint i=0; i<n; ++i )
        {
            if( x1[i]<=0 )
                mexErrMsgTxt("Input vector should be >0.");
            x[i] = (unsigned int) (x1[i]-1);
        }
    }
    else
    {
        nc = n1; 
        // retrieve from last argument
        n = (int) *mxGetPr(prhs[3]);
        x = new ulong[n];
        code=cvector(0,nc);
        // copy code in uchar format
        for( uint i=0; i<nc; ++i )
            code[i+1] = (uchar) x1[i];
    }

    arithcode acode;
    const unsigned int nradd = 256; // radix for coding
    // initialize the coder 
    acode.ilob      = (ulong*) lvector( 1, NWK ); // Allocate space within acode.
    acode.iupb      = (ulong*) lvector( 1, NWK );
    acode.ncumfq    = (ulong*) lvector( 1, MC+2 );
    arcmak( nfreq, nchh, nradd, &acode );

    /* message initialization */
    ulong lc = 1;
    ulong zero = 0;
    arcode( &zero,&code,&n,&lc,0,&acode );
    if( dir==+1 )
    {
        /* here we arithmetically encode mess(1:n) */
        for( uint j=0;j<n;++j ) 
        {
            // code x[j]
            arcode(&x[j],&code,&nc,&lc,+1,&acode);
        }
        /* message termination, send nch */
        arcode(&nchh,&code,&nc,&lc,1,&acode);
    }
    else
    {
        /* here we decode the message, hopefully to get the original back */
        for( uint j=0; j<n; ++j ) 
        {
            ulong tmp;
            arcode(&tmp,&code,&nc,&lc,-1,&acode);
            if( tmp==nchh ) 
                break; // end of message
            else 
                x[j] = tmp+1;
        }
    }

    // output result
    if( dir==+1 )
    {
        plhs[0] = mxCreateNumericMatrix( lc-1, 1, mxDOUBLE_CLASS, mxREAL );
        double* y = mxGetPr( plhs[0] );
        for( uint i=0; i<=lc-2; ++i )
            y[i] = (double) code[i+1];
    }
    else
    {
        plhs[0] = mxCreateNumericMatrix( n, 1, mxDOUBLE_CLASS, mxREAL );
        double* y = mxGetPr( plhs[0] );
        for( uint i=0; i<n; ++i )
            y[i] = (double) x[i];
    }


    delete [] x;
    delete [] nfreq;
    free_cvector(code,0,MAXLINE);
    free_lvector(acode.ncumfq,1,MC+2);
    free_lvector(acode.iupb,1,NWK);
    free_lvector(acode.ilob,1,NWK);
}
