/*=================================================================
% perform_arithmetic_coding - perform adaptive arithmetic coding
%
%   y = perform_arithmetic_coding_mex(x,dir,known_size,known_bounds);
%
%   This is a wrapper to the code of Christophe Bernard.
%   
%   Copyright (c) 2005 Gabriel Peyré
*=================================================================*/

#include "mex.h"
#include "ac.h"

const char* filename = "tmp.tmp";

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 input arguments are required."); 
    if( nlhs!=1 ) 
        mexErrMsgTxt("1 output arguments are required."); 

	// first argument : input array
	int n = mxGetM(prhs[0]); 
	int p = mxGetN(prhs[0]);
    if( p>1 )
        mexErrMsgTxt("Works only for vector arrays."); 
    int dir = (int) *mxGetPr(prhs[1]);
    double* x = mxGetPr(prhs[0]);

    if( mxGetM(prhs[3])!=2 && mxGetN(prhs[3])!=2 )
        mexErrMsgTxt("known_bounds must be of size 2."); 
    int known_size = (int) *mxGetPr(prhs[2]);
    // lower and upper bounds on the int to code
    int known_a = (int) mxGetPr(prhs[3])[0];
    int known_b = (int) mxGetPr(prhs[3])[1];
    bool coding_size = true;
    if( known_size>0 )          // the decoder will know the size and provide it
        coding_size = false;
    bool coding_bounds = false;
    if( known_a==known_b && known_b==-1 )
        coding_bounds = true;  // the decoder will not know the bounds and provide it    

    if( dir==1 )
    {
        ac_encoder ace;
        ac_model acm;
        // compute range of the data
        int a = 1<<20, b = -(1<<20); 
        for( int i=0; i<n; ++i )
        {
            if( x[i]<a )
                a = (int) x[i];
            if( x[i]>b )
                b = (int) x[i];
        }
        int subrange = b-a+1;

        // init
        ac_encoder_init(&ace, filename);
        ac_model_init(&acm, subrange, NULL, 1);

        // code the size of the image and range
        if( coding_size )
            ac_encode_bytes(&ace, 4, &n);
        else if( known_size!=n )
            mexErrMsgTxt("The provided size does not match the real size."); 
        if( coding_bounds )
            ac_encode_bytes(&ace, 4, &a);
        else if( known_a>a )
            mexErrMsgTxt("The provided bound does not match the real size."); 
        if( coding_bounds )
            ac_encode_bytes(&ace, 4, &b);
        else if( known_b<b )
            mexErrMsgTxt("The provided bound does not match the real size."); 

        // perform coding
        for( int i=0; i<n; ++i ) 
            ac_encode_symbol(&ace, &acm, (int) (x[i]-a), 0);
        // end coding
        ac_encoder_done(&ace);


        // reopen output file
        FILE* fin = fopen( filename, "rb" );
        if( fin==NULL )
            mexErrMsgTxt("Cannot open file."); 
        int nbr_bytes = 0;

        // compute number of bytes
        while( getc(fin)!=EOF )
            nbr_bytes++;
        fclose(fin);

        // retrieve results in byte
        fin = fopen( filename, "rb" );
        if( fin==NULL )
            mexErrMsgTxt("Cannot open file."); 
        plhs[0] = mxCreateNumericMatrix( nbr_bytes, 1, mxDOUBLE_CLASS, mxREAL );
        double* y = mxGetPr( plhs[0] );
        for( int i=0; i<nbr_bytes; ++i )
        {
            y[i] = (double) getc(fin);
        }
        fclose(fin);
    }
    else
    {
        // write data to a file byte by byte
        FILE* fin = fopen( filename, "wb" );
        if( fin==NULL )
            mexErrMsgTxt("Cannot open file."); 
        for( int i=0; i<n; ++i )
        {
            char c = (char) x[i];
            fwrite( &c, sizeof(char), 1, fin );
        }
        fclose(fin);
        // initialize coder
        ac_decoder acd;
        ac_model acm;
        ac_decoder_init(&acd, filename);

        // retrieve size of the data
        int n, a, b;
        if( coding_size )
            ac_decode_bytes(&acd, 4, &n);
        else
            n = known_size;
        if( coding_bounds )
            ac_decode_bytes(&acd, 4, &a);
        else
            a = known_a;
        if( coding_bounds )
            ac_decode_bytes(&acd, 4, &b);
        else
            b = known_b;
        int subrange = b-a+1;
        ac_model_init(&acm, subrange, NULL, 1);

        // retrieve the data
        plhs[0] = mxCreateNumericMatrix( n, 1, mxDOUBLE_CLASS, mxREAL );
        double* y = mxGetPr( plhs[0] );
        for( int i=0; i<n; ++i ) 
            y[i] = ac_decode_symbol(&acd, &acm, 0)+a;
    }

}
