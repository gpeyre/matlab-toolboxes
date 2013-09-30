/*=================================================================
% perform_arithmetic_coding - perform adaptive arithmetic coding.
%
%   This is a wrapper to the arithmetic coder of the 
%       Baseline Wavelet Transform Coder Construction Kit
%       http://www.cs.dartmouth.edu/~gdavis
%
%   y = perform_arithmetic_coding_escape(x,dir,known_size,known_bounds);
%   
%   Copyright (c) 2006 Gabriel Peyré
*=================================================================*/

#define HISTO_CAPACITY 1024*16

#include <iostream>
#include <fstream>
#include <math.h>
#include "coder/global.hh"
#include "coder/entropy.hh"
#include "mex.h"

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
        coding_bounds = true;  // the decoder will not know the bounds and pr

    // create encoderovide it    
    int capacity = HISTO_CAPACITY;       // capacity of histogram for arithmetic coder  
    EscapeCoder* entropy = new EscapeCoder(capacity);

    if( dir==1 )
    {
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
		if( subrange>=HISTO_CAPACITY )
			mexErrMsgTxt( "Too much symbols." );

        // Open output file
        ofstream outfile(filename, ios::out | ios::trunc | ios::binary);
        if( !outfile )
            mexErrMsgTxt("Cannot open file."); 
        // Create I/O interface object for arithmetic coder
        Encoder* encoder = new Encoder( outfile );

        // set the number of symbols
        entropy->setNSym( subrange );

        // code the size of the image and range
        if( coding_size )
            encoder->writePositive(n);
        else if( known_size!=n )
            mexErrMsgTxt("The provided size does not match the real size."); 
        if( coding_bounds )
            encoder->writeInt(a);
        else if( known_a>a )
            mexErrMsgTxt("The provided bound does not match the real size."); 
        if( coding_bounds )
            encoder->writeInt(b);
        else if( known_b<b )
            mexErrMsgTxt("The provided bound does not match the real size."); 

        // perform coding
        for( int i=0; i<n; i++ ) 
        {
            // rescale to positive integers
            int symbol = (int) (x[i]-a);
            assert( symbol>=0 && symbol<subrange );
            entropy->write(encoder, symbol, TRUE);
        }
        // finish encoding
        encoder->flush();
        delete encoder;
        outfile.close();

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

        // open compressed image file
        ifstream infile( filename, ios::in | ios::nocreate | ios::binary);
        if( !infile )
            error( "Unable to open file %s", filename );
        // Create I/O interface object for arithmetic decoder
        Decoder *decoder = new Decoder( infile );

        // retrieve size of the data
        int a, b;
        if( coding_size )
            n = decoder->readPositive();
        else
            n = known_size;
        if( coding_bounds )
            a = decoder->readInt();
        else
            a = known_a;
        if( coding_bounds )
            b = decoder->readInt();
        else
            b = known_b;
        int subrange = b-a+1;

        // set the number of symbols
        entropy->setNSym( subrange );

        // retrieve the data
        plhs[0] = mxCreateNumericMatrix( n, 1, mxDOUBLE_CLASS, mxREAL );
        double* y = mxGetPr( plhs[0] );

        for( int i=0; i<n; ++i ) 
        {
            int symbol = entropy->read( decoder, TRUE );
            assert( symbol<subrange && symbol >= 0);
            y[i] = a + symbol;
        }
    }
}
