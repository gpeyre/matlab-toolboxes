
/*------------------------------------------------------------------------------*/
/** 
 *  \file   perform_haar_transform.h
 *  \brief  Definition of class \c perform_haar_transform
 *  \author Gabriel Peyré
 *  \date   5-21-2005
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _PERFORM_HAAR_TRANSFORM_H_
#define _PERFORM_HAAR_TRANSFORM_H_

#include <math.h>

void perform_haar_transform_1d( const double* x, double* y, int n, int Jmin, int dir=1 )
{
    const double a = 1.0/sqrt(2.0);
    // copy to perform inplace computation
    memcpy( y,x,n*sizeof(double) );

    if( dir==1 )
    {
        int k = 1, t;
        double tmp;
        /* forward transform */
        while( k<=n-1 )
        {
            t = 0;
            while( t+k<=n-1 )
            {
                // perform a transform for y[t] and y[t+k]
                tmp = a*(y[t]-y[t+k]);
                y[t] = a*(y[t]+y[t+k]);
                y[t+k] = tmp;
                t += 2*k;
            }
            k *= 2;
        }
    }
    else
    {
        int l = log2(n-1);
        int k = 1 << l;		// 2^l
        int t;
        double tmp;
        /* backward transform */
        while( k>0 )
        {
            t = 0;
            while( t+k<=n-1 )
            {
                // perform a transform for y[t] and y[t+k]
                tmp = a*(y[t]-y[t+k]);
                y[t] = a*(y[t]+y[t+k]);
                y[t+k] = tmp;
                t += 2*k;
            }
            if( k==1 )
                return;
            k = k>>1;
        }
    }
}

#define X(i,j) x[i+(j)*n]
#define Y(i,j) y[i+(j)*n]
void perform_haar_transform_2d( const double* x, double* y, int n, int Jmin, int dir=1 )
{
    const double a = 1.0/sqrt(2.0);
    double tmp;
    int k = 1, tx, ty;
    // copy to perform inplace computation
    memcpy( y,x,n*n*sizeof(double) );

    if( dir==1 )
    {
        /* forward transform */
        while( k<=n-1 )
        {
            // transform on X
            for( ty=0; ty<n; ty+=k )
            {
                for( tx=0; tx+k<n; tx+=2*k )
                {
                    // perform a transform for y[t] and y[t+k]
                    tmp			= a*( Y(tx,ty) - Y(tx+k,ty) );
                    Y(tx,ty)	= a*( Y(tx,ty) + Y(tx+k,ty) );
                    Y(tx+k,ty)	= tmp;
                }
            }
            // transform on Y
            for( tx=0; tx<n; tx+=k )
            {
                for( ty=0; ty+k<n; ty+=2*k )
                {
                    // perform a transform for y[t] and y[t+k]
                    tmp			= a*( Y(tx,ty) - Y(tx,ty+k) );
                    Y(tx,ty)	= a*( Y(tx,ty) + Y(tx,ty+k) );
                    Y(tx,ty+k)	= tmp;
                }
            }
            k = k<<1;	// k*=2;
        }
    }
    else
    {
        int l = log2(n-1);
        k = 1 << l;		// 2^l
        /* backward transform */
        while( k>0 )
        {
            // transform on Y
            for( tx=0; tx<n; tx+=k )
            {
                for( ty=0; ty+k<n; ty+=2*k )
                {
                    // perform a transform for y[t] and y[t+k]
                    tmp			= a*( Y(tx,ty) - Y(tx,ty+k) );
                    Y(tx,ty)	= a*( Y(tx,ty) + Y(tx,ty+k) );
                    Y(tx,ty+k)	= tmp;
                }
            }
            // transform on X
            for( ty=0; ty<n; ty+=k )
            {
                for( tx=0; tx+k<n; tx+=2*k )
                {
                    // perform a transform for y[t] and y[t+k]
                    tmp			= a*( Y(tx,ty) - Y(tx+k,ty) );
                    Y(tx,ty)	= a*( Y(tx,ty) + Y(tx+k,ty) );
                    Y(tx+k,ty)	= tmp;
                }
            }
            if( k==1 )
                return;
            k = k>>1;
        }
    }
}


#endif // _PERFORM_HAAR_TRANSFORM_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
