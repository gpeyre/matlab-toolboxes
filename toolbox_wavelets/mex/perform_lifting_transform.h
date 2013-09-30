/*------------------------------------------------------------------------------*/
/** 
 *  \file   perform_lifting_transform.h
 *  \author Gabriel Peyré
 *  \date   5-21-2005
 */ 
/*------------------------------------------------------------------------------*/

#ifndef _PERFORM_LIFTING_TRANSFORM_H_
#define _PERFORM_LIFTING_TRANSFORM_H_

#include <math.h>

#define STEP_PREDICT 0
#define STEP_UPDATE 1
#define STEP_SCALE 2


// coarse is x(0:h:n) and details is x(h/2:h:n)
void step_predict( double* x, int n, int h, double alpha )
{
    int h2 = h>>1;  // h/2
    double x2;
    for( int k=h2; k<n; k+=h )
    {
        // d1[k] = d[k] + alpha*(s[k]+s[k+1])
        if( k+h2<n )
            x2 = x[k+h2];
        else
            x2 = x[k-h2]; // use symmetry
        x[k] = x[k] + alpha * ( x[k-h2]+x2 );
    }
}

void step_update( double* x, int n, int h, double beta )
{
    int h2 = h>>1;  // h/2
    double x1, x2;
    for( int k=0; k<n; k+=h )
    {
        // s1[k] = s[k] + beta*(d[k]+d[k-1]).
        if( k-h2>0 )
            x1 = x[k-h2];
        else
            x1 = x[k+h2]; // use symmetry
        if( k+h2<n )
            x2 = x[k+h2];
        else
            x2 = x[k-h2]; // use symmetry
        x[k] = x[k] + beta * ( x1+x2 );
    }
}

void step_scale( double* x, int n, int h, double zeta )
{
    int h2 = h>>1;  // h/2
    for( int k=0; k<n; k+=h )
    {
        // s1=s*zeta; d1=d/zeta;
        x[k] *= zeta;
        if( k+h2<n )
            x[k+h2] /= zeta;
    }
}

void perform_lifting_transform_1d_fwd( double* x, int n, 
                                  int nbr_steps, const int* step_type, const double* step_param, int Jmin )
{
    int Jmax = log2( n )-1;
    int h = 1;
    for( int j=1; j<=Jmax-Jmin+1; ++j )
    {
        h *= 2;
        // coarse is x(0:h:n) and details is x(h/2:h:n)
        for( int s=0; s<nbr_steps; ++s ) 
        {
            switch( step_type[s] )
            {
            case STEP_PREDICT:
                step_predict( x, n, h, step_param[s] );
                break;
            case STEP_UPDATE:
                step_update( x, n, h, step_param[s] );
                break;
            case STEP_SCALE:
                step_scale( x, n, h, step_param[s] );
                break;
            }
        }
    }
}


void perform_lifting_transform_1d_bwd( double* x, int n, 
                                      int nbr_steps, const int* step_type, const double* step_param, int Jmin )
{   
    int Jmax = log2( n )-1;
    int h = 1<<(Jmax-Jmin+2);
    for( int j=Jmax-Jmin+1; j>0; --j )
    {
        h /= 2;
        // coarse is x(0:h:n) and details is x(h/2:h:n)
        for( int s=nbr_steps-1; s>=0; --s ) 
        {
            switch( step_type[s] )
            {
            case STEP_PREDICT:
                step_predict( x, n, h, -step_param[s] );
                break;
            case STEP_UPDATE:
                step_update( x, n, h, -step_param[s] );
                break;
            case STEP_SCALE:
                step_scale( x, n, h, 1.0/step_param[s] );
                break;
            }
        }
    }
}

void perform_lifting_transform_1d( double* x, int n, 
                                  int nbr_steps, const int* step_type, const double* step_param, 
                                  int Jmin = 0, int dir=1 )
{
    if( dir==1 )
        perform_lifting_transform_1d_fwd( x, n, nbr_steps, step_type, step_param, Jmin );
    else
        perform_lifting_transform_1d_bwd( x, n, nbr_steps, step_type, step_param, Jmin );
}

void perform_lifting_transform_2d( double* x, int n, int p,
                                  int nbr_steps, const int* step_type, const double* step_param, 
                                  int Jmin = 0, int dir=1 )
{

}

#endif // _PERFORM_LIFTING_TRANSFORM_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////
