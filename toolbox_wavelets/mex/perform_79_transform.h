/*------------------------------------------------------------------------------*/
/** 
*  \file   perform_79_transform.h
*  \author Gabriel Peyré
*  \date   5-21-2005
*/ 
/*------------------------------------------------------------------------------*/

#ifndef _PERFORM_79_TRANSFORM_H_
#define _PERFORM_79_TRANSFORM_H_

#include "perform_lifting_transform.h"

void perform_79_transform_1d( double* x, int n, int Jmin, int dir=1 )
{
    // perform the transform
    const int step_type[] = {0,1,0,1,2};
    const double alpha = -1.586134342;
    const double beta  = -0.05298011854;
    const double gamma = 0.8829110762;
    const double delta = 0.4435068522;
    const double zeta  = 1.149604398;
    const double step_param[] = {alpha,beta,gamma,delta,zeta};
    perform_lifting_transform_1d( x, n, 5, step_type, step_param, Jmin, dir );
}

void perform_79_transform_2d( double* x, int n, int Jmin, int dir=1 )
{
    // not yet implemented
}


#endif // _PERFORM_79_TRANSFORM_H_


///////////////////////////////////////////////////////////////////////////////
//  Copyright (c) Gabriel Peyré
///////////////////////////////////////////////////////////////////////////////
//                               END OF FILE                                 //
///////////////////////////////////////////////////////////////////////////////