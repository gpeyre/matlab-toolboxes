/*=================================================================
% perform_moment_transform - compute a wavelet haar transform
%
% w = perform_moment_transform(v,pos, monomials, dir, part)
%
%	monomials is of size 2x(2*k2)
%	pos is of size 2xn  where n is the number of points
%	part is of size nb_pts_max x P  where P is the number of partition.
%   
%   Copyright (c) 2004 Gabriel Peyré
*=================================================================*/

#include "perform_moment_transform.h"


/* Global variables */
double* v = NULL;
double* w = NULL;
double* pos = NULL;
double* monomials = NULL;
double** part = NULL;
int* G;
int* si;

