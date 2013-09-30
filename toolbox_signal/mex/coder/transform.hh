/*---------------------------------------------------------------------------*/
// Baseline Wavelet Transform Coder Construction Kit
//
// Geoff Davis
// gdavis@cs.dartmouth.edu
// http://www.cs.dartmouth.edu/~gdavis
//
// Copyright 1996 Geoff Davis 9/11/96
//
// Permission is granted to use this software for research purposes as
// long as this notice stays attached to this software.
//
/*---------------------------------------------------------------------------*/
#ifndef _TRANSFORM_
#define _TRANSFORM_
#include "wavelet.hh"
#include "image.hh"
/*---------------------------------------------------------------------------*/

class WaveletTransform 
{
public:
  WaveletTransform  (Wavelet *wavelet, Image *image,
		     int nsteps, int symmetric = -1);
  WaveletTransform  (Wavelet *wavelet, int hsize, int vsize, int
		     nsteps = 0, int symmetric = -1);

  // copy constructor
  WaveletTransform  (const WaveletTransform &WaveletTransform);
  ~WaveletTransform ();
  
  WaveletTransform& operator= (const WaveletTransform &WaveletTransform);
  
  void transform (Image *image, Wavelet *wavelet = NULL, 
		  int nsteps = 0, int symmetric = -1);
  void invert    (Image *invertedImage);

  void mallatToLinear (Real *mallat);
  void linearToMallat (Real *mallat);
  
  int hsize, vsize;     // size of transformed image
  Real *value;          // array of transform coefficients
  Wavelet *wavelet;     // wavelet used for transform
  int nsteps;           // # of transform iterations to perform
  int symmetric;        // TRUE for symmetric extensions, FALSE for periodic

  Real *subband (int n) { return subbandPtr[n]; };
  Real *subband (int scale, int orientation)
                        { return subbandPtr[3*(scale-1)+orientation+1]; };
  
  // pointers to start of subbands
  int nSubbands;     // # of subbands in transform
  int *subbandSize;    // # of coeffs in each subband
  int *subbandHsize;   // horizontal size of subband
  int *subbandVsize;   // vertical size of subband
  Real **subbandPtr; // pointer to start of each subband

protected:
  void init ();
  void freeAll ();
};


#endif
/*---------------------------------------------------------------------------*/
