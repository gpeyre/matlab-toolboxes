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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transform.hh"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

WaveletTransform::WaveletTransform  (Wavelet *wavelet, Image *image, 
				     int nsteps, int symmetric) : 
  wavelet (wavelet), nsteps(nsteps), symmetric(symmetric)
{
  value = NULL;

  if (image != NULL) {
    hsize = image->hsize;
    vsize = image->vsize;
    transform (image, wavelet, nsteps, symmetric);
  } else {
    hsize = vsize = 0;
  }
}

/*---------------------------------------------------------------------------*/

WaveletTransform::WaveletTransform (Wavelet *wavelet, int hsize, int
				    vsize, int nsteps, int symmetric):
                         hsize(hsize), vsize(vsize), wavelet
			 (wavelet), nsteps(nsteps), symmetric(symmetric)
{
  nsteps = 0;
  symmetric = -1;
  init ();
  for (int i = 0; i < hsize*vsize; i++)
    value[i] = 0;
}

/*---------------------------------------------------------------------------*/

WaveletTransform::WaveletTransform (const WaveletTransform &t)
{
  wavelet = t.wavelet;
  hsize = t.hsize;
  vsize = t.vsize;
  nsteps = t.nsteps;
  symmetric = t.symmetric;

  if (t.value == NULL) {
    value = NULL;
  } else {
    init ();
    for (int i = 0; i < hsize*vsize; i++)
      value[i] = t.value[i];
  }
}

/*---------------------------------------------------------------------------*/

WaveletTransform::~WaveletTransform ()
{
   freeAll ();
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::init () 
{
   int i;

   value = new Real [hsize*vsize];

   nSubbands = 3 * nsteps + 1;
   subbandSize = new int [nSubbands];
   subbandHsize = new int [nSubbands];
   subbandVsize = new int [nSubbands];
   subbandPtr = new Real* [nSubbands];

   int *lowHsize = new int [nsteps];
   int *lowVsize = new int [nsteps];
   int *highHsize = new int [nsteps];
   int *highVsize = new int [nsteps];

   lowHsize[nsteps-1] = (hsize+1)/2;
   lowVsize[nsteps-1] = (vsize+1)/2;
   highHsize[nsteps-1] = hsize/2; 
   highVsize[nsteps-1] = vsize/2;
  
   for (i = nsteps-2; i >= 0; i--) {
     lowHsize[i] = (lowHsize[i+1]+1)/2;
     lowVsize[i] = (lowVsize[i+1]+1)/2;
     highHsize[i] = lowHsize[i+1]/2;
     highVsize[i] = lowVsize[i+1]/2; 
   }

   subbandPtr[0] = value;
   subbandHsize[0] = lowHsize[0];
   subbandVsize[0] = lowVsize[0];
   subbandSize[0] = subbandHsize[0]*subbandVsize[0];

   for (i = 0; i < nsteps; i++) {
     subbandHsize[3*i+1] = highHsize[i];
     subbandVsize[3*i+1] = lowVsize[i];
     subbandHsize[3*i+2] = lowHsize[i];
     subbandVsize[3*i+2] = highVsize[i];
     subbandHsize[3*i+3] = highHsize[i];
     subbandVsize[3*i+3] = highVsize[i];
   }

   for (i = 1; i < nSubbands; i++) {
     subbandSize[i] = subbandHsize[i]*subbandVsize[i];
     subbandPtr[i] = subbandPtr[i-1] + subbandSize[i-1];
   }

   delete [] lowHsize;
   delete [] lowVsize;
   delete [] highHsize;
   delete [] highVsize;
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::freeAll ()
{
  if (value != NULL) {
    delete [] value;
    delete [] subbandSize;
    delete [] subbandHsize;
    delete [] subbandVsize;
    delete [] subbandPtr;
  }
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::transform (Image *image, Wavelet *newWavelet, 
				  int steps, int isSymmetric)
{
  // clear out old info and set up subband pointers
  freeAll ();
  hsize = image->hsize;
  vsize = image->vsize;
  wavelet = newWavelet;
  nsteps = steps;
  symmetric = isSymmetric;
  init ();

  Real *temp = new Real [hsize*vsize];
  wavelet->transform2d (image->value, temp, hsize, vsize, nsteps,
			 symmetric);

  // linearize data
  mallatToLinear (temp);
  delete [] temp;
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::invert (Image *invertedImage)
{
  Real *temp = new Real [hsize*vsize];

  // put data in Mallat format
  linearToMallat (temp);

  wavelet->invert2d (temp, invertedImage->value, hsize, vsize,
		      nsteps, symmetric);

  delete [] temp;
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::mallatToLinear (Real *mallat)
{
  int i, j, k;
  
  int *lowHsize = new int [nsteps];
  int *lowVsize = new int [nsteps];
  
  lowHsize[nsteps-1] = (hsize+1)/2;
  lowVsize[nsteps-1] = (vsize+1)/2;
  
  for (i = nsteps-2; i >= 0; i--) {
    lowHsize[i] = (lowHsize[i+1]+1)/2;
    lowVsize[i] = (lowVsize[i+1]+1)/2;
  }
  
  // move transformed image (in Mallat order) into linear array structure
  // special case for LL subband
  for (j = 0; j < subbandVsize[0]; j++)
    for (i = 0; i < subbandHsize[0]; i++)
      subbandPtr[0][j*subbandHsize[0]+i] = 
	mallat[j*hsize+i];
  
  for (k = 0; k < nsteps; k++) {
    for (j = 0; j < subbandVsize[k*3+1]; j++)
      for (i = 0; i < subbandHsize[k*3+1]; i++)
	subbandPtr[k*3+1][j*subbandHsize[k*3+1]+i] = 
	  mallat[j*hsize+(lowHsize[k]+i)];

    for (j = 0; j < subbandVsize[k*3+2]; j++)
      for (i = 0; i < subbandHsize[k*3+2]; i++)
	subbandPtr[k*3+2][j*subbandHsize[k*3+2]+i] = 
	  mallat[(lowVsize[k]+j)*hsize+i];

    for (j = 0; j < subbandVsize[k*3+3]; j++)
      for (i = 0; i < subbandHsize[k*3+3]; i++)
	subbandPtr[k*3+3][j*subbandHsize[k*3+3]+i] = 
	  mallat[(lowVsize[k]+j)*hsize+(lowHsize[k]+i)];
  }

  delete [] lowHsize;
  delete [] lowVsize;
}

/*---------------------------------------------------------------------------*/

void WaveletTransform::linearToMallat (Real *mallat)
{
  int i, j, k;

  int *lowHsize = new int [nsteps];
  int *lowVsize = new int [nsteps];
  
  lowHsize[nsteps-1] = (hsize+1)/2;
  lowVsize[nsteps-1] = (vsize+1)/2;
  
  for (i = nsteps-2; i >= 0; i--) {
    lowHsize[i] = (lowHsize[i+1]+1)/2;
    lowVsize[i] = (lowVsize[i+1]+1)/2;
  }
  
  // put linearized image in Mallat format
  // special case for LL subband
  for (j = 0; j < subbandVsize[0]; j++)
    for (i = 0; i < subbandHsize[0]; i++)
      mallat[j*hsize+i] = subbandPtr[0][j*subbandHsize[0]+i];

  for (k = 0; k < nsteps; k++) {
    for (j = 0; j < subbandVsize[k*3+1]; j++)
      for (i = 0; i < subbandHsize[k*3+1]; i++)
	mallat[j*hsize+(lowHsize[k]+i)] = 
	  subbandPtr[k*3+1][j*subbandHsize[k*3+1]+i];

    for (j = 0; j < subbandVsize[k*3+2]; j++)
      for (i = 0; i < subbandHsize[k*3+2]; i++)
	mallat[(lowVsize[k]+j)*hsize+i] = 
	  subbandPtr[k*3+2][j*subbandHsize[k*3+2]+i];

    for (j = 0; j < subbandVsize[k*3+3]; j++)
      for (i = 0; i < subbandHsize[k*3+3]; i++)
	mallat[(lowVsize[k]+j)*hsize+(lowHsize[k]+i)] = 
	  subbandPtr[k*3+3][j*subbandHsize[k*3+3]+i];
  }
  
  delete [] lowHsize;
  delete [] lowVsize;
}

/*---------------------------------------------------------------------------*/
