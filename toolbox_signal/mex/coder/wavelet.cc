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
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "global.hh"
#include "image.hh"
#include "wavelet.hh"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Wavelet::Wavelet (FilterSet *filterset)
{
  analysisLow = filterset->analysisLow;
  analysisHigh = filterset->analysisHigh;
  synthesisLow = filterset->synthesisLow;
  synthesisHigh = filterset->synthesisHigh;
  symmetric = filterset->symmetric;

  // amount of space to leave for padding vectors for symmetric extensions
  npad = max(analysisLow->size, analysisHigh->size);
}

/*---------------------------------------------------------------------------*/

Wavelet::~Wavelet ()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Wavelet::transform1d (Real *input, Real *output, int size,
			    int nsteps, int sym_ext)
{
   int i;
   int currentIndex = 0;
   Real *data[2];
   int lowSize = size, highSize;

   // If form of extension unspecified, default to symmetric
   // extensions for symmetrical filters and periodic extensions for
   // asymmetrical filters
   if (sym_ext == -1)
     sym_ext = symmetric;

   // data[0] and data[1] are padded with npad entries on each end
   data [0] = new Real [2*npad+size];
   data [1] = new Real [2*npad+size];

   for (i = 0; i < size; i++)
     data[currentIndex][npad+i] = input[i];

   while (nsteps--)  {
     if (lowSize <= 2 && symmetric == 1) {
       warning ("Reduce # of transform steps or increase signal size");
       warning ("  or switch to periodic extension");
       error ("Low pass subband is too small");
     }

     // Transform
     printf ("transforming, size = %d\n", lowSize);
     transform_step (data[currentIndex], data[1-currentIndex], 
		     lowSize, sym_ext);
     
     highSize = lowSize/2;
     lowSize = (lowSize+1)/2;
     
     // Copy high-pass data to output signal
     copy (data[1-currentIndex] + npad + lowSize, output +
	   lowSize, highSize);

     for (i = 0; i < lowSize+highSize; i++)
       printf ("%5.2f ", data[1-currentIndex][npad+i]);
     printf ("\n\n");
     
     // Now pass low-pass data (first 1/2 of signal) back to
     // transform routine
     currentIndex = 1 - currentIndex;
   }
   
   // Copy low-pass data to output signal
   copy (data[currentIndex] + npad, output, lowSize);
   
   delete [] data [1];
   delete [] data [0];
}

/*---------------------------------------------------------------------------*/

void Wavelet::invert1d (Real *input, Real *output, int size, 
			 int nsteps, int sym_ext)
{
   int i;
   int currentIndex = 0;
   Real *data[2];

   // If form of extension unspecified, default to symmetric
   // extensions for symmetrical filters and periodic extensions for
   // asymmetrical filters
   if (sym_ext == -1)
     sym_ext = symmetric;

   int *lowSize = new int [nsteps];
   int *highSize = new int [nsteps];

   lowSize[0] = (size+1)/2;
   highSize[0] = size/2;

   for (i = 1; i < nsteps; i++) {
     lowSize[i] = (lowSize[i-1]+1)/2;
     highSize[i] = lowSize[i-1]/2;
   }

   data [0] = new Real [2*npad+size];
   data [1] = new Real [2*npad+size];

   copy (input, data[currentIndex]+npad, lowSize[nsteps-1]);

   while (nsteps--)  {
     
     // grab the next high-pass component
     copy (input + lowSize[nsteps], 
	   data[currentIndex]+npad+lowSize[nsteps], highSize[nsteps]);
     
     // Combine low-pass data (first 1/2^n of signal) with high-pass
     // data (next 1/2^n of signal) to get higher resolution low-pass data
     invert_step (data[currentIndex], data[1-currentIndex], 
		  lowSize[nsteps]+highSize[nsteps], sym_ext);
     
     // Now pass low-pass data (first 1/2 of signal) back to
     // transform routine
     currentIndex = 1 - currentIndex;
   }

   // Copy inverted signal to output signal
   copy (data[currentIndex]+npad, output, size);

   delete [] highSize;
   delete [] lowSize;

   delete [] data [1];
   delete [] data [0];
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Wavelet::transform2d (Real *input, Real *output, int hsize, int vsize,
			    int nsteps, int sym_ext)
{
   int j;
   int hLowSize = hsize, hHighSize;
   int vLowSize = vsize, vHighSize;

   // If form of extension unspecified, default to symmetric
   // extensions for symmetrical filters and periodic extensions for
   // asymmetrical filters
   if (sym_ext == -1)
     sym_ext = symmetric;

   Real *temp_in = new Real [2*npad+max(hsize,vsize)];
   Real *temp_out = new Real [2*npad+max(hsize,vsize)];

   copy (input, output, hsize*vsize);

   while (nsteps--)  {
      if ((hLowSize <= 2 || vLowSize <= 2) && sym_ext == 1) {
	warning ("Reduce # of transform steps or increase signal size");
	warning ("  or switch to periodic extension");
	error ("Low pass subband is too small");
      }

      // Do a convolution on the low pass portion of each row
      for (j = 0; j < vLowSize; j++)  {
	 // Copy row j to data array
	 copy (output+(j*hsize), temp_in+npad, hLowSize);
	 
	 // Convolve with low and high pass filters
	 transform_step (temp_in, temp_out, hLowSize, sym_ext);

	 // Copy back to image
	 copy (temp_out+npad, output+(j*hsize), hLowSize);
      }

      // Now do a convolution on the low pass portion of  each column
      for (j = 0; j < hLowSize; j++)  {
	 // Copy column j to data array
	 copy (output+j, hsize, temp_in+npad, vLowSize);
	 
	 // Convolve with low and high pass filters
	 transform_step (temp_in, temp_out, vLowSize, sym_ext);

	 // Copy back to image
	 copy (temp_out+npad, output+j, hsize, vLowSize);
      }

      // Now convolve low-pass portion again
      hHighSize = hLowSize/2;
      hLowSize = (hLowSize+1)/2;
      vHighSize = vLowSize/2;
      vLowSize = (vLowSize+1)/2;
   }

   delete [] temp_out;
   delete [] temp_in;
}

/*---------------------------------------------------------------------------*/

void Wavelet::invert2d (Real *input, Real *output, int hsize, int vsize,
			 int nsteps, int sym_ext)
{
   int i, j;

   // If form of extension unspecified, default to symmetric
   // extensions for symmetrical filters and periodic extensions for
   // asymmetrical filters
   if (sym_ext == -1)
     sym_ext = symmetric;

   int *hLowSize = new int [nsteps],
       *hHighSize = new int [nsteps];
   int *vLowSize = new int [nsteps],
       *vHighSize = new int [nsteps];

   hLowSize[0] = (hsize+1)/2;
   hHighSize[0] = hsize/2;
   vLowSize[0] = (vsize+1)/2;
   vHighSize[0] = vsize/2;

   for (i = 1; i < nsteps; i++) {
     hLowSize[i] = (hLowSize[i-1]+1)/2;
     hHighSize[i] = hLowSize[i-1]/2;
     vLowSize[i] = (vLowSize[i-1]+1)/2;
     vHighSize[i] = vLowSize[i-1]/2;
   }

   Real *temp_in = new Real [2*npad+max(hsize,vsize)];
   Real *temp_out = new Real [2*npad+max(hsize,vsize)];

   copy (input, output, hsize*vsize);

   while (nsteps--)  {
      // Do a reconstruction for each of the columns
      for (j = 0; j < hLowSize[nsteps]+hHighSize[nsteps]; j++)  {
	 // Copy column j to data array
	 copy (output+j, hsize, temp_in+npad, 
	       vLowSize[nsteps]+vHighSize[nsteps]);
	 
	 // Combine low-pass data (first 1/2^n of signal) with high-pass
	 // data (next 1/2^n of signal) to get higher resolution low-pass data
	 invert_step (temp_in, temp_out,
		      vLowSize[nsteps]+vHighSize[nsteps], sym_ext);

	 // Copy back to image
	 copy (temp_out+npad, output+j, hsize,
	       vLowSize[nsteps]+vHighSize[nsteps]);
      }

      // Now do a reconstruction pass for each row
      for (j = 0; j < vLowSize[nsteps]+vHighSize[nsteps]; j++)  {
	 // Copy row j to data array
	 copy (output + (j*hsize), temp_in+npad, 
	       hLowSize[nsteps]+hHighSize[nsteps]);

	 // Combine low-pass data (first 1/2^n of signal) with high-pass
	 // data (next 1/2^n of signal) to get higher resolution low-pass data
	 invert_step (temp_in, temp_out,
		      hLowSize[nsteps]+hHighSize[nsteps], sym_ext);
	 
	 // Copy back to image
	 copy (temp_out+npad, output + (j*hsize), 
	       hLowSize[nsteps]+hHighSize[nsteps]);
      }
   }

   delete [] hLowSize;
   delete [] hHighSize;
   delete [] vLowSize;
   delete [] vHighSize;

   delete [] temp_in;
   delete [] temp_out;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// input and output are padded with npad values at the beginning and
// at the end

void Wavelet::transform_step (Real *input, Real *output, int size, 
			      int sym_ext)
{
  int i, j;

  int lowSize = (size+1)/2;
  int left_ext, right_ext;

  if (analysisLow->size %2) {
    // odd filter length
    left_ext = right_ext = 1;
  } else {
    left_ext = right_ext = 2;
  }

  if (sym_ext)
    symmetric_extension (input, size, left_ext, right_ext, 1);
  else
    periodic_extension (input, size);
    
  //                      coarse  detail
  // xxxxxxxxxxxxxxxx --> HHHHHHHHGGGGGGGG
  for (i = 0; i < lowSize; i++)  {
    output[npad+i] = 0.0;
    for (j = 0; j < analysisLow->size; j++)  {
      output [npad+i] += 
	input[npad + 2*i + analysisLow->firstIndex + j] *
	analysisLow->coeff[j];
    }
  }
  
  for (i = lowSize; i < size; i++)  {
    output[npad+i] = 0.0;
    for (j = 0; j < analysisHigh->size; j++)  {
      output [npad+i] += 
	input[npad + 2*(i-lowSize) + analysisHigh->firstIndex + j] * 
	analysisHigh->coeff[j];
    }
  }
}

/*---------------------------------------------------------------------------*/

void Wavelet::invert_step (Real *input, Real *output, int size, int sym_ext)
{
   int i, j;
   int left_ext, right_ext, symmetry;
   // amount of low and high pass -- if odd # of values, extra will be
   //   low pass
   int lowSize = (size+1)/2, highSize = size/2;

   symmetry = 1;
   if (analysisLow->size % 2 == 0) {
     // even length filter -- do (2, X) extension
     left_ext = 2;
   } else {
     // odd length filter -- do (1, X) extension
     left_ext = 1;
   }

   if (size % 2 == 0) {
     // even length signal -- do (X, 2) extension
     right_ext = 2;
   } else {
     // odd length signal -- do (X, 1) extension
     right_ext = 1;
   }

   Real *temp = new Real [2*npad+lowSize];
   for (i = 0; i < lowSize; i++) {
     temp[npad+i] = input[npad+i];
   }

   if (sym_ext)
     symmetric_extension (temp, lowSize, left_ext, right_ext, symmetry);
   else
     periodic_extension (temp, lowSize);

   // coarse  detail
   // HHHHHHHHGGGGGGGG --> xxxxxxxxxxxxxxxx
   for (i = 0; i < 2*npad+size; i++)
     output[i] = 0.0;

   int firstIndex = synthesisLow->firstIndex;
   int lastIndex = synthesisLow->size - 1 + firstIndex;

   for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)  {
      for (j = 0; j < synthesisLow->size; j++)  {
	output[npad + 2*i + firstIndex + j] +=
	  temp[npad+i] * synthesisLow->coeff[j];
      }
   }

   left_ext = 2;

   if (analysisLow->size % 2 == 0) {
     // even length filters
     right_ext = (size % 2 == 0) ? 2 : 1;
     symmetry = -1;
   } else {
     // odd length filters
     right_ext = (size % 2 == 0) ? 1 : 2;
     symmetry = 1;
   }

   for (i = 0; i < highSize; i++) {
     temp[npad+i] = input[npad+lowSize+i];
   }
   if (sym_ext)
     symmetric_extension (temp, highSize, left_ext, right_ext,
			  symmetry);
   else
     periodic_extension (temp, highSize);


   firstIndex = synthesisHigh->firstIndex;
   lastIndex = synthesisHigh->size - 1 + firstIndex;

   for (i = -lastIndex/2; i <= (size-1-firstIndex)/2; i++)  {
      for (j = 0; j < synthesisHigh->size; j++)  {
	output[npad + 2*i + firstIndex + j] +=
	  temp[npad+i] * synthesisHigh->coeff[j];
      }
   }

   delete [] temp;
}

/*---------------------------------------------------------------------------*/
// Do symmetric extension of data using prescribed symmetries
//   Original values are in output[npad] through output[npad+size-1]
//   New values will be placed in output[0] through output[npad] and in
//      output[npad+size] through output[2*npad+size-1] (note: end values may
//      not be filled in)
//   left_ext = 1 -> extension at left bdry is   ... 3 2 1 | 0 1 2 3 ...
//   left_ext = 2 -> extension at left bdry is ... 3 2 1 0 | 0 1 2 3 ...
//   right_ext = 1 or 2 has similar effects at the right boundary
//
//   symmetry = 1  -> extend symmetrically
//   symmetry = -1 -> extend antisymmetrically

void Wavelet::symmetric_extension (Real *output, int size, int left_ext, int
				   right_ext, int symmetry)
{
  int i;
  int first = npad, last = npad + size-1;

  if (symmetry == -1) {
    if (left_ext == 1)
      output[--first] = 0;
    if (right_ext == 1)
      output[++last] = 0;
  }
  int originalFirst = first;
  int originalLast = last;
  int originalSize = originalLast-originalFirst+1;

  int period = 2 * (last - first + 1) - (left_ext == 1) - (right_ext == 1);

  if (left_ext == 2)
    output[--first] = symmetry*output[originalFirst];
  if (right_ext == 2)
    output[++last] = symmetry*output[originalLast];

  // extend left end
  int nextend = min (originalSize-2, first);
  for (i = 0; i < nextend; i++) {
    output[--first] = symmetry*output[originalFirst+1+i];
  }

  // should have full period now -- extend periodically
  while (first > 0) {
    first--;
    output[first] = output[first+period];
  }

  // extend right end
  nextend = min (originalSize-2, 2*npad+size-1 - last);
  for (i = 0; i < nextend; i++) {
    output[++last] = symmetry*output[originalLast-1-i];
  }

  // should have full period now -- extend periodically
  while (last < 2*npad+size-1) {
    last++;
    output[last] = output[last-period];
  }
}

/*---------------------------------------------------------------------------*/
// Do periodic extension of data using prescribed symmetries
//   Original values are in output[npad] through output[npad+size-1]
//   New values will be placed in output[0] through output[npad] and in
//      output[npad+size] through output[2*npad+size-1] (note: end values may
//      not be filled in)

void Wavelet::periodic_extension (Real *output, int size)
{
  int first = npad, last = npad + size-1;

  // extend left periodically
  while (first > 0) {
    first--;
    output[first] = output[first+size];
  }

  // extend right periodically
  while (last < 2*npad+size-1) {
    last++;
    output[last] = output[last-size];
  }
}

/*---------------------------------------------------------------------------*/
