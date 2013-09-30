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
#include <iostream>
#include <fstream.h>
#include "transform.hh"
#include "coeffset.hh"
#include "allocator.hh"
#include "quantizer.hh"
/*---------------------------------------------------------------------------*/
void decompress     (Image **image, Wavelet *wavelet, int nStages, 
		     int capacity, int paramPrecision, char *filename, 
		     int nQuant, int monolayer);
/*---------------------------------------------------------------------------*/

int main (int argc, char **argv)
{
  char *program = argv[0];

  if (argc < 3)  {
    fprintf (stderr, "Usage: %s [encoded image][decoded image]\n",
	     program);
    exit(0);
  }

  char *infile_name = argv[1];
  char *outfile_name = argv[2];
  printf ("Reading compressed image %s, writing %s\n\n", 
	  infile_name, outfile_name);

  // Create a new wavelet from the 7/9 Antonini filterset
  Wavelet *wavelet = new Wavelet (&Antonini);

  // Coding parameters
  int nStages = 5;          // # of stages in the wavelet transform
  int capacity = 512;       // capacity of histogram for arithmetic coder
  int paramPrecision = 4;   // precision for stored quantizer parameters
  int nQuant = 10;          // # of quantizers to examine for allocation
  int monolayer = FALSE;    // TRUE for non-embedded uniform quantizer
                            // FALSE for multilayer quantizer

  Image *reconstruct;
  decompress (&reconstruct, wavelet, nStages, capacity, paramPrecision,
	      infile_name, nQuant, monolayer);
#ifdef PGM
  reconstruct->savePGM (outfile_name);
#else
  reconstruct->saveRaw (outfile_name);
#endif
  delete reconstruct;
  delete wavelet;

  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Decompress an image
//
// Image **image         Place for decompressed image to be stored
// Wavelet *wavelet      Wavelet to use for transform
// int nStages           # of stages to use in transform
// int capacity          Capacity of histograms for arithmetic coder
// int paramPrecision    Precision for storing quantizer parameters --
//                          precision = n means params are stored with
//                          accuracy 2^(-n)
// int maxQuant          Maximum # of different quantizer resolutions
//                          used on any subband.  
// char *filename        Name of compressed file
// int monolayer         TRUE for non-embedded quantizer, FALSE for embedded

/*---------------------------------------------------------------------------*/
void decompress (Image **image, Wavelet *wavelet, int nStages, 
		 int capacity, int paramPrecision, char *filename, 
		 int maxQuant, int monolayer)
{
  int i;

  // open compressed image file
  ifstream infile (filename, ios::in | ios::nocreate | ios::binary);
  if (!infile) {
    error ("Unable to open file %s", filename);
  }

  // Create I/O interface object for arithmetic decoder
  Decoder *decoder = new Decoder (infile);

  // Read image dimensions from file
  int hsize = decoder->readPositive ();
  int vsize = decoder->readPositive ();
  //  printf ("hsize = %d, vsize = %d\n", hsize, vsize);

  // Create an empty transform of the appropriate size -- fill it in
  //    as coefficients are decoded
  WaveletTransform *transform = 
    new WaveletTransform (wavelet, hsize, vsize, nStages);

  // For each subband allocate a CoeffSet, an EntropyCoder, and a
  //    Quantizer (don't need to know anything about errors here)
  int nSets = transform->nSubbands;
  CoeffSet **coeff = new CoeffSet* [nSets];
  EntropyCoder **entropy = new EntropyCoder* [nSets];
  Quantizer **quant = new Quantizer* [nSets];
  // Quantizer precision for each subband
  int *precision = new int [nSets];

  for (i = 0; i < nSets; i++) {
    if (monolayer) {
      // Use uniform quantizer and single layer escape coder for each
      //   subband
      entropy[i] = new EscapeCoder (capacity);
      // Assume all subbands have pdf's centered around 0 except the
      //   low pass subband 0
      quant[i]   = new UniformQuant 
	((MonoLayerCoder *)entropy[i], paramPrecision, i != 0);
    } else {
      // Use a layered quantizer with dead zones and a layered entropy
      //   coder for each subband
      entropy [i] = new LayerCoder (maxQuant, i != 0, capacity);
      // Assume all subbands have pdf's centered around 0 except the
      //   low pass subband 0
      quant [i] = new LayerQuant 
	((MultiLayerCoder *)entropy[i], paramPrecision, i != 0, maxQuant);
    }
    // Indicate that each set of coefficients to be read corresponds
    //    to a subband
    coeff[i]   = new CoeffSet (transform->subband(i),
			       transform->subbandSize[i], quant[i]);
  }

  for (i = 0; i < nSets; i++) {
    // Read quantizer parameters for each subband
    coeff[i]->readHeader (decoder, precision[i]);
    /*
    if (monolayer) {
      printf ("coeff[%d]->qmin = %g  qmax = %g  qmean = %g  precision=%d\n", i,
	      ((UniformQuant *)(coeff[i]->quant))->qmin, 
	      ((UniformQuant *)(coeff[i]->quant))->qmax,
	      ((UniformQuant *)(coeff[i]->quant))->qmean,
	      precision[i]);
    } else {
      printf ("coeff[%d]->qmin = %g  qmax = %g  qmean = %g precision=%d\n", i,
	      ((LayerQuant *)(coeff[i]->quant))->qmin, 
	      ((LayerQuant *)(coeff[i]->quant))->qmax,
	      ((LayerQuant *)(coeff[i]->quant))->qmean,
	      precision[i]);
    }
    */
  }
  for (i = 0; i < nSets; i++) {
    // Read, decode, and dequantize coefficients for each subband
    coeff[i]->decode (decoder, precision[i]);
  }

  // Close file
  delete decoder;
  infile.close ();

  // Clean up
  for (i = 0; i < nSets; i++) {
    delete entropy[i];
    delete quant[i];
    delete coeff[i];
  }
  delete [] entropy;
  delete [] quant;
  delete [] coeff;
  delete [] precision;

  // Allocate image and invert transform
  *image = new Image (hsize, vsize);
  transform->invert (*image);
  delete transform;
}
  
/*---------------------------------------------------------------------------*/
