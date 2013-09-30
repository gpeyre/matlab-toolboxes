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
#include <iostream.h>
#include <fstream.h>
#include "transform.hh"
#include "coeffset.hh"
#include "allocator.hh"
#include "quantizer.hh"
/*---------------------------------------------------------------------------*/
void compress       (Image *image, Wavelet *wavelet, int nStages, 
		     int capacity, Real p, Real *weight,
		     int paramPrecision, int budget, int nQuant, 
		     Real minStepSize, char *filename, int monolayer);
/*---------------------------------------------------------------------------*/

int main (int argc, char **argv)
{
  char *program = argv[0];

#ifdef PGM
  if (argc < 4)  {
    fprintf (stderr, 
	     "Usage: %s [image][output][ratio]\n",
	     program);
    fprintf (stderr, 
	     "image: image to be compressed (in PGM format)\n");
    fprintf (stderr, 
	     "output: name of compressed image\n");
    fprintf (stderr, 
	     "ratio: compression ratio\n");
    exit(0);
  }
  char *infile_name = argv[1];
  char *outfile_name = argv[2];
  Real ratio = atof(argv[3]);

  // Load the image to be coded
  Image *image = new Image (infile_name);

  int budget = (int)((Real)(image->hsize*image->vsize)/ratio);  // (assumes 8 bit pixels)
  printf ("Reading %d x %d image %s, writing %s\nCompression ratio %g:1\n", 
	  image->hsize, image->vsize, infile_name, outfile_name, ratio);

#else
  if (argc < 6)  {
    fprintf (stderr, 
	     "Usage: %s [image][width][height][output][ratio]\n",
	     program);
    fprintf (stderr, 
	     "image: image to be compressed (in RAW format)\n");
    fprintf (stderr, 
	     "width, height: width and height of image to be compressed\n");
    fprintf (stderr, 
	     "output: name of compressed image\n");
    fprintf (stderr, 
	     "ratio: compression ratio\n");
    exit(0);
  }
  char *infile_name = argv[1];
  int hsize = atoi(argv[2]);
  int vsize = atoi(argv[3]);
  char *outfile_name = argv[4];
  Real ratio = atof(argv[5]);
  int budget = (int)((Real)(hsize*vsize)/ratio);  // (assumes 8 bit pixels)
  printf ("Reading %d x %d image %s, writing %s\nCompression ratio %g:1\n", 
	  hsize, vsize, infile_name, outfile_name, ratio);

  // Load the image to be coded
  Image *image = new Image (infile_name, hsize, vsize);
#endif

  // Create a new wavelet from the 7/9 Antonini filterset
  Wavelet *wavelet = new Wavelet (&Antonini);

  // Coding parameters
  Real p = 2.0;             // exponent for L^p error metric
  int nStages = 5;          // # of stages in the wavelet transform
  int capacity = 512;       // capacity of histogram for arithmetic coder  
  int paramPrecision = 4;   // precision for stored quantizer parameters
  int nQuant = 10;          // # of quantizers to examine for allocation
  Real minStepSize = 0.05;  // smallest quantizer step to consider
  int monolayer = FALSE;    // TRUE for non-embedded uniform quantizer
                            // FALSE for multilayer quantizer

  Real *weight =            // perceptual weights for coefficient sets
    new Real[3*nStages+1];
  // for now give all sets equal weight
  for (int i = 0; i < 3*nStages+1; i++)
    weight[i] = 1.0;

  compress (image, wavelet, nStages, capacity, p, weight, 
	    paramPrecision, budget, nQuant, minStepSize, outfile_name, 
	    monolayer);

  delete [] weight;
  delete image;
  delete wavelet;
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Compress an image and save to a file
//
// Image *image          Image to be compressed
// Wavelet *wavelet      Wavelet to use for transform
// int nStages           # of stages to use in transform
// int capacity          Capacity of histograms for arithmetic coder
// Real p                Exponent for L^p error metric
// Real *weight          Perceptual weights for subbands
// int paramPrecision    Precision for storing quantizer parameters --
//                          precision = n means params are stored with
//                          accuracy 2^(-n)
// int budget            Total # of bytes for compressed image
// int nQuant            # of different quantizer resolutions to
//                          consider for each subband.  Step size for
//                          quantizer k is roughly (max coeff - min coeff)/2^k
// Real minStepSize      # minimum quantizer step size to consider --
//                          prevents too much effort from being
//                          expended on fine-scale subbands
// char *filename        name for compressed file
// int monolayer         TRUE for non-embedded quantizer, FALSE for embedded

/*---------------------------------------------------------------------------*/
void compress (Image *image, Wavelet *wavelet, int nStages, 
	       int capacity, Real p, Real *weight,
	       int paramPrecision, int budget, int nQuant, 
	       Real minStepSize, char *filename, int monolayer)
{
  int i;
  // Compute the wavelet transform of the given image
  WaveletTransform *transform = 
    new WaveletTransform (wavelet, image, nStages); 

  // For each subband allocate a CoeffSet, an error metric, an
  // EntropyCoder, and a Quantizer
  int nSets = transform->nSubbands;
  CoeffSet **coeff = new CoeffSet* [nSets];
  ErrorMetric **err = new ErrorMetric* [nSets];
  EntropyCoder **entropy = new EntropyCoder* [nSets];
  Quantizer **quant = new Quantizer* [nSets];

  for (i = 0; i < nSets; i++) {
    // Use an L^p error metric for each subband
    err[i] = new LpError (p);

    if (monolayer) {
      // Use uniform quantizer and single layer escape coder for each
      //   subband
      entropy[i] = new EscapeCoder  (capacity);
      // Assume all subbands have pdf's centered around 0 except the
      //   low pass subband 0
      quant[i]   = new UniformQuant ((MonoLayerCoder *)entropy[i],
				     paramPrecision, i != 0, err[i]);
    } else {
      // Use a layered quantizer with dead zones and a layered entropy
      //   coder for each subband
      entropy [i] = new LayerCoder (nQuant, i != 0, capacity);
      // Assume all subbands have pdf's centered around 0 except the
      //   low pass subband 0
      quant [i] = new LayerQuant ((MultiLayerCoder *)entropy[i],
				  paramPrecision, i != 0, nQuant, err[i]);
    }
    // Partition the wavelet transformed coefficients into subbands --
    //   each subband will have a different quantizer
    coeff[i]   = new CoeffSet (transform->subband(i),
			       transform->subbandSize[i], quant[i]);

    // For each subband determine the rate and distortion for each of
    //    the possible nQuant quantizers
    coeff[i]->getRateDist (nQuant, minStepSize);
  }
  
  Allocator *allocator = new Allocator ();
  // Use rate/distortion information for each subband to find bit
  //    allocation that minimizes total (weighted) distortion subject
  //    to a byte budget
  budget -= nSets * 4;  // subtract off approximate size of header info
  allocator->optimalAllocate (coeff, nSets, budget, TRUE, weight);
  printf ("Target rate = %d bytes\n", budget);
  // Display the resulting allocation
  allocator->print (coeff, nSets);

  // Open output file
  ofstream outfile (filename, ios::out | ios::trunc | ios::binary);
  if (!outfile) {
    error ("Unable to open file %s", filename);
  }
  // Create I/O interface object for arithmetic coder
  Encoder *encoder = new Encoder (outfile);

  // Write image size to output file
  encoder->writePositive (image->hsize);
  encoder->writePositive (image->vsize);
  //  printf ("hsize = %d, vsize = %d\n", image->hsize, image->vsize);

  for (i = 0; i < nSets; i++) {
    // Write quantizer parameters for each subband to file
    coeff[i]->writeHeader (encoder, allocator->precision[i]);
    /*
    if (monolayer) {
      printf ("coeff[%d]->qmin = %g  qmax = %g  qmean = %g  precision=%d\n", i,
    	      ((UniformQuant *)(coeff[i]->quant))->qmin, 
	      ((UniformQuant *)(coeff[i]->quant))->qmax,
	      ((UniformQuant *)(coeff[i]->quant))->qmean,
	      allocator->precision[i]);
    } else {
      printf ("coeff[%d]->qmin = %g  qmax = %g  qmean = %g precision=%d\n", i,
	      ((LayerQuant *)(coeff[i]->quant))->qmin, 
	      ((LayerQuant *)(coeff[i]->quant))->qmax,
	      ((LayerQuant *)(coeff[i]->quant))->qmean,
	      allocator->precision[i]);
    }
    */
  }
  for (i = 0; i < nSets; i++) {
    // Quantize and write entropy coded coefficients for each subband
    coeff[i]->encode (encoder, allocator->precision[i]);
  }

  // Flush bits from arithmetic coder and close file
  encoder->flush ();
  delete encoder;
  outfile.close ();

  // Clean up
  for (i = 0; i < nSets; i++) {
    delete err[i];
    delete entropy[i];
    delete quant[i];
    delete coeff[i];
  }
  delete [] err;
  delete [] entropy;
  delete [] quant;
  delete [] coeff;
  delete allocator;
  delete transform;
}

/*---------------------------------------------------------------------------*/
