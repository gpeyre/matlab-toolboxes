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
#include <iostream.h>
#include <math.h>
#include "global.hh"
#include "entropy.hh"
#include "quantizer.hh"
#include "coder.hh"
#include "allocator.hh"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Allocator::Allocator ()
{
  precision = NULL;
}

/*---------------------------------------------------------------------------*/

Allocator::~Allocator ()
{
  if (precision != NULL)
    delete [] precision;
}

/*---------------------------------------------------------------------------*/

void Allocator::resetPrecision (int nSets)
{
  if (precision != NULL)
    delete [] precision;

  precision = new int [nSets];
  for (int i = 0; i < nSets; i++)
    precision[i] = 0;
}

/*---------------------------------------------------------------------------*/

void Allocator::optimalAllocate (CoeffSet **coeff, int nSets, 
				 int budget, int augment, Real *weight)
{
  Real bitBudget = 8*budget;

  Real lambda, lambdaLow, lambdaHigh;
  Real rateLow, rateHigh, currentRate;
  Real distLow, distHigh, currentDist;

  resetPrecision (nSets);
  
  lambdaLow = 0.0;
  allocateLambda (coeff, nSets, lambdaLow, rateLow, distLow, weight);
  if (rateLow < bitBudget) // this uses the largest possible # of bits
    return;                //   -- if this is within the budget, do it
  
  lambdaHigh = 1000000.0;
  Real lastRateHigh = -1;
  do {
    // try to use the smallest possible # of bits
    allocateLambda (coeff, nSets, lambdaHigh, rateHigh, distHigh, weight);

    // if this is still > bitBudget, try again w/ larger lambda
    if (rateHigh > bitBudget && lastRateHigh != rateHigh) {
      lambdaLow = lambdaHigh;
      rateLow = rateHigh;
      distLow = distHigh;
      lambdaHigh *= 10.0;
    }
  } while (rateHigh > bitBudget && lastRateHigh != rateHigh);

  // give up when changing lambda has no effect on things
  if (lastRateHigh == rateHigh)
    return;

  // Note rateLow will be > rateHigh
  if (rateLow < bitBudget) 
    error ("Failed to bracket bit budget = %d: rateLow = %g rateHigh = %g\n", 
	   budget, rateLow, rateHigh);
  
  while (lambdaHigh - lambdaLow > 0.01)  {
    lambda = (lambdaLow + lambdaHigh)/2.0;
    
    allocateLambda (coeff, nSets, lambda, currentRate, currentDist, weight);
    
    if (currentRate > bitBudget)
      lambdaLow = lambda;
    else
      lambdaHigh = lambda;
  }
  
  if (currentRate > bitBudget)  {
    lambda = lambdaHigh;
    allocateLambda (coeff, nSets, lambda, currentRate, currentDist, weight);
  }

  if (augment)
    greedyAugment (coeff, nSets, bitBudget-currentRate, weight);
}

/*---------------------------------------------------------------------------*/

void Allocator::allocateLambda  (CoeffSet **coeff, int nSets, 
				 Real lambda, Real &optimalRate, 
				 Real &optimalDist, Real *weight)
{
   int i, j;
   Real G, minG, minRate, minDist;

   optimalRate = optimalDist = 0.0;

   // want to minimize G = distortion + lambda * rate
   
   // loop through all rate-distortion curves
   for (i = 0; i < nSets; i++)  {
     minG = minRate = minDist = MaxReal;
     
     for (j = 0; j < coeff[i]->nQuant; j++) {
       G = weight[i]*coeff[i]->dist[j] + lambda * coeff[i]->rate[j];
       if (G < minG)  {
	 minG = G;
	 minRate = coeff[i]->rate[j];
	 minDist = weight[i]*coeff[i]->dist[j];
	 precision[i] = j;
       }
     }
     
     optimalRate += minRate;
     optimalDist += minDist;
   }
   //   printf ("lambda = %g  optimal rate = %g, optimal dist = %g\n",
   //	   lambda, optimalRate, optimalDist); 
}

/*---------------------------------------------------------------------------*/

void Allocator::greedyAugment (CoeffSet **coeff, int nSets, 
			       Real bitsLeft, Real *weight)
{
  int bestSet, newPrecision = -1;
  Real delta, maxDelta, bestDeltaDist, bestDeltaRate = 0;

  do {
    bestSet = -1;
    maxDelta = 0;
	    
    // Find best coeff set to augment 
    for (int i = 0; i < nSets; i++) {
      for (int j = precision[i]+1; j < coeff[i]->nQuant; j++) {
	Real deltaRate = coeff[i]->rate[j] -
	  coeff[i]->rate[precision[i]];
	Real deltaDist = -weight[i]*(coeff[i]->dist[j] -
	  coeff[i]->dist[precision[i]]);
	
	if (deltaRate != 0 && deltaRate <= bitsLeft) {
	  delta = deltaDist / deltaRate;
	  
	  if (delta > maxDelta) {
	    maxDelta = delta;
	    bestDeltaRate = deltaRate;
	    bestDeltaDist = deltaDist;
	    bestSet = i;
	    newPrecision = j;
	  }
	}
      }
    }
    
    if (bestSet != -1) {
      precision[bestSet] = newPrecision;
      bitsLeft -= bestDeltaRate;
    }
  } while (bestSet != -1);
}

/*---------------------------------------------------------------------------*/

void Allocator::print (CoeffSet **coeff, int nSets)
{
  Real totalRate = 0, totalDist = 0;
  int totalData = 0;

  printf ("Set  Precision     Rate                Distortion\n");
  for (int i = 0; i < nSets; i++) {
    printf ("%2d:     %2d       %9.2f   %5.2f   %11.2f   %7.2f\n", i, precision[i],
	    coeff[i]->rate[precision[i]], 
	    coeff[i]->rate[precision[i]]/(Real)coeff[i]->nData, 
	    coeff[i]->dist[precision[i]], 
	    coeff[i]->dist[precision[i]]/(Real)coeff[i]->nData);
    totalRate += coeff[i]->rate[precision[i]];
    totalDist += coeff[i]->dist[precision[i]];
    totalData += coeff[i]->nData;
  }
  Real rms = sqrt(totalDist/(Real)totalData);
  Real psnr = 20.0 * log(255.0/rms)/log(10.0);

  printf ("\n");
  printf ("total rate = %g\n", totalRate/8.0);
  printf ("total dist = %g\n", totalDist);
  printf ("total coeffs = %d\n", totalData);
  printf ("RMS error = %g\n", rms);
  printf ("PSNR (transform domain) = %g\n", psnr);
}

/*---------------------------------------------------------------------------*/
