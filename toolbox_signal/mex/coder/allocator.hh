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
// allocator.hh
//
// Given rate/distortion curves for nSets collections of transform
// coefficients (contained in CoeffSet objects), performs a
// constrained optimization of quantizer resolutions.  An array of
// quantizer precisions, precision[i], is found so that the sum (over
// i) of weight[i]*distortion[i][precision[i]] is minimized subject to
// the constraint that the sum (over i) of cost[i][precision[i]] is
// less than or equal to the given budget.
//
// Functions:
// ----------
// optimalAllocate     Does bit allocation using an algorithm described
//                     in Y. Shoham and A. Gersho, "Efficient bit
//                     allocation for an arbitrary set of quantizers,"
//                     IEEE Transactions on Acoustics, Speech, and
//                     Signal Processing, Vol. 36, No. 9,
//                     pp. 1445-1453, Sept 1988.
//
// greedyAugment       The Shoham & Gersho algorithm doesn't yield
//                     optimal allocations for all possible budgets.
//                     The optimalAllocate routine returns the best
//                     allocation that doesn't exceed the given
//                     budget.  GreedyAugment uses marginal analysis
//                     to greedily increase individual quantizer
//                     precisions until we reach the budget.
//                     Allocations will still be a little under budget
//                     but shouldn't be by much.  Note that the header
//                     info is not included in the overall budget.
//
// print               Prints out the current allocation
//
/*---------------------------------------------------------------------------*/
#ifndef _ALLOCATOR_
#define _ALLOCATOR_
#include "coeffset.hh"
/*---------------------------------------------------------------------------*/

class Allocator {
public:
  Allocator ();
  ~Allocator ();

  void optimalAllocate (CoeffSet **coeff, int nSets, 
			int budget, int augment, Real *weight);

  void greedyAugment   (CoeffSet **coeff, int nSets, Real bitsLeft, 
			Real *weight);

  void print           (CoeffSet **coeff, int nSets);

  void allocateLambda  (CoeffSet **coeff, int nSets, 
			Real lambda, Real &optimalRate, 
			Real &optimalDist, Real *weight);
  void resetPrecision  (int nSets);

  int nSets, *precision;
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
