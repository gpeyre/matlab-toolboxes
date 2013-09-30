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
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include "global.hh"
#include "quantizer.hh"
#include "coeffset.hh"
/*---------------------------------------------------------------------------*/

CoeffSet::CoeffSet (Real *data, int nData, Quantizer *quant) :
  nData (nData), data (data), quant (quant)
{
  rate = dist = NULL;
}

/*---------------------------------------------------------------------------*/
CoeffSet::~CoeffSet ()
{
  if (rate != NULL) {
    delete [] rate;
    delete [] dist;
  }
}

/*---------------------------------------------------------------------------*/

void CoeffSet::getRateDist (int newNQuant, Real minStepSize)
{
  nQuant = newNQuant;
  rate = new Real [nQuant];
  dist = new Real [nQuant];

  quant->setDataEncode (data, nData);

  for (int i = 0 ;i < nQuant; i++) {
    quant->getRateDist (i, minStepSize, rate[i], dist[i]);
  }
}

/*---------------------------------------------------------------------------*/
