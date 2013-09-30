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
#include <iostream.h>
#include <stdio.h>
#include <math.h>
#include "global.hh"
#include "coder.hh"
/*---------------------------------------------------------------------------*/

Encoder::Encoder (ostream &out, ostream &log)
{
   bitout   = new BitOut (out, log);
   arith    = new ArithEncoder (*bitout);
   intcoder = new CdeltaEncode (bitout);
}

/*---------------------------------------------------------------------------*/

Encoder::~Encoder ()
{
   delete intcoder;
   delete arith;
   delete bitout;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Decoder::Decoder (istream &in, ostream &log)
{
  bitin    = new BitIn (in, log);
  intcoder = new CdeltaDecode (bitin);
  arith    = NULL;
}

/*---------------------------------------------------------------------------*/

Decoder::~Decoder ()
{
   delete intcoder;
   if (arith != NULL) delete arith;
   delete bitin;
}

/*---------------------------------------------------------------------------*/
