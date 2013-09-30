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
// entropy.hh
//
// The EntropyCoder class provides a front-end to the arithmetic
// coding and histogram maintainence routines.  A big problem in
// developing efficient entropy coders is the initialization of the
// pdf's.  Rather than initializing our coders with, say, a
// generalized Gaussian, the coders included below are designed to
// adapt quickly to whatever coefficient distribution is found.  We
// provide two basic types of entropy coders, MonoLayerCoders and
// MultiLayerCoders.
//
// The MonoLayerCoder base class is a set of coders designed for use
// with single layer quantizers such as the UniformQuantizer.  The
// EscapeCoder class is a derived class that works well with highly
// concentrated pdf's such as the generalized Gaussians found in
// subband coding.  The entropy coder maintains two different
// frequency tables, each containing counts for each possible
// quantized value.  The first histogram also contains an extra
// "escape" symbol.  In the first histogram all symbol frequencies are
// set to 0 except for the escape symbol, which is set to 1.  In the
// second histogram, the frequency of all symbols is set to 1.  When a
// symbol is encountered for the first time, an escape is sent using
// the first histogram and then the symbol is output using the second
// histogram.  Symbol that have already been seen are coded using the
// first histogram.  This procedure prevents rare symbols from
// increasing the overall cost of coding symbols early on.  The
// procedure is described in detail in Chapter 6 of _Text Compression_
// by Bell, Cleary, and Witten.
// 
// The MultiLayerCoder takes a different approach to fast adaptation.
// Quantizer bins are hierarchical, and quantization consists of
// specifying an initial bin and a sequence of refinements.  All
// histograms are very small and adapt quickly.  This entropy coder is
// more complex than the single layer coder, since the histogram used
// depends on the level of refinement of the current symbol as well as
// previous refinements for that symbol.  Full details may be found in
// the description of layerd PCM quantization in D. Taubman and
// A. Zakhor, "Multirate 3-D subband coding of video", IEEE
// Transactions on Image Processing, Vol 3, No. 5, Sept, 1994.
//
/*---------------------------------------------------------------------------*/
#ifndef _ENTROPY_
#define _ENTROPY_
#include "iHisto.h"
#include "coder.hh"
/*---------------------------------------------------------------------------*/

class EntropyCoder {
public:
  EntropyCoder  (int histoCapacity);
  virtual ~EntropyCoder () {};

  // write symbol to encoder, update histogram if update flag set, return cost
  virtual Real write (Encoder *encoder, int symbol, char update, 
			int context1 = -1, int context2 = -1) = 0;

  // read symbol from decoder, update histogram if update flag set,
  //   return symbol
  virtual int read     (Decoder *decoder, char update, 
			int context1 = -1, int context2 = -1) = 0;

  // get cost of symbol, update histogram if update flag set
  virtual Real cost  (int symbol, char update, 
		        int context1 = -1, int context2 = -1) = 0;

  // update histogram, return symbol cost
  Real updateCost    (int symbol, int context1 = -1, int context2 = -1)
                       { return cost (symbol, TRUE, context1, context2); }

  // update histogram, don't compute symbol cost
  void update          (int symbol, int context1 = -1, int context2 = -1)
                       { cost (symbol, TRUE, context1, context2); }

  // reset histogram(s)
  virtual void reset  () = 0;

  inline Real entropy (Real p) { return -log(p)/Log2; };

  int histoCapacity;
};

/*---------------------------------------------------------------------------*/

class MonoLayerCoder : public EntropyCoder 
{
public:
  MonoLayerCoder (int histoCapacity);
  virtual ~MonoLayerCoder () {};

  // reset # of symbols in histogram
  virtual void setNSym (int nSym) = 0;

  int nSym;
};

/*---------------------------------------------------------------------------*/

class MultiLayerCoder : public EntropyCoder 
{
public:
  MultiLayerCoder (int histoCapacity);
  virtual ~MultiLayerCoder () {};
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class EscapeCoder : public MonoLayerCoder {
public:
  EscapeCoder (int histoCapacity);
  EscapeCoder (int nSym, int histoCapacity);
  ~EscapeCoder ();

  // write symbol to encoder, update histogram if update flag set, return cost
  Real write (Encoder *encoder, int symbol, char update, 
		int context1 = -1, int context2 = -1);

  // read symbol from decoder, update histogram if update flag set,
  //    return symbol
  int read  (Decoder *decoder, char update, 
	     int context1 = -1, int context2 = -1);

  // get cost of symbol, update histogram if update flag set
  Real cost (int symbol, char update, int context1 = -1, int context2 = -1);

  // reset histograms
  void   reset ();

  // reset # of symbols
  void setNSym (int newNSym);

  iHistogram *freq, *uniform;
  char *seen;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class LayerCoder : public MultiLayerCoder {
public:
  LayerCoder (int nLayers, int signedSym, int capacity);
  ~LayerCoder ();

  // write symbol to encoder, update histogram if update flag set, return cost
  Real write (Encoder *encoder, int symbol, char update, 
		int context1 = -1, int context2 = -1);
  
  // read symbol from decoder, update histogram if update flag set,
  //    return symbol 
  int read     (Decoder *decoder, char update, int context1 = -1, 
		int context2 = -1);

  // get cost of symbol, update histogram if update flag set
  Real cost (int symbol, char update, int layer, int context);

  // reset histograms
  void   reset ();

  int signedSym;        // TRUE if symbols signed, FALSE if unsigned
  int nLayers;          // # of bitplane layers
  int *nFreq;           // # of frequency counts for each layer
  iHistogram *** freq;  // frequency counts for each context
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
