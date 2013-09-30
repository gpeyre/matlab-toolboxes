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
// coeffset.hh
//
// Each CoeffSet object corresponds to a subset of transform
// coefficients.  the object contains rate/distortion curves for its
// associated set of coefficients which it obtains from the quantizer
// passed in during construction.  Also included are functions for
// encoding and decoding the set given a quantizer precision level and
// functions for reading and writing the quantizer parameters.
//
// Functions:
// ----------
// getRateDist      Get rate/distortion curves for the set.
// encode           Encode coefficients with the specified precision
//                  and write them using the given encoder.
// decode           Read quantized coefficients from the given decoder
//                  and dequantize them.
// writeHeader      Write quantizer parameters using the given encoder
// readHeader       Read quantizer parameters using the given decoder
//
/*---------------------------------------------------------------------------*/
#ifndef _COEFFSET_
#define _COEFFSET_
#include "coder.hh"
#include "quantizer.hh"
/*---------------------------------------------------------------------------*/

class CoeffSet {
public:
  CoeffSet (Real *data, int ndata, Quantizer *quant);
  ~CoeffSet ();

  void getRateDist (int nQuant, Real minStepSize);

  void writeHeader (Encoder *encoder, int precision)
        { quant->writeHeader (encoder, precision); };
  void readHeader  (Decoder *decoder, int &precision)
        { quant->readHeader  (decoder, precision); };

  void encode (Encoder *encoder, int precision)
        { quant->quantize (encoder, precision); };
  void decode (Decoder *decoder, int precision)
        { quant->setDataDecode (data, nData);
	  quant->dequantize (decoder, precision); };

  int  nData;
  Real *data;
  Quantizer *quant;

  // The data in each CoeffSet is quantized using quantizers with
  // different levels of precision.  From this we generate a
  // rate/distortion curve which is used to optimize bit allocation.
  int  nQuant;    // # of quantizers used
  Real *rate;     // cost (in bits) to perform each quantization
  Real *dist;     // distortion incurred using each quantizer 
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
