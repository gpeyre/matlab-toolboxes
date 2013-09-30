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
// quantizer.hh
//
// Routines for quantizing coefficients using various levels of
// precision.  Two basic styles of quantizer are included.  
//
// The UniformQuantizer class is a single-layer, non-embedded
// quantizer. It maps each coefficient to a single symbol which is
// later entropy coded.  The fineness of the quantization is
// controlled by the precision parameter.  The UniformQuantizer class
// codes symbols using a MonoLayer entropy coder.
//
// The LayerQuantizer class is a multi-layer, embedded quantizer.
// It maps each coefficient to a sequence of symbols that correspond
// to an initial coarse quantization and subsequent refinements.  The
// precision parameter controls the number of refinements used.  By
// conditioning the probabilities of refinements on coarser-scale
// quantizations of the same symbol we can achieve the equal or greater
// efficiency than we can with a single layer coder.  Greater
// efficiency is usually possible because the coarse-scale histograms
// are very small (2 or 3 symbols) and adapt quickly.  The
// LayerQuantizer class entropy codes symbols using a MultiLayer
// entropy coder that takes into account the use of multiple symbols
// for each coefficient and conditions the coding appropriately.
// 
// Functions:
// ----------
// setDataEncode    Set data to be encoded.  This must be called
//                  before invoking quantize or getRateDist because
//                  some internals need to be initialized.
// setDataDecode    Set data to be decoded.  This must be called
//                  before invoking dequantize because some
//                  internals need to be initialized.
// getRateDist      Generate a rate/distortion curve for a set of
//                  coefficients.
// quantize         Quantize a set of coefficients, entropy code them,
//                  and write the result using the given encoder
// dequantize       Read a set of coefficients from the given decoder,
//                  decode them and dequantize them.
// writeHeader      Write quantizer parameters using the given encoder.
// readHeader       Set quantizer parameters using a stored header
//                  read using the given decoder.  
/*---------------------------------------------------------------------------*/
#ifndef _QUANTIZER_
#define _QUANTIZER_
#include "entropy.hh"
#include "metric.hh"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Quantizer {
public:
  Quantizer  (ErrorMetric *err);
  virtual ~Quantizer () {};

  virtual void setDataEncode  (Real *data, int nData) = 0;
  virtual void setDataDecode  (Real *data, int nData, int imax = -1,
			       int imin = 1, int imean = -1) = 0;

  virtual void getRateDist  (int precision, Real minStepSize, 
			     Real &rate, Real &dist) = 0; 

  virtual void quantize     (Encoder *encoder, int precision) = 0;
  virtual void dequantize   (Decoder *decoder, int precision) = 0;

  virtual void writeHeader (Encoder *encoder, int precision) = 0;
  virtual void readHeader  (Decoder *decoder, int &precision) = 0;

  void getStats ();

  int realToInt (Real x, int precision) 
       { return (int)(x * exp (precision * Log2) + 0.5); };
  Real intToReal (int n, int precision)
       { return ((Real)n / exp (precision * Log2)); };

  ErrorMetric *err;
  Real *data;
  int nData;
  Real max, min, mean, var, sum, sumSq;
  Real initialDist;
};

/*---------------------------------------------------------------------------*/

class UniformQuant : public Quantizer {
public:
  UniformQuant  (MonoLayerCoder *entropy, int paramPrecision, 
		 int zeroCenter, ErrorMetric *err = NULL); 
  ~UniformQuant () {};

  void setDataEncode (Real *newData, int newNData);
  void setDataDecode (Real *newData, int newNData, 
		      int imax = -1, int imin = 1, int imean = -1);
  void getRateDist  (int precision, Real minStepSize, Real &rate, Real &dist);

  void quantize     (Encoder *encoder, int precision);
  void dequantize   (Decoder *decoder, int precision);

  void writeHeader (Encoder *encoder, int precision);
  void readHeader  (Decoder *decoder, int &precision);

  void setParams (int paramPrecision, Real max, Real min, Real mean);

  int imin, imax, imean;
  Real qmin, qmax, qmean;
  MonoLayerCoder *entropy;
  int paramPrecision, zeroCenter;
};

/*---------------------------------------------------------------------------*/

class LayerQuant : public Quantizer {
public:
  LayerQuant  (MultiLayerCoder *entropy, int paramPrecision, 
	       int signedSym, int nLayers, ErrorMetric *err = NULL);
  ~LayerQuant ();

  void setDataEncode  (Real *data, int nData);
  void setDataDecode  (Real *data, int nData, int imax = -1,
		       int imin = 1, int imean = -1);

  void getRateDist  (int precision, Real minStepSize, 
		     Real &rate, Real &dist); 

  void quantize     (Encoder *encoder, int precision);
  void dequantize   (Decoder *decoder, int precision);

  void writeHeader (Encoder *encoder, int precision);
  void readHeader  (Decoder *decoder, int &precision);

  void quantizeLayer   (Encoder *encoder);
  void dequantizeLayer (Decoder *decoder);
  void resetLayer      ();

  void setParams (int paramPrecision, Real max, Real min, Real mean);

  MultiLayerCoder *entropy;
  int paramPrecision, signedSym, nLayers;
  int imin, imax, imean;
  Real qmin, qmax, qmean;

  int currentLayer, *context;
  Real *layerRate, *layerDist, *residual;
  Real threshold;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
