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
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include "global.hh"
#include "quantizer.hh"
/*---------------------------------------------------------------------------*/
Quantizer::Quantizer (ErrorMetric *err) : err (err)
{
  data = NULL;
  nData = 0;
  max = min = sum = sumSq = mean = var = 0;
  initialDist = 0;
}

/*---------------------------------------------------------------------------*/
void Quantizer::getStats ()
{
  max = -MaxReal;
  min = MaxReal;
  sum = sumSq = 0;

  for (int i = 0; i < nData; i++) {
    if (data[i] < min)
      min = data[i];
    if (data[i] > max)
      max = data[i];
    sum += data[i];
    sumSq += square(data[i]);
  }
  mean = sum / (Real)nData;
  var = sumSq / (Real)nData - square (mean);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
UniformQuant::UniformQuant (MonoLayerCoder *entropy, int paramPrecision,
			    int zeroCenter, ErrorMetric *err) : 
  Quantizer (err), entropy (entropy), paramPrecision (paramPrecision), 
  zeroCenter (zeroCenter)
{
}

/*---------------------------------------------------------------------------*/
void UniformQuant::setDataEncode (Real *newData, int newNData)
{
  data = newData;
  nData = newNData;

  getStats ();
  
  if (zeroCenter) {
    max = fabs(max) > fabs(min) ? fabs(max) : fabs(min);
    min = -max;
  }
  imax = realToInt (max, paramPrecision);
  qmax = intToReal (imax, paramPrecision);
  // Make sure qmax >= max and -qmax <= min
  while (qmax < max)
    qmax = intToReal (++imax, paramPrecision);
  
  if (zeroCenter) {
    imin = -imax;
    qmin = -qmax;
  } else {
    imin = realToInt (min, paramPrecision);
    qmin = intToReal (imin, paramPrecision);
    // Make sure qmin <= min
    while (qmin > min)
      qmin = intToReal (--imin, paramPrecision);
  }

  imean = realToInt (mean, paramPrecision);
  qmean = intToReal (imean, paramPrecision);

  initialDist = 0;
  if (zeroCenter)
    for (int i = 0; i < nData; i++)
      initialDist += (*err)(data[i]);
  else
    for (int i = 0; i < nData; i++)
      initialDist += (*err)(data[i] - qmean);
}

/*---------------------------------------------------------------------------*/
void UniformQuant::setDataDecode (Real *newData, int newNData, 
				  int imax, int imin, int imean)
{
  data = newData;
  nData = newNData;

  if (imin < imax) {
    qmax = intToReal (imin, paramPrecision);
    
    if (zeroCenter) {
      qmin = -qmax;
    } else {
      qmin = intToReal (imin, paramPrecision);
    }
    qmean = intToReal (imean, paramPrecision);
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::getRateDist (int precision, Real minStepSize, 
				Real &rate, Real &dist) 
{
  if (precision > 0) {
    const int  nSteps = (1<<precision)-1;
    const Real stepSize = (qmax-qmin)/(Real)nSteps;
    const Real recipStepSize = 1.0/stepSize;
    
    if (stepSize < minStepSize) {
      rate = MaxReal;
      dist = MaxReal;
      return;
  }
    
    entropy->setNSym (nSteps);
    rate = dist = 0;
    
    for (int i = 0; i < nData; i++) {
      int symbol = (int)((data[i]-qmin)*recipStepSize);
      assert (symbol < nSteps && symbol >= 0);
      rate += entropy->cost (symbol, TRUE);
      Real reconstruct = qmin + ((Real)symbol + 0.5) * stepSize;
      dist += (*err) (data[i]-reconstruct);
    }
  } else {
    rate = dist = 0;
    if (zeroCenter) {
      for (int i = 0; i < nData; i++) {
	dist += (*err) (data[i]);
      }
    } else {
      for (int i = 0; i < nData; i++) {
	dist += (*err) (data[i] - qmean);
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::quantize (Encoder *encoder, int precision)
{
  if (precision > 0) {
    const int  nSteps = (1<<precision)-1;
    const Real stepSize = (qmax-qmin)/(Real)nSteps;
    const Real recipStepSize = 1.0/stepSize;
    
    entropy->setNSym (nSteps);
    
    for (int i = 0; i < nData; i++) {
      int symbol = (int)((data[i]-qmin)*recipStepSize);
      assert (symbol < nSteps && symbol >= 0);
      entropy->write (encoder, symbol, TRUE);
    }
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::dequantize (Decoder *decoder, int precision)
{
  if (precision > 0) {
    const int  nSteps = (1<<precision)-1;
    const Real stepSize = (qmax-qmin)/(Real)nSteps;
    int symbol;
    
    entropy->setNSym (nSteps);
    
    for (int i = 0; i < nData; i++) {
      symbol = entropy->read (decoder, TRUE);
      assert (symbol < nSteps && symbol >= 0);
      data[i] = qmin + ((Real)symbol + 0.5) * stepSize;
    }
  } else {
    for (int i = 0; i < nData; i++) {
      data[i] = qmean;
    }
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::writeHeader (Encoder *encoder, int precision)
{
  encoder->writeNonneg (precision);

  if (precision > 0) {
    encoder->writeInt (imax);
    if (!zeroCenter)
      encoder->writeInt (imin);
  } else {
    if (!zeroCenter)
      encoder->writeInt (imean);
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::readHeader  (Decoder *decoder, int &precision)
{
  precision = decoder->readNonneg ();

  if (precision > 0) {
    imax = decoder->readInt ();
    qmax = intToReal (imax, paramPrecision);
    
    if (zeroCenter) {
      qmin = -qmax;
    } else {
      imin = decoder->readInt ();
      qmin = intToReal (imin, paramPrecision);
    }
    qmean = 0;
  } else {
    if (!zeroCenter) {
      imean = decoder->readInt ();
      qmean = intToReal (imean, paramPrecision);
    } else {
      qmean = 0;
      imean = realToInt (qmean, paramPrecision);
    }
    qmax = qmin = qmean;
  }
}

/*---------------------------------------------------------------------------*/
void UniformQuant::setParams (int newParamPrecision, Real newMax, 
			      Real newMin, Real newMean)
{
  paramPrecision = newParamPrecision;
  qmax = newMax;
  qmin = newMin;
  qmean = newMean;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
LayerQuant::LayerQuant  (MultiLayerCoder *entropy, int paramPrecision,
			 int signedSym, int nLayers, ErrorMetric *err) :
  Quantizer (err), entropy (entropy), paramPrecision (paramPrecision), 
  signedSym (signedSym), nLayers (nLayers)
{
  currentLayer = -1;
  layerRate = new Real [nLayers];
  layerDist = new Real [nLayers];
  context = NULL;
  residual = NULL;
}

/*---------------------------------------------------------------------------*/
LayerQuant::~LayerQuant ()
{
  delete [] layerRate;
  delete [] layerDist;
  if (context != NULL)
    delete [] context;
  if (residual != NULL)
    delete [] residual;
}

/*---------------------------------------------------------------------------*/
void LayerQuant::setDataEncode  (Real *newData, int newNData)
{
  data = newData;
  nData = newNData;

  getStats ();
  
  if (signedSym) {
    max = fabs(max) > fabs(min) ? fabs(max) : fabs(min);
    min = -max;
  }
  imax = realToInt (max, paramPrecision);
  qmax = intToReal (imax, paramPrecision);
  // Make sure qmax >= max and -qmax <= min
  while (qmax < max)
    qmax = intToReal (++imax, paramPrecision);
  
  if (signedSym) {
    imin = -imax;
    qmin = -qmax;
  } else {
    imin = realToInt (min, paramPrecision);
    qmin = intToReal (imin, paramPrecision);
    // Make sure qmin <= min
    while (qmin > min)
      qmin = intToReal (--imin, paramPrecision);
  }
  if (signedSym)
    threshold = qmax/2.0;
  else
    threshold = (qmax - qmin)/2.0;
  
  currentLayer = -1;
  if (context != NULL)
    delete [] context;
  if (residual != NULL)
    delete [] residual;
  context = new int [nData];
  residual = new Real [nData];

  resetLayer ();

  initialDist = 0;
  for (int i = 0; i < nData; i++)
    initialDist += (*err)(residual[i]);

  entropy->reset ();
}

/*---------------------------------------------------------------------------*/
void LayerQuant::setDataDecode  (Real *newData, int newNData, int imax, 
				 int imin, int imean)
{
  data = newData;
  nData = newNData;

  if (imin < imax) {
    qmax = intToReal (imin, paramPrecision);
    qmean = intToReal (imean, paramPrecision);
    
    if (signedSym) {
      qmin = -qmax;
      threshold = qmax/2.0;
    } else {
      qmin = intToReal (imin, paramPrecision);
      threshold = (qmax - qmin)/2.0;
    }
    
  }

  currentLayer = -1;
  if (context != NULL)
    delete [] context;
  context = new int [nData];
  if (signedSym) {
    for (int i = 0; i < nData; i++) {
      data[i] = 0;
      context[i] = 0;
    }
  } else {
    for (int i = 0; i < nData; i++) {
      data[i] = qmean;
      context[i] = 0;
    }
  }

  resetLayer ();
  entropy->reset ();
}

/*---------------------------------------------------------------------------*/
void LayerQuant::getRateDist (int precision, Real minStepSize, 
			      Real &rate, Real &dist)
{
  assert (precision <= nLayers);

  Real currentRate = 0;
  Real currentDist = initialDist;
  
  for (int i = 0; i < precision; i++) {
    // rates & distortions have been computed for layers up to currentLayer
    if (i <= currentLayer) {
      currentRate += layerRate[i];
      currentDist += layerDist[i];
    } else {
      if (threshold > minStepSize) {
	quantizeLayer (NULL);
	currentRate += layerRate[i];
	currentDist += layerDist[i];
      } else {
	layerRate[i] = MaxReal;
	layerDist[i] = -MaxReal;
	currentRate = MaxReal;
	currentDist = MaxReal;
      }
    }
  }
  rate = currentRate;
  dist = currentDist;
}

/*---------------------------------------------------------------------------*/
void LayerQuant::quantize (Encoder *encoder, int precision)
{
  resetLayer ();
  entropy->reset ();

  for (int i = 0; i < precision; i++) {
    quantizeLayer (encoder);
  }
}

/*---------------------------------------------------------------------------*/
void LayerQuant::resetLayer ()
{
  if (residual != NULL) {
    if (signedSym) {
      for (int i = 0 ;i < nData; i++) {
	residual[i] = data[i];
	context[i] = 0;
      }
    } else {
      // on the first layer we remove the mean
      for (int i = 0 ;i < nData; i++) {
	residual[i] = data[i] - 0.5*(qmax+qmin);
	context[i] = 0;
      }
    }
  } else {
    if (signedSym) {
      for (int i = 0 ;i < nData; i++) {
	context[i] = 0;
      }
    } else {
      // on the first layer we remove the mean
      for (int i = 0 ;i < nData; i++) {
	context[i] = 0;
      }
    }
  }
  currentLayer = -1;

  if (signedSym)
    threshold = qmax/2.0;
  else
    threshold = (qmax - qmin)/2.0;
}

/*---------------------------------------------------------------------------*/
// deltaRate = bits required to code current layer
// deltaDist = reduction in distortion from current layer

void LayerQuant::quantizeLayer (Encoder *encoder)
{
  const Real halfThreshold = 0.5 * threshold;
  const Real threeHalvesThreshold = 1.5 * threshold;
  Real deltaRate = 0, deltaDist = 0;
  int symbol;

  currentLayer++;

  //  printf ("current layer = %d, threshold = %g\n", currentLayer, threshold);
  if (signedSym) {
    for (int i = 0; i < nData; i++) {
      deltaDist -= (*err)(residual[i]);  // subtract off old error
      //      Real oldResid = residual[i];
      if (context[i] == 0) {
	if (residual[i] > threshold) {
	  symbol = 1;
	  residual[i] -= threeHalvesThreshold;
	} else if (residual[i] < -threshold) {
	  symbol = -1;
	  residual[i] += threeHalvesThreshold;
	} else {
	  symbol = 0;
	}
      } else {
	if (residual[i] > 0) {
	  symbol = 1;
	  residual[i] -= halfThreshold;
	} else {
	  symbol = 0;
	  residual[i] += halfThreshold;
	}
      }

      //      printf ("%d (%d): %g -> %g\n", symbol, context[i], oldResid,
      //	      residual[i]);
      deltaDist += (*err)(residual[i]);  // add in new error
      deltaRate += entropy->write (encoder, symbol, TRUE, currentLayer,
				   context[i]);
      context[i] = 2*context[i] + symbol;
    }

  } else {
    for (int i = 0; i < nData; i++) {
      deltaDist -= (*err)(residual[i]);  // subtract off old error

      if (residual[i] > 0) {
	symbol = 1;
	residual[i] -= halfThreshold;
      } else {
	symbol = 0;
	residual[i] += halfThreshold;
      }

      deltaDist += (*err)(residual[i]);  // add in new error
      deltaRate += entropy->write (encoder, symbol, TRUE, currentLayer,
				   context[i]);
      context[i] = 2*context[i] + symbol;   
    }
  }

  //  printf ("layer= %d  cost = %g\n", currentLayer, deltaRate);
  threshold *= 0.5;

  layerRate[currentLayer] = deltaRate;
  layerDist[currentLayer] = deltaDist;
}

/*---------------------------------------------------------------------------*/
void LayerQuant::dequantize   (Decoder *decoder, int precision)
{
  resetLayer ();

  for (int i = 0; i < precision; i++) {
    dequantizeLayer (decoder);
  }
}

/*---------------------------------------------------------------------------*/
void LayerQuant::dequantizeLayer (Decoder *decoder)
{
  int symbol;

  currentLayer++;

  //  printf ("current layer = %d, threshold = %g\n", currentLayer, threshold);
  if (signedSym) {
    for (int i = 0; i < nData; i++) {
      symbol = entropy->read (decoder, TRUE, currentLayer,
			      context[i]);

      //      int oldData = data[i];
      if (context[i] == 0)
	data[i] += 1.5*threshold * symbol;
      else 
	data[i] += (symbol - 0.5) * threshold;
      //      printf ("%d (%d):  %g -> %g  (%g)\n", symbol, context[i],
      //	      oldData, data[i], threshold);
	  
      context[i] = 2*context[i] + symbol;
    }
  } else {
    for (int i = 0; i < nData; i++) {
      symbol = entropy->read (decoder, TRUE, currentLayer, context[i]);

      data[i] += (symbol - 0.5) * threshold;
      context[i] = 2*context[i] + symbol;
    }
  }

  threshold *= 0.5;
}

/*---------------------------------------------------------------------------*/
void LayerQuant::writeHeader (Encoder *encoder, int precision)
{
  encoder->writeNonneg (precision);

  if (precision > 0) {
    encoder->writeInt (imax);
    if (!signedSym)
      encoder->writeInt (imin);
  } else {
    if (!signedSym)
      encoder->writeInt (imean);
  }

  //  printf ("precision = %d\n", precision);
  //  printf ("qmax = %g  qmin = %g  qmean = %g\n", qmax, qmin, qmean);
}

/*---------------------------------------------------------------------------*/
void LayerQuant::readHeader  (Decoder *decoder, int &precision)
{
  precision = decoder->readNonneg ();

  //  printf ("precision = %d\n", precision);

  if (precision > 0) {
    imax = decoder->readInt ();
    qmax = intToReal (imax, paramPrecision);
    
    if (signedSym) {
      qmin = -qmax;
      qmean = 0;
    } else {
      imin = decoder->readInt ();
      qmin = intToReal (imin, paramPrecision);
      qmean = 0.5 * (qmax + qmin);
    }
  } else {
    if (!signedSym) {
      imean = decoder->readInt ();
      qmean = intToReal (imean, paramPrecision);
    } else {
      qmean = 0;
      imean = realToInt (qmean, paramPrecision);
    }
    qmax = qmin = qmean;
  }
  //  printf ("qmax = %g  qmin = %g  qmean = %g\n", qmax, qmin, qmean);
}

/*---------------------------------------------------------------------------*/
void LayerQuant::setParams (int newParamPrecision, Real newMax, 
			    Real newMin, Real newMean)
{
  paramPrecision = newParamPrecision;
  qmax = newMax;
  qmin = newMin;
  qmean = newMean;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
