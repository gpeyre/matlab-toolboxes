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
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
EntropyCoder::EntropyCoder (int histoCapacity) :
                            histoCapacity (histoCapacity)
{
}

/*---------------------------------------------------------------------------*/
MonoLayerCoder::MonoLayerCoder (int histoCapacity) :
                                EntropyCoder (histoCapacity)
{
}

/*---------------------------------------------------------------------------*/
MultiLayerCoder::MultiLayerCoder (int histoCapacity) :
                                  EntropyCoder (histoCapacity)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
EscapeCoder::EscapeCoder (int histoCapacity) : 
                          MonoLayerCoder (histoCapacity)
{
  freq = uniform = NULL;
  seen = NULL;
}

/*---------------------------------------------------------------------------*/

EscapeCoder::EscapeCoder (int histoCapacity, int nSym) : 
                          MonoLayerCoder (histoCapacity)
{
  seen = NULL;
  setNSym (nSym);
}

/*---------------------------------------------------------------------------*/
EscapeCoder::~EscapeCoder ()
{
  if (seen != NULL) {
    delete [] seen;
    delete freq;
    delete uniform;
  }
}

/*---------------------------------------------------------------------------*/

void EscapeCoder::setNSym (int newNSym)
{
  if (seen != NULL) {
    delete [] seen;
    delete freq;
    delete uniform;
  }

  nSym = newNSym;

  freq = new iHistogram (nSym+1, histoCapacity);
  uniform = new iHistogram (nSym, histoCapacity);

  seen = new char [nSym];
  for (int i = 0; i < nSym; i++)
    seen[i] = FALSE;

  reset ();
}

/*---------------------------------------------------------------------------*/
// write symbol to encoder, update histogram if update flag set, return cost
Real EscapeCoder::write (Encoder *encoder, int symbol, char update, 
			   int context1, int context2)
{
  const int Escape = nSym;
  Real bits;

  context1 = context2 = -1;  // prevents warning from compiler

  if (seen[symbol]) {
    bits = freq->Entropy(symbol);
    if (encoder != NULL)
      encoder->writeSymbol (symbol, freq);
  } else {
    if (encoder != NULL) {
      encoder->writeSymbol (Escape, freq);
      encoder->writeSymbol (symbol, uniform);
    }
    bits = freq->Entropy(Escape) + uniform->Entropy(symbol);

    if (update)
      seen[symbol] = TRUE;
  }
  if (update)
    freq->IncCount(symbol);

  return bits;
}

/*---------------------------------------------------------------------------*/
// read symbol from decoder, update histogram if update flag set
int EscapeCoder::read (Decoder *decoder, char update, 
		       int context1, int context2)
{
  const int Escape = nSym;
  int symbol;

  context1 = context2 = -1;  // prevents warning from compiler

  symbol = decoder->readSymbol (freq);
  assert (symbol >= 0 && symbol <= nSym);
  if (symbol == Escape) {
    symbol = decoder->readSymbol (uniform);
    assert (symbol >= 0 && symbol < nSym);
  }

  if (update)
    freq->IncCount(symbol);

  return symbol;
}

/*---------------------------------------------------------------------------*/
Real EscapeCoder::cost (int symbol, char update, int context1, int context2)
{
  const int Escape = nSym;
  Real bits;

  context1 = context2 = -1;  // prevents warning from compiler

  if (seen[symbol]) {
    bits = freq->Entropy(symbol);

  } else {
    bits = freq->Entropy(Escape) + uniform->Entropy(symbol);

    if (update)
      seen[symbol] = TRUE;
  }
  if (update)
    freq->IncCount(symbol);

  return bits;
}

/*---------------------------------------------------------------------------*/
void EscapeCoder::reset ()
{
  int *temp = new int [nSym+1];
  int i;

  for (i = 0; i < nSym; i++) // no symbols observed -- set all counts to 0
    temp[i] = 0;
  temp[nSym] = 1;                // except for escape symbol
  freq->InitCounts (temp);

  for (i = 0; i < nSym; i++) // uniform histogram -- set all counts to 1
    temp[i] = 1;
  uniform->InitCounts (temp);

  delete [] temp;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
LayerCoder::LayerCoder (int nLayers, int signedSym, int capacity) :
  MultiLayerCoder (capacity), signedSym(signedSym), nLayers(nLayers)

{
  int i, j;
  nFreq = new int [nLayers];
  freq = new iHistogram** [nLayers];

  if (signedSym)
    nFreq[0] = 3;
  else
    nFreq[0] = 2;
  freq[0] = new iHistogram* [nFreq[0]];

  for (i = 1; i < nLayers; i++) {
    if (signedSym)
      nFreq[i] = 2 * nFreq[i-1] + 1;
    else 
      nFreq[i] = 2 * nFreq[i-1];
    freq[i] = new iHistogram* [nFreq[i]];
  }

  for (i = 0; i < nLayers; i++) {
    for (j = 0; j < nFreq[i]; j++) {
      freq[i][j] = new iHistogram (3, capacity);
    }
  }

  reset ();
}

/*---------------------------------------------------------------------------*/
LayerCoder::~LayerCoder ()
{
  for (int i = 0; i < nLayers; i++) {
    for (int j = 0; j < nFreq[i]; j++) {
      delete freq[i][j];
    }
    delete [] freq[i];
  }
  delete [] freq;
  delete [] nFreq;
}

/*---------------------------------------------------------------------------*/
// for signed symbol coder, symbols will be
//   -1, 0, 1  for context = 0
//    0, 1     for context > 0
//    -1, 0    for context < 0
//   1 will be added to all incoming symbols
// for unsigned coder, symbols will be
//    0, 1  -- 1 will be added to all incoming symbols

void LayerCoder::reset ()
{
  int plusCounts[3] =  {0, 1, 1};
  int zeroCounts[3] =  {1, 1, 1};
  int i, j;

  if (signedSym) {
    for (i = 0; i < nLayers; i++) {
      for (j = -nFreq[i]/2; j < 0; j++)
	freq[i][j+nFreq[i]/2]->InitCounts (plusCounts);
	//	freq[i][j+nFreq[i]/2]->InitCounts (minus_counts);
      freq[i][nFreq[i]/2]->InitCounts (zeroCounts);
      for (j = 1; j <= nFreq[i]/2; j++)
	freq[i][j+nFreq[i]/2]->InitCounts (plusCounts);
    }
  } else {
    for (i = 0; i < nLayers; i++) {
      for (j = 0; j < nFreq[i]; j++) {
	freq[i][j]->InitCounts (plusCounts);
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
// write symbol to encoder, update histogram if update flag set, return cost
Real LayerCoder::write (Encoder *encoder, int symbol, char update, 
			  int layer, int context)
{
  symbol++;
  if (signedSym) {
    context += nFreq[layer]/2;
  }
  if (encoder != NULL)
    encoder->writeSymbol (symbol, freq[layer][context]);
  Real bits = freq[layer][context]->Entropy(symbol);
  
  if (update)
    freq[layer][context]->IncCount(symbol);
  return bits;
}

/*---------------------------------------------------------------------------*/
// read symbol from decoder, update histogram if update flag set
int LayerCoder::read  (Decoder *decoder, char update, int layer, int context)
{
  int symbol;

  if (signedSym) {
    context += nFreq[layer]/2;
  }

  symbol = decoder->readSymbol (freq[layer][context]);
  if (update)
    freq[layer][context]->IncCount(symbol);

  symbol--;
  return symbol;
}

/*---------------------------------------------------------------------------*/
Real LayerCoder::cost (int symbol, char update, int layer, int context)
{
  symbol++;
  if (signedSym) {
    context += nFreq[layer]/2;
  }
  Real bits = freq[layer][context]->Entropy(symbol);
  
  if (update)
    freq[layer][context]->IncCount(symbol);
  return bits;
}

/*---------------------------------------------------------------------------*/
