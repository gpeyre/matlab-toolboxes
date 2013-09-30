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
// coder.hh
//
// Front end to I/O routines for arithmetic coder.  Also has functions
// for reading and writing arbitrarily sized integers and individual
// bits.
//
// Note: in the current implementation there is no problem switching
// from reading/writing integers/bits to using the arithmetic coder.
// However, switching from the arithmetic coder to reading/writing
// integers/bits is problematic, since the arithmetic coder keeps a
// number of bits in an internal buffer.  This results in a loss of
// synchronization and will produce errors.  For this reason we write
// all header information requiring the output of integers/bits first
// before coding any coefficients.
//
/*---------------------------------------------------------------------------*/
#ifndef _CODER_
#define _CODER_
/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <iostream.h>
#include <iomanip.h>
#include "iHisto.h"
#include "BitIO.h"
#include "Arith.h"
#include "IntCoding.hh"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class Encoder {
public:
  BitOut *bitout;

  ArithEncoder *arith;
  CdeltaEncode *intcoder;

  Encoder  (ostream &out, ostream &log=cerr);
  ~Encoder ();

  // write symbol
  void writeSymbol (int symbol, iHistogram *h)
              { arith->Encode(h->Count(symbol), h->LeftCount(symbol),
			      h->TotalCount()); };
  
  void writePositive (int i)     { intcoder->EncodePositive (i); };
  void writeNonneg (int i)       { intcoder->EncodeNonNegative (i); };
  void writeInt (int i)          { intcoder->Encode (i); };
  void writeBit (int bit)        { bitout->output_bit (bit); }
  void writeNBits (int n, int i) { while(n--){writeBit((i&(1 << n))!=0);} };

  void flush () { arith->flush (); };  
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class Decoder {
public:
  BitIn *bitin;
  ArithDecoder *arith;
  CdeltaDecode *intcoder;
  
  Decoder  (istream &in, ostream &log=cerr);
  ~Decoder ();
  
  int readSymbol (iHistogram *h) 
                        { if (arith==NULL) arith = new ArithDecoder(*bitin); 
                          return arith->Decode(*h); };
  int readPositive ()   { return intcoder->DecodePositive (); };
  int readNonneg ()     { return intcoder->DecodeNonNegative (); };
  int readInt ()        { return intcoder->Decode (); };
  int readBit ()        { return bitin->input_bit (); }
  int readNBits (int n) { int i = 0; while (n--) 
                           {i |= readBit()*(1<<n);} return i; }
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
