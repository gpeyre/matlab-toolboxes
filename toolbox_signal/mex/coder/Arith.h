// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//

// These objects implement arithmetic coding as described in
// Bell, Cleary, and Witten "Text Compression" Prentice Hall.

#ifndef _ARITH_CODER_
#define _ARITH_CODER_

#include <iostream.h>
#include "iHisto.h"
#include "BitIO.h"

class CodingValues {
public:
    long low, high, value;
    static const int CodeValueBits;
    static const long MaxFreq;
    static const long One;
    static const long Qtr;
    static const long Half;
    static const long ThreeQtr;
};		

class ArithEncoder : public CodingValues {    
    BitOut &output;
    inline void bpf(int bit)
    {
	output.output_bit(bit);
	for (int i = 0; i < bitsToFollow; i++) {
	    output.output_bit(1 - bit);
	}
	bitsToFollow = 0;
    }
public:
    ArithEncoder(BitOut &bo);
    void Encode(int count, int countLeft, int countTot);
    void flush(void);
    int bitsToFollow;
};

class ArithDecoder : public CodingValues {
    BitIn &input;
public:
    ArithDecoder(BitIn &bi);
    int Decode(iHistogram &h);
};

#endif

