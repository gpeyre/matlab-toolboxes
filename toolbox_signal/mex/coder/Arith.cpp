// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//

// These objects implement arithmetic coding as described in
// Bell, Cleary, and Witten "Text Compression" Prentice Hall.

#include <iostream.h>
#include <iomanip.h>
#include <assert.h>
#include <math.h>
#include "global.hh"
#include "BitIO.h"
#include "iHisto.h"
#include "Arith.h"

const int CodingValues::CodeValueBits = 16;
const long CodingValues::MaxFreq =  ((long)1 << (CodeValueBits - 2)) - 1;
const long CodingValues::One = ((long)1 << CodeValueBits) - 1;
const long CodingValues::Qtr = One / 4 + 1;
const long CodingValues::Half = 2 * Qtr;
const long CodingValues::ThreeQtr = 3 * Qtr;


ArithEncoder::ArithEncoder(BitOut &bo) : output(bo)
{
    low = 0;
    high = One;
    bitsToFollow = 0;
}


void
ArithEncoder::flush(void)
{
    for (int i = 0; i < CodeValueBits; i++) {
	if (low >= Half) {
	    bpf(1);
	    low -= Half;
	} else {
	    bpf(0);
	}
	low = 2 * low;
    }
    output.flush();
}

void ArithEncoder::Encode(int count, int countLeft, int countTot)
{
    //    cerr << "Encode(" << count << ", " << countLeft << ", " << countTot << ")\n";


    assert(count);
	
    long range = high - low + 1;
    high = low + (range * (countLeft + count)) / countTot - 1;
    low = low + (range * countLeft) / countTot;

    while (1) {
	if (high < Half) {
	    bpf(0);
	} else if (low >= Half) {
	    bpf(1);
	    low -= Half;
	    high -= Half;
	} else if (low >= Qtr && high < ThreeQtr) {
	    bitsToFollow++;
	    low -= Qtr;
	    high -= Qtr;
	} else {
	    break;
	}
	low = 2 * low;
	high = 2 * high + 1;
    }
}

ArithDecoder::ArithDecoder(BitIn &bi) : input(bi)
{
    low = 0;
    high = One;
    value = 0;
    for (int i = 0; i < CodeValueBits; i++) {
        value = (value << 1) + input.input_bit();
    }

}

int
ArithDecoder::Decode(iHistogram &h)
{
    long range; 
    int cum;
    int answer;
    int ct, ctLeft, ctTotal;

    ctTotal = h.TotalCount();
    range = high - low + 1;
    cum = (((long)(value - low) + 1) * ctTotal - 1) / range;
    answer = h.Symbol(cum);
    
    ct =  h.Count(answer); ctLeft = h.LeftCount(answer);

    //    cerr << "Decoder : cum= " << cum <<  "-> " << answer << "\n";
    //    cerr << "        ct = " << ct << " ctLeft = " << ctLeft << " ctTotal = " << ctTotal << "\n";

    high = low + (range * (ctLeft + ct)) / ctTotal - 1;
    low = low + (range * ctLeft) / ctTotal;
    while (1) {                 
        if (high < Half) {
	    ;
        } else if (low >= Half) { 
            value -= Half;
            low -= Half; 
            high -= Half;
        }
        else if (low >= Qtr   
		 && high < ThreeQtr) {
            value -= Qtr;
            low -= Qtr;    
            high -= Qtr;
        }
        else break;       
        low = 2 * low;
        high = 2 * high + 1;  
        value = 2 * value + input.input_bit(); 
    }
    return answer;
}
