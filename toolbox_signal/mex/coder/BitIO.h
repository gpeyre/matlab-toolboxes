// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//


#ifndef _BITIO_
#define _BITIO_

#include <stdlib.h>
#include <iostream.h>
#include <assert.h>

#define PRINTBITS 0


// The BitIn object allows you to read bits out of an istream with
//optional logging in ascii to an ostream.
//
// Remember that istreams can be istrstreams, which read out of in
// core buffers.
//
// Casting a BitIn to void * tells you whether there 
// is anything in it.

class BitIn {
public:
    BitIn(istream &is, ostream &log=cerr);
    inline int BitCount() { return bitCount; }
    inline int input_bit(void)
    {
	int bit;
	if (!bitsInBuf) {
	    assert(!input.eof()); /* more bits to read */
	    bitsInBuf = 8;
	    input.get(bitBuf);
	}
		    
	bit = (bitBuf & 128) >> 7;
	bitBuf = bitBuf << 1;
	bitsInBuf--;
#if PRINTBITS
	bitLog << "bit(" << bit << ")\n";
#endif
	bitCount++;
	return bit;
    }
    inline BitIn& operator>>(int &i) { i = input_bit(); return *this; }
    inline BitIn& operator>>(char &c) { c = input_bit(); return *this; }
    operator void *() {
	if (bitsInBuf) return this;
	if (input.get(bitBuf)) {
	    bitsInBuf = 8;
	}
	return input;
    }
private:
    unsigned char bitBuf;
    int bitsInBuf;
    int bitCount;
    istream &input;    
    ostream &bitLog;
};
	

//
// BitOut allows you to write bits to an ostream.
// Remember that an ostream can be an in core buffer if
// you want (ostrstream).	
		
class BitOut {
public:
    BitOut(ostream &out, ostream &log=cerr);
    inline int  BitCount() { return bitsInBuf + 8 * nBytes; }

    void output_bit(int bit)
    {
#if PRINTBITS
	bitLog << "bit(" << bit << ")\n";
#endif
	assert(bitsInBuf < 8);
	assert((bit == 0) || (bit == 1));
		    
	bitBuf = (bitBuf << 1) | bit;
	bitsInBuf++;

	if (bitsInBuf == 8) {
	    output.put(bitBuf);
	    nBytes++;
	    bitsInBuf = 0;
	}
		    
    }
    inline BitOut& operator<<(int i) { output_bit(i); return *this; }
    void flush(void);
private:
    unsigned char bitBuf;
    int bitsInBuf;
    int nBytes;
    ostream &output;
    ostream &bitLog;
};

#endif
