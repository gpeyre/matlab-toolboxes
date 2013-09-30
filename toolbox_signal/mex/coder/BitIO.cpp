// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//

#include <stdlib.h>
#include <iostream.h>
#include <assert.h>
#include <math.h>
#include "global.hh"
#include "BitIO.h"


BitIn::BitIn(istream &is, ostream &log) : input(is), bitLog(log)
{
    bitsInBuf = 0;
    bitCount = 0;
}

BitOut::BitOut(ostream &out, ostream &log) : output(out), bitLog(log)
{
    bitsInBuf = 0;
    nBytes = 0;
}

void BitOut::flush(void)
{
    while (bitsInBuf) { output_bit(0); }
}


