// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//

/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <iostream.h>
#include <iomanip.h>

#include "BitIO.h"
#include "IntCoding.hh"
/*---------------------------------------------------------------------------*/

#define PrintBits 0

/*---------------------------------------------------------------------------*/

#if PrintBits
#define PB(x) cerr << (x);
#define PNL   cerr << "\n";
#define PCol  cerr << ":";
#else
#define PB(x)
#define PNL
#define PCol
#endif

#define PRINTCALLS 0

/*---------------------------------------------------------------------------*/
void CdeltaEncode::EncodePositive(int i)
{
  int digits;
  assert(i > 0 && i < (1 << 30));
  
#if PRINTCALLS
  cerr << "CdeltaEncode::EncodePositive(" << i << ")\n";
#endif
  
  for (digits = 30; digits > 0; digits--) {
    if (i >= (1 << digits)) {
      PB(0);	    
      output->output_bit(0);
    }
  }
  PB(1);
  output->output_bit(1);
  PCol;
  for (digits = 29; digits >= 0; digits--) {
    if (i >= (2 << digits)) {
      if (i & (1 << digits)) {
	PB(1);		
	output->output_bit(1);
      } else {
	PB(0);
	output->output_bit(0);
      }
    }
  }
  PNL;
}

/*---------------------------------------------------------------------------*/

void CdeltaEncode::EncodeNonNegative(int i)
{
#if PRINTCALLS
  cerr << "CdeltaEncode::EncodeNonNegative(" << i << ")\n";
#endif
  EncodePositive(i + 1);
}

/*---------------------------------------------------------------------------*/
void CdeltaEncode::Encode(int i)
{
#if PRINTCALLS
  cerr << "CdeltaEncode::Encode(" << i << ")\n";
#endif
  if (i < 0) {
    output->output_bit(1);
    EncodePositive(-i);
  } else {
    output->output_bit(0);
    EncodeNonNegative(i);
  }
}

/*---------------------------------------------------------------------------*/

int CdeltaDecode::DecodePositive(void)
{
  int digits = 0;
  int value = 1;
  
  while(input->input_bit() == 0) {

    PB(0);
    digits++;
  }
  PB(1);
  PCol;
  for (int i = 0; i < digits; i++) {
    value = 2 * value + input->input_bit();
    PB(value & 1);
  }
  PNL;
  
#if PRINTCALLS
  cerr << "CdeltaDecode::DecodePositive()->" << value << "\n";
#endif
  
  return value;
}

/*---------------------------------------------------------------------------*/
int CdeltaDecode::DecodeNonNegative(void)
{
  int value = DecodePositive() - 1;
#if PRINTCALLS
  cerr << "CdeltaDecode::DecodeNonNegative()->" << value << "\n";
#endif
  return value;
}

/*---------------------------------------------------------------------------*/

int CdeltaDecode::Decode(void)
{
  int value;

  if (input->input_bit()) {
    value =  -DecodePositive();
  } else {
    value =  DecodeNonNegative();
  }
#if PRINTCALLS
  cerr << "CdeltaDecode::Decode()->" << value << "\n";
#endif
  return value;
}

/*---------------------------------------------------------------------------*/
#ifdef UNIT_TEST

void aout(void *c, void *data, int n)
{
  BitIn *bi = new BitIn((unsigned char *)data, n);
  CdeltaDecode *cde = new CdeltaDecode(bi);
  int i = cde->DecodePositive();
  cout << "count = " << i << "\n";
  for (int j = 0; j < i; j++) {
    cout << cde->Decode() << "\n";
  }
}

/*---------------------------------------------------------------------------*/

void main()
{
  int i;
  FILTER *f = new FILTER(0, (FilterOutput *)aout);
  BitOut *bo = new BitOut(f);
  CdeltaEncode *cde = new CdeltaEncode(bo);
  
  cde->EncodePositive(11);
  for (i = -5; i <= 5; i++) {
    cde->Encode(i);
  }
  bo->flush();
}

#endif
/*---------------------------------------------------------------------------*/




