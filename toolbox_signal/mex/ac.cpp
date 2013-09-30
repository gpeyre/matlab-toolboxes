/*--------------------------------------------------------------*/
/* Prototype de compression d'images N&B                        */
/* (c) 2002-2004 Christophe Bernard, Let It Wave                */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "ac.h"
#include <math.h>

#define Code_value_bits 16

#define Top_value (((long)1<<Code_value_bits)-1)
#define First_qtr (Top_value/4+1)
#define Half	  (2*First_qtr)
#define Third_qtr (3*First_qtr)
#define Max_frequency 16383

static void output_bit (ac_encoder *, int);
static void bit_plus_follow (ac_encoder *, int);
static int input_bit (ac_decoder *);
static void update_model (ac_model *, int, int);

static void output_interval(ac_encoder *);
static void input_interval(ac_decoder *);

/*static void flatten_model (ac_model *);*/

int ac_return_status=0;

#define error(m) fprintf(stderr, m)

#define check(b,m) if( !(b) ) error(m)

void fallintotrap() {}

static void
output_bit (ac_encoder *ace, int bit)
{
  ace->buffer >>= 1;
  if (bit)
    ace->buffer |= 0x80;
  ace->bits_to_go -= 1;
  ace->total_bits += 1;
  if (ace->bits_to_go==0)  {
    if (ace->fp)
      putc (ace->buffer, ace->fp);
    ace->bits_to_go = 8;
  }

  return;
}

static void
bit_plus_follow (ac_encoder *ace, int bit)
{
  output_bit (ace, bit);
  while (ace->fbits > 0)  {
    output_bit (ace, !bit);
    ace->fbits -= 1;
  }

  return;
}

static int
input_bit (ac_decoder *acd)
{
  int t;

  if (acd->bits_to_go==0)  {
    acd->buffer = getc(acd->fp);
    if (acd->buffer==EOF)  {
      acd->garbage_bits += 1;
      if (acd->garbage_bits>Code_value_bits-2)
        error ("arithmetic decoder bad input file");
    }
    acd->bits_to_go = 8;
  }

  t = acd->buffer&1;
  acd->buffer >>= 1;
  acd->bits_to_go -= 1;

  return t;
}

size_t
ac_rawput(ac_encoder *ace, char *data, size_t size)
{
  if(ace->fp) return fwrite(data, size, 1, ace->fp);
  else return size;
}

size_t
ac_rawget(ac_decoder *acd, char *data, size_t size)
{
  return fread(data, size, 1, acd->fp);
}

static void
update_model (ac_model *acm, int sym, int flatten_mod)
{
  int i;

  if ((acm->cfreq[0]==Max_frequency)||(flatten_mod))  {
    int cum = 0;
    acm->cfreq[acm->nsym] = 0;
    for (i = acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = (acm->freq[i] + 1) / 2;
      cum += acm->freq[i];
      acm->cfreq[i] = cum;
    }
  }

if ((flatten_mod))  {
    int cum = 0;
    acm->cfreq[acm->nsym] = 0;
    for (i = acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = (acm->freq[i] + 1) / 2;
      cum += acm->freq[i];
      acm->cfreq[i] = cum;
    }
  }

if ((flatten_mod))  {
    int cum = 0;
    acm->cfreq[acm->nsym] = 0;
    for (i = acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = (acm->freq[i] + 1) / 2;
      cum += acm->freq[i];
      acm->cfreq[i] = cum;
    }
  }

  if ((flatten_mod))  {
    int cum = 0;
    acm->cfreq[acm->nsym] = 0;
    for (i = acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = (acm->freq[i] + 1) / 2;
      cum += acm->freq[i];
      acm->cfreq[i] = cum;
    }
  }

  acm->freq[sym] += 1;
  for (i=sym; i>=0; i--)
    acm->cfreq[i] += 1;

  return;
}

/*
static void
flatten_model (ac_model *acm)
{
  int i;
  int cum = 0;

    acm->cfreq[acm->nsym] = 0;
    for (i = acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = (acm->freq[i] + 1) / 2;
      cum += acm->freq[i];
      acm->cfreq[i] = cum;
      }
   return;
}*/

void
ac_encoder_init (ac_encoder *ace, const char *fn)
{

  if (fn)  {
    ace->fp = fopen (fn, "wb"); /* open in binary mode */
    check (!ace->fp, "arithmetic encoder could not open file");
  }  else  {
    ace->fp = NULL;
  }

  ace->bits_to_go = 8;

  ace->low = 0;
  ace->high = Top_value;
  ace->fbits = 0;
  ace->buffer = 0;

  ace->total_bits = 0;

  return;
}

void
ac_encoder_done (ac_encoder *ace)
{
  ace->fbits += 1;
  if (ace->low < First_qtr)
    bit_plus_follow (ace, 0);
  else
    bit_plus_follow (ace, 1);
  if (ace->fp)
    putc (ace->buffer >> ace->bits_to_go, ace->fp);

  if (ace->fp)
    fclose (ace->fp);

  return;
}

void
ac_decoder_init (ac_decoder *acd, const char *fn)
{
  int i;

  acd->fp = fopen (fn, "rb"); /* open in binary mode */
  check (!acd->fp, "arithmetic decoder could not open file");

  acd->bits_to_go = 0;
  acd->garbage_bits = 0;

  acd->value = 0;
  for (i=1; i<=Code_value_bits; i++)  {
    acd->value = 2*acd->value + input_bit(acd);
  }
  acd->low = 0;
  acd->high = Top_value;

  return;
}

void
ac_decoder_done (ac_decoder *acd)
{
  fclose (acd->fp);

  return;
}

void
ac_model_init (ac_model *acm, int nsym, int *ifreq, int adapt)
{
  int i;

  acm->nsym = nsym;
  acm->freq = (int *) (void *) calloc (nsym, sizeof (int));
  check (!acm->freq, "arithmetic coder model allocation failure");
  acm->cfreq = (int *) (void *) calloc (nsym+1, sizeof (int));
  check (!acm->cfreq, "arithmetic coder model allocation failure");
  acm->adapt = adapt;

  if(ifreq)  
  {
    acm->cfreq[acm->nsym] = 0;
    for (i=acm->nsym-1; i>=0; i--)  {
      acm->freq[i] = ifreq[i];
      acm->cfreq[i] = acm->cfreq[i+1] + acm->freq[i];
    }
    if (acm->cfreq[0] > Max_frequency)
      error ("arithmetic coder model max frequency exceeded");
  }  
  else  
  {
    for (i=0; i<acm->nsym; i++) {
      acm->freq[i] = 1;
      acm->cfreq[i] = acm->nsym - i;
    }
    acm->cfreq[acm->nsym] = 0;
  }

  return;
}

void
ac_model_done (ac_model *acm)
{
  acm->nsym = 0;
  free (acm->freq);
  acm->freq = NULL;
  free (acm->cfreq);
  acm->cfreq = NULL;

  return;
}


double
ac_encoder_bits (ac_encoder *ace)
{
  return ace->total_bits+ace->fbits+
    log(Top_value/(double)(ace->high-ace->low))/log(2);
}

void
ac_encode_symbol (ac_encoder *ace, ac_model *acm, int sym, int flatten_mod)
{
  long range;

  check( sym<0||sym>=acm->nsym, "symbol out of range");

  range = (long)(ace->high-ace->low)+1;
  ace->high = ace->low + (range*acm->cfreq[sym])/acm->cfreq[0]-1;
  ace->low = ace->low + (range*acm->cfreq[sym+1])/acm->cfreq[0];

  for (;;)  {
    if (ace->high<Half)  {
      bit_plus_follow (ace, 0);
    }  else if (ace->low>=Half)  {
      bit_plus_follow (ace, 1);
      ace->low -= Half;
      ace->high -= Half;
    }  else if (ace->low>=First_qtr && ace->high<Third_qtr)  {
      ace->fbits += 1;
      ace->low -= First_qtr;
      ace->high -= First_qtr;
    }  else
      break;
    ace->low = 2*ace->low;
    ace->high = 2*ace->high+1;
  }

  if (acm->adapt)
    update_model (acm, sym, flatten_mod);

  return;
}

int
ac_decode_symbol (ac_decoder *acd, ac_model *acm,  int flatten_mod)
{
  long range;
  int cum;
  int sym;


  range = (long)(acd->high-acd->low)+1;
  cum = (((long)(acd->value-acd->low)+1)*acm->cfreq[0]-1)/range;

  for (sym = 0; acm->cfreq[sym+1]>cum; sym++)
    /* do nothing */ ;

  check (sym<0||sym>=acm->nsym, "symbol out of range");

  acd->high = acd->low + (range*acm->cfreq[sym])/acm->cfreq[0]-1;
  acd->low = acd->low +  (range*acm->cfreq[sym+1])/acm->cfreq[0];

  for (;;)  {
    if (acd->high<Half)  {
      /* do nothing */
    }  else if (acd->low>=Half)  {
      acd->value -= Half;
      acd->low -= Half;
      acd->high -= Half;
    }  else if (acd->low>=First_qtr && acd->high<Third_qtr)  {
      acd->value -= First_qtr;
      acd->low -= First_qtr;
      acd->high -= First_qtr;
    }  else
      break;
    acd->low = 2*acd->low;
    acd->high = 2*acd->high+1;
    acd->value = 2*acd->value + input_bit(acd);
  }

  if (acm->adapt)
    update_model (acm, sym, flatten_mod);

  return sym;
}

void
ac_encode_uniform(ac_encoder *ace, int nsym, int sym)
{
  int range;
  
  check(sym<0 || sym>=nsym, "symbol out of range");

  range = (int)(ace->high-ace->low)+1;
  ace->high = ace->low + (range*(nsym-sym))/nsym-1;
  ace->low = ace->low + (range*(nsym-1-sym))/nsym;

  output_interval(ace);
}

int
ac_decode_uniform(ac_decoder *acd, int nsym)
{
  int range;
  int cum;
  int sym;

  range = (int)(acd->high-acd->low)+1;
  cum = (((int)(acd->value-acd->low)+1)*nsym-1)/range;

  sym= nsym-cum-1; 

  check( sym<0 || sym>=nsym , "symbol out of range");

  acd->high = acd->low + (range*(nsym-sym))/nsym-1;
  acd->low = acd->low +  (range*(nsym-1-sym))/nsym;

  input_interval(acd);
  return sym;
}

void 
ac_encode_bytes(ac_encoder *ace, int size, void *data)
{
  unsigned char* c= (unsigned char*) data;
  int i;

  for (i= 0; i<size; i++){
    ac_encode_uniform(ace, 256, c[i]);
  }
}

void 
ac_decode_bytes(ac_decoder *acd, int size, void *data)
{
  unsigned char *c= (unsigned char*) data;
  int i;

  for (i= 0; i<size; i++){
    c[i]= ac_decode_uniform(acd, 256);
  }
}
/*--------------------------------------------------------------*/
/* Manipulation d'intervalles                                   */
/*--------------------------------------------------------------*/

static void
input_interval(ac_decoder *acd)
{
  check(acd->high<=acd->low, "decoder state error");

  for (;;)  {
    if (acd->high<Half)  {
      /* do nothing */
    }  else if (acd->low>=Half)  {
      acd->value -= Half;
      acd->low -= Half;
      acd->high -= Half;
    }  else if (acd->low>=First_qtr && acd->high<Third_qtr)  {
      acd->value -= First_qtr;
      acd->low -= First_qtr;
      acd->high -= First_qtr;
    }  else
      break;
    acd->low = 2*acd->low;
    acd->high = 2*acd->high+1;
    acd->value = 2*acd->value + input_bit(acd);
  }
}

void output_interval(ac_encoder *ace)
{
  check(ace->high<=ace->low, "encoder state error");
  check(ace->high>65536, "encoder state error");
  check(ace->low<0, "encoder state error");

  for(;;)  
  {
    if (ace->high<Half)  
    {
      bit_plus_follow (ace, 0);
    }  
    else if(ace->low>=Half)  
    {
      bit_plus_follow (ace, 1);
      ace->low -= Half;
      ace->high -= Half;
    }  
    else if(ace->low>=First_qtr && ace->high<Third_qtr)  
    {
      ace->fbits += 1;
      ace->low -= First_qtr;
      ace->high -= First_qtr;
    }  
    else
      break;
    ace->low = 2*ace->low;
    ace->high = 2*ace->high+1;
  }
}
