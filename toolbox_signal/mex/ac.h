/*--------------------------------------------------------------*/
/* Prototype de compression d'images N&B                        */
/* (c) 2002-2004 Christophe Bernard, Let It Wave                */
/*--------------------------------------------------------------*/

#ifndef AC_H
#define AC_H

#include <stdio.h>

typedef struct {
  FILE *fp;
  long low;
  long high;
  long fbits;
  int buffer;
  int bits_to_go;
  long total_bits;
} ac_encoder;

typedef struct {
  FILE *fp;
  long value;
  long low;
  long high;
  int buffer;
  int bits_to_go;
  int garbage_bits;
} ac_decoder;

typedef struct {
  int nsym;
  int *freq;
  int *cfreq;
  int adapt;
} ac_model;

void ac_encoder_init (ac_encoder *, const char *);
void ac_encoder_done (ac_encoder *);
void ac_decoder_init (ac_decoder *, const char *);
void ac_decoder_done (ac_decoder *);
void ac_model_init (ac_model *, int, int *, int);
void ac_model_done (ac_model *);
double ac_encoder_bits (ac_encoder *);
void ac_encode_symbol (ac_encoder *, ac_model *, int, int);
int ac_decode_symbol (ac_decoder *, ac_model *, int);
size_t ac_rawput(ac_encoder *ace, char *data, size_t size);
size_t ac_rawget(ac_decoder *acd, char *data, size_t size);
void ac_encode_uniform(ac_encoder *ace, int nsym, int sym);
int ac_decode_uniform(ac_decoder *acd, int nsym);
void ac_encode_bytes(ac_encoder *ace, int size, void *data);
void ac_decode_bytes(ac_decoder *acd, int size, void *data);

#endif
