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
#ifndef _GLOBAL_
#define _GLOBAL_
/*---------------------------------------------------------------------------*/
// global parameters
/*---------------------------------------------------------------------------*/
//#define DEBUG
// Use PGM images as default (comment this line out to use Raw images)
#define PGM
//#define DOUBLE_REAL
#define FLOAT_REAL

#ifdef DOUBLE_REAL
typedef double Real;
#endif
#ifdef FLOAT_REAL
typedef float Real;
#endif
/*---------------------------------------------------------------------------*/
// standard #defines
/*---------------------------------------------------------------------------*/
#define TRUE  1
#define FALSE 0

#define BACKSPACE 8
#define BS        8
#define ESC       27
/*---------------------------------------------------------------------------*/
// useful constants
/*---------------------------------------------------------------------------*/
#ifdef DOUBLE_REAL
const Real eps = 1.e-15;
const double MaxReal = 1.0e+100;
#endif
#ifdef FLOAT_REAL
const Real eps = 1.e-8;
const double MaxReal = 1.0e+100;
#endif

const Real Pi = (Real)3.14159265358979;
const Real TwoPi = 2.0 * Pi;
const Real Sqrt2 = (Real)sqrt(2.0);
const Real Log2 = (Real)log(2.0);

/*---------------------------------------------------------------------------*/
// helpful inline functions
/*---------------------------------------------------------------------------*/
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))

/*---------------------------------------------------------------------------*/
inline Real mod (Real x, Real N) {
   Real xmodN = x - N*((int)(x/N));
   if (xmodN < 0) xmodN += N;
   return xmodN;
}

/*---------------------------------------------------------------------------*/
inline Real square (Real x) { return (x*x); }
/*---------------------------------------------------------------------------*/
inline int  isquare (int x) { return (x*x); }
/*---------------------------------------------------------------------------*/
inline int  sign (Real x)   { return (x > 0 ? 1 : x < 0 ? -1 : 0); }
/*---------------------------------------------------------------------------*/
inline int log2 (int x) {
   int count = 0;

   while (x > 1)  {
      x >>= 1;
      count++;
   }
   return count;
}

/*---------------------------------------------------------------------------*/
// functions in global.cc
/*---------------------------------------------------------------------------*/
void init ();
void shut_down ();
volatile void error (char *format, ...);
void warning (char *format, ...);
void no_more_memory ();
/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
