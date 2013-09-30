/*
    From NR book: 

    The routine arcmak constructs the cumulative frequency distribution table used
    to partition the interval at each stage. In the principal routine arcode, when an
    interval of size jdif is to be partitioned in the proportions of some n to some ntot,
    say, then we must compute (n*jdif)/ntot. With integer arithmetic, the numerator
    is likely to overflow; and, unfortunately, an expression like jdif/(ntot/n) is not
    equivalent. In the implementation below, we resort to double precision floating
    arithmetic for this calculation. Not only is this inefficient, but different roundoff
    errors can (albeit very rarely) make different machines encode differently, though any
    one type of machine will decode exactly what it encoded, since identical roundoff
    errors occur in the two processes. For serious use, one needs to replace this floating
    calculation with an integer computation in a double register (not available to the
    C programmer).

    The internally set variable minint, which is the minimum allowed number
    of discrete steps between the upper and lower bounds, determines when new lowsignificance
    digits are added. minint must be large enough to provide resolution of
    all the input characters. That is, we must have pi × minint > 1 for all i. A value
    of 100Nch, or 1.1/ min pi, whichever is larger, is generally adequate. However, for
    safety, the routine below takes minint to be as large as possible, with the product
    minint*nradd just smaller than overflow. This results in some time inefficiency,
    and in a few unnecessary characters being output at the end of a message. You can
    decrease minint if you want to live closer to the edge.

    A final safety feature in arcmak is its refusal to believe zero values in the table
    nfreq; a 0 is treated as if it were a 1. If this were not done, the occurrence in a
    message of a single character whose nfreq entry is zero would result in scrambling
    the entire rest of the message. If you want to live dangerously, with a very slightly
    more efficient coding, you can delete the IMAX( ,1) operation.
*/

#define NRANSI
#include "nrutil.h"
#include <limits.h>
#define MC 512
#ifdef ULONG_MAX
#define MAXINT (ULONG_MAX >> 1)
#else
#define MAXINT 2147483647
#endif

typedef struct {
    unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;

/*
    Given a table nfreq[1..nchh] of the frequency of occurrence of nchh symbols, and given
    a desired output radix nradd, initialize the cumulative frequency table and other variables for
    arithmetic compression in the structure acode.
*/

void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
            arithcode *acode)
{
    unsigned long j;

    if (nchh > MC) nrerror("input radix may not exceed MC in arcmak.");
    if (nradd > 256) nrerror("output radix may not exceed 256 in arcmak.");

    acode->minint=MAXINT/nradd;
    acode->nch=nchh;
    acode->nrad=nradd;
    acode->ncumfq[1]=0;
    for (j=2;j<=acode->nch+1;j++)
        acode->ncumfq[j]=acode->ncumfq[j-1]+IMAX(nfreq[j-1],1);
    acode->ncum=acode->ncumfq[acode->nch+2]=acode->ncumfq[acode->nch+1]+1;
}
#undef MC
#undef MAXINT
#undef NRANSI
