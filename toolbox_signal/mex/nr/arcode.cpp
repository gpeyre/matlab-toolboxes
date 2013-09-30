#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include "nrutil.h"
#define NWK 20
#define JTRY(j,k,m) ((long)((((double)(k))*((double)(j)))/((double)(m))))

typedef unsigned char uchar;
typedef unsigned long ulong;

typedef struct {
    ulong *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;

/*

    Compress (isign = 1) or decompress (isign = -1) the single character ich into or out of
    the character array *codep[1..lcode], starting with byte *codep[lcd] and (if necessary)
    incrementing lcd so that, on return, lcd points to the first unused byte in *codep. Note
    that the structure acode contains both information on the code, and also state information on
    the particular output being written into the array *codep. An initializing call with isign=0
    is required before beginning any *codep array, whether for encoding or decoding. This is in
    addition to the initializing call to arcmak that is required to initialize the code itself. A call
    with ich=nch (as set in arcmak) has the reserved meaning “end of message.”

*/

void arcode(ulong *ich, uchar **codep, ulong *lcode,
            ulong *lcd, int isign, arithcode *acode)
{
    void arcsum(ulong iin[], ulong iout[], ulong ja,
        int nwk, ulong nrad, ulong nc);
    void nrerror(char error_text[]);
    int j,k;
    ulong ihi,ja,jh,jl,m;

    if( !isign ) 
    {
        acode->jdif=acode->nrad-1;
        for (j=NWK;j>=1;j--) 
        {
            acode->iupb[j]=acode->nrad-1;
            acode->ilob[j]=0;
            acode->nc=j;
            if( acode->jdif > acode->minint ) 
                return;
            acode->jdif=(acode->jdif+1)*acode->nrad-1;
        }
        nrerror("NWK too small in arcode.");
    } 
    else 
    {
        if (isign > 0) 
        {
            if( *ich>acode->nch ) 
                nrerror("bad ich in arcode.");
        }
        else 
        {
            ja=(*codep)[*lcd]-acode->ilob[acode->nc];
            for( j=acode->nc+1;j<=NWK;j++ ) 
            {
                ja *= acode->nrad;
                ja += ((*codep)[*lcd+j-acode->nc]-acode->ilob[j]);
            }
            ihi=acode->nch+1;
            *ich=0;
            while (ihi-(*ich) > 1) 
            {
                m=(*ich+ihi)>>1;
                if( ja>=JTRY(acode->jdif,acode->ncumfq[m+1],acode->ncum) )
                    *ich=m;
                else 
                    ihi=m;
            }
            if (*ich == acode->nch) return;
        }
        jh=JTRY(acode->jdif,acode->ncumfq[*ich+2],acode->ncum);
        jl=JTRY(acode->jdif,acode->ncumfq[*ich+1],acode->ncum);
        acode->jdif=jh-jl;
        arcsum(acode->ilob,acode->iupb,jh,NWK,acode->nrad,acode->nc);
        arcsum(acode->ilob,acode->ilob,jl,NWK,acode->nrad,acode->nc);
        for (j=acode->nc;j<=NWK;j++) {
            if (*ich != acode->nch && acode->iupb[j] != acode->ilob[j]) break;
            if (*lcd > *lcode) {
                fprintf(stderr,"Reached the end of the 'code' array.\n");
                fprintf(stderr,"Attempting to expand its size.\n");
                ulong old_size = *lcode;
                *lcode += *lcode/2;
                // uchar* new_array = new uchar[*lcode];
                uchar* new_array = cvector(0,*lcode);
                memcpy( new_array, *codep, old_size*sizeof(uchar)+1 );
                // delete [] (*codep+1);
                free_cvector(*codep,0,2048);
                *codep = new_array;
#if 0
                *codep = (uchar *) realloc(*codep, (unsigned)(*lcode*sizeof(uchar)));
#endif
                if( *codep==NULL ) 
                {
                        nrerror("Size expansion failed");
                }
            }
            if( isign>0 ) 
                (*codep)[*lcd]=(uchar)acode->ilob[j];
            ++(*lcd);
        }
        if( j>NWK ) 
            return;
        acode->nc=j;
        for( j=0;acode->jdif<acode->minint;j++ )
            acode->jdif *= acode->nrad;
        if (acode->nc-j < 1) 
            nrerror("NWK too small in arcode.");
        if( j ) 
        {
            for( k=acode->nc;k<=NWK;k++ ) 
            {
                acode->iupb[k-j]=acode->iupb[k];
                acode->ilob[k-j]=acode->ilob[k];
            }
        }
        acode->nc -= j;
        for (k=NWK-j+1;k<=NWK;k++) 
            acode->iupb[k]=acode->ilob[k]=0;
    }
    return;
}
#undef NWK
#undef JTRY
