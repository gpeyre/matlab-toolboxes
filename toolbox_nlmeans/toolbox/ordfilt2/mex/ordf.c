/* Copyright 1993-2003 The MathWorks, Inc. */

/*
 * ORDF.MEX
 *
 * B = ORDFA(A, ORDER, OFFSETS, S, [STARTROW STARTCOL], [MB NB])
 * A        --- input matrix, can be 2-D real numeric, or logical
 * ORDER    --- order-statistic to compute
 * OFFSETS  --- linear index offsets that form neighborhood
 * S        --- additive offsets; this is an optional input; 
 *              A must be double if this input is specified
 * STARTROW --- which pixel in A to start from (one-based)
 * STARTCOL --- which pixel in A to start from (one-based)
 * MB       --- row dimension of output
 * NB       --- column dimension of output
 *
 * A must be padded; no boundary handling is included.  That's why
 * STARTROW and STARTCOL must be provided and why B is not the same
 * size as A.
 *
 * See ORDFILT2.M
 */

#include "mex.h"

static char rcsid[] = "$Revision: 1.2.6.4 $";


void ValidateInputs(int nlhs, 
                    mxArray *plhs[], 
                    int nrhs, 
                    const mxArray *prhs[], 
                    const mxArray **A, 
                    int *startRow, 
                    int *startCol,
                    int *order, 
                    int **offsets, 
                    int *numOffsets,
		    double **add,
                    int *Mb, 
                    int *Nb,
                    int *Md,
                    int *Nd)
{
    const mxArray *orderArray;
    const mxArray *offsetsArray;
    const mxArray *startsArray;
    const mxArray *outsizeArray;
    const mxArray *domainsizeArray;
    const mxArray *sArray;
    double *p;
    int k;
    int idx;
    int lastIdx;
    int testIdx;
    int Ma;
    


    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("Images:ordf:tooManyOutputs",
                          "%s",
                          "Too many outputs.");
    }

    if (nrhs < 6)
    {
        mexErrMsgIdAndTxt("Images:ordf:toofewInputs",
                          "%s",
                          "Too few inputs.");
    }
    if (nrhs > 7)
    {
        mexErrMsgIdAndTxt("Images:ordf:tooManyInputs",
                          "%s",
                          "Too many inputs.");
    }

    *A = prhs[0];
    orderArray      = prhs[1];
    offsetsArray    = prhs[2];
    startsArray     = prhs[3];
    outsizeArray    = prhs[4];
    domainsizeArray = prhs[5];


    if ( nrhs == 6)
    {
	/* Additive offsets are not used in this case */
	*add = NULL;
    }
    else /* additive offsets were specified */
    {
	sArray = prhs[6];

	*numOffsets = mxGetNumberOfElements(offsetsArray);

        /* Get Additive offsets */
	if (mxGetNumberOfElements(sArray) != *numOffsets)
	{
            mexErrMsgIdAndTxt("Images:ordf:numElementsNotEqualForOFFSETSAndS",
                              "%s",
                              "S must have the same"
                              " number of elements as OFFSETS.");
	}
	*add = mxGetPr(sArray);
    }

    /* Get order */
    *order = (int) mxGetScalar(orderArray);
    *order -= 1;  /* change from one-based to zero-based */
    
    /* Get offsets */
    p = mxGetPr(offsetsArray);
    *numOffsets = mxGetNumberOfElements(offsetsArray);
    *offsets = (int *) mxCalloc(*numOffsets, sizeof(int));
    for (k = 0; k < *numOffsets; k++)
    {
        (*offsets)[k] = (int) p[k];
    }
        
    /* Get starts */
    if (!mxIsDouble(startsArray) || (mxGetNumberOfElements(startsArray) != 2))
    {
        mexErrMsgIdAndTxt("Images:ordf:STARTMustBe2elemDoubleVector",
                          "%s",
                          "START must be a 2-element double vector.");
    }
    p = mxGetPr(startsArray);
    *startRow = (int) p[0] - 1; /* convert from one-based to zero-based */
    *startCol = (int) p[1] - 1;
    if ((((double) *startRow) != (p[0] - 1)) ||
        (((double) *startCol) != (p[1] - 1)))
    {
        mexErrMsgIdAndTxt("Images:ordf:STARTMustContainIntegers",
                          "%s",
                          "START must contain integers.");
    }
    
    /* Get output size */
    if (!mxIsDouble(outsizeArray) || 
        (mxGetNumberOfElements(outsizeArray) != 2))
    {
        mexErrMsgIdAndTxt("Images:ordf:OUTSIZEMustBe2elemDoubleVector",
                          "%s",
                          "OUTSIZE must be a 2-element double vector.");
    }
    p = mxGetPr(outsizeArray);
    *Mb = (int) p[0];
    *Nb = (int) p[1];
    if ((((double) *Mb) != p[0]) ||
        (((double) *Nb) != p[1]))
    {
        mexErrMsgIdAndTxt("Images:ordf:OUTSIZEMustContainIntegers",
                          "%s",
                          "OUTSIZE must contain integers.");
    }

    /* Get domain size */
    if (!mxIsDouble(domainsizeArray) || 
        (mxGetNumberOfElements(domainsizeArray) != 2))
    {
        mexErrMsgIdAndTxt("Images:ordf:DOMAINSIZEMustBe2elemDoubleVector",
                          "%s",
                          "DOMAIN must be a 2-element double vector.");
    }
    p = mxGetPr(domainsizeArray);
    *Md = (int) p[0];
    *Nd = (int) p[1];
    if ((((double) *Md) != p[0]) ||
        (((double) *Nd) != p[1]))
    {
        mexErrMsgIdAndTxt("Images:ordf:DOMAINSIZEMustContainIntegers",
                          "%s",
                          "DOMAINSIZE must contain integers.");
    }


    /* Consistency checks */
    if ((*startRow < 0) ||
        (*startCol < 0) ||
        (*Mb < 0) ||
        (*Nb < 0) ||
        (*startRow >= mxGetM(*A)) ||
        (*startCol >= mxGetN(*A)) ||
        (*startRow + *Mb - 1 >= mxGetM(*A)) ||
        (*startCol + *Nb - 1 >= mxGetN(*A)))
    {
        mexErrMsgIdAndTxt("Images:ordf:invalidInput",
                          "%s",
                          "Invalid input.");
    }
    
    /* Check for seg-v bait */
    Ma = mxGetM(*A);
    idx = *startCol * Ma + *startRow;
    lastIdx = Ma * mxGetN(*A) - 1;
    for (k = 0; k < *numOffsets; k++)
    {
        testIdx = idx + (*offsets)[k];
        if ((testIdx < 0) || (testIdx > lastIdx))
        {
            mexErrMsgIdAndTxt("Images:ordf:outOfRangeOFFSETS",
                              "%s",
                              "Out of range value in OFFSETS.");
        }
    }
}
/* select_uint8
 *
 * Selection algorithm using QuickSort-style partitioning. 
 * It is based on work presented in Sedgewick, Algorithms, 2d ed, 
 * 1988, pp. 115-130. This function modifies the values in a[].
 */

#define TYPE uint8_T
TYPE select_uint8
#include "ordf_select.h"

#define TYPE uint16_T
TYPE select_uint16
#include "ordf_select.h"

#define TYPE uint32_T
TYPE select_uint32
#include "ordf_select.h"

#define TYPE int8_T
TYPE select_int8
#include "ordf_select.h"

#define TYPE int16_T
TYPE select_int16
#include "ordf_select.h"

#define TYPE int32_T
TYPE select_int32
#include "ordf_select.h"

/* There is no need to sort the values in case of logical
 * input, since all the values are either 0 or 1 
 */
mxLogical select_logical(mxLogical *a, int N, int order)
{
    int i;
    int numZeros = N;
    /* Note that order is really order-1 (i.e. it was 
       converted to an array index for use in the other
       select_ routines */
    int trueOrder = order+1;
    
    /* Count the number of zeros in this vector */
    for(i=0; i<N; i++)
    {
	numZeros-=(int)a[i];
    }

    if(trueOrder <= numZeros)
	return false;
    else
	return true;
}

#define TYPE double
#define CHECK_NANS
TYPE select_double
#include "ordf_select.h"

#define TYPE float
#define CHECK_NANS
TYPE select_single
#include "ordf_select.h"

/* ordfilt2_uint8
 *
 * This method produces the output array using select_uint8 function.
 */ 

#define TYPE uint8_T
#define SELECT select_uint8
void ordfilt2_uint8
#include "ordf_ordfilt2.h"

#define TYPE uint16_T
#define SELECT select_uint16
void ordfilt2_uint16
#include "ordf_ordfilt2.h"

#define TYPE uint32_T
#define SELECT select_uint32
void ordfilt2_uint32
#include "ordf_ordfilt2.h"

#define TYPE int8_T
#define SELECT select_int8
void ordfilt2_int8
#include "ordf_ordfilt2.h"

#define TYPE int16_T
#define SELECT select_int16
void ordfilt2_int16
#include "ordf_ordfilt2.h"

#define TYPE int32_T
#define SELECT select_int32
void ordfilt2_int32
#include "ordf_ordfilt2.h"

#define TYPE mxLogical
#define SELECT select_logical
void ordfilt2_logical
#include "ordf_ordfilt2.h"

#define TYPE double
#define SELECT select_double
void ordfilt2_double
#include "ordf_ordfilt2.h"

#define TYPE float
#define SELECT select_single
void ordfilt2_single
#include "ordf_ordfilt2.h"

/* method which handles additive offsets */
#define TYPE double
#define SELECT select_double
#define ADD_OFFSET
void ordfilt2_add_offset
#include "ordf_ordfilt2.h"

/*
 * ordfilt_hist_uint8
 *
 * Implements histogram approach to computing order statistic.
 */
#define TYPE uint8_T
#define MAX_VAL_TYPE MAX_uint8_T
#define MIN_VAL_TYPE MIN_uint8_T
void ordfilt_hist_uint8
#include "ordf_hist.h"

#define TYPE int8_T
#define MAX_VAL_TYPE MAX_int8_T
#define MIN_VAL_TYPE MIN_int8_T
#define SIGNED
void ordfilt_hist_int8
#include "ordf_hist.h"

#define TYPE uint16_T
#define MAX_VAL_TYPE MAX_uint16_T
#define MIN_VAL_TYPE MIN_uint16_T
void ordfilt_hist_uint16
#include "ordf_hist.h"

#define TYPE int16_T
#define MAX_VAL_TYPE MAX_int16_T
#define MIN_VAL_TYPE MIN_int16_T
#define SIGNED
void ordfilt_hist_int16
#include "ordf_hist.h"


/* These defines specify the size of the domain matrix beyond which
   the histogram approach (faster for a larger domain) of calculating
   the order-statistic is used instead of quicksort style method.
   These values were determined experimentally.
*/
#define HIST8_REQ_ROWS 7
#define HIST8_REQ_COLS 1

#define HIST16_REQ_ROWS     3
#define HIST16_REQ_ELEMENTS 520

/* This method determines if it's more efficient to use the 
   histogram based filtering in case of INT8 and UINT8
   data type.
*/
bool useHistMethodForINT8(int numOffsets, int Md, int Nd)
{
    /* Invoke histogram based method only if a full domain 
       matrix was specified and it contained more than a specified
       number of rows and columns.
    */
    if( (numOffsets == Md*Nd) && 
        (Md >= HIST8_REQ_ROWS) && 
        (Nd >= HIST8_REQ_COLS) )
    {    
        return true;
    }
    else
    {
        return false;
    }
} 

/* This method determines if it's more efficient to use the 
   histogram based filtering in case of INT16 and UINT16
   data type.
*/
bool useHistMethodForINT16(int numOffsets, int Md, int Nd)
{
    /* Invoke histogram based method only if a full domain 
       matrix was specified and it contained more than a specified
       number of rows and more than a specified number of elements.
    */
    if( (numOffsets == Md*Nd) && 
        (Md >= HIST16_REQ_ROWS) && 
        (numOffsets >= HIST16_REQ_ELEMENTS) )
    {
        return true;
    }
    else
    {
        return false;
    }
} 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int startRow;  /* zero-based */
    int startCol;  /* zero-based */
    int order;     /* zero-based */
    int Mb;
    int Nb;
    int Md;
    int Nd;
    int Ma;
    mxClassID classA;
    int sizeB[2];
    int *offsets;
    int numOffsets;
    const mxArray *A;
    mxArray *B;
    double *add;

    ValidateInputs(nlhs, plhs, nrhs, prhs, &A, &startRow, &startCol,
                   &order, &offsets, &numOffsets, &add, &Mb, &Nb,
                   &Md, &Nd);

    /* Create output */
    classA = mxGetClassID(A);
    sizeB[0] = Mb;
    sizeB[1] = Nb;

    if (mxIsLogical(A))
    {
        B = mxCreateLogicalArray(2, sizeB);
    }
    else
    {
        B = mxCreateNumericArray(2, sizeB, classA, mxREAL);
    }
    
    Ma = mxGetM(A);

    switch (classA)
    {
      case mxUINT8_CLASS:
        if( useHistMethodForINT8(numOffsets, Md, Nd) )
        {
            ordfilt_hist_uint8((uint8_T *)mxGetData(A), 
                               (uint8_T *)mxGetData(B),
                               startRow, startCol, Mb, Nb,
                               Ma, order, offsets,
                               numOffsets, Md, Nd);
        }
        else
        {
            ordfilt2_uint8((uint8_T *) mxGetData(A), (uint8_T *) mxGetData(B),
                           startRow, startCol, Mb, Nb, Ma, order, offsets,
                           numOffsets);
        }
        break;
        
      case mxUINT16_CLASS:
        if( useHistMethodForINT16(numOffsets, Md, Nd) )
        {
            ordfilt_hist_uint16((uint16_T *)mxGetData(A), 
                                (uint16_T *)mxGetData(B),
                                startRow, startCol, Mb, Nb,
                                Ma, order, offsets,
                                numOffsets, Md, Nd);
        }
        else
        {
            ordfilt2_uint16((uint16_T *) mxGetData(A), 
                            (uint16_T *) mxGetData(B),
                            startRow, startCol, Mb, Nb, Ma, 
                            order, offsets, numOffsets);
        }
        
        break;

      case mxUINT32_CLASS:
        ordfilt2_uint32((uint32_T *) mxGetData(A), (uint32_T *) mxGetData(B),
                        startRow, startCol, Mb, Nb, Ma, order, offsets,
                        numOffsets);
        break;

      case mxINT8_CLASS:
        if( useHistMethodForINT8(numOffsets, Md, Nd) )
        {
            ordfilt_hist_int8((int8_T *)mxGetData(A), 
                              (int8_T *)mxGetData(B),
                              startRow, startCol, Mb, Nb,
                              Ma, order, offsets,
                              numOffsets, Md, Nd);
        }
        else
        {
            ordfilt2_int8((int8_T *) mxGetData(A), (int8_T *) mxGetData(B),
                          startRow, startCol, Mb, Nb, Ma, order, offsets,
                          numOffsets);
        }
        break;
        
      case mxINT16_CLASS:
        if( useHistMethodForINT16(numOffsets, Md, Nd) )
        {
            ordfilt_hist_int16((int16_T *)mxGetData(A), 
                               (int16_T *)mxGetData(B),
                               startRow, startCol, Mb, Nb,
                               Ma, order, offsets,
                               numOffsets, Md, Nd);
        }
        else
        {
            ordfilt2_int16((int16_T *) mxGetData(A), 
                           (int16_T *) mxGetData(B),
                           startRow, startCol, Mb, Nb, Ma, 
                           order, offsets, numOffsets);
        }

        break;

      case mxINT32_CLASS:
        ordfilt2_int32((int32_T *) mxGetData(A), (int32_T *) mxGetData(B),
                       startRow, startCol, Mb, Nb, Ma, order, offsets,
                       numOffsets);
        break;

      case mxLOGICAL_CLASS:
        ordfilt2_logical(mxGetLogicals(A), mxGetLogicals(B),
                         startRow, startCol, Mb, Nb, Ma, order, offsets,
                         numOffsets);
        break;
      
      case mxDOUBLE_CLASS:
        if(add==NULL)
	{
	    ordfilt2_double((double *) mxGetData(A), (double *) mxGetData(B),
			    startRow, startCol, Mb, Nb, Ma, order, offsets,
			    numOffsets);
	}
	else /* additive offsets were specified */
	{   
	    ordfilt2_add_offset((double *) mxGetData(A), 
                                (double *) mxGetData(B),
				startRow, startCol, Mb, Nb, 
                                Ma, order, offsets,
				numOffsets, add);
	}
	break;
        
      case mxSINGLE_CLASS:
        ordfilt2_single((float *) mxGetData(A), (float *) mxGetData(B),
                        startRow, startCol, Mb, Nb, Ma, order, offsets,
                        numOffsets);
        break;
  
      default:
        /* Should have errored in ValidateInputs() */
        mexErrMsgIdAndTxt("Images:ordf:invalidInputClass",
                          "%s",
                          "Invalid input class.");
    }

    mxFree((void *) offsets);

    plhs[0] = B;
}
