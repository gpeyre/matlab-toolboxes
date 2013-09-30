/*
 * Copyright 1993-2003 The MathWorks, Inc.
 * $Revision: 1.1.6.3 $  $Date: 2003/01/26 06:00:45 $
 *
 * This file contains a function body for order filtering
 * using histogram approach. In the case of large filter 
 * kernel, it is faster than the quick sort function.
 *
 * The origin of the idea for this type of approach comes
 * from the following paper:
 * A fast two-dimensional median filtering algorithm
 * T.S.Huang, G.J.Yang, G.Y.Tang
 * IEEE transactions on acoustics, speech and signal 
 * processing Vol ASSP 27 No 1, February 1979.
 *
 * NOTE:  This method is designed only for full domain 
 *        matrices.
 */

#define VAL_RANGE (MAX_VAL_TYPE-MIN_VAL_TYPE+1)

(TYPE *pA, TYPE *pB,
 int startRow, int startCol, 
 int Mb, int Nb, int Ma, int order, 
 int *offsets, int numOffsets,
 int Md, int Nd)
{
    int    *upperRowOffsets, *lowerRowOffsets;
    int     col, row;
    int     k;
    double *hist;
    double  sum;
    double  numPixBelowTrackedIndex, numPixAboveTrackedIndex;
    TYPE   *p;
    int     pixIndex, pixIndexTracked;

    /* This method should not be used at all if there are
       fewer than 3 rows of data in the domain matrix */
    if(Md < 3)
    {
        mexErrMsgIdAndTxt("Images:ordf:tooFewRowsInDomainForHistMethod",
                          "%s",
                          "There should be at least three rows in the "
                          "domain matrix in order to run histogram based "
                          "algorithm.");
    }

    order++; /* convert to 1 based index */

    /* Allocate space for the histogram */
    hist = mxCalloc(VAL_RANGE, sizeof(double));

    /* Allocate space for first and last rows */
    upperRowOffsets = mxCalloc(Nd, sizeof(*upperRowOffsets));
    lowerRowOffsets = mxCalloc(Nd, sizeof(*lowerRowOffsets));

    /* Store the offsets */
    for (k=0; k < Nd; k++)
    {
        upperRowOffsets[k] = offsets[k*Md];
        lowerRowOffsets[k] = offsets[(k+1)*Md-1];
    }

    for (col = 0; col < Nb; col++)
    {
        p = pA + (startCol+col)*Ma + startRow;

        /* Compute the histogram for first location of 
         * every column
         */
        memset((void *)hist, 0, VAL_RANGE*sizeof(double));
        for (k = 0; k < numOffsets; k++)
        {
            pixIndex = *(p + offsets[k]);
#ifdef SIGNED  /* switch from pixel values to 0 based index */
            pixIndex -= MIN_VAL_TYPE;
#endif
            hist[pixIndex]++;
        }

        /* traverse the histogram looking for the result */
        for (pixIndex=0, sum=0.0; pixIndex < VAL_RANGE; pixIndex++)
        {
            numPixBelowTrackedIndex = sum;
            sum += hist[pixIndex];
            if(sum >= order) { break; }
        }

        for (row = 0; row < Mb; row++)
        {
            /* place result in the output array and move on */
#ifdef SIGNED
            *pB++ = (TYPE)(pixIndex+MIN_VAL_TYPE);
#else
            *pB++ = (TYPE)pixIndex;
#endif

            /* store the current index of the pixel value in 
               the position specified by order
            */
            pixIndexTracked = pixIndex;

            /* Remove upper row results from the histogram */
            for (k = 0; k < Nd; k++)
            {
                pixIndex = *(p + upperRowOffsets[k]);
#ifdef SIGNED
                pixIndex -= MIN_VAL_TYPE;
#endif
                hist[pixIndex]--;

                /* Maintain the count of pixels below the 
                   pixel specified by the order */
                if(pixIndex < pixIndexTracked)
                {
                    numPixBelowTrackedIndex--;
                }
            }

            /* point to the next portion of the image */
            p++;

            /* Add lower row to the histogram */
            for (k = 0; k < Nd; k++)
            {
                pixIndex = *(p + lowerRowOffsets[k]);
#ifdef SIGNED
                pixIndex -= MIN_VAL_TYPE;
#endif
                hist[pixIndex]++;

                if(pixIndex < pixIndexTracked)
                {
                    numPixBelowTrackedIndex++;
                }
            }

            /* Compute the next result, but don't traverse the
               entire histogram.  Instead, use the information
               gained by counting the number of pixels below 
               the pixel specified by the order */
            numPixAboveTrackedIndex = numPixBelowTrackedIndex+
                hist[pixIndexTracked];

            if(order <= numPixBelowTrackedIndex)
            {
                /* travel to the left since:
                 * order <= numPixBelowTrackedIndex
                 */ 

                /* traverse the histogram looking for the result */
                for (pixIndex=pixIndexTracked-1; pixIndex >= 0; pixIndex--)
                {
                    numPixBelowTrackedIndex -= hist[pixIndex];
                    if(numPixBelowTrackedIndex < order) { break; }
                }
            }
            else if(order > numPixAboveTrackedIndex )
            {
                /* travel to the right since:
                 * order > numPixBelowTrackedIndex+hist[pixIndexTracked] 
                 */
                for(pixIndex=pixIndexTracked+1,sum=numPixAboveTrackedIndex;
                    pixIndex < VAL_RANGE; pixIndex++)
                {
                    numPixBelowTrackedIndex = sum;
                    sum += hist[pixIndex];
                    if(sum >=order) { break; }
                }
            }
            else
            {
                /* nothing needs to be done when:
                 * (numPixBelowTrackedIndex < order) &&
                 * (order <= (numPixBelowTrackedIndex+hist[pixIndexTracked]))
                 */
                pixIndex = pixIndexTracked;
            }
        }
    }

    /* Clean up! */
    mxFree((void *) upperRowOffsets);
    mxFree((void *) lowerRowOffsets);
    mxFree((void *) hist);
}

#undef TYPE
#undef MAX_VAL_TYPE
#undef MIN_VAL_TYPE
#undef SIGNED
#undef VAL_RANGE
