/*
 * Copyright 1993-2003 The MathWorks, Inc. 
 * $Revision: 1.2.4.1 $  $Date: 2003/01/26 06:00:46 $
 *
 * This file contains a function body for 2-D order-statistic
 * filtering.
 */

(TYPE *pA, TYPE *pB,
 int startRow, int startCol, 
 int Mb, int Nb, int Ma, int order, 
 int *offsets, int numOffsets
#ifdef ADD_OFFSET
 , double *add
#endif /* ADD_OFFSET */
)
{
    TYPE *vector;
    TYPE *p;
    int col;
    int row;
    int k;

    vector = mxCalloc(numOffsets, sizeof(*vector));
    
    for (col = 0; col < Nb; col++)
    {
        p = pA + (startCol+col)*Ma + startRow;
        for (row = 0; row < Mb; row++)
        {
            for (k = 0; k < numOffsets; k++)
            {
#ifdef ADD_OFFSET
		vector[k] = *(p + offsets[k]) + add[k];
#else
                vector[k] = *(p + offsets[k]);
#endif

            }
            *pB++ = 
		SELECT
		(vector, numOffsets, order);
            p++;
        }
    }

    mxFree((void *) vector);
}

#undef TYPE
#undef SELECT
#undef ADD_OFFSET

