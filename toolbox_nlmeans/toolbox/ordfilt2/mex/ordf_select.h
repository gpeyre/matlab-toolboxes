/*
 * Copyright 1993-2003 The MathWorks, Inc. 
 * $Revision: 1.2.4.2 $  $Date: 2003/01/26 06:00:47 $
 *
 * This file contains a function body for selection algorithm using 
 * QuickSort-style partitioning. It is based on work presented 
 * in Sedgewick, Algorithms, 2d ed, 1988, pp. 115-130 
 * This function modifies the values in a[].  
 */

(TYPE *a, int N, int order)
{
    int i;
    int j;
    int left;
    int right;
    int mid;
    TYPE key;
    TYPE tmp;


#ifdef CHECK_NANS
    int k;

    /* If a NaN is present, return NaN. */
    for (k = 0; k < N; k++)
    {
        if (mxIsNaN(a[k]))
        {
            return((TYPE)mxGetNaN());
        }
    }
#endif /* CHECK_NANS */


    /* The initial partition is the entire array. */
    left = 0;
    right = N - 1;

    /* Perform the scanning only if the current partition has at least */
    /* three elements. */
    while (right > (left+1)) 
    {
        /* Use median-of-three method to select the partitioning key. */
        /* This method makes worst-case behavior much less likely and */
        /* avoids the needs for sentinel values outside the array. */

        /* First we sort the first, middle, and last element of the */
        /* partition. */
        mid = (left + right) >> 1;
        if (a[left] > a[mid])
        {
            tmp = a[left];
            a[left] = a[mid];
            a[mid] = tmp;
        }
        if (a[left] > a[right])
        {
            tmp = a[left];
            a[left] = a[right];
            a[right] = tmp;
        }
        if (a[mid] > a[right])
        {
            tmp = a[mid];
            a[mid] = a[right];
            a[right] = tmp;
        }
        
        /* Now swap the middle value with the next-to-last value. */
        tmp = a[mid];
        a[mid] = a[right - 1];
        a[right - 1] = tmp;
        
        /* Use the median of the three values as the partitioning key */
        key = a[right - 1];
        
        /* Start the partitioning scan.  Note that because we sorted */
        /* the left and right values, we can start the comparisons */
        /* with (left+1) and (right-2).  This is how we avoid the */
        /* need for sentinel values. */
        i = left;
        j = right - 1;
        for (;;)
        {
            while (a[++i] < key) ;      mxAssert(i <= (N-1), "");
            while (a[--j] > key) ;      mxAssert(j >= 0, "");
            if (i >= j) break;     /* pointers crossed; end scan */
            /* swap values at end of current interval */
            tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
        /* One last swap needed at end of scan */
        tmp = a[i];
        a[i] = a[right-1];
        a[right-1] = tmp;

        /* Select the left or right subpartition depending on */
        /* the value of order. */
        if (i >= order) right = i-1;
        if (i <= order) left = i+1;
    }

    if ((right - left) == 1)
    {
        /* Last partition has two elements that may not be sorted. */
        if (a[left] > a[right])
        {
            tmp = a[left];
            a[left] = a[right];
            a[right] = tmp;
        }
    }
            
    return(a[order]);
}

#undef TYPE
#undef CHECK_NANS


