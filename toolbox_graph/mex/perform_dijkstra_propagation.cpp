/*

	[D,S] = perform_dijkstra_propagation(W,start_verts,end_verts,nb_iter_max,H);
    D is the distance to starting points.
    S is the state : dead=-1, open=0, far=1.
    
    Copyright (c) 2005 Gabriel Peyré

TODO : 
* ajouter un tableau bool pour eventuellement restreindre la propagation (mettre null si pas restriction).
Utiliser cette methode pour speed up le farthest seeding.
* faire de meme dans le FM.
* faire un code pour FM sur mesh via une structure de triangulation locale / 1ring (eviter le loading du fichier off), 
ou alors utiliser GW pour faire les calculs.

*/

#include <math.h>
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "fheap/fib.h"
#include "fheap/fibpriv.h"
#include "mex.h"


#define kDead -1
#define kOpen 0
#define kFar 1

/* Global variables */
int n;			// number of vertices
double* D = NULL;
double* S = NULL;
double* W = NULL;
double* start_points = NULL;
double* end_points = NULL;
double* H = NULL;
int nb_iter_max = 100000;
int nb_start_points = 0;
int nb_end_points = 0;
fibheap_el** heap_pool = NULL;
// sparse array
int* irs = NULL; // returns a pointer to the row indices
int* jcs = NULL;

typedef bool (*T_callback_insert_node)(int i, int ii);

inline 
bool end_points_reached(const int i)
{
	for( int k=0; k<nb_end_points; ++k )
	{
		if( i==((int)end_points[k]) )
			return true;
	}
	return false;
}

inline 
int compare_points(void *x, void *y)
{
	int a = (int) x;
	int b = (int) y;
	if( H==NULL )
		return cmp( D[a], D[b] );
	else
		return cmp( D[a]+H[a], D[b]+H[b] );
}

// test the heap validity
void check_heap( int i )
{
	for( int x=0; x<n; ++x )
		{
			if( heap_pool[x]!=NULL )
			{
				int j = (int) heap_pool[x]->fhe_data;
				if( H==NULL )
				{
					if( D[i]>D[j] )
						mexErrMsgTxt("Problem with heap.\n");
				}
				else
				{
					if( D[i]+H[i]>D[j]+H[j] )
						mexErrMsgTxt("Problem with heap.\n");
				}
			}
		}
}

// select to test or not to test (debug purpose)
// #define CHECK_HEAP check_heap(i,j);
#define CHECK_HEAP


void perform_dijkstra_propagation(T_callback_insert_node callback_insert_node = NULL)
{
	// create the Fibonacci heap
	struct fibheap* open_heap = fh_makeheap();
	fh_setcmp(open_heap, compare_points);

	// initialize points
	for( int i=0; i<n; ++i )
	{
		D[i] = GW_INFINITE;
		S[i] = kFar;
	}

	// record all the points
	heap_pool = new fibheap_el*[n]; 
	memset( heap_pool, NULL, n*sizeof(fibheap_el*) );

	// inialize open list
	for( int k=0; k<nb_start_points; ++k )
	{
		int i = (int) start_points[k];

		if( D[i]==0 )
			mexErrMsgTxt("start_points should not contain duplicates.");

		heap_pool[i] = fh_insert( open_heap, (void*) i );			// add to heap
		D[i] = 0;
		S[i] = kOpen;
	}

	// perform the front propagation
	int num_iter = 0;
	bool stop_iteration = GW_False;
	while( !fh_isempty(open_heap) && num_iter<nb_iter_max && !stop_iteration )
	{
		num_iter++;

		// current point
		int i = (int) fh_extractmin( open_heap );
		heap_pool[i] = NULL;
		S[i] = kDead;
		stop_iteration = end_points_reached(i);

		CHECK_HEAP;

        // index in   irs[jcs[k]]...irs[jcs[k+1]-1] are the node connected to node k
        // values are W[jcs[k]]...W[jcs[k]-1]
		// recurse on each neighbor of i
		for( int k=jcs[i]; k<jcs[i+1]; ++k )
		{
			int ii = irs[k];
            double P = W[k]; // graph weight

			bool bInsert = true;
			if( callback_insert_node!=NULL )
				bInsert = callback_insert_node(ii,ii);
			if( ii>=0 && ii<n && bInsert )
			{
				// compute its neighboring values
				double a1 = D[i] + P;

				if( ((int) S[ii]) == kDead )
				{
					// check if action has change. Should not happen for Dijkstra
					if( a1<D[ii] )
					    mexWarnMsgTxt("The update is not monotone");
#if 1
					if( a1<D[ii] )	// should not happen for FM
						D[ii] = a1;
#endif
				}
				else if( ((int) S[ii]) == kOpen )
				{
					// check if action has change.
					if( a1<D[ii] )
					{
						D[ii] = a1;
						// Modify the value in the heap
						fibheap_el* cur_el = heap_pool[ii];
						if( cur_el!=NULL )
							fh_replacedata( open_heap, cur_el, cur_el->fhe_data );	// use same data for update
						else
							mexErrMsgTxt("Error in heap pool allocation."); 
					}
				}
				else if( ((int) S[ii]) == kFar )
				{
					if( D[ii]!=GW_INFINITE )
						mexErrMsgTxt("Distance must be initialized to Inf");  
					S[ii] = kOpen;
					// distance must have change.
					D[ii] = a1;
					// add to open list
					heap_pool[ii] = fh_insert( open_heap, (void*) ii );			// add to heap	
				}
				else 
					mexErrMsgTxt("Unkwnown state."); 
			}	// end switch
		}		// end for
	}			// end while


	// free heap
	fh_deleteheap(open_heap);
	// free fibheap pool
	GW_DELETEARRAY(heap_pool);
}


void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{ 
	/* retrive arguments */
	if( nrhs<4 ) 
		mexErrMsgTxt("4 or 5 input arguments are required."); 
	if( nlhs<1 ) 
		mexErrMsgTxt("1 or 2 output arguments are required."); 


     /* dealing with sparse array */
    W   = mxGetPr(prhs[0]); // returns a pointer to the numerical values
 	irs     = mxGetIr(prhs[0]); // returns a pointer to the row indices
    jcs     = mxGetJc(prhs[0]); // returns a pointer to the column pointer array
    n = mxGetM(prhs[0]);
    

	// second argument : start_points
	start_points = mxGetPr(prhs[1]);
	nb_start_points = mxGetN(prhs[1]);
	// third argument : end_points
	end_points = mxGetPr(prhs[2]);
	nb_end_points = mxGetN(prhs[2]);
	// third argument : nb_iter_max
	nb_iter_max = (int) *mxGetPr(prhs[3]);
	// second argument : heuristic
	if( nrhs==5 )
	{
		H = mxGetPr(prhs[4]);
		if( H!=NULL && (mxGetM(prhs[4])!=n || mxGetN(prhs[4])!=1) )
			mexErrMsgTxt("H must be of size n x 1."); 
	}
	else
		H = NULL;
	// first ouput : distance
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); 
	D = mxGetPr(plhs[0]);
	// second output : state
	if( nlhs>=2 )
	{
		plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL); 
		S = mxGetPr(plhs[1]);
	}
	else
	{
		S = new double[n];
	}

	// launch the propagation
	perform_dijkstra_propagation();

	if( nlhs<2 )
		GW_DELETEARRAY(S);
	return;
}
