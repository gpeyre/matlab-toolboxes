// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//
#ifndef iHistInclude

#define iHistInclude

#include <math.h>

#define iHistDebug 0

#if iHistDebug
#    define iHistDBbegin(h) if (h->debug) {
#    define iHistDBend      }
#else
#    define iHistDBbegin(h) if (0) {
#    define iHistDBend	    }
#endif

// iHistogram:
// provides histogram functionality needed for arithmetic coding.
// Symbols are non-negative integers in the range [0,nsyms-1].
// The initial count for all symbols is 0.
//
// This object implements Alistair Moffat's linear time algorithm for
// adaptive arithmetic coding.
//


// this is a function for filling in the counts in an application
// specific way

typedef void iHistBinFun(void *closure, int *counts);

class iHistogram {
protected:
    int totalCount;
    int maxCt;
    int nsyms;
    int *count;			// pos to ct
    int *treeCount;		// pos to ctLeft
    int *symToPos;		// symbol to position
    int *posToSym;		// position to symbol
    void BuildTreeCount(void);	// rebuild treecount.
    // we do this after making some non-incrmental
    // change to the counts, like scaling.

    void ScaleCounts(void);	// Divides counts by two, except that
    //  if a count was originally non-zero,
    //  doesn't let count drop below 1.
    // DOES NOT REBUILD treeCount!

    inline int Parent(int child) { return (child - 1) >> 1; }
    inline int Lchild(int parent) { return parent + parent + 1; }
    inline int Rchild(int parent) { return parent + parent + 2; }


    Real onelog2;
    void Reorganise(int sym);	// move sym as far to left as possible
    
    int  sorted(int from, int to); //  predicate for an assert
    
    void swapSyms(int pos_a, int pos_b); // swap syms (no tree update) given pos-es.
    
    void pc(int from, int to);
    void qs(int from, int to); // pivot is count[from]
    
    void qsortSyms(void);	// sorts syms, keeping symToPos & posToSym
    // up to date.

public:
    int debug;

    iHistogram(int nsyms, int maxct); 
    ~iHistogram(); 
    
    // following two Init fcns update the count array.
    // rebuilding internal structures afterward.
    // symbol 0 goes in position 0, etc.
    void InitCounts(int *cnts);
    void InitCounts(iHistBinFun *f, void *closure); // can look at old cts
    // if desired.

    int TotalCount();		// total frequency count (scaled)
    int LeftCount(int sym);	// count to left of symbol
    
    // freq of symbol
    inline int Count(int sym) {	return count[symToPos[sym]];}

    void IncCount(int sym);	// increase freq of symbol
    int Symbol(int val);	// get symbol from value

    // Entropy of symbol w.r.t. histogram
    Real  Entropy(int sym)
    { return -log(Count(sym)/(Real)totalCount) * onelog2; }

    void print(ostream &strm);   // print the object
    void printerr(void);        // print the object to stderr
private:
};


#endif
