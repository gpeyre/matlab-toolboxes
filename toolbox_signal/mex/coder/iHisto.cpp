// Copyright 1996 John M. Danskin 7/30/96
//
// Permission is granted to use this software for any purpose as
// long as this notice stays attached to this software.
//
#include <math.h>
#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <assert.h>
#include <math.h>
#include "global.hh"
#include "iHisto.h"

iHistogram::iHistogram(int n, int maxct)
{
#if iHistDebug
    debug = 1;
#endif
    onelog2 = 1/log(2.0);
    nsyms = n;
    maxCt = maxct;
    totalCount = 0;
    count = new int[nsyms];
    treeCount = new int[nsyms];
    symToPos = new int[nsyms];
    posToSym = new int[nsyms];

    assert(count && treeCount && symToPos && posToSym);

    for (int i = 0; i < nsyms; i++) {
	count[i] = 0;
	treeCount[i] = 0;
	symToPos[i] = i;
	posToSym[i] = i;
    }    
}

iHistogram::~iHistogram()
{
    delete [] count;
    delete [] treeCount;
    delete [] symToPos;
    delete [] posToSym;
}

#if 0
Real iHistogram::Entropy(int sym)
{
    return -log(Count(sym)/(Real)totalCount) * onelog2;
}
#endif

void iHistogram::BuildTreeCount(void)
{
    int i;
    for (i = 0; i < nsyms; i++) { treeCount[i] = 0; }

    for (i = nsyms - 1; i > 0; i--) {
	int pos = i;
	while (pos > 0) {
	    int parent = Parent(pos);
	    if (pos == Lchild(parent)) {
		treeCount[parent] += count[i];
	    }
	    pos = parent;
	}
    }
}


void iHistogram::swapSyms(int pos_a, int pos_b)
{
    int sym_a = posToSym[pos_a];
    int sym_b = posToSym[pos_b];
    assert(symToPos[sym_a] == pos_a);
    assert(symToPos[sym_b] == pos_b);
    
    symToPos[sym_a] = pos_b;
    symToPos[sym_b] = pos_a;
    
    posToSym[pos_a] = sym_b;
    posToSym[pos_b] = sym_a;
    
    int t = count[pos_a];
    count[pos_a] = count[pos_b];
    count[pos_b] = t;
    
    assert(symToPos[sym_a] == pos_b);
    assert(symToPos[sym_b] == pos_a);
    assert(posToSym[pos_a] == sym_b);
    assert(posToSym[pos_b] == sym_a);
    
}

int iHistogram::sorted(int from, int to)
{
    iHistDBbegin(this);
    cerr << "sorted(" << from << ", " << to << ")->";
    iHistDBend;
    for (int i = from; i < to - 1; i++) {
	if (count[i] < count[i + 1]) {
	    iHistDBbegin(this);
	    cerr << 0 << "\n";
	    iHistDBend;
	    return 0;
	}
    }
    iHistDBbegin(this);
    cerr << 1 << "\n";
    iHistDBend;
    return 1;
}

void iHistogram::pc(int from, int to)
{
    cerr << "<" << from << ", " << to << ">";
    for (int i = from; i < to; i++) {
	cerr << setw(4) << count[i];
    }
    cerr << "\n";
}

void iHistogram::qs(int from, int to)
{
    int pivot = from;
    int ub = to;
    
    iHistDBbegin(this);
    cerr << "qs(" << from << ", " << to << ")\n";
    cerr << "pivoting on " << count[pivot] << "\n";
    pc(from, to);
    iHistDBend;
    
    if (sorted(from, to)) return;

    swapSyms(pivot, (from + to) / 2); // pivot on middle
    
    while (pivot + 1 < ub) {
	if (count[pivot + 1] < count[pivot]) {
	    ub--;
	    iHistDBbegin(this);
	    cerr << "swapping count[" << pivot + 1 << "]=" << count[pivot + 1] <<
		" and count[" << ub << "]=" << count[ub] << "\n";
	    iHistDBend;
	    swapSyms(pivot + 1, ub);
	} else {
	    iHistDBbegin(this);
	    cerr << "swapping count[" << pivot + 1 << "]=" << count[pivot + 1] <<
		" and count[" << pivot << "]=" << count[pivot] << "\n";
	    iHistDBend;
	    swapSyms(pivot, pivot + 1);	
	    pivot++;
	}
	iHistDBbegin(this);
	pc(from, to);
	iHistDBend;
    }
    qs(from, pivot);
    qs(pivot + 1, to);
    
#if iHistDebug
    assert(sorted(from, to));
#endif

    iHistDBbegin(this);
    pc(from, to);
    iHistDBend;
    
}

void iHistogram::qsortSyms(void)
{
    qs(0, nsyms);
}


void iHistogram::ScaleCounts(void)
{
    totalCount = 0;
    for (int i = 0; i < nsyms; i++) {
	if (count[i] == 0) {
	    ;			// skip this one
	} else {
	    count[i] = count[i] / 2;
	    if (!count[i]) {
		count[i] = 1;
	    }
	    totalCount += count[i];
	}
    }
}

void iHistogram::InitCounts(int *cnts)
{
    totalCount = 0;
    for (int i = 0; i < nsyms; i++) {
	symToPos[i] = i;
	posToSym[i] = i;
	totalCount += cnts[i];
	count[i] = cnts[i];
    }
    while (totalCount >= maxCt) { ScaleCounts(); }
    qsortSyms();
    BuildTreeCount();
}

void iHistogram::InitCounts(iHistBinFun *f, void *closure)\
{
    totalCount = 0;
    f(closure, count);
    for (int i = 0; i < nsyms; i++) {
	symToPos[i] = i;
	posToSym[i] = i;
	totalCount += count[i];
    }
    while (totalCount >= maxCt) { ScaleCounts(); }
    qsortSyms();
    BuildTreeCount();
}


int iHistogram::TotalCount() { return totalCount; }

int iHistogram::LeftCount(int sym)
{
    int lc = 0;
    int pos = symToPos[sym];
    assert((sym >= 0) && (sym < nsyms));

    lc = treeCount[pos];

    while (pos) {
	int parent = Parent(pos);
	if (Rchild(parent) == pos) {
	    lc += treeCount[parent] + count[parent];
	}
	pos = parent;
    }

    return lc;
}

#if 0
int iHistogram::Count(int sym)
{
    assert((sym >= 0) && (sym < nsyms));

    return count[symToPos[sym]];
}
#endif

void iHistogram::Reorganise(int sym)
{
    int pos = symToPos[sym];

    if (!pos || count[pos - 1] > count[pos]) {
	return;
    }

    // find leftmost count == count[pos], then swap values

    int right = pos;
    int left  = 0;

    // make sure it isn't under left at start. We'll never put it 
    // there later.
    if (count[left] == count[pos]) { right = left; }

    while (left != right) {
	int middle = (left + right) / 2;
	if (count[middle] > count[pos]) {
	    //middle is to left of leftmost count == count[pos]
	    left = middle + 1;
	} else {
	    //leftmost count == count[pos] is either under middle, or to left
	    //of middle.
	    right = middle;
	}
    }

    assert(count[right] == count[pos]);
    assert(!right || count[right - 1] > count[pos]);

    int rsym = posToSym[right];

    symToPos[sym] = right;
    posToSym[pos] = posToSym[right];
    symToPos[rsym] = pos;
    posToSym[right] = sym;

    
}

void iHistogram::IncCount(int sym)
{
    iHistDBbegin(this) {
	cerr << "IncCount(" << sym << ")\n";
    } iHistDBend;
    
    assert((sym >= 0) && (sym < nsyms));

    Reorganise(sym);

    int pos = symToPos[sym];

    count[pos]++;
    totalCount++;

    while(pos > 0) {
	int parent = Parent(pos);
	if (pos == Lchild(parent)) {
	    treeCount[parent]++;
	}
	pos = parent;
    }

    if (totalCount >= maxCt) {
	ScaleCounts(); 
	BuildTreeCount();
    }
    assert(totalCount < maxCt);
}

int iHistogram::Symbol(int val)
{
    assert(val >= 0 && val < totalCount);

    iHistDBbegin(this) {
	cerr << "Symbol(" << val << ")\n";
    } iHistDBend;

    int pos = 0;    
    int leftSum = 0;
    
    while(1) {
	iHistDBbegin(this) {
	    cerr << "Symbol: pos= " << pos << " leftSum= " << leftSum << " ";
	    cerr << "treeCount[" << pos << "]=" << treeCount[pos] << " ";
	    cerr << "count[" << pos << "]=" << count[pos] << "\n";
	} iHistDBend;

	if (leftSum + treeCount[pos] > val) {
	    // look at left subtree
	    pos = Lchild(pos);
	} else {
	    if (val >= leftSum + treeCount[pos] + count[pos]) {
		// look at right subtree
		leftSum += treeCount[pos] + count[pos];
		pos = Rchild(pos);
	    } else {
		assert(val >= leftSum &&
		       val < leftSum + treeCount[pos] + count[pos]);
		iHistDBbegin(this) {
		    cerr << "Symbol: posToSym[" << pos << "]=" << posToSym[pos]
			 << "\n";
		} iHistDBend;
		return posToSym[pos];
	    }
	}
	assert(pos > 0 && pos < nsyms);
    }
}

void iHistogram::print(ostream &strm)
{
    strm << "iHistogram " << (void *)this << "\n";

    strm << "\ttotalCount=" << totalCount << "\n";

    strm << "\tmaxCt=" << maxCt << "\n";

    strm << "\tnsyms=" << nsyms << "\n";

    int i;
    strm << "symToPos  ";
    for (i = 0; i < nsyms; i++) { strm << setw(6) << symToPos[i] << " "; }
    strm << "\n";
    strm << "count     ";
    for (i = 0; i < nsyms; i++) { strm << setw(6) << count[i] << " "; }
    strm << "\n";
    strm << "treeCount ";
    for (i = 0; i < nsyms; i++) { strm << setw(6) << treeCount[i] << " "; }
    strm << "\n";
    strm << "posToSym  ";
    for (i = 0; i < nsyms; i++) { strm << setw(6) << posToSym[i] << " "; }
    strm << "\n";
}

void iHistogram::printerr() { print(cerr); }

