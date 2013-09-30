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
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "wavelet.hh"
/*---------------------------------------------------------------------------*/
Real HaarCoeffs [] = { 1.0/Sqrt2, 1.0/Sqrt2 };

// A few Daubechies filters

Real Daub4Coeffs [] = { 0.4829629131445341,  0.8365163037378077,
		          0.2241438680420134, -0.1294095225512603 };

Real Daub6Coeffs [] = { 0.3326705529500825,  0.8068915093110924,
		          0.4598775021184914, -0.1350110200102546,
		         -0.0854412738820267,  0.0352262918857095 };

Real Daub8Coeffs [] = { 0.2303778133088964,  0.7148465705529154,
			  0.6308807679398587, -0.0279837694168599,
			 -0.1870348117190931,  0.0308413818355607,
			  0.0328830116668852, -0.0105974017850690 };

// Filter from Eero Simoncelli's PhD thesis -- used in Edward Adelson's EPIC wavelet coder
// These are probably the filter coefficients used in Shapiro's EZW paper
Real AdelsonCoeffs[] = { 0.028220367, -0.060394127, -0.07388188, 
			 0.41394752,   0.7984298,   0.41394752, 
			 -0.07388188, -0.060394127, 0.028220367 };

// 7/9 Filter from M. Antonini, M. Barlaud, P. Mathieu, and
// I. Daubechies, "Image coding using wavelet transform", IEEE
// Transactions on Image Processing", Vol. pp. 205-220, 1992.

Real AntoniniSynthesis [] = { -6.453888262893856e-02,
			      -4.068941760955867e-02,
			       4.180922732222124e-01,
			       7.884856164056651e-01,
			       4.180922732222124e-01,
			      -4.068941760955867e-02,
			      -6.453888262893856e-02 };
Real AntoniniAnalysis[] =  {  3.782845550699535e-02,
			     -2.384946501937986e-02,
			     -1.106244044184226e-01,
			      3.774028556126536e-01,
			      8.526986790094022e-01,
			      3.774028556126537e-01,
			     -1.106244044184226e-01,
			     -2.384946501937986e-02,
			      3.782845550699535e-02 };

// Unpublished 18/10 filter from Villasenor's group

Real Villa1810Synthesis [] = { 9.544158682436510e-04,
			  -2.727196296995984e-06,
			  -9.452462998353147e-03,
			  -2.528037293949898e-03,
			   3.083373438534281e-02,
			  -1.376513483818621e-02,
			  -8.566118833165798e-02,
			   1.633685405569902e-01,
			   6.233596410344172e-01,
			   6.233596410344158e-01,
			   1.633685405569888e-01,
			  -8.566118833165885e-02,
			  -1.376513483818652e-02,
			   3.083373438534267e-02,
			  -2.528037293949898e-03,
			  -9.452462998353147e-03,
			  -2.727196296995984e-06,
			   9.544158682436510e-04};
Real Villa1810Analysis [] = {  2.885256501123136e-02,
			   8.244478227504624e-05,
			  -1.575264469076351e-01,
			   7.679048884691438e-02,
			   7.589077294537618e-01,
			   7.589077294537619e-01,
			   7.679048884691436e-02,
			  -1.575264469076351e-01,
			   8.244478227504624e-05,
			   2.885256501123136e-02};

// Filters from Chris Brislawn's tutorial code

Real BrislawnAnalysis [] = {  0.037828455506995,  -0.023849465019380, 
			       -0.110624404418423,   0.377402855612654,
				0.852698679009403, 
				0.377402855612654,  -0.110624404418423,
			       -0.023849465019380,   0.037828455506995 };
Real BrislawnSynthesis [] = { -0.064538882628938, -0.040689417609558,
				 0.418092273222212, 
				 0.788485616405664,
				 0.418092273222212, 
				-0.040689417609558, -0.064538882628938};

Real Brislawn2Analysis [] = {  0.026913419, -0.032303352,
			        -0.241109818,  0.054100420,
			         0.899506092,  0.899506092, 
			         0.054100420, -0.241109818, 
			        -0.032303352,  0.026913419};
Real Brislawn2Synthesis [] = { 0.019843545,  0.023817599, 
				-0.023257840,  0.145570740,  
				 0.541132748,  0.541132748, 
				 0.145570740, -0.023257840, 
				 0.023817599, 0.019843545 };

// Filters from J. Villasenor, B. Belzer, J. Liao, "Wavelet Filter
// Evaluation for Image Compression." IEEE Transactions on Image
// Processing, Vol. 2, pp. 1053-1060, August 1995.

Real Villa1Analysis [] = {
  3.782845550699535e-02,
  -2.384946501937986e-02,
  -1.106244044184226e-01,
  3.774028556126536e-01,
  8.526986790094022e-01,
  3.774028556126537e-01,
  -1.106244044184226e-01,
  -2.384946501937986e-02,
  3.782845550699535e-02
};

Real Villa1Synthesis [] = {
  -6.453888262893856e-02,
  -4.068941760955867e-02,
  4.180922732222124e-01,
  7.884856164056651e-01,
  4.180922732222124e-01,
  -4.068941760955867e-02,
  -6.453888262893856e-02
};

Real Villa2Analysis [] = {
  -8.472827741318157e-03,
  3.759210316686883e-03,
  4.728175282882753e-02,
  -3.347508104780150e-02,
  -6.887811419061032e-02,
  3.832692613243884e-01,
  7.672451593927493e-01,
  3.832692613243889e-01,
  -6.887811419061045e-02,
  -3.347508104780156e-02,
  4.728175282882753e-02,
  3.759210316686883e-03,
  -8.472827741318157e-03
};
Real Villa2Synthesis [] = {
  1.418215589126359e-02,
  6.292315666859828e-03,
  -1.087373652243805e-01,
  -6.916271012030040e-02,
  4.481085999263908e-01,
  8.328475700934288e-01,
  4.481085999263908e-01,
  -6.916271012030040e-02,
  -1.087373652243805e-01,
  6.292315666859828e-03,
  1.418215589126359e-02
};

Real Villa3Analysis [] = {
  -1.290777652578771e-01,
  4.769893003875977e-02,
  7.884856164056651e-01,
  7.884856164056651e-01,
  4.769893003875977e-02,
  -1.290777652578771e-01
};
Real Villa3Synthesis [] = {
  1.891422775349768e-02,
  6.989495243807747e-03,
  -6.723693471890128e-02,
  1.333892255971154e-01,
  6.150507673110278e-01,
  6.150507673110278e-01,
  1.333892255971154e-01,
  -6.723693471890128e-02,
  6.989495243807747e-03,
  1.891422775349768e-02
};

Real Villa4Analysis [] = {
  -1.767766952966369e-01,
  3.535533905932738e-01,
  1.060660171779821e+00,
  3.535533905932738e-01,
  -1.767766952966369e-01
};
Real Villa4Synthesis [] = {
  3.535533905932738e-01,
  7.071067811865476e-01,
  3.535533905932738e-01
};

Real Villa5Analysis [] = {
  7.071067811865476e-01,
  7.071067811865476e-01
};
Real Villa5Synthesis [] = {
  -8.838834764831845e-02,
  8.838834764831845e-02,
  7.071067811865476e-01,
  7.071067811865476e-01,
  8.838834764831845e-02,
  -8.838834764831845e-02
};

Real Villa6Analysis [] = {
  3.314563036811943e-02,
  -6.629126073623885e-02,
  -1.767766952966369e-01,
  4.198446513295127e-01,
  9.943689110435828e-01,
  4.198446513295127e-01,
  -1.767766952966369e-01,
  -6.629126073623885e-02,
  3.314563036811943e-02

};
Real Villa6Synthesis [] = {
  3.535533905932738e-01,
  7.071067811865476e-01,
  3.535533905932738e-01
};

// Filter

Real OdegardAnalysis[] = {
   5.2865768532960523e-02,
  -3.3418473279346828e-02,
  -9.3069263703582719e-02,
   3.8697186387262039e-01,
   7.8751377152779212e-01,
   3.8697186387262039e-01,
  -9.3069263703582719e-02,
  -3.3418473279346828e-02,
   5.2865768532960523e-02
};
Real OdegardSynthesis[] = {
  -8.6748316131711606e-02,
  -5.4836926902779436e-02,
   4.4030170672498536e-01,
   8.1678063499210640e-01,
   4.4030170672498536e-01,
  -5.4836926902779436e-02,
  -8.6748316131711606e-02
};

FilterSet Haar      (FALSE, HaarCoeffs,          2, 0);
FilterSet Daub4     (FALSE, Daub4Coeffs,         4, 0);
FilterSet Daub6     (FALSE, Daub6Coeffs,         6, 0);
FilterSet Daub8     (FALSE, Daub8Coeffs,         8, 0);
FilterSet Antonini  (TRUE,  AntoniniAnalysis,    9, -4, 
		            AntoniniSynthesis,   7, -3);
FilterSet Villa1810 (TRUE,  Villa1810Analysis,  10, -4,
		            Villa1810Synthesis, 18, -8);
FilterSet Adelson   (TRUE,  AdelsonCoeffs,       9, -4);
FilterSet Brislawn  (TRUE,  BrislawnAnalysis,    9, -4, 
		            BrislawnSynthesis,   7, -3);
FilterSet Brislawn2 (TRUE,  Brislawn2Analysis,  10, -4,
		            Brislawn2Synthesis, 10, -4);

FilterSet Villa1 (TRUE, Villa1Analysis,  9, -4, Villa1Synthesis,  7, -3);
FilterSet Villa2 (TRUE, Villa2Analysis, 13, -6, Villa2Synthesis, 11, -5);
FilterSet Villa3 (TRUE, Villa3Analysis,  6, -2, Villa3Synthesis, 10, -4);
FilterSet Villa4 (TRUE, Villa4Analysis,  5, -2, Villa4Synthesis,  3, -1);
FilterSet Villa5 (TRUE, Villa5Analysis,  2,  0, Villa5Synthesis,  6, -2);
FilterSet Villa6 (TRUE, Villa6Analysis,  9, -4, Villa6Synthesis,  3, -1);


FilterSet Odegard (TRUE, OdegardAnalysis, 9, -4, OdegardSynthesis, 7, -3);
/*---------------------------------------------------------------------------*/
// Destructor

Filter::~Filter ()
{
  if (coeff != NULL)
    delete [] coeff;
}

/*---------------------------------------------------------------------------*/

void Filter::init (int filterSize, int filterFirst, Real *data)
{
   size = filterSize;
   firstIndex = filterFirst;
   center = -firstIndex;

   coeff = new Real [size];
   if (data != NULL) {
     for (int i = 0; i < size; i++)
       coeff[i] = data[i];
   } else {
     for (int i = 0; i < size; i++)
       coeff[i] = 0;
   }
}

/*---------------------------------------------------------------------------*/

void Filter::copy (const Filter& filter)
{
  if (coeff != NULL)
    delete [] coeff;
  init (filter.size, filter.firstIndex, filter.coeff);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FilterSet::FilterSet (int symmetric, 
	     Real *anLow, int anLowSize, int anLowFirst,  
	     Real *synLow, int synLowSize, int synLowFirst) : 
             symmetric(symmetric) 
{
  int i, sign;

  analysisLow = new Filter (anLowSize, anLowFirst, anLow);

  // If no synthesis coeffs are given, assume wavelet is orthogonal
  if (synLow == NULL)  {
    synthesisLow = new Filter (*analysisLow);

    // For orthogonal wavelets, compute the high pass filter using
    // the relation g_n = (-1)^n h_{1-n}^*
    // (or equivalently g_{1-n} = (-1)^{1-n} h_n^*)

    analysisHigh = new Filter (analysisLow->size, 2 - analysisLow->size -
				analysisLow->firstIndex);
      
    // Compute (-1)^(1-n) for first n
    if (analysisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;
    
    for (i = 0; i < analysisLow->size; i++)  {
      analysisHigh->coeff[1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex] = 
	sign * analysisLow->coeff[i];
      assert (1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex >= 0);
      assert (1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex < analysisHigh->size);
      sign *= -1;
    }

    // Copy the high pass analysis filter to the synthesis filter
    synthesisHigh = new Filter (*analysisHigh);    
    
  } else {
    // If separate synthesis coeffs given, assume biorthogonal
    
    synthesisLow = new Filter (synLowSize, synLowFirst, synLow);

    // For orthogonal wavelets, compute the high frequency filter using
    // the relation g_n = (-1)^n complement (h~_{1-n}) and
    //              g~_n = (-1)^n complement (h_{1-n})
    // (or equivalently g_{1-n} = (-1)^{1-n} complement (h~_n))
    
    analysisHigh = new Filter (synthesisLow->size, 2 - synthesisLow->size -
			 synthesisLow->firstIndex);
    
    // Compute (-1)^(1-n) for first n
    if (synthesisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;
    
    for (i = 0; i < synthesisLow->size; i++)  {
      analysisHigh->coeff[1 - i - synthesisLow->firstIndex -
			  analysisHigh->firstIndex] = 
	sign * synthesisLow->coeff[i];
      assert (1 - i - synthesisLow->firstIndex - 
		           analysisHigh->firstIndex >= 0);
      assert (1 - i - synthesisLow->firstIndex - 
		           analysisHigh->firstIndex < analysisHigh->size);
      sign *= -1;
    }

    synthesisHigh = new Filter 
                           (analysisLow->size, 2 - analysisLow->size -
			    analysisLow->firstIndex); 
    
    // Compute (-1)^(1-n) for first n
    if (analysisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;

    for (i = 0; i < analysisLow->size; i++)  {
      synthesisHigh->coeff[1 - i - analysisLow->firstIndex -
			   synthesisHigh->firstIndex] = 
	sign * analysisLow->coeff[i];
      assert (1 - i - analysisLow->firstIndex - 
		           synthesisHigh->firstIndex >= 0);
      assert (1 - i - analysisLow->firstIndex - 
		           synthesisHigh->firstIndex < synthesisHigh->size);
      sign *= -1;
    }
  }
}

/*---------------------------------------------------------------------------*/

FilterSet::FilterSet (const FilterSet& filterset)
{
  copy (filterset);
}

/*---------------------------------------------------------------------------*/

FilterSet::~FilterSet ()
{
  delete analysisLow;
  delete analysisHigh;
  delete synthesisLow;
  delete synthesisHigh;
}

/*---------------------------------------------------------------------------*/

FilterSet& FilterSet::operator= (const FilterSet filterset)
{
  delete analysisLow;
  delete analysisHigh;
  delete synthesisLow;
  delete synthesisHigh;
  copy (filterset);
  return *this;
}

/*---------------------------------------------------------------------------*/

void FilterSet::copy (const FilterSet& filterset)
{
  symmetric = filterset.symmetric;
  analysisLow = new Filter (*(filterset.analysisLow));
  analysisHigh = new Filter (*(filterset.analysisHigh));
  synthesisLow = new Filter (*(filterset.synthesisLow));
  synthesisHigh = new Filter (*(filterset.synthesisHigh));
}

/*---------------------------------------------------------------------------*/
