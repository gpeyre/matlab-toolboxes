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
#ifndef _IMAGE_
#define _IMAGE_
#include "global.hh"
/*---------------------------------------------------------------------------*/
// Basic image class

// Improvements that would be useful:
//    It would be nice if one could load and save in formats other
//    than raw byte format

class Image {
public:
  // Create a blank image with width=hsize, height=vsize
  //    If hsize is unspecified, creates an image of width=0, height=0
  //    If vsize is unspecified, creates a square image of width=height=hsize
  Image  (int hsize = -1, int vsize = -1);

  // Copy constructor
  Image  (const Image& image);

  // Loads an image of size hsize by vsize from the specified file
  //    The file is assumed to contain an image in raw byte format
  Image  (const char *filename, int hsize, int vsize = -1);

  // Loads an image of size hsize by vsize from the specified file
  //    The file is assumed to contain an image in PGM format
  Image  (const char *filename);

  // Destructor
  ~Image ();
  
  // Assignment operator
  Image& operator= (const Image &image);

  // Alternative means of addressing pixels
  Real& operator() (int i) { return value[i]; }
  // Second method of addressing pixels
  Real& operator() (int x, int y) { return value[y*hsize+x]; }
  
  // Loads an image from the specified file.  The file is assumed to
  //   contain an image in raw byte format of size hsize by vsize
  void loadRaw   (const char *filename);

  // Saves an image to the specified file.  The image is written in
  //   raw byte format.
  void saveRaw   (const char *filename);

  // Load an image from a PGM file.
  void loadPGM (const char *filename);
  // Save an image to a PGM file.
  void savePGM (const char *filename);

  // Compares one image to another and returns the PSNR
  Real compare_psnr (const Image *im2);

  // Compares one image to another and returns the mean squared error
  Real compare_mse (const Image *im2);

  // Resize an image -- crop original image or pad with 0's to make fit
  //   new boundaries
  void resize (int new_hsize, int new_vsize);

  unsigned char realToChar (Real x) 
       { return (x>255?255:(x<0?0:(unsigned char)(x+0.5))); };

  Real *value;
  int hsize, vsize;

protected:
  void PGMSkipComments (FILE* infile, unsigned char* ch);
  unsigned int PGMGetVal (FILE* infile);
  void PGMLoadData (FILE *infile, const char *filename);
  void PGMSaveData (FILE* outfile);
};

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
