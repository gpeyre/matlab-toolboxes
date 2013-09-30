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
#include <stdlib.h>
#include <math.h>
#include "global.hh"
#include "image.hh"
/*---------------------------------------------------------------------------*/
// Create a blank image with width=hsize, height=vsize
//    If hsize is unspecified, creates an image with width=0, height=0
//    If vsize is unspecified, creates a square image with width =
//       height = hsize 

Image::Image  (int new_hsize, int new_vsize) : hsize(new_hsize),
    vsize(new_vsize)
{
  if (hsize == -1)
    hsize = vsize = 0;
  if (vsize == -1)
    vsize = hsize;
    
  value = new Real [hsize*vsize];
  if (value == NULL)
    error ("Can't allocate memory for image of size %d by %d\n",
	   hsize, vsize);
}

/*---------------------------------------------------------------------------*/
// Copy constructor

Image::Image (const Image& image)
{
   int i;

   hsize = image.hsize;
   vsize = image.vsize;
   value = new Real [hsize*vsize];
   if (value == NULL)
     error ("Can't allocate memory for image of size %d by %d\n",
	    hsize, vsize);

   for (i = 0; i < hsize*vsize; i++)
     value[i] = image.value[i];
}

/*---------------------------------------------------------------------------*/
// Loads a raw image of size hsize by vsize from the specified file
//    The file is assumed to contain an image in raw byte format

Image::Image (const char *filename, int new_hsize, int new_vsize) : 
  hsize(new_hsize), vsize(new_vsize)
{
  if (hsize == -1)
    hsize = vsize = 0;
  if (vsize == -1)
    vsize = hsize;
    
  value = new Real [hsize*vsize];
  if (value == NULL)
    error ("Can't allocate memory for image of size %d by %d\n",
	   hsize, vsize);

   loadRaw (filename);
}

/*---------------------------------------------------------------------------*/
// Loads a PGM image from a file. Sets hsize, vsize.

Image::Image (const char *filename)
{
  vsize = hsize = 0;
  value = NULL;
  loadPGM (filename);
}

/*---------------------------------------------------------------------------*/
// Destructor

Image::~Image ()
{
   hsize = vsize = -1;
   delete [] value;
}

/*---------------------------------------------------------------------------*/
// Assignment operator

Image &Image::operator= (const Image& image)
{
  delete [] value;
  hsize = image.hsize;
  vsize = image.vsize;
  value = new Real [hsize*vsize];

  for (int i = 0; i < hsize*vsize; i++)
    value[i] = image.value[i];

  return *this;
}

/*---------------------------------------------------------------------------*/
// Loads an image from the specified file.  The file is assumed to
//   contain an image in raw byte format of size hsize by vsize

void Image::loadRaw (const char *filename)
{
   FILE *infile;
   unsigned char *buffer;
   int i;

   infile = fopen (filename, "rb");
   if (infile == NULL)
      error ("Unable to open file %s\n", filename);

   buffer = new unsigned char [hsize * vsize];

   if (fread (buffer, hsize*vsize, sizeof(unsigned char), infile) != 1)
     error ("Read < %d chars when loading file %s\n", hsize*vsize, filename);

   for (i = 0; i < hsize*vsize; i++)
     value[i] = (Real)buffer[i];

   delete [] buffer;
   fclose (infile);
}

/*---------------------------------------------------------------------------*/
// Saves an image to the specified file.  The image is written in
//   raw byte format.

void Image::saveRaw (const char *filename)
{
  FILE *outfile;
  unsigned char *buffer;
  int i;
  
  outfile = fopen (filename, "wb+");
  if (outfile == NULL)
    error ("Unable to open file %s\n", filename);
  
  buffer = new unsigned char [hsize*vsize];
  
  for (i = 0; i < hsize*vsize; i++)
    buffer[i] = realToChar(value[i]);
	 
  fwrite (buffer, hsize*vsize, 1, outfile);
  
  delete [] buffer;
  fclose (outfile);
}

/*---------------------------------------------------------------------------*/
// Private stuff to load pgms/ppms

void Image::PGMSkipComments (FILE* infile, unsigned char* ch)
{
  while ((*ch == '#')) {
    while (*ch != '\n') { *ch = fgetc(infile); }
    while (*ch <  ' ' ) { *ch = fgetc(infile); }
  }
} // Comment(s)

/*---------------------------------------------------------------------------*/
// Get a number from a pgm file header, skipping comments etc.

unsigned int Image::PGMGetVal (FILE* infile)
{
  unsigned int tmp;
  unsigned char ch;
  do { ch = fgetc(infile); } while ((ch <= ' ') && (ch != '#'));
  PGMSkipComments(infile, &ch);
  ungetc(ch, infile);
  if (fscanf(infile,"%u",&tmp) != 1) {
    printf("%s\n","Error parsing file!");
    exit(1);
  }
  return(tmp);
}

/*---------------------------------------------------------------------------*/
// Loads a binary (P5) PGM image from the specified file. Sets hsize and
// vsize to the correct values for the file.

void Image::loadPGM (const char *filename)
{
  FILE* infile;
  unsigned char ch = ' ';

  infile = fopen (filename, "rb");
  if (infile == NULL)
    error ("Unable to open file %s\n", filename);

  // Look for type indicator
  while ((ch != 'P') && (ch != '#')) { ch = fgetc(infile); }
  PGMSkipComments(infile, &ch);
  char ftype = fgetc(infile); // get type, 5 or 6

  // Look for x size, y size, max grey level
  int xsize = (int)PGMGetVal(infile);
  int ysize = (int)PGMGetVal(infile);
  int maxg = (int)PGMGetVal(infile);

  // Do some consistency checks

  if ( (hsize <= 0) && (vsize <= 0) ) {
    resize (xsize, ysize);
    if (value == NULL)
      error ("Can't allocate memory for image of size %d by %d\n",
  	     hsize, vsize);
  } else {
    if ((xsize != hsize) || (ysize != vsize)) {
      error ("File dimensions conflict with image settings\n");
    }
  }

  if (ftype == '5') {
    printf("File %s is of type PGM, is %d x %d with max gray level %d\n",
           filename, hsize, vsize, maxg);
    PGMLoadData(infile, filename);
  }
  if (ftype == '6') {
    printf("File %s is of type PPM, is %d x %d with max gray level %d\n",
           filename, hsize, vsize, maxg);
    error("Attempt to load a PPM as a PGM\n");
  }

  fclose(infile);
}

/*---------------------------------------------------------------------------*/
// Loads the data segment of a PGM image from the specified file.

void Image::PGMLoadData (FILE *infile, const char *filename)
{
   unsigned char *buffer;
   int i;

   buffer = new unsigned char [hsize * vsize];

   long fp = -1*hsize*vsize;
   fseek(infile, fp, SEEK_END);

   if (fread (buffer, hsize*vsize, sizeof(unsigned char), infile) != 1)
     error ("Read < %d chars when loading file %s\n", hsize*vsize, filename);

   for (i = 0; i < hsize*vsize; i++)
     value[i] = (Real)buffer[i];

   delete [] buffer;
}

/*---------------------------------------------------------------------------*/
// Saves an image to the specified file.  The image is written in
// (P5) PGM byte format

void Image::savePGM (const char *filename)
{
  FILE *outfile;
  
  outfile = fopen (filename, "wb+");
  if (outfile == NULL)
    error ("Unable to open file %s\n", filename);

  fprintf(outfile, "P5\n#%s\n%d %d\n255\n", filename, hsize, vsize);
  PGMSaveData(outfile);
  fclose(outfile);
}

/*---------------------------------------------------------------------------*/
// Saves image data to the specified file.  The image is written in
// (P5) PGM raw byte format

void Image::PGMSaveData (FILE* outfile)
{
  unsigned char *buffer;
  int i;
  buffer = new unsigned char [hsize*vsize];
  
  for (i = 0; i < hsize*vsize; i++)
    buffer[i] = realToChar(value[i]);
	 
  fwrite (buffer, hsize*vsize, 1, outfile);
  
  delete [] buffer;
}

/*---------------------------------------------------------------------------*/
// Compares one image to another and returns the PSNR

Real Image::compare_psnr (const Image *im2)
{
  if (im2->hsize != hsize || im2->vsize != vsize)
    error ("Cannot compare images of different sizes (%dx%d and %dx%d)\n",
	   hsize, vsize, im2->hsize, im2->vsize);
  
  Real total_error = 0.0;
  
  for (int i = 0; i < hsize*vsize; i++)
    total_error += square (value[i] - im2->value[i]);

  Real RMS = sqrt(total_error/(Real)(hsize*vsize));
  Real PSNR = 20.0 * log(255.0/RMS)/log(10.0);

  return PSNR;
}

/*---------------------------------------------------------------------------*/
// Compares one image to another and returns the mean squared error

Real Image::compare_mse (const Image *im2)
{
  if (im2->hsize != hsize || im2->vsize != vsize)
    error ("Cannot compare images of different sizes (%dx%d and %dx%d)\n",
	   hsize, vsize, im2->hsize, im2->vsize);
  
  Real total_error = 0.0;
  
  for (int i = 0; i < hsize*vsize; i++)
    total_error += square (value[i] - im2->value[i]);

  Real MSE = total_error/(Real)(hsize*vsize);
  
  return MSE;
}

/*---------------------------------------------------------------------------*/
// Resize an image -- crop original image or pad with 0's to make fit
//   new boundaries

void Image::resize (int new_hsize, int new_vsize)
{
  Real *new_value = new Real [new_hsize*new_vsize];

  for (int i = 0; i < new_hsize*new_vsize; i++)
    new_value[i] = 0;

  for (int j = 0; j < min(vsize, new_vsize); j++)
    for (int i = 0; i < min(hsize, new_hsize); i++)
      new_value[j*new_hsize+i] = value[j*hsize+i];

  hsize = new_hsize;
  vsize = new_vsize;

  if (value != NULL)
    delete [] value;
  value = new_value;
}

/*---------------------------------------------------------------------------*/
