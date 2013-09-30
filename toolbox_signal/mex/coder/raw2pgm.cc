/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include "global.hh"
#include "image.hh"
/*---------------------------------------------------------------------------*/

int main (int argc, char **argv)
{
   Image *image;

   char *program = argv[0];

   if (argc != 5)  {
      fprintf (stderr, 
	       "Convert an image in pbm/pgm format to raw pixel format\n");
      fprintf (stderr, 
	       "Usage: %s [raw image name][height][width][pgm image name]\n", 
	       program);
      
      return 1;
   }

   int hsize = atoi(argv[3]);
   int vsize = atoi(argv[2]);
   image = new Image (argv[1], hsize, vsize);
   
   image->savePGM (argv[4]);
}

/*---------------------------------------------------------------------------*/

