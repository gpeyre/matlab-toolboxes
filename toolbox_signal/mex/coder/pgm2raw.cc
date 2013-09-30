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

   if (argc != 3)  {
      fprintf (stderr, 
	       "Convert an image in pbm/pgm format to raw pixel format\n");
      fprintf (stderr, 
	       "Usage: %s [pgm image name][raw image name]\n", program);
      
      return 1;
   }

   image = new Image (argv[1]);
   image->saveRaw (argv[2]);
}

/*---------------------------------------------------------------------------*/

