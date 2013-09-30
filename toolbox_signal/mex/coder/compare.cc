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
   Image *image1, *image2;

   char *program = argv[0];

#ifdef PGM
   if (argc < 2 || argc > 3)  {
      fprintf (stderr, 
	       "Compare 2 raw byte images.  Gives RMS error and PSNR\n");
      fprintf (stderr, 
	       "Usage: %s [image 1][image 2]\n", program);
      
      return 1;
   }

   image1 = new Image (argv[1]);
   image2 = new Image (argv[2]);

   if ((image1->hsize != image2->hsize) || (image1->vsize != image2->vsize))
     error ("Unequal image sizes:  %s is size %d x %d, %s is size %d x %d\n", 
	    argv[1], argv[2], image1->hsize, image1->vsize, image2->hsize,
	    image2->vsize);

   printf ("comparing %s to %s (%d x %d images)\n", argv[1], argv[2],
	   image1->hsize, image1->vsize);
#else
   int hsize, vsize;

   if (argc < 4 || argc > 5)  {
      fprintf (stderr, 
	       "Compare 2 raw byte images.  Gives RMS error and PSNR\n");
      fprintf (stderr, 
	       "Usage: %s [image 1][image 2][width][height]\n", program);
      
      return 1;
   }

   hsize = atoi(argv[3]);
   if (argc == 5)
     vsize = atoi(argv[4]);
   else
     vsize = hsize;
   printf ("comparing %s to %s (%d x %d images)\n", argv[1], argv[2],
	   hsize, vsize);
   image1 = new Image (argv[1], hsize, vsize);
   image2 = new Image (argv[2], hsize, vsize);
#endif

   Real mse = image1->compare_mse (image2);
   printf ("mean square error = %g\n", mse);
   printf ("root mean square error = %g\n", sqrt(mse));
   if (mse > 0)
     printf ("PSNR = %g dB\n", 20.0 * log(255.0/sqrt(mse))/log(10.0));
   else
     printf ("PSNR = Infinite (images are identical)\n");

   return 0;
}

/*---------------------------------------------------------------------------*/

