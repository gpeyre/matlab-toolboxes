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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <new.h>
#include <math.h>
#include <assert.h>
#include "global.hh"
/*---------------------------------------------------------------------------*/
#ifdef DEBUG
static FILE *debug_file;
static int debug_file_open = FALSE;
#endif

/*---------------------------------------------------------------------------*/
// function called when out of memory -- put a debugger breakpoint here
//   if trying to locate cause of Out of memory error

void no_more_memory ()
{
   error ("Out of memory");
}

/*---------------------------------------------------------------------------*/
// Initialize system-level things

void init()
{
  // Call no_more_memory when unable to malloc
  set_new_handler (no_more_memory);
  
#ifdef DEBUG
  debug_file = fopen ("debug.log", "w+");
  debug_file_open = (debug_file != NULL);
#endif
}

/*---------------------------------------------------------------------------*/
// Close down system-level stuff

void shut_down ()
{
#ifdef DEBUG
   fclose (debug_file);
   debug_file_open = FALSE;
#endif
}

/*---------------------------------------------------------------------------*/

volatile void error (char *format, ...)
{
   va_list list;

   va_start (list, format);

   printf ("Error: ");
   vprintf (format, list);
   va_end (list);
   printf ("\n");

#ifdef DEBUG
   if (debug_file_open) {
     fprintf (debug_file, "Error: ");
     vfprintf (debug_file, format, list);
     fprintf (debug_file, "\n");
     fflush (debug_file);
   }
#endif

   assert(0);
}

/*---------------------------------------------------------------------------*/

void warning (char *format, ...)
{
   va_list list;

   va_start (list, format);

#ifdef DEBUG
   if (debug_file_open) {
     fprintf (debug_file, "Warning: ");
     vfprintf (debug_file, format, list);
     fprintf (debug_file, "\n");
     fflush (debug_file);
   }
#endif

   printf ("Warning: ");
   vprintf (format, list);
   va_end (list);
   printf ("\n");
}

/*---------------------------------------------------------------------------*/
