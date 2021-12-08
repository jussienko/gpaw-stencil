/*  Copyright (C) 2015       CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "bmgs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double mytime(void);

void unpack(bmgsstencil *stencil, const double *in, double* buf, int range);

int main(int argc, char** argv)
{
   
   bmgsstencil stencil;
   double h[3] = {0.2, 0.2, 0.2};
   long gpts[3] = {64, 64, 64};
   int nbands = 256;
   int order = 7;

   int repeats = 10;

   int nthreads = 1;

   if (argc > 1)
     gpts[0] = gpts[1] = gpts[2] = atoi(argv[1]);
   if (argc == 3)
     nbands = atoi(argv[2]);

   printf("Grid dimension %i, number of bands %i\n", gpts[0], nbands);
#ifdef _OPENMP
#pragma omp parallel 
{
   #pragma omp master
   nthreads = omp_get_num_threads();
}
#endif
   printf("Number of threads %i\n", nthreads);

   int ngpts = gpts[0] * gpts[1] * gpts[2];
     
   double *input = (double *) malloc(ngpts * nbands * sizeof(double));
   double *output = (double *) malloc(ngpts * nbands * sizeof(double));

   #pragma omp parallel
   for (int n=0; n < nbands; n++)
     {
       #pragma omp for schedule(static) collapse(2)
       for (int ix=0; ix < gpts[0]; ix++)
         for (int iy=0; iy < gpts[1]; iy++)
           for (int iz=0; iz < gpts[2]; iz++)
             {
               int ind = n * ngpts + iz + iy * gpts[2] + ix * gpts[2] * gpts[1];
               input[ind] = h[0]*h[0] * (ix*ix + iy*iy + iz*iz);
               output[ind] = 0.0;
             }
     }
				     
   stencil = bmgs_laplace(order, 0.5, h, gpts);

   int ngpts2 = (2*order + gpts[0]) * (2*order + gpts[1]) * (2*order + gpts[2]);
   double *buffer = (double *) malloc(ngpts2 * sizeof(double));

   /* Warm-up */
   for (int r = 0; r < 0; r++) {
     #pragma omp parallel
     for (int n = 0; n < nbands; n++)
       {
         const double *in = input + n * ngpts;
         double *out = output + n * ngpts;
         unpack(&stencil, in, buffer, order);
         bmgs_fd(&stencil, buffer, out);
       }
     }
   
   double t0 = mytime();
   for (int r = 0; r < repeats; r++) {
     #pragma omp parallel
     for (int n = 0; n < nbands; n++)
       {
	     const double *in = input + n * ngpts;
         double *out = output + n * ngpts;
         unpack(&stencil, in, buffer, order);
         bmgs_fd(&stencil, buffer, out);
       }
   }
   
   double t1 = mytime();

   double gflop = repeats*nbands*ngpts * 1.0e-9; 
   gflop *= stencil.ncoefs*2;

   printf("Time %6.3f s, %6.3f GFLOPS\n", t1 - t0, 
	  gflop / (t1 - t0));

   double s = 0.0;
   for (int n = 0; n < nbands * ngpts; n++)
     s += output[n];
   printf("Sum:    %f\n", s);

}

double mytime(void)
{
  struct timeval t;
  gettimeofday(&t, (struct timezone *)NULL);
  return (double)t.tv_sec + t.tv_usec*0.000001;
}

void unpack(bmgsstencil *s, const double *in, double *buf, 
	    int range)
{

  const int startb[3] = {range, range, range};

  int sizea[3];
  int sizeb[3];

  int ng2 = 1;
  for (int i=0; i < 3; i++)
    {
      sizea[i] = s->n[i];
      sizeb[i] = s->n[i] + 2*range;
      ng2 *= sizeb[i];
    }

  // Zero the boundaries
 #pragma omp single
 memset(buf, 0, ng2 * sizeof(double));

  // Copy the central part
  bmgs_paste(in, sizea, buf, sizeb, startb);
}

