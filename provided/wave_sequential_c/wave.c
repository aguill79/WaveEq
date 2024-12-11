// Example physical wave program
//   --a medium is bounded by a square box that act as reflective boundary
//   --conditions. A pulse is made in the middle of this box
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//-------------------------------------------------------------------------
// Included Portable timer functions 
//-------------------------------------------------------------------------
#include "timer.h"


// gnu_plot compatable output
void output_wave_gnu(int N, float *u, char *file) {
   int i,j;
   char fname[30];
   FILE *fp;
   sprintf(fname,"%s.gdt",file); // create gnuplot data file
   if ((fp = fopen(fname,"w")) == NULL) {
      perror("Unable to open gnuplot output file");
      exit(1);
   };

   for (i=1; i<N-1; i++) {
      for (j=1; j<N-1; j++) {
         fprintf(fp,"%d %d %f\n",i,j,u[i*N+j]);
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

// pdmdump compatable output
void pbmdump(int N, float *u, char * file) {

   float min, max, range;
   int i,j;
   FILE  *fp;

   min=max=u[0];
   for (i=0; i<=N; i++) 
      for (j=0; j<=N; j++) {
         if (min > u[i*N+j]) min=u[i*N+j];
         if (max < u[i*N+j]) max=u[i*N+j];
      }
   range = max-min;

   if ((fp = fopen(file,"wb")) == NULL) {
      perror("Unable to open PBM dump file");
   } else {
      fprintf(fp,"P5\n%d %d\n255\n",N,N);
      for (i=0; i<N; i++) 
         for (j=0; j<N; j++) 
            fputc((char) (((u[i*N+j]-min)/range)*255.0), fp);
      fclose(fp);
   }
}


int main(int argc, char *argv[]) {

   char  file[32] = "dump";
   float *u=NULL,*unew=NULL,*uold=NULL,*tmp_ptr=NULL;
   float rho,rho2,h,dt;
   int i,j,t,size;
   int N,steps;
   int iN,iS,jE,jW;
   double timer_val;


   if (argc < 3) {
      perror("Usage: seq_wave <size> <steps> [<output file>] ");
      exit(1);
   }
   sscanf(argv[1],"%d",&N);

   sscanf(argv[2],"%d",&steps);
   printf ("N = %d steps = %d \n",N,steps);

   if (argc == 4) {
      sscanf(argv[3],"%s",file);
   }
   printf("Output will be dumped to <%s> \n\n",file);

   size = N*N; // size = N^2

   if ((u = (float *) malloc(size*sizeof(float)))==NULL) {
      printf("Malloc Error: Memory Allocation Problem for size N=%d\n",N);
      exit(1);
   }
   if ((uold = (float *) malloc(size*sizeof(float)))==NULL) {
      printf("Malloc Error: Memory Allocation Problem for size N=%d\n",N);
      exit(1);
   }
   if ((unew = (float *) malloc(size*sizeof(float)))==NULL) {
      printf("Malloc Error: Memory Allocation Problem for size N=%d\n",N);
      exit(1);
   }

   // Start timer
   startTimer(&timer_val);

   h = 1/(float) N;
   dt =  h/(float) sqrt((double)2.0);
   rho = dt/h;

   // create a disturbance in the middle of the wave plane
   for (i=0; i<N-1; i++) {
      for (j=0; j<N-1; j++) { // set background to
         u[i*N+j]=0.0;        // ref level 0.0
      }
   }

   for (i=N/2-1;i<=N/2;i++) {
      for (j=N/2-1;j<=N/2;j++) { // depress middle 4
         u[i*N+j] = -1.0;        // elements to -1
      }
   }

   u[(N/2-2)*N+N/2-2] = 1.2; // raise up a 
   u[(N/2-2)*N+N/2-1] = 1.4; // volumetrically equivalent
   u[(N/2-2)*N+N/2]   = 1.4; // crest along the edges
   u[(N/2-2)*N+N/2+1] = 1.2; // of the depressed region

   u[(N/2-1)*N+N/2-2] = 1.4;
   u[(N/2-1)*N+N/2+1] = 1.4;
   u[(N/2)*N+N/2-2]   = 1.4;
   u[(N/2)*N+N/2+1]   = 1.4;

   u[(N/2+1)*N+N/2-2] = 1.2;
   u[(N/2+1)*N+N/2-1] = 1.4;
   u[(N/2+1)*N+N/2]   = 1.4;
   u[(N/2+1)*N+N/2+1] = 1.2;

   // set up same conditions one time step earlier -- deriviative 0
   for (i=0; i<N-1; i++) {
      for (j=0; j<N-1; j++) {
         uold[i*N+j]=u[i*N+j];
      }
   }

   rho2=rho*rho;
   for (t=0; t<steps; t++) {
      for (i=0; i<N; i++) {
         iS = (i<N-1 ? (i+1):(N-3));
         iN = (i>0 ? (i-1):(2));
         for (j=0; j<N; j++) {
            jE = (j<N-1 ? (j+1):(N-3));
            jW = (j>0 ? (j-1):(2));
            unew[i*N+j] = 2.0*(1-2.0*rho2)*u[i*N+j] + rho2*(u[iS*N+j] +
                        u[iN*N+j] + u[i*N+jW] + u[i*N+jE]) - uold[i*N+j];
         }
      }
      // swap pointers instead of moving data
      tmp_ptr = uold;
      uold    = u;
      u       = unew;
      unew    = tmp_ptr;
   }

   // End timer
   timer_val=stopNreadTimer(&timer_val);
   printf("Processing time : %f (ms)\n",timer_val*1000.);

   output_wave_gnu(N,u,file);

}

