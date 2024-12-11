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

#define TILE_WIDTH 16 // block x and y dimensions

// iceil macro
// returns an integer ceil value where integer numerator is first parameter
// and integer denominator is the second parameter. iceil is the rounded
// up value of numerator/denominator when there is a remainder
#define iceil(num,den) (num/den+(num%den>0))


__global__ void WaveKernel(float *unew, float *u, float *uold, int NN,
                           float rho2) {

   // Calculate the row index of the Pd element and M
   int Row = blockIdx.y*TILE_WIDTH+threadIdx.y + 1;
   int Col = blockIdx.x*TILE_WIDTH+threadIdx.x + 1;

   if ((Row < NN-1) && (Col < NN-1)) {
     unew[Row*NN+Col]=2.0*(1-2.0*rho2)*u[Row*NN+Col] + rho2*(u[(Row+1)*NN+Col]+ 
                      u[(Row-1)*NN+Col] + u[Row*NN+(Col-1)] + u[Row*NN+(Col+1)])
                     - uold[Row*NN+Col];
   }
}

// gnu_plot compatable output
void output_wave_gnu(int NN, float *u, char *file) {
   int i,j;
   char fname[30];
   FILE *fp;
   sprintf(fname,"%s.gdt",file); // create gnuplot data file
   if ((fp = fopen(fname,"w")) == NULL) {
      perror("Unable to open gnuplot output file");
      exit(1);
   };

   for (i=1; i<NN-1; i++) {
      for (j=1; j<NN-1; j++) {
         fprintf(fp,"%d %d %f\n",i,j,u[i*NN+j]);
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
   float *u=NULL,*uold=NULL,*tmp_ptr=NULL;
   float rho,rho2,h,dt;
   int i,j,t,size;
   int N,NN,steps;
   cudaError_t error_id;
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

   NN = N+2; // size of expanded matrix that includes boundary conditions
   size = NN*NN; // size = (N+2)^2 leaves room for boundary conditions

   if ((u = (float *) malloc(size*sizeof(float)))==NULL) {
      printf("Malloc Error: Memory Allocation Problem for size N=%d\n",N);
      exit(1);
   }
   if ((uold = (float *) malloc(size*sizeof(float)))==NULL) {
      printf("Malloc Error: Memory Allocation Problem for size N=%d\n",N);
      exit(1);
   }

   // Start timer
   startTimer(&timer_val);

   h = 1/(float) N;
   dt =  h/(float) sqrt((double)2.0);
   rho = dt/h;

   // create a disturbance in the middle of the wave plane
   for (i=0; i<NN; i++) {
      for (j=0; j<NN; j++) { // set background to
         u[i*NN+j]=0.0;        // ref level 0.0
      }
   }

   for (i=NN/2-1;i<=NN/2;i++) {
      for (j=NN/2-1;j<=NN/2;j++) { // depress middle 4
         u[i*NN+j] = -1.0;        // elements to -1
      }
   }

   // excite the region of interest
   u[(NN/2-2)*NN+NN/2-2] = 1.2; // raise up a 
   u[(NN/2-2)*NN+NN/2-1] = 1.4; // volumetrically equivalent
   u[(NN/2-2)*NN+NN/2]   = 1.4; // crest along the edges
   u[(NN/2-2)*NN+NN/2+1] = 1.2; // of the depressed region

   u[(NN/2+1)*NN+NN/2-2] = 1.4;
   u[(NN/2+1)*NN+NN/2+1] = 1.4;
   u[(NN/2)*NN+NN/2-2]   = 1.4;
   u[(NN/2)*NN+NN/2+1]   = 1.4;

   u[(NN/2+1)*NN+NN/2-2] = 1.2;
   u[(NN/2+1)*NN+NN/2+1] = 1.4;
   u[(NN/2+1)*NN+NN/2-2]   = 1.4;
   u[(NN/2+1)*NN+NN/2+1] = 1.2;

   // set up same conditions one time step earlier -- deriviative 0
   for (i=0; i<NN; i++) {
      for (j=0; j<NN; j++) {
         uold[i*NN+j]=u[i*NN+j];
      }
   }

   rho2=rho*rho;

   float *uoldd,*ud,*unewd;

   // Allocate device memory and Transfer host arrays u and uold 
   error_id=cudaMalloc((void **) &uoldd, size*sizeof(float));
   if (error_id != cudaSuccess) {
      printf( "Device Memory allocation for uoldd failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaMemcpy(uoldd, uold, size*sizeof(float), cudaMemcpyHostToDevice);
   if (error_id != cudaSuccess) {
      printf( "Memory Copy for uoldd failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaMalloc((void **) &ud, size*sizeof(float));
   if (error_id != cudaSuccess) {
      printf( "Device Memory allocation for ud failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaMemcpy(ud, u, size*sizeof(float), cudaMemcpyHostToDevice);
   if (error_id != cudaSuccess) {
      printf( "Memory Copy for ud failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   // Allocate device memory of unew array for results
   // (note memory is copied because boundary conditions will be needed
   //  later when pointers are swapped)
   error_id=cudaMalloc((void **) &unewd, size*sizeof(float));
   if (error_id != cudaSuccess) {
      printf( "Device Memory allocation for unewd failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaMemcpy(unewd, u, size*sizeof(float), cudaMemcpyHostToDevice);
   if (error_id != cudaSuccess) {
      printf( "Memory Copy for unewd failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }


   // Setup the kernel execution configuration parameters
   dim3 dimGrid;
   dimGrid.x = iceil(N,TILE_WIDTH);
   dimGrid.y = iceil(N,TILE_WIDTH);
   dim3 dimBlock(TILE_WIDTH,TILE_WIDTH);

   for (t=0; t<steps; t++) {

      // Launch the kernel!!!
      WaveKernel<<<dimGrid, dimBlock>>>(unewd, ud, uoldd, NN, rho2);

      error_id=cudaGetLastError();
      if (error_id != cudaSuccess) {
         printf( "Attempted Launch of WaveKernel returned %d\n-> %s\n",
             (int)error_id, cudaGetErrorString(error_id) );
         exit(EXIT_FAILURE);
      }

      // swap device pointers instead of moving data
      tmp_ptr = uoldd;
      uoldd    = ud;
      ud       = unewd;
      unewd    = tmp_ptr;
   }

   // Transfer P from device to host
   error_id=cudaMemcpy(u,unewd,size*sizeof(float) ,cudaMemcpyDeviceToHost);
   if (error_id != cudaSuccess) {
      printf( "Memory Copy back to host for u failed--returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }


   // Free device matrices
   error_id=cudaFree(uoldd);
   if (error_id != cudaSuccess) {
      printf( "Cuda could not free memory uoldd -- returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaFree(ud);
   if (error_id != cudaSuccess) {
      printf( "Cuda could not free memory ud -- returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }

   error_id=cudaFree(unewd);
   if (error_id != cudaSuccess) {
      printf( "Cuda could not free memory unewd -- returned %d\n-> %s\n",
          (int)error_id, cudaGetErrorString(error_id) );
      exit(EXIT_FAILURE);
   }


   // End timer
   timer_val=stopNreadTimer(&timer_val);
   printf("Processing time : %f (ms)\n",timer_val*1000.);

   output_wave_gnu(NN,u,file);

}

