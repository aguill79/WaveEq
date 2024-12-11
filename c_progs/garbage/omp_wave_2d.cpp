// This is a modified version of homework assignment 6 for the 2D heat 
// transfer problem. This program will be modified to accomodate the 
// calculation of a 2D wave equation
// ************************************************************************
//
// Old compile instructs
/*
To compile on the ASA-X
   GNU Compiler
      g++ heat_2d_serial.cpp -o heat_2d_serial -Ofast

To execute on the ASA-X
   GNU Compiler
      run_script heat_2d_serial_gnu.sh
      where heat_2d_serial_gnu.sh is a script file that contains
         #!/bin/bash
         ./heat_2d_serial 10000 5 S 
         # execute a 10000 x 10000 point 2d-heat transfer problem 
         # for 5 iternations and suppress its output 
*/

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <cmath>
#include <omp.h>

#define TIMER_CLEAR  (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_ELAPSED   (double) (tv2.tv_usec- tv1.tv_usec)/1000000.0+(tv2.tv_sec-tv1.tv_sec)
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
struct timeval tv1,tv2;

// Global Constants

// These are the two input parameters. 
int n;                // number of non boundary condition rows in problem 
int num_iterations;   // number of successive iterations before terminating
int total_rows;       // total number of rows including 
                      // including boundary condition rows
                      
double v = 1;         // Wave velocity
double L = 10;         // Physical length and width of grid
                      // This is needed to calculate the step size
double T = 25;        // Physical time duration
                      // Used to calculate time step size

int total_cols;       // total number of columns including
                      // boundary condition columns
                      
double **wave;         // pointer to shared temperature
                      // array row-ordered storage

double *wave_buf;     // next iteration pointer to shared temperature
                      // array row-ordered storage


// I kind of like this macro style
// Old style macro to give the illusion of 2D memory
#define Wave(x,y,t) wave[t][(x)*total_cols+y] 
#define Wave_buf(x,y) wave_buf[(x)*total_cols+y] // *(temp_buf+x*total_cols+y)

double gaussian(double x, double y){

    double precision = 0.0001;

    // Place the gaussian at the center of the gridu
    double x0 = n/2;
    double y0 = n/2;
    double k = 1.0;

    double u;

    u = exp(-k*(pow(x-x0,2) + pow(y-y0,2)));

    //if(abs(u) < precision){
    //    u = 0;
    //}

    return u;
}

// This function represents the boundary conditions for rate of change at 
// t=0. For this purpose of the project, dU(x,y,0)/dx=dU(x,y,0)/dy=0
double g(double x, double y){
    return 0;
}

// routine to initialize the wave vector and the wave at the boundary
void init_wave(void) {
    //const int fireplace_start = 0.3 * (double) n;
    //const int fireplace_end = 0.7 * (double) n;

    // Set leftmost boundary condition on process 0
    for (int row=0;row < total_rows; row++) {
        for (int col=0;col<total_cols;col++) {
            Wave(row,col,0) = gaussian(row,col);

            //if (row == 0) {
            //    if (col<=fireplace_start || col > fireplace_end) {
            //        Temp(row,col) = ROOM_TEMP; // temp[row*total_cols+col];
            //    }
            //    else {
            //        Temp(row,col) = FIREPLACE_TEMP;
            //    }
            //}
            //else {
            //    Temp(row,col) = ROOM_TEMP;
            //}
        }
    }
}



// This function needs to be modified to use a single parallel region that 
// launches p threads
//
// Define the parallel region using pragma parallel and then define within 
// in it pragma for sections
void compute_wave() {

    if(num_iterations <=1) return;

    // I'm trying to figure out how to incorporate different step sizes
    // For now, we will just consider the length of the grid to equal n
    // so the step size will be 1
    //double dx = L/n;
    double dx = 1;
    double dt = T/num_iterations;

    double c = pow(v*dt/dx, 2);

    cout << "Some factors: " << endl;
    cout << "dx " << dx << endl;
    cout << "dt " << dt << endl;
    cout << "c " << c << endl;

    for (int j=1;j<=n;j++) {
        for (int k=1;k<=n;k++) {
            Wave(j,k,1)=0.5*(c*(Wave(j-1,k,1)+Wave(j+1,k,1)+
                    Wave(j,k-1,1)+Wave(j,k+1,1)) + 2*(1-2*c)
                    *Wave(j,k,1)) + dt*g(j,k);

            //Wave_buf(j,k,1)=0.5*(c*(Wave(j-1,k)+Wave(j+1,k)+
            //        Wave(j,k-1)+Wave(j,k+1)) + 2*(1-2*c)
            //        *Wave(j,k)) + dt*g(j,k);
        }
    }

    //for (int j=1;j<=n;j++) {
    //    for (int k=1;k<=n;k++) { 
    //        Wave(j,k)=Wave_buf(j,k);
    //    }
    //}

    #pragma omp parallel
    {
        int my_temp;
        
        for (int i=1;i<num_iterations;i++) {

            #pragma omp for
            for (int j=1;j<=n;j++) {
                my_temp = omp_get_thread_num();
                //cout << "Thread number:" << my_temp << endl;

                for (int k=1;k<=n;k++) {
                    Wave_buf(j,k)=0.5*(c*(Wave(j-1,k,i)+Wave(j+1,k,i)+
                            Wave(j,k-1,i)+Wave(j,k+1,i)) + 2*(1-2*c)
                            *Wave(j,k,i)) - Wave(j,k,i-1);

                    //Wave_buf(j,k)=0.5*(c*(Wave(j-1,k)+Wave(j+1,k)+
                    //        Wave(j,k-1)+Wave(j,k+1)) + 2*(1-2*c)
                    //        *Wave(j,k));
                }
            }
            
            #pragma omp for
            for (int j=1;j<=n;j++) {
                for (int k=1;k<=n;k++) { 
                    Wave(j,k,i)=Wave_buf(j,k);
                }
            }

        }
    }

}
// routine to display temperature values at each point including the 
// boundary points
void print_wave(void) {
    cout << "Wave Including Boundary Points" << endl;

    int print = num_iterations;
    for (int row=0;row<total_rows;row++) {
        for (int col=0;col<total_cols;col++) {
            cout << setw(5) << setprecision(3) << Wave(row,col,num_iterations-1) << " ";
        }
        cout << endl << endl << flush;
    }
    cout << endl;
    
}
// Routine that performs a simple 64 integer checksum
// of the binary contents of the final Temp array
// This is used to perform a quick comparison of the
// results to insure that modifications to the original
// program did not affect the accuracy of the computation
unsigned long long int checksum(void) {
    void *ptr;
    unsigned long long int sum = 0;
    for (int row=0;row<total_rows;row++) {
        for (int col=0;col<total_cols;col++) {
            ptr=(void *) &Wave(row,col,num_iterations-1);
            sum += *(unsigned long long int *) ptr;
        }
    }
    return sum;
} 

void print_usage_instructions(char *command) {
    cout << "Usage: " << command << 
        " [Dim n = num row/cols square matrices] [Number Iterations] " << 
        "<x>" << endl <<
        "   where optional argument x = " << endl <<
        "       H suppress output -- print n,runtime and" << endl <<
        "                            checksum in human readable format" << endl <<
        "       S suppress output -- print checkSum" << endl <<
        "       C suppress output -- print n,runtime" << endl <<
        "                            in CSV format" << endl <<
        "       G suppress output -- print n,runtime" << endl <<
        "                            in gnuplot format" << endl;

        exit(1);
}

int main (int argc, char *argv[]){

    if (argc!=3 && argc!=4) {
        print_usage_instructions(argv[0]);
    }

    // get total number of points not counting boundary points
    // from first command line argument 
    // Warning No Error Checking 
    n = atoi(argv[1]);

    // get total number of iterations to run simulation
    // Warning No Error Checking
    num_iterations = atoi(argv[2]);

    // set total columns plus boundary points
    total_rows = n+2; // total rows plus boundary points

    // set total columns plus boundary points
    total_cols = n+2; // total columns plus boundary points

    // dynamically allocate shared memory to
    // temp and temp_buf arrays
    //wave = new double [total_rows*total_cols]; 
    wave = new double *[num_iterations];
    for(int t=0; t<num_iterations; t++){
        wave[t] = new double[total_rows*total_cols]; 
    }

    wave_buf = new double [total_rows*total_cols];

    // initialize temperature matrix
    init_wave();

    // begin timer
    TIMER_CLEAR;
    TIMER_START;

    // compute new temps
    compute_wave(); 

    // stop timer
    TIMER_STOP;

    // print out the results if there is no suppress output argument
    if (argc==3) {
        print_wave(); // print out the results
        // print time in normal human readable format
        cout << "Execution Time = " << TIMER_ELAPSED << " Seconds"
             << endl;
        cout << "64 bit Checksum = " << checksum() << endl;
    }
    // if there exists a 4th argument, then suppress the output
    else {
        // print time in gnuplot format
        if (*argv[3]=='G') {
            cout << n << " " << TIMER_ELAPSED << endl;
        }
        // print time in CSV format 
        else if (*argv[3]=='C') {
            cout << n << "," << TIMER_ELAPSED << endl;
        }
        // print 64 bit checkSum
        else if (*argv[3]=='S') {
            cout << "64 bit Checksum = " << checksum() << endl;
        }
        else if (*argv[3]=='H') {
           // print time and Checksum in normal human readable format
           cout << "Number of active data points =" << n << endl;
           cout << "Execution Time = " << TIMER_ELAPSED << " Seconds"
                << endl;
            cout << "64 bit Checksum = " << checksum() << endl;
        }
    }

    for(int t=0; t<num_iterations; t++){
        delete [] wave[t];
    }

    delete [] wave;
    delete wave_buf;
}
