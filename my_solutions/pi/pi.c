/*

   This program will numerically compute the integral of

   4/(1+x*x) 

   from 0 to 1.  The value of this integral is pi -- which 
   is great since it gives us an easy way to check the answer.

   The is the original sequential program.  It uses the timer
   from the OpenMP runtime library

History: Written by Tim Mattson, 11/99.

*/
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <omp.h>
static long num_steps = 100000000;
#define NUM_THREADS 4
#define NUM_TESTS 3
#define PAD (sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(double))
//(64/sizeof(double) - 1) // because i have 2 cores with 2 threads each, i thought that each core should have 2 values in its cache line. see info.txt for more information
double step;

double block_dist_without_padding(double);
double block_dist(double, int);
double cyclic_dist_without_padding(double);
double cyclic_dist(double, int);
double critical_section_cyclic_dist(double);

double get_minimum(double a, double b) {
  return (a < b ? a : b);
}
int main ()
{
//  double pi, sum = 0.0;
  double start_time, run_time;
  double pi;
  //int pad[] = {PAD-2, PAD-1, PAD, 1+PAD};
  int pad[] = {PAD};
  step = 1.0/(double) num_steps;

  printf("default padding value is: %ld", PAD);

  //omp_set_num_threads(NUM_THREADS);

  run_time = 0.0;
  start_time = omp_get_wtime();
  for (int j = 0; j < NUM_TESTS; j++) {
    start_time = omp_get_wtime();
    pi = block_dist_without_padding(step);
    run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
  }
  printf("\n pi(block dist without padding) with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);

  run_time = 0.0;
  for (int i = 0; i < sizeof(pad)/sizeof(int); i++) {
    for (int j = 0; j < NUM_TESTS; j++) {
      start_time = omp_get_wtime();
      pi = block_dist(step, pad[i]);
      run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
    }
    printf("\n pi(block dist with padding: %d) with %ld steps is %lf in %lf seconds\n ",pad[i],num_steps,pi,run_time);
  } 

  run_time = 0.0;
  start_time = omp_get_wtime();
  for (int j = 0; j < NUM_TESTS; j++) {
    start_time = omp_get_wtime();
    pi = cyclic_dist_without_padding(step);
    run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
  }
  printf("\n pi(cyclic dist without padding) with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);

  //start_time = omp_get_wtime();
  //pi = cyclic_dist_without_padding(step);
  //run_time = omp_get_wtime() - start_time;
  run_time = 0.0;
  for (int i = 0; i < sizeof(pad)/sizeof(int); i++) {
    for (int j = 0; j < NUM_TESTS; j++) {
      start_time = omp_get_wtime();
      pi = cyclic_dist(step, pad[i]);
      run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
    }
    printf("\n pi(cyclic dist with padding: %d) with %ld steps is %lf in %lf seconds\n ",pad[i],num_steps,pi,run_time);
  }
  run_time = 0.0;
  start_time = omp_get_wtime();
  for (int j = 0; j < NUM_TESTS; j++) {
    start_time = omp_get_wtime();
    pi = critical_section_cyclic_dist(step);
    run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
  }
  printf("\n pi(cyclic dist with critical section) with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
}	  
double block_dist_without_padding(double step) {
  double sums[NUM_THREADS] = {0.0};
  //double sums[NUM_THREADS][PAD];
  double pi, sum = 0.0;
//  bool flags[NUM_THREADS] = {false};
  int nthreads;
//  for (int i = 0; i < NUM_THREADS; i++) { sums[i] = 0.0; flags[i] = 0; }
#pragma omp parallel
  {
    double x = 0.0;
    int num_threads = omp_get_num_threads();
    int thread_i = omp_get_thread_num();
    if (thread_i == 0) nthreads = num_threads;

    int nsteps = num_steps / num_threads;
    int rsteps = num_steps % num_threads;
    nsteps += (thread_i < rsteps ? 1 : 0);

    int offset = nsteps * thread_i;
    offset += (thread_i < rsteps ? thread_i : rsteps);

    double lsum = 0.0;
    for (int i=offset;i< offset + nsteps; i++){
      x = (i-0.5)*step;
      lsum += 4.0/(1.0+x*x);
    }
    sums[thread_i] = lsum;
    //sums[thread_i][0] = lsum;
//    flags[thread_i] = true;
  }
//  for (int i = 0; i < NUM_THREADS; i++) if (flags[i] == 1) sum += sums[i];
  for (int i = 0; i < nthreads; i++) sum += sums[i];
  //for (int i = 0; i < nthreads; i++) sum += sums[i][0];
  pi = step * sum;
  return pi;
}

double block_dist(double step, int pad) {
  //double sums[NUM_THREADS] = {0.0};
  double sums[NUM_THREADS][pad]; // 2D arrays in C are stored in row-major order
  double pi, sum = 0.0;
//  bool flags[NUM_THREADS] = {false};
  int nthreads;
//  for (int i = 0; i < NUM_THREADS; i++) { sums[i] = 0.0; flags[i] = 0; }
#pragma omp parallel
  {
    double x = 0.0;
    int num_threads = omp_get_num_threads();
    int thread_i = omp_get_thread_num();
    if (thread_i == 0) nthreads = num_threads;

    int nsteps = num_steps / num_threads;
    int rsteps = num_steps % num_threads;
    nsteps += (thread_i < rsteps ? 1 : 0);

    int offset = nsteps * thread_i;
    offset += (thread_i < rsteps ? thread_i : rsteps);

    double lsum = 0.0;
    for (int i=offset;i< offset + nsteps; i++){
      x = (i-0.5)*step;
      lsum += 4.0/(1.0+x*x);
    }
    //sums[thread_i] = lsum;
    sums[thread_i][0] = lsum;
//    flags[thread_i] = true;
  }
//  for (int i = 0; i < NUM_THREADS; i++) if (flags[i] == 1) sum += sums[i];
  //for (int i = 0; i < nthreads; i++) sum += sums[i];
  for (int i = 0; i < nthreads; i++) sum += sums[i][0];
  pi = step * sum;
  return pi;
}
double cyclic_dist_without_padding(double step) {
  double sums[NUM_THREADS] = {0.0};
//  bool flags[NUM_THREADS] = {false};
  double pi;
  int nthreads;
  #pragma omp parallel
  {
    double x = 0.0;
    int thread_i = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    if (thread_i == 0) nthreads = num_threads;

//    sums[thread_i][0] = 0.0;
    for (int i = thread_i; i < num_steps; i+= num_threads) {
      x = (i-0.0) * step;
      sums[thread_i] += 4.0/(1.0 + x*x);
    }
 //   flags[thread_i] = true;
  }
  pi = 0.0;
//  for (int i = 0; i < NUM_THREADS; i++) if (flags[i] == true) pi += sums[i];
  for (int i = 0; i < nthreads; i++) pi += sums[i];
 
  pi *= step;
  return pi;
}

double cyclic_dist(double step, int pad) {
  //double sums[NUM_THREADS] = {0.0};
  double sums[NUM_THREADS][pad];
//  bool flags[NUM_THREADS] = {false};
  double pi;
  int nthreads;
  #pragma omp parallel
  {
    double x = 0.0;
    int thread_i = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    if (thread_i == 0) nthreads = num_threads;

    sums[thread_i][0] = 0.0;
    for (int i = thread_i; i < num_steps; i+= num_threads) {
      x = (i-0.0) * step;
      //sums[thread_i] += 4.0/(1.0 + x*x);
      sums[thread_i][0] += 4.0/(1.0 + x*x);
    }
 //   flags[thread_i] = true;
  }
  pi = 0.0;
//  for (int i = 0; i < NUM_THREADS; i++) if (flags[i] == true) pi += sums[i];
  //for (int i = 0; i < nthreads; i++) pi += sums[i];
  for (int i = 0; i < nthreads; i++) pi += sums[i][0];
 
  pi *= step;
  return pi;
}
double critical_section_cyclic_dist(double step) {
  double pi = 0.0;
  #pragma omp parallel
  {
    double x = 0.0;
    double psum = 0.0;
    int thread_i = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    for (int i = thread_i; i < num_steps; i += num_threads) {
      x = (i-0.0) * step;
      psum += 4.0 / (1.0 + x*x);
    }
  #pragma omp critical
    pi += psum;
  }
  pi *= step;
  return pi;
}
/*
 
#include <stdio.h>
#include <omp.h>
static long num_steps = 100000000;
double step;
int main ()
{
	  int i;
	  double x, pi, sum = 0.0;
	  double start_time, run_time;

	  step = 1.0/(double) num_steps;

        	 
	  start_time = omp_get_wtime();
	  for (i=1;i<= num_steps; i++){
		  x = (i-0.5)*step;
		  sum = sum + 4.0/(1.0+x*x);
	  }
  
	  pi = step * sum;
	  run_time = omp_get_wtime() - start_time;
	  printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
}	  





*/



