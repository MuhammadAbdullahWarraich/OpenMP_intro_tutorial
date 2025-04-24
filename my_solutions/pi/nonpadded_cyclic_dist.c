#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <omp.h>
static long num_steps = 100000000;
#define NUM_THREADS 4
#define NUM_TESTS 3
double step;

double cyclic_dist(double);
double get_minimum(double a, double b) {
  return (a < b ? a : b);
}

int main () {
  double start_time, run_time;
  double pi;
  step = 1.0/(double) num_steps;

  omp_set_num_threads(NUM_THREADS);

  run_time = 0.0;
  for (int j = 0; j < NUM_TESTS; j++) {
    start_time = omp_get_wtime();
    pi = cyclic_dist(step);
    run_time = (j == 0 ? omp_get_wtime() - start_time : get_minimum(run_time, omp_get_wtime() - start_time));
  }
  printf("\n pi(cyclic dist without padding) with %ld steps is %lf in %lf seconds\n ", num_steps, pi, run_time);

  return 0;
}	  

double cyclic_dist(double step) {
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

