#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>

double get_time()
{
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}


#include <math.h>
#include "puma560.c"

#define RUNS 1000000

double T[4*4];
double J[6*NDOF];
double tau[NDOF];
double M[NDOF*NDOF];
double c[NDOF];
double g[NDOF];
double parms[NDYNPARM];
double q[NDOF];
double dq[NDOF];
double ddq[NDOF];

int main(){


  int runs_init = 1;
  int runs_mult = 10;
  int runs_max = RUNS;

  double t0;
  double tf;
  int i;
  int runs;

  // Warm-up
  t0 = get_time();
  for(i=0;i<10;++i){
      func( T, J, M, c, g, parms, q, dq );
  }
  tf = get_time();

  // Benchmark
  for(runs=runs_init; runs<=runs_max; runs*=runs_mult){
    t0 = get_time();
    for(i=0;i<runs;++i){
        func( T, J, M, c, g, parms, q, dq );
    }
    tf = get_time();
    printf("%10d %f\n", runs, (tf-t0));
  }


  return 0;
}
