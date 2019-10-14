// Timing a function
#include "parallel_setup.h"
#include <ctime>

/*
double now()
{
  return (double) clock() / CLOCKS_PER_SEC;
}
*/


#ifdef ENABLE_PARALLEL
  double now()
  {
    return omp_get_wtime();
  }
#else
  double now()
  {
    return (double) clock() / CLOCKS_PER_SEC;
  }
#endif


using CG_solver_f = Vec<double> (*)(Ellpack<double, 7>, Vec<double>, long, double);

double timeit(CG_solver_f f, Ellpack<double, 7> ellp, Vec<double> vec, long MAX_ITER, double r_tol)
{
  auto start = now();
  (*f)(ellp, vec, MAX_ITER, r_tol);
  auto time = now() - start;
  return time;
}
