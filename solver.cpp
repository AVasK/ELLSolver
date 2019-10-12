// Solver

#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include "ELLPACK.hpp"
#include "vec.hpp"

using namespace ELLPACK;

template <typename T>
Vec<T> CG_solver(Ellpack<T,7> A, Vec<T> b, long MAX_ITER, T r_tol)
{
  size_t N = A.getN();
  Vec<T> x_0 (N, 0.0);
  auto r_prev = b - A * x_0;

  auto convergence = false;
  auto k = 1;
  Vec<T> p_k (N);
  Vec<T> p_prev (N);
  T rho_k (N);
  T rho_prev (N);
  Vec<T> z_k (N);
  Vec<T> r_k (N);
  Vec<T> q_k (N);
  Vec<T> x_k (N);
  Vec<T> x_prev (N);

  do
  {
    z_k = DiagInv<T>(A) * r_prev;
    rho_prev = rho_k;
    rho_k = r_prev * z_k;

    p_prev = p_k;

    if ( k == 1 )
    {
      p_k = z_k;
    }
    else
    {
      auto beta_k = rho_k / rho_prev;
      p_k = z_k + beta_k * p_prev;
    }

    q_k = A * p_k;
    auto alpha_k = rho_k / (p_k * q_k);

    x_prev = x_k;
    x_k = x_prev + alpha_k * p_k;

    r_prev = r_k;
    r_k = r_prev - alpha_k * q_k;

    if (( rho_k < r_tol ) or ( k <= MAX_ITER))
    {
      convergence = true;
    }
    else {
      k += 1;
    }
  }
  while ( convergence != true );
  return p_k;
}

#endif
