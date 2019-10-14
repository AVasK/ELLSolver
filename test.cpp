#include <iostream>
#include "ELLPACK.hpp"
#include "parallel_setup.h"
#include "gen.hpp"
#include "solver.cpp"
#include "timing.hpp"

using namespace ELLPACK;

int main()
{
  #ifdef ENABLE_PARALLEL
  omp_set_num_threads(4);
  std::cout << "[ OMP ENABLED ]\n";
  #endif

  /*
  auto ellp = Ellpack<double>(4); // 4x4 matrix

  ellp[0][0] = 1;
  ellp[0][1] = 2;
  ellp[1][1] = 1;
  ellp[1][2] = 1;
  ellp[2][2] = 1;
  ellp[3][3] = 1;
  std::cout << ellp << "\n * \n";
  Vec<double> vec = {1,1,1,1};
  std::cout << vec;
  std::cout << "\n = \n" << ellp * vec << "\n";


  auto d = DiagInv<double>(ellp);
  std::cout << "\n" << d;
  std::cout << d * vec << "\n";

  std::cout << "SOLVING...";
  Vec<double> b = {3,2,1,1};

  */

  auto N = 200;
  auto ellp = makerELL<double>(N,N,N);
  auto b = Vec<double>(N*N*N, 100.3);
  //auto c = Vec<double>(N*N*N, 2.3);

  auto t = timeit( &CG_solver<double>, ellp, b, 10000000, 0.0000001 );
  //std::cout << "\n" << v << "\n";
  std::cout << "Took " << t << " s.\n";
  std::cout << "Done!";
}
