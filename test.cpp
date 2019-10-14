#include <iostream>
#include "ELLPACK.hpp"
#include "parallel_setup.h"
#include "gen.hpp"
#include "solver.cpp"
#include "timing.hpp"

using namespace ELLPACK;
using namespace std;

auto NX = 100;
auto NY = 100;
auto NZ = 100;
auto TOL=0.001;
auto MAXIT = 1000;
auto NTHREADS = 4;
auto TEST_BASE_OPS = false;


pair<int,double> compare(string arg, vector<string> opt)
{
  string str, s_num;
  int idx = -1;
  for (int i = 0; i < arg.size(); i++)
  {
    if (idx != -1) s_num += arg[i];

    if (arg[i] == '=')
    {
      idx = i;
    }
    if (idx == -1) str += arg[i];

  }

  if (idx == -1)
  {
    return pair<int, double>(-1, 0);
  }

  for (auto i = 0; i < opt.size(); i++)
  {
    if (str == opt[i])
    {
      auto num = stof(s_num, nullptr);
      return pair<int, double>(i, num);
    }
  }

  return pair<int, double>(-1, 0);
}

int main(int argc, char * argv[])
{

  std::vector<string> options = {
    "nx",
    "ny",
    "nz",
    "tol",
    "maxit",
    "nt",
    "qa"
  };
  if (argc == 1)
  {
    std::cout << "Defaults used. Provide arguments as follows: \n"
              << "nx=<int> ny=<int> nz=<int> размер решетки топологии для генерации\n"
              << "tol=<double> невязка относительно нормы правой части\n"
              << "maxit=<int>  максимальное число итераций\n"
              << "nt=<int>     число нитей"
              << "qa           флаг тестирования базовых операций\n";
  }

  for (int arg = 1; arg < argc; arg++)
  {
    auto s_arg = string(argv[arg]);
    if (s_arg == "qa")
    {
      TEST_BASE_OPS = true;
    }
    else
    {
      auto opt = compare(s_arg, options);
      int opt_n; double num;
      std::tie(opt_n, num) = opt;

      std::cout << "Set " << options[opt_n] << " = " << num << "\n";

      switch(opt_n)
      {
          case -1:
            std::cout << "F'd up\n";
            exit(1);
          case 0:
            NX = num;
          case 1:
            NY = num;
          case 2:
            NZ = num;
          case 3:
            TOL = num;
          case 4:
            MAXIT = num;
          case 5:
            NTHREADS = num;
      }
    }
  }


  #ifdef ENABLE_PARALLEL
  omp_set_num_threads(NTHREADS);
  std::cout << "[ OMP ENABLED ]\n";
  #endif


  auto ellp = makerELL<double>(NX, NY, NZ);
  //std::cout << ellp << "\n";
  auto b = GenCos<double>(NX * NY * NZ);
  //auto c = Vec<double>(N*N*N, 2.3);

  auto t = timeit( &CG_solver<double>, ellp, b, MAXIT, TOL );
  auto v = CG_solver(ellp, b, MAXIT, TOL);
  std::cout << "\n" << L2_norm(v) << "\n";
  std::cout << "Took " << t << " s.\n";
  std::cout << "Done!";


  /* TESTING WITH MANUALLY PICKED VALUES
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

}
