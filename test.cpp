#include <iostream>
#include "ELLPACK.hpp"
#include "vec.hpp"
#include "solver.cpp"

using namespace ELLPACK;

int main()
{
  auto ellp = Ellpack<double>(4); // 4x4 matrix
  ellp[0][0] = 1;
  ellp[1][1] = 2;
  ellp[2][2] = 3;
  ellp[3][3] = 4;
  std::cout << ellp << "\n * \n";
  Vec<double> vec = {2,2,2,2};
  std::cout << vec;
  std::cout << "\n = \n" << ellp * vec << "\n";

  auto d = DiagInv<double>(ellp);
  std::cout << "\n" << d;
  std::cout << d * vec << "\n";

  std::cout << "SOLVING...";
  Vec<double> b = {2,4,6,8};
  auto v = CG_solver(ellp, b, 100, 0.001);
  std::cout << "\n" << v << "\n";
}
