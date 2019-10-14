#include <iostream>
#include "ELLPACK.hpp"
#include "gen.hpp"

using namespace ELLPACK;

int main()
{
  /*
  auto mat = EllMatrix<int> (3);
  std::cout << mat;
  */

  /*
  //Ellpack<int> ellmat (5);
  auto ellmat = Ellpack<int>(2000);
  ellmat[1][1] = 1;
  ellmat[1][2] = 2;
  //std::cout << ellmat[1][1] << '\n';
  std::cout << ellmat;
  std::cout << ellmat.getCoeffMat() << "\n\n";
  std::cout << ellmat.getCoordMat() << "\n";
  */

  //auto ellp = makerELL<double>(2, 2, 2);
  //std::cout << ellp << "\n---\n" << ellp.getCoeffMat() << "\n";
  //ellp.num_init();
  //std::cout << "\n";
  //std::cout << "S = " << ellp.sum() << "\n";
  //std::cout << ellp.submatrix(0, 0, 7);

  auto ellp = Ellpack<int>(4); // 4x4 matrix
  ellp[0][0] = 1;
  ellp[0][3] = 1;
  ellp[1][1] = 2;
  ellp[1][2] = 2;
  ellp[2][0] = 2;
  std::cout << ellp << "\n * \n";
  Vec<int> vec = {1,2,2,1};
  std::cout << vec;
  std::cout << "\n = \n" << ellp * vec << "\n";

  /*
  ellp[3][0] = 7;
  ellp[3][3] = 7;
  ellp[3][0] = 3;
  std::cout << ellp;
  */

}
