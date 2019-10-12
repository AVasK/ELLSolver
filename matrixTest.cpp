#include <iostream>
#include "ELLPACK.hpp"

using namespace ELLPACK;

template <typename T>
Ellpack<T> makerELL(int Nx, int Ny, int Nz)
{
  long N = Nx*Ny*Nz;
  Ellpack<T> mat(N);
  long i; // poitions in matrix
  long int Nx_Ny = Nx * Ny; // saving result of potentially heavy computation to avoid in-loop computation
  // POSSIBLE PROBLEM: CACHE MISS, THEN OPTIMIZATION WILL MAKE THINGS WORSE.

  for (int X = 0; X < Nx; X++)
  {
      for (int Y = 0; Y < Ny; Y++)
      {
          for (int Z = 0; Z < Nz; Z++)
          {
              i = Z*(Nx_Ny) + Y*Nx + X; // potential for loop unrolling & fma?
              mat[i][i] = 1;
              // check Z coord
              if (Z > 0)    { mat[i][i-Nx_Ny] = 1; }
              if (Z < Nz-1) { mat[i][i+Nx_Ny] = 1; }
              // check Y coord
              if (Y > 0)    { mat[i][i-Nx] = 1; }
              if (Y < Ny-1) { mat[i][i+Nx] = 1; }
              // check X coord
              if (X > 0)    { mat[i][i-1] = 1; }
              if (X < Nx-1) { mat[i][i+1] = 1; }
          }
      }
  }
  return mat;
}



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
