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
  mat.num_init();
  return mat;
}
