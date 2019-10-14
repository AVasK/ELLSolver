#include "parallel_setup.h"

 int main()
 {
   Vec<int> a = {1,1,2,1};
   auto b = Vec<int>(4);
   Vec<int> c (4, 2);
   auto d = c;
   d[2] = 7;

   std::cout << a << "\n" << b << "\n" << c << "\n" << d << "\n";

   std::cout << "a * 2 = " << a * 2 << "\n";

   axpy(2, a, c);
   std::cout << "axpy(2, x, y) = " << c << "\n";
 }
