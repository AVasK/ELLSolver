#include <omp.h>

#include <iostream>

int main()
{
  omp_set_num_threads(8);
  long double sum = 0;
  auto start = omp_get_wtime();
  #pragma omp parallel shared(sum)
  {
    #pragma omp for _schedule_ reduction(+ : sum)
    for (long long i = 0; i < 50000000; i++)
    {
      sum += i;
    }
  }
  auto time = omp_get_wtime() - start;
  std::cout << "sum = " << sum << " time = " << time << "\n";
}
