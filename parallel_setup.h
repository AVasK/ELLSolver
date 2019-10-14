#ifdef OMP_H_
  #define ENABLE_PARALLEL
#elif defined _OPENMP
  #define ENABLE_PARALLEL
#endif


#ifdef ENABLE_PARALLEL
  #include "vec_p.hpp"
  #include <omp.h>
#else
  #include "vec_p.hpp"
#endif
