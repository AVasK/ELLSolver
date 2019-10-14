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

#define _schedule_ schedule(guided)
// static, 100 : 3.88395 s.
// static : 2.70023 s.
// dynamic: 2.90254 s.
// dyn, 100: 1.95408 s.
// dyn, 1000:  1.79846 s.
// guided :  1.77173 s.
// sequential : 3.75919 s.
