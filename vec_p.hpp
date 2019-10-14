#ifndef VECP_H_
#define VECP_H_

#include "parallel_setup.h"
#include <iostream>
#include <vector>
#include <cstring>

template <typename T>
class Vec
{
private:
  std::vector<T> values;
  size_t _N;
public:
  Vec(std::initializer_list<T> list) : values (list) { _N = values.size(); }

  explicit
  Vec(size_t N) : values (N), _N(N) {}

  Vec(size_t N, T set_value) : values (N, set_value), _N(N) {}

  size_t size() { return _N; }

  T& operator[] (size_t idx);
  T  operator[] (size_t idx) const;

  T operator* (const Vec<T> & other);
  Vec operator* (T value);
  Vec operator+ (const Vec & other);
  Vec operator- (const Vec & other);

  friend void axpy(T a, const Vec & x, const Vec & y, Vec & out)
  {
    #pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = 0; i < x._N; i++)
    {
      out[i] = a*x[i] + y[i];
    }
  }


  friend void copy(const Vec & from, Vec & to)
  {
    std::memcpy(&from.values[0], &to.values[0], sizeof(T)*from.size());
  }


  friend void scalar(T alpha, const Vec & v, Vec & out)
  {
    #pragma omp parallel for
    for (size_t i = 0; i < v._N; i++)
    {
      out[i] = v[i] * alpha;
    }
  }


  friend void axpby(T a, const Vec & x, T b, const Vec & y, Vec & out)
  {
    #pragma omp parallel for
    for (size_t i = 0; i < x._N; i++)
    {
      out[i] = a*x[i] + b*y[i];
    }
  }
};

template <typename T>
Vec<T> operator*(T a, Vec<T> vec) { return vec * a; }

template <typename T>
std::ostream & operator<< (std::ostream& os, Vec<T> vector);

// Including template-class definitions:
#include "vec_p.tpp"


#endif
