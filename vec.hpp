#ifndef VEC_H_
#define VEC_H_

#include <iostream>
//#include <array>
#include <vector>
#include <cmath>

template <typename T>
class Vec
{
private:
  std::vector<T> values;
public:
  Vec(std::initializer_list<T> list) : values (list) {}

  explicit
  Vec(size_t N) : values (N) {}

  Vec(size_t N, int set_value) : values (N, set_value) {}

  size_t size() { return values.size(); }

  T& operator[] (size_t idx);
  T  operator[] (size_t idx) const;

  T operator* (Vec<T> other);
  Vec operator* (T value);
  Vec operator+ (Vec other);
  Vec operator- (Vec other);
  Vec axpy(T a, Vec x, Vec & y);
};



template <typename T>
Vec<T> operator*(T a, Vec<T> vec) { return vec * a; }

template <typename T>
std::ostream & operator<< (std::ostream& os, Vec<T> vector);

// Including template-class definitions:
#include "vec.tpp"

#endif
