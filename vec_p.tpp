// parallel_vec

template <typename T>
T& Vec<T>::operator[] (size_t idx)
{
  return values[idx];
}

template <typename T>
T Vec<T>::operator[] (size_t idx) const
{
  return values[idx];
}

template <typename T>
std::ostream & operator<< (std::ostream& os, Vec<T> vector)
{
  size_t N = vector.size();

  for (size_t i = 0; i < N-1; i++)
  {
    os << vector[i] << ", ";
  }

  os << vector[N-1];
  return os;
}

// OPERATIONS ARE PARALLEL (BASELINE):

// DOT (Vector product)
template <typename T>
T Vec<T>::operator* (const Vec<T> & other)
{
  size_t N = values.size();
  T sum;

  #pragma omp parallel shared(sum)
  #pragma omp for _schedule_ reduction(+:sum)
  for (size_t i = 0; i < N; i++)
  {
    sum += (*this)[i] * other[i];
  }

  return sum;
}


// Multiply by number
template <typename T>
Vec<T> Vec<T>::operator* (T value)
{
  size_t N = values.size();
  auto temp = Vec<T>(N);

  #pragma omp parallel for _schedule_
  for (size_t i = 0; i < N; i++)
  {
    temp[i] = (*this)[i] * value;
  }
  return temp;
}

// Vector addition
template <typename T>
Vec<T> Vec<T>::operator+ (const Vec<T> & other)
{
  size_t N = values.size();
  auto temp = Vec<T>(N);

  #pragma omp parallel for _schedule_
  for (size_t i = 0; i < N; i++)
  {
    temp[i] = values[i] + other[i];
  }
  return temp;
}

template <typename T>
Vec<T> v_abs(const Vec<T> & v)
{
  Vec<T> temp;
  # pragma omp parallel for _schedule_
  for (size_t i = 0; i < v.size(); i++)
  {
    temp[i] = abs(v[i]);
  }
  return temp;
}

// Vector subtraction
template <typename T>
Vec<T> Vec<T>::operator- (const Vec<T> & other)
{
  size_t N = values.size();
  auto temp = Vec<T>(N);

  #pragma omp parallel for _schedule_
  for (size_t i = 0; i < N; i++)
  {
    temp[i] = values[i] - other[i];
  }
  return temp;
}

template <typename T>
T L2_norm(Vec<T> v)
{
  return sqrt(v * v);
}

template <typename T>
T C_norm(Vec<T> v)
{
  auto v_ = v_abs(v);
  return sqrt(v_ * v_);
}
