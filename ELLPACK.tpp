// ELLPACK classmethod implementation
#include <stdexcept>



using namespace ELLPACK;

// SLICE:

template <typename T, size_t n_nonzero>
T& Slice<T, n_nonzero>::operator[] (unsigned j)
{
  if (j < n_nonzero && j >= 0)
  {
    return slice[j];
  }
  throw std::out_of_range("Out of bounds for slice");
}

template <typename T, size_t n_nonzero>
const T& Slice<T, n_nonzero>::operator[] (unsigned j) const
{
  if (j < n_nonzero && j >= 0)
  {
    return slice[j];
  }
  throw std::out_of_range("Out of bounds for slice");
}

// returns an index at slice, where slice[idx] == j, -1 otherwise
template <typename T, size_t n_nonzero>
unsigned Slice<T, n_nonzero>::has (unsigned j) const
{
  for (unsigned idx = 0; idx < n_nonzero; idx++)
  {
    if ( slice[idx] == j )
    {
      return idx;
    }
  }
  return -1;
}


template <typename T, size_t n_nonzero>
unsigned Slice<T, n_nonzero>::findEmpty() const
{
  return this->has(0);
}

// EllMatrix:

template <typename T, size_t n_nonzero>
Slice<T, n_nonzero> EllMatrix<T, n_nonzero>::operator[] (unsigned i)
{
  return Slice<T, n_nonzero>( &matrix[ n_nonzero*i ] );
}

template <typename T, size_t n_nonzero>
const Slice<T, n_nonzero> EllMatrix<T, n_nonzero>::operator[] (unsigned i) const
{
  return Slice<T, n_nonzero>( &matrix[ n_nonzero*i ] );
}

template <typename T, size_t n_nonzero>
EllMatrix<T, n_nonzero>::EllMatrix(const EllMatrix& old)
{
  _N = old._N;
  matrix = new T[n_nonzero * _N];
  /*std::*/memcpy(matrix, old.matrix, old._N * n_nonzero * sizeof(T));
}


// Ellpack:

template <typename T, size_t n_nonzero>
void Ellpack<T, n_nonzero>::num_init() // initialize non-zeros as stated in the task.
{
  auto diag = std::vector<std::tuple<int, int>>();
  T sum = 0; // running sum of off-diagonal elements

  for (unsigned i = 0; i < _N; i++)
  {
    for (unsigned j = 0; j < n_nonzero; j++)
    {
      unsigned j_real;
      if ( (j_real = coord[i][j]) != i )
      {
        if ( coeff[i][j] != 0 )
        // init off-diagonal
          sum += std::fabs( coeff[i][j] = T(std::cos(i*j_real + 3.14)) ); // what is 3.14 needed for?
      }
      else
      {
        // push (i,j)
        diag.push_back(std::tuple<int,int>(i, j));
      }
    }
  }

  for (auto coord : diag)
  {
    int i, j;
    std::tie(i, j) = coord;

    coeff[i][j] = T(1.5 * sum);
  }
}

template <typename T, size_t n_nonzero>
Ellpack<T, n_nonzero> Ellpack<T, n_nonzero>::submatrix(int row, int col, size_t size)
{
  auto res = Ellpack<T, n_nonzero>(size);

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      res[i][j] = (T) (*this)[i+row][j+col];
    }
  }
  return res;
}


template <typename T, size_t n_nonzero>
T Ellpack<T, n_nonzero>::sum() const
{
  T s = 0;
  for (int row = 0; row < _N; row++)
  {
    for (int col = 0; col < n_nonzero; col++)
    {
      if (coeff[row][col] > 0)
      {
        s += coeff[row][col];
      }
    }
  }
  return s;
}

template <typename T, size_t n_nonzero>
std::ostream & operator<< (std::ostream & os, Ellpack<T, n_nonzero> m)
{
  // If the matrix is huge, don't print it all,
  // just the upperleft corner.
  if ( m.getN() > 40 )
  {
    os << "_HEAD_____\n";
    os << m.submatrix(0,0,40);
    return os;
  }

  for (unsigned i = 0; i < m.getN(); i++)
  {
    for (unsigned j = 0; j < m.getN(); j++)
    {
      os.precision(3); // DEBUG
      os << m[i][j] << " "; // DEBUG: Replace with ' ';
    }
    os << "\n";
  }
  return os;
}

// SpMV
template <typename T, size_t n_nonzero>
Vec<T> Ellpack<T, n_nonzero>::operator* (const Vec<T> & v_in)
{
  /*
  if ( _N != v_in.size() )
  {
    throw std::logic_error("sizes of matrix & vector don't match");
  }
  */

  Vec<T> res(_N); // uninitialized result vec.

  #pragma omp parallel for
  for ( int row = 0; row < _N; row++ )
  {
    T sum = 0; // running sum for rows

    for ( int col = 0; col < n_nonzero; col++ )
    {
      auto val = coeff[row][col];
      auto idx = coord[row][col];

      sum += v_in[idx] * val;
    }
    // row(th) element of vector is = to sum.
    res[row] = sum;
  }
  return res;
}


template <typename T, size_t n_nonzero>
Ellpack<T, n_nonzero>::EllRow::RowProxy::operator T() const
{
  //std::cout << "Getting an item @ " << i << ", " << j << '\n';
  unsigned j_idx;
  if ( (j_idx = obj.coord[i].has(j)) != -1 )
  {
    // non-zero entry
    return obj.coeff[i][j_idx];
  }
  else
  {
    // zero
    return 0;
  }
}


template < typename T, size_t n_nonzero >
typename Ellpack<T, n_nonzero>::EllRow::RowProxy &
Ellpack<T, n_nonzero>::EllRow::RowProxy::operator= (T value)
{
  //std::cout << "Changing value @" << i << ", " << j << " to " << value << "?\n";
  // finding empty place in Slice to append a new value
  int j_idx;

  if ( (j_idx = obj.coord[i].has(j)) != -1 )
  {
    // Such an index already exists.
    obj.coeff[i][j_idx] = value;
  }
  else if ( (j_idx = obj.coeff[i].findEmpty()) != -1 )
  {
    // Inserting new
    obj.coord[i][j_idx] = j;
    obj.coeff[i][j_idx] = value;
  }
  else
  {
    throw std::out_of_range("Too many elements per row");
  }
  return *this;
}
