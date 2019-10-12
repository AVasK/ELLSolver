// Ellpack-Itpack matrix

#ifndef ELLPACK_H_
#define ELLPACK_H_

#include <iostream>
#include <cmath>
#include <vector>
#include "vec.hpp" // for SpMV

namespace ELLPACK {

  template <typename T, size_t n_nonzero>
  class Slice
  {
  private:
    T * slice;

  public:
    Slice(T * addr) : slice( addr ) {}
    Slice(const Slice &) = default;

    unsigned has(unsigned j) const;

    unsigned findEmpty() const;

    inline
    T& operator[] (unsigned);

    const T& operator[] (unsigned) const;
  };


  template <typename T, size_t n_nonzero>
  class EllMatrix
  {
  private:
    T * matrix;
    //shared_ptr<T> matrix;
    size_t _N;

  public:
    EllMatrix(int N) : matrix ( new T[ N * n_nonzero ]()/*0-init*/ ), _N( N ){}

    ~EllMatrix() { delete[] matrix; }

    EllMatrix(const EllMatrix&);

    EllMatrix operator=(EllMatrix) = delete;

    Slice<T, n_nonzero> operator[] (unsigned);
    const Slice<T, n_nonzero> operator[] (unsigned) const;

    inline
    T& operator() (unsigned i, unsigned j) { return matrix[ n_nonzero*i + j ]; }

    friend
    std::ostream & operator<< (std::ostream & os, EllMatrix mat)
    {
      for (unsigned i = 0; i < mat._N; i++)
      {
        for (unsigned j = 0; j < n_nonzero; j++)
        {
          os << mat(i, j) << ' ';
        }
        os << '\n';
      }
      return os;
    }

  };


  template <typename T, size_t n_nonzero=7>
  class Ellpack
  {
  private:
    EllMatrix<T, n_nonzero> coeff;
    EllMatrix<T, n_nonzero> coord;
    size_t _N;

    // Proxy classes:

    // Row representation
    class EllRow
    {
    private:
      Ellpack<T, n_nonzero> & obj;  // refering back to Ellpack object
      unsigned i;                 // keeping the number of row

      // Proxy to distinguish between l-value & r-value uses of []
      class RowProxy
      {
      private:
        Ellpack & obj;
        unsigned i;
        unsigned j; // keeping the column number for access or modification

      public:
        RowProxy(Ellpack & _obj_, unsigned _i_, unsigned _j_) : obj(_obj_), i(_i_), j(_j_) {}
        RowProxy & operator= (T value);
        operator T() const;
      };

    public:
      EllRow(Ellpack & _obj_, unsigned _i_) : obj(_obj_), i (_i_) {}

      inline
      RowProxy operator[](unsigned j) { return RowProxy(obj, i, j); }
    };

  public:
    Ellpack(int N) : coeff (N), coord (N), _N (N) {}

    inline
    size_t getN() { return _N; }

    T sum() const;

    void num_init();  // initialize non-zeros as stated in the task.

    inline
    EllRow operator[] (unsigned i) { return EllRow(*this, i); }

    Ellpack submatrix(int row, int col, size_t size);

    inline const EllMatrix<T, n_nonzero> & getCoordMat() { return coord; }

    inline const EllMatrix<T, n_nonzero> & getCoeffMat() { return coeff; }

    // SpMV for solver
    Vec<T> operator* (Vec<T> v);
  };


  template <typename T>
  class DiagInv
  {
  private:
    std::vector<T> diag;
    size_t _N;
  public:

    size_t getN() { return _N; }

    DiagInv(size_t N) : _N ( N ), diag( std::vector<T> (N) ) {}

    template <typename T2, size_t n_nonzero>
    DiagInv(Ellpack<T2, n_nonzero> & ellp)
    {
      _N = ellp.getN();
      diag = std::vector<T> (_N);

      for (size_t i = 0; i < _N; i++)
      {
        diag[i] = 1 / ellp[i][i];
      }
    }

    friend DiagInv operator* (Vec<T> vec, DiagInv<T> diag)
    {
      return diag * vec;
    }

    Vec<T> operator* (Vec<T> vec)
    {
      auto temp = Vec<T>(_N);
      for (size_t i = 0; i < _N; i++)
      {
        temp[i] = vec[i] * diag[i];
      }
      return temp;
    }

    friend std::ostream & operator<< (std::ostream & os, DiagInv<T> self)
    {
      size_t N = self.getN();

      for ( size_t i = 0; i < N; i++ )
      {
        for ( size_t j = 0; j < N; j++ )
        {
          if ( i == j )
          {
            os << self.diag[i] << " ";
          }
          else {
            os << " " << " ";
          }
        }
        os << "\n";
      }
      return os;
    }

  };

  // operator<< for Ellpack
  template <typename TX, size_t n>
  std::ostream & operator<< (std::ostream & os, Ellpack<TX, n>);

  #include "ELLPACK.tpp"
}
#endif
