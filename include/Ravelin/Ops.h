/// Class for general operations
class OPS
{
  public:
    static bool rel_equal(REAL x, REAL y);
    static bool rel_equal(REAL x, REAL y, REAL tol);

    /// Does a axpy operation
    template <class V1, class V2>
    static V2& axpy(REAL alpha, const V1& v1, V2& v2)
    {
      if (v1.size() != v2.size())
        throw MissizeException();

      CBLAS::axpy(v1.size(), v1.data(), v1.inc(), v2.data(), v2.inc()); 
      return v2;
    }

    /// Does a transposition operation
    template <class M1, class M2>
    static M2& transpose(const M1& m1, M2& m2)
    {
      // resize m2
      m2.resize(m1.columns(), m1.rows());

      // setup necessary constants
      const unsigned LD1 = m1.leading_dim();
      const unsigned LD2 = m2.leading_dim();
      const REAL* d1 = m1.data(); 
      REAL* d2 = m2.data(); 

      // do the transposition
      for (unsigned i=0, i1=0; i< m1.rows(); i++)
      {
        CBLAS::copy(m1.columns(), d1+i, LD1, d2+i1, 1); 
        i1 += LD2;
      }

      return m2;
    }

    /// Multiples a matrix by the transpose of a matrix or vector
    template <class T, class U, class V>
    static V& mult_transpose(const T& x, const U& y, V& z)
    {
      // verify that we can multiply these
      if (x.columns() != y.columns())
        throw MissizeException();

      // resize the result
      z.resize(x.rows(), y.rows());

      // do the multiplication
      CBLAS::gemm(CblasNoTrans, CblasTrans, x.rows(), y.rows(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

      return z;
    }

    /// Multiples a matrix by a matrix or vector
    template <class T, class U, class V>
    static V& mult(const T& x, const U& y, V& z)
    {
      // verify that we can multiply these
      if (x.columns() != y.rows())
        throw MissizeException();

      // resize the result
      z.resize(x.rows(), y.columns());

      // do the multiplication
      if (y.columns() > 1)
        // matrix
        CBLAS::gemm(CblasNoTrans, CblasNoTrans, x.rows(), y.columns(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);
      else
        // vector
        CBLAS::gemv(CblasColMajor, CblasNoTrans, x.rows(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.inc(), 0.0, z.data(), z.inc());

      return z;
    }

    /// Multiples the transpose of a matrix or vector by a matrix or vector
    template <class T, class U, class V>
    static V& transpose_mult(const T& x, const U& y, V& z)
    {
      // verify that we can multiply these
      if (x.rows() != y.rows())
        throw MissizeException();

      // resize the result
      z.resize(x.columns(), y.columns());

      // multiply
      CBLAS::gemm(CblasTrans, CblasNoTrans, x.columns(), y.columns(), x.rows(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

      return z;
    }

    /// Multiples a matrix by a matrix or vector
    template <class T, class U, class V>
    static V& transpose_mult_transpose(const T& x, const U& y, V& z)
    {
      // verify that we can multiply these
      if (x.rows() != y.columns())
        throw MissizeException();

      // resize the result
      z.resize(x.columns(), y.rows());

      // do the multiplication
      CBLAS::gemm(CblasTrans, CblasTrans, x.columns(), y.rows(), x.rows(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

      return z;
    }

    template <class U, class V, class M>
    static M& outer_prod(const U& x, const V& y, M& z)
    {
      // make sure both are vectors
      if (x.columns() != 1 || y.columns() != 1)
        throw MissizeException();

      // get n and m
      unsigned m = x.rows();
      unsigned n = y.rows();
    
      // resize the matrix and set to zero
      z.set_zero(m,n);

      // call BLAS outer product routine
      if (m > 0 && n > 0)
        CBLAS::ger(CblasColMajor, m, n, x.data(), x.inc(), y.data(), y.inc(), z.data(), z.leading_dim());
    
      return z;
    }

}; // end class 

