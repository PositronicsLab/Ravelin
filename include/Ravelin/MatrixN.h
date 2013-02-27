/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef MATRIXN
#error This class is not to be included by the user directly. Use MatrixNd.h or MatrixNf.h instead.
#endif

/// A generic, possibly non-square matrix.  
/**
 * The underlying data is stored in column-major format (e.g., the element at row 1, column 0 is the second element in the flat array).
 */
class MATRIXN
{
  friend class SHAREDMATRIXN;

  public:
    MATRIXN();
    MATRIXN(unsigned rows, unsigned columns);
    MATRIXN(const MATRIXN& source);
    MATRIXN(unsigned rows, unsigned columns, const REAL* array);
    MATRIXN(const VECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN(const MATRIX3& m);
//    MATRIXN(const POSE& m);
    MATRIXN& set_identity();
    MATRIXN& set_identity(unsigned sz);
    bool is_symmetric(REAL tolerance = -1.0) const;
    static MATRIXN identity(unsigned dim);
    REAL norm_inf() const;
    static MATRIXN construct_variable(unsigned rows, unsigned cols, ...);
    virtual ~MATRIXN() { }
    MATRIXN& zero_upper_triangle();
    MATRIXN& zero_lower_triangle();
    static VECTORN& diag_mult(const VECTORN& d, const VECTORN& v, VECTORN& result);
    static MATRIXN& diag_mult(const VECTORN& d, const MATRIXN& m, MATRIXN& result);
    static MATRIXN& diag_mult_transpose(const VECTORN& d, const MATRIXN& m, MATRIXN& result);
    MATRIXN& select_square(const std::vector<bool>& indices, MATRIXN& result) const;
    MATRIXN& set(const VECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN& operator=(const MATRIX3& m);
//    MATRIXN& set(const POSE& m);
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    virtual MATRIXN& resize(unsigned rows, unsigned columns, bool preserve = false);
    MATRIXN& remove_row(unsigned i);
    MATRIXN& remove_column(unsigned i);
    
    static MATRIXN zero(unsigned rows, unsigned columns);
    MATRIXN& negate();
    MATRIXN& set_zero();
    virtual MATRIXN& transpose();
    static MATRIXN& transpose(const MATRIXN& m, MATRIXN& result);
    virtual MATRIXN& operator=(const MATRIXN& source);
    static MATRIXN& mult(const MATRIXN& m1, const MATRIXN& m2, MATRIXN& result) { return m1.mult(m2, result); }
    MATRIXN& operator+=(const MATRIXN& m);
    MATRIXN& operator-=(const MATRIXN& m);
    MATRIXN& operator/=(REAL scalar);
    MATRIXN& operator*=(REAL scalar);
    REAL* data() { return _data.get(); }
    const REAL* data() const { return _data.get(); }    
    void compress();
    unsigned leading_dim() const { return _rows; }
    unsigned inc() const { return 1; }

    /// Sets this to a m x n sized zero matrix
    MATRIXN& set_zero(unsigned m, unsigned n) { return resize(m,n).set_zero(); }

    REAL operator()(unsigned i, unsigned j) const;
    REAL& operator()(const unsigned i, const unsigned j);

    #define MATRIXX MATRIXN
    #include "MatrixCommon.inl"
    #undef MATRIXX
    #define XMATRIXN MATRIXN
    #include "XMatrixN.inl"
    #undef XMATRIXN

  protected:
    SAFESTATIC FastThreadable<MATRIXN> _n;
    SAFESTATIC FastThreadable<VECTORN> _workv;

    boost::shared_array<REAL> _data;
    unsigned _rows;
    unsigned _columns;
    unsigned _capacity;
}; // end class

std::ostream& operator<<(std::ostream& out, const MATRIXN& m);
std::istream& operator>>(std::istream& in, MATRIXN& m);

