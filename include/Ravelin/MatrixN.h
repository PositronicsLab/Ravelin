/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef MATRIXN
#error This class is not to be included by the user directly. Use MatrixNd.h or MatrixNf.h instead.
#endif

class MATRIX2;
class MATRIX3;

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
    MATRIXN(unsigned rows, unsigned columns, const REAL* array);
    MATRIXN(const VECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN(const SHAREDVECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN(const CONST_SHAREDVECTORN& v, Transposition trans = eNoTranspose);
//    MATRIXN(const POSE& m);
    MATRIXN(const MATRIX2& m);
    MATRIXN(const MATRIX3& m);
    MATRIXN(const MATRIXN& m);
    MATRIXN(const SHAREDMATRIXN& m);
    MATRIXN(const CONST_SHAREDMATRIXN& m);
    MATRIXN& set_identity();
    MATRIXN& set_identity(unsigned sz);
    bool is_symmetric(REAL tolerance = -1.0) const;
    static MATRIXN identity(unsigned dim);
    REAL norm_inf() const;
    static MATRIXN construct_variable(unsigned rows, unsigned cols, ...);
    virtual ~MATRIXN() { }
    MATRIXN& zero_upper_triangle();
    MATRIXN& zero_lower_triangle();
    MATRIXN& set(const VECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN& set(const SHAREDVECTORN& v, Transposition trans = eNoTranspose);
    MATRIXN& set(const CONST_SHAREDVECTORN& v, Transposition trans = eNoTranspose);
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
    virtual MATRIXN& operator=(const MATRIXN& source);
    MATRIXN& operator=(const MATRIX2& source);
    MATRIXN& operator=(const MATRIX3& source);
    MATRIXN& operator=(const SHAREDMATRIXN& source);
    MATRIXN& operator=(const CONST_SHAREDMATRIXN& source);
    MATRIXN& operator=(const VECTORN& v) { return set(v, eNoTranspose); }
    MATRIXN& operator=(const SHAREDVECTORN& v) { return set(v, eNoTranspose); }
    MATRIXN& operator=(const CONST_SHAREDVECTORN& v) { return set(v, eNoTranspose); }
    static MATRIXN& mult(const MATRIXN& m1, const MATRIXN& m2, MATRIXN& result) { return m1.mult(m2, result); }
    MATRIXN& operator+=(const MATRIXN& m);
    MATRIXN& operator-=(const MATRIXN& m);
    MATRIXN& operator/=(REAL scalar);
    MATRIXN& operator*=(REAL scalar);
    REAL* data() { return _data.get(); }
    const REAL* data() const { return _data.get(); }    
    void free_memory() { resize(0,0); compress(); }
    void compress() { _data.compress(); }
    unsigned leading_dim() const { return _rows; }
    unsigned inc() const { return 1; }

    /// Gets the total number of elements in this matrix
    unsigned size() const { return _rows * _columns; }

    const REAL& operator()(unsigned i, unsigned j) const;
    REAL& operator()(const unsigned i, const unsigned j);

    #define MATRIXX MATRIXN
    #include "MatrixCommon.inl"
    #undef MATRIXX
    #define XMATRIXN MATRIXN
    #include "XMatrixN.inl"
    #undef XMATRIXN

  protected:
    #ifndef REENTRANT
    static FastThreadable<MATRIXN> _n;
    static FastThreadable<VECTORN> _workv;
    #endif

    SharedResizable<REAL> _data;
    unsigned _rows;
    unsigned _columns;
}; // end class

std::ostream& operator<<(std::ostream& out, const MATRIXN& m);
std::istream& operator>>(std::istream& in, MATRIXN& m);

