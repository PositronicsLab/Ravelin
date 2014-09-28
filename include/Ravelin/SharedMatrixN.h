/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SHAREDMATRIXN
#error This class is not to be included by the user directly. Use SharedMatrixNd.h or SharedMatrixNf.h instead.
#endif

class MATRIXN;
class CONST_SHAREDMATRIXN;
class SHAREDMATRIXN;

/// A generic, possibly non-square matrix using constant shared data 
/**
 * The underlying data is stored in column-major format (e.g., the element at row 1, column 0 is the second element in the flat array).
 */
class CONST_SHAREDMATRIXN
{
  friend class MATRIXN;
  friend class SHAREDMATRIXN;

  public:
    CONST_SHAREDMATRIXN();
    CONST_SHAREDMATRIXN(const SHAREDMATRIXN& source);
    CONST_SHAREDMATRIXN(const CONST_SHAREDMATRIXN& source);
    CONST_SHAREDMATRIXN(unsigned rows, unsigned cols, unsigned leading_dim, unsigned start, SharedResizable<REAL> data);
    const SHAREDMATRIXN get() const;
    void reset_from(const CONST_SHAREDMATRIXN& source);
    void reset_from(const SHAREDMATRIXN& source);
    bool is_symmetric(REAL tolerance = (REAL) -1.0) const;
    REAL norm_inf() const;
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    unsigned leading_dim() const { return _ld; }
    CONST_SHAREDMATRIXN& resize(unsigned rows, unsigned columns, bool preserve = false);
    const REAL& operator()(unsigned i, unsigned j) const;
    const REAL* data() const { return _data.get()+_start; }    

    /// Resets the shared matrix
    void reset() { _data.reset(); _rows = _start = _columns = _ld = 0; }

    #include "ConstMatrixCommon.inl"
    #include "ConstMatrixN.inl"

  protected:
    SharedResizable<REAL> _data;
    unsigned _rows;
    unsigned _start;
    unsigned _columns;
    unsigned _ld;       // the leading dimension of the matrix
}; // end class

/// A generic, possibly non-square matrix using shared data 
/**
 * The underlying data is stored in column-major format (e.g., the element at row 1, column 0 is the second element in the flat array).
 */
class SHAREDMATRIXN
{
  friend class MATRIXN;
  friend class CONST_SHAREDMATRIXN;

  public:
    SHAREDMATRIXN();
    SHAREDMATRIXN(const SHAREDMATRIXN& source);
    SHAREDMATRIXN(unsigned rows, unsigned cols, unsigned leading_dim, unsigned start, SharedResizable<REAL> data);
    void reset_from(const SHAREDMATRIXN& source);
    SHAREDMATRIXN& set_identity();
    bool is_symmetric(REAL tolerance = (REAL) -1.0) const;
    REAL norm_inf() const;
    SHAREDMATRIXN& zero_upper_triangle();
    SHAREDMATRIXN& zero_lower_triangle();
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    unsigned leading_dim() const { return _ld; }
    SHAREDMATRIXN& resize(unsigned rows, unsigned columns, bool preserve = false);
    SHAREDMATRIXN& negate();
    SHAREDMATRIXN& set_zero();
    SHAREDMATRIXN& operator/=(REAL scalar);
    SHAREDMATRIXN& operator*=(REAL scalar);
    SHAREDMATRIXN& operator=(const MATRIX3& m);
    SHAREDMATRIXN& operator=(const MATRIXN& source);
    SHAREDMATRIXN& operator=(const SHAREDMATRIXN& source);
    SHAREDMATRIXN& operator=(const CONST_SHAREDMATRIXN& source);
    REAL& operator()(unsigned i, unsigned j);
    const REAL& operator()(unsigned i, unsigned j) const;
    REAL* data() { return _data.get()+_start; }
    const REAL* data() const { return _data.get()+_start; }    
 
    /// Resets the shared matrix
    void reset() { _data.reset(); _rows = _start = _columns = _ld = 0; }

    #define XMATRIXN SHAREDMATRIXN
    #include "XMatrixN.inl"
    #undef XMATRIXN
    #define MATRIXX SHAREDMATRIXN
    #include "MatrixCommon.inl"
    #undef MATRIXX

  protected:
    SharedResizable<REAL> _data;
    unsigned _rows;
    unsigned _start;
    unsigned _columns;
    unsigned _ld;       // the leading dimension of the matrix
}; // end class


std::ostream& operator<<(std::ostream& out, const SHAREDMATRIXN& m);
std::ostream& operator<<(std::ostream& out, const CONST_SHAREDMATRIXN& m);
std::istream& operator>>(std::istream& in, SHAREDMATRIXN& m);

