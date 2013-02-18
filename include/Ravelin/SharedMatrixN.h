/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SHAREDMATRIXN
#error This class is not to be included by the user directly. Use SharedMatrixNd.h or SharedMatrixNf.h instead.
#endif

/// A generic, possibly non-square matrix.  
/**
 * The underlying data is stored in column-major format (e.g., the element at row 1, column 0 is the second element in the flat array).
 */
class SHAREDMATRIXN
{
  friend class MATRIXN;

  public:
    SHAREDMATRIXN();
    SHAREDMATRIXN(const SHAREDMATRIXN& source);
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
    const REAL& operator()(unsigned i, unsigned j) const;
    REAL& operator()(const unsigned i, const unsigned j);
    REAL* data() { return _data.get()+_start; }
    const REAL* data() const { return _data.get()+_start; }    

    #define MATRIXX SHAREDMATRIXN
    #include "MatrixCommon.inl"
    #undef MATRIXX
    #define XMATRIXN SHAREDMATRIXN
    #include "XMatrixN.inl"
    #undef XMATRIXN

  protected:
    boost::shared_array<REAL> _data;
    unsigned _rows;
    unsigned _start;
    unsigned _columns;
    unsigned _ld;       // the leading dimension of the matrix
}; // end class

std::ostream& operator<<(std::ostream& out, const SHAREDMATRIXN& m);
std::istream& operator>>(std::istream& in, SHAREDMATRIXN& m);

