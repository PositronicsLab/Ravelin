/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef MATRIX2
#error This class is not to be included by the user directly. Use Matrix2d.h or Matrix2f.h instead.
#endif

class MATRIXN;
class SHAREDMATRIXN;
class CONST_SHAREDMATRIXN;

/// A general 2x2 matrix
class MATRIX2
{
  public:
    MATRIX2() {}
    MATRIX2(const REAL* array);
    MATRIX2(REAL m00, REAL m01, REAL m10, REAL m11);
    MATRIX2(const MATRIX2& source) { operator=(source); }
    MATRIX2(const MATRIXN& m) { operator=(m); }
    MATRIX2(const SHAREDMATRIXN& m) { operator=(m); }
    MATRIX2(const CONST_SHAREDMATRIXN& m) { operator=(m); }
    unsigned size() const { return 4; }
    unsigned rows() const { return 2; }
    unsigned columns() const { return 2; }
    bool is_symmetric(REAL tolerance) const;
    bool orthonormalize();
    bool is_orthonormal() const;    
    REAL det() const;
    MATRIX2& invert();
    static MATRIX2 invert(const MATRIX2& m);
    MATRIX2 inverse() const { MATRIX2 m = *this; return m.invert(); }
    void set_rot_Z(REAL angle);
    static MATRIX2 rot_Z(REAL angle);
    static MATRIX2 transpose(const MATRIX2& m);
    void transpose();
    static bool valid_rotation(const MATRIX2& R);
    static MATRIX2 identity() { MATRIX2 m; m.set_identity(); return m; }
    static MATRIX2 zero() { MATRIX2 m; m.set_zero(); return m; }
    ORIGIN2 get_row(unsigned i) const;
    ORIGIN2 get_column(unsigned i) const;
    void set_identity();
    void set_zero();
    ORIGIN2 mult(const ORIGIN2& v) const;
    ORIGIN2 transpose_mult(const ORIGIN2& v) const;
    MATRIX2 mult(const MATRIX2& m) const;
    MATRIX2 transpose_mult(const MATRIX2& m) const;
    MATRIX2 transpose_mult_transpose(const MATRIX2& m) const;
    MATRIX2 mult_transpose(const MATRIX2& m) const;
    MATRIX2& operator=(const MATRIX2& source);
    MATRIX2& operator=(const MATRIXN& source);
    MATRIX2& operator=(const SHAREDMATRIXN& source);
    MATRIX2& operator=(const CONST_SHAREDMATRIXN& source);
    MATRIX2& operator+=(const MATRIX2& m);
    MATRIX2& operator-=(const MATRIX2& m);
    MATRIX2& operator*=(const MATRIX2& m) { return *this = *this * m; }
    MATRIX2& operator*=(REAL scalar);
    MATRIX2& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    ORIGIN2 operator*(const ORIGIN2& o) const { return mult(o); } 
    MATRIX2 operator+(const MATRIX2& m) const { MATRIX2 n = *this; n += m; return n; }
    MATRIX2 operator-(const MATRIX2& m) const { MATRIX2 n = *this; n -= m; return n; }
    MATRIX2 operator*(const MATRIX2& m) const { return mult(m); }
    MATRIX2 operator*(REAL scalar) const { MATRIX2 m = *this; m *= scalar; return m; }
    MATRIX2 operator/(REAL scalar) const { return operator*(1.0/scalar); }
    MATRIX2 operator-() const; 
    unsigned leading_dim() const { return 2; }
    unsigned inc() const { return 1; }
    const REAL& xx() const { return _data[0]; }
    const REAL& yx() const { return _data[1]; }
    const REAL& xy() const { return _data[2]; }
    const REAL& yy() const { return _data[3]; }
    REAL& xx() { return _data[0]; }
    REAL& yx() { return _data[1]; }
    REAL& xy() { return _data[2]; }
    REAL& yy() { return _data[3]; }
    REAL& operator()(unsigned i, unsigned j);
    const REAL& operator()(unsigned i, unsigned j) const;
    MATRIX2& resize(unsigned rows, unsigned columns, bool preserve = false);
    const REAL* data(unsigned i) const;
    REAL* data(unsigned i);

    /// Gets a constant pointer to the beginning of the matrix array 
    const REAL* data() const { return _data; }

    /// Gets a pointer to the beginning of the matrix array
    REAL* data() { return _data; }

    #define MATRIXX MATRIX2
    #include "MatrixCommon.inl"
    #undef MATRIXX

  private:
    bool orthonormalize(VECTOR2& a, VECTOR2& b);
    REAL _data[4];
}; // end class

std::ostream& operator<<(std::ostream& out, const MATRIX2& m);

