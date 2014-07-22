/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef MATRIX3
#error This class is not to be included by the user directly. Use Matrix3d.h or Matrix3f.h instead.
#endif

class QUAT;
class AANGLE;
class MATRIXN;
class SHAREDMATRIXN;
class SHAREDVECTORN;
class CONST_SHAREDMATRIXN;
class CONST_SHAREDVECTORN;

/// A 3x3 matrix that may be used for orientation, inertia tensors, etc.
class MATRIX3
{
  public:
    MATRIX3() { }
    MATRIX3(const MATRIX3& source) { operator=(source); }
    MATRIX3(const REAL* array);
    MATRIX3(const QUAT& q) { operator=(q); }
    MATRIX3(REAL, REAL, REAL, REAL, REAL, REAL, REAL, REAL, REAL);
    MATRIX3& operator=(const QUAT& q);
    MATRIX3& operator=(const AANGLE& a);
    MATRIX3(const MATRIXN& m) { operator=(m); }
    MATRIX3(const SHAREDMATRIXN& m) { operator=(m); }
    MATRIX3(const CONST_SHAREDMATRIXN& m) { operator=(m); }
    REAL norm_inf() const;
    unsigned rows() const { return 3; }
    unsigned columns() const { return 3; }
    bool is_symmetric(REAL tolerance) const;
    bool orthonormalize();
    bool is_orthonormal() const;    
    REAL det() const;
    MATRIX3& invert();
    const REAL& xx() const { return _data[0]; }
    REAL& xx() { return _data[0]; }
    const REAL& xy() const { return _data[3]; }
    REAL& xy() { return _data[3]; }
    const REAL& xz() const { return _data[6]; }
    REAL& xz() { return _data[6]; }
    const REAL& yx() const { return _data[1]; }
    REAL& yx() { return _data[1]; }
    const REAL& yy() const { return _data[4]; }
    REAL& yy() { return _data[4]; }
    const REAL& yz() const { return _data[7]; }
    REAL& yz() { return _data[7]; }
    const REAL& zx() const { return _data[2]; }
    REAL& zx() { return _data[2]; }
    const REAL& zy() const { return _data[5]; }
    REAL& zy() { return _data[5]; }
    const REAL& zz() const { return _data[8]; }
    REAL& zz() { return _data[8]; }
    static VECTOR3 calc_differential(const MATRIX3& R1, const MATRIX3& R2);
    static MATRIX3 invert(const MATRIX3& m);
    MATRIX3 inverse() const { MATRIX3 m = *this; return m.invert(); }
    void set_rot_X(REAL angle);
    static MATRIX3 rot_X(REAL angle);
    void set_rot_Y(REAL angle);
    static MATRIX3 rot_Y(REAL angle);
    void set_rot_Z(REAL angle);
    static MATRIX3 rot_Z(REAL angle);
    static MATRIX3 skew_symmetric(REAL a, REAL b, REAL c);
    static MATRIX3 skew_symmetric(const ORIGIN3& o);
    static MATRIX3 skew_symmetric(const VECTOR3& v);
    static VECTOR3 inverse_skew_symmetric(const MATRIX3& R);
    static MATRIX3 transpose(const MATRIX3& m);
    void transpose();
    static bool valid_rotation(const MATRIX3& R);
    static bool valid_rotation_scale(const MATRIX3& R);
    static MATRIX3 identity() { MATRIX3 m; m.set_identity(); return m; }
    MATRIX3& set_zero(unsigned m, unsigned n) { resize(m,n); set_zero(); return *this; }
    static MATRIX3 zero() { MATRIX3 m; m.set_zero(); return m; }
    MATRIX3& set_identity();
    MATRIX3& set_zero() { std::fill_n(_data, 9, (REAL) 0.0); return *this; }
    ORIGIN3 transpose_mult(const ORIGIN3& v) const;
    ORIGIN3 mult(const ORIGIN3& v) const;
    MATRIX3 mult(const MATRIX3& m) const;
    MATRIX3 transpose_mult(const MATRIX3& m) const;
    MATRIX3 mult_transpose(const MATRIX3& m) const;
    MATRIX3 transpose_mult_transpose(const MATRIX3& m) const;
    MATRIX3& operator=(const MATRIX3& source);
    MATRIX3& operator=(const MATRIXN& source);
    MATRIX3& operator=(const SHAREDMATRIXN& source);
    MATRIX3& operator=(const CONST_SHAREDMATRIXN& source);
    MATRIX3& operator+=(const MATRIX3& m);
    MATRIX3& operator-=(const MATRIX3& m);
    MATRIX3& operator*=(REAL scalar);
    MATRIX3& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    MATRIX3 operator+(const MATRIX3& m) const { MATRIX3 n = *this; n += m; return n; }
    MATRIX3 operator-(const MATRIX3& m) const { MATRIX3 n = *this; n -= m; return n; }
    MATRIX3 operator*(REAL scalar) const { MATRIX3 m = *this; m *= scalar; return m; }
    MATRIX3 operator/(REAL scalar) const { return operator*(1.0/scalar); }
    MATRIX3 operator-() const; 
    MATRIX3 operator*(const MATRIX3& m) const { MATRIX3 result; mult(m, result); return result; }
    ORIGIN3 operator*(const ORIGIN3& v) const { ORIGIN3 result; mult(v, result); return result; }
    unsigned leading_dim() const { return 3; }
    unsigned inc() const { return 1; }
    ORIGIN3 get_column(unsigned i) const;
    ORIGIN3 get_row(unsigned i) const;
    REAL& operator()(unsigned i, unsigned j);
    const REAL& operator()(unsigned i, unsigned j) const;
    MATRIX3& resize(unsigned rows, unsigned columns, bool preserve = false);
    const REAL* data(unsigned i) const;
    REAL* data(unsigned i);

    /// Gets the total number of elements in this matrix
    unsigned size() const { return 9; }

    /// Gets constant pointer to the beginning of the matrix array
    const REAL* data() const { return _data; }

    /// Gets pointer to the beginning of the matrix array
    REAL* data() { return _data; }

    #define MATRIXX MATRIX3
    #include "MatrixCommon.inl"
    #undef MATRIXX    

  private:
    bool orthonormalize(VECTOR3& a, VECTOR3& b, VECTOR3& c);
    REAL _data[9];
}; // end class

std::ostream& operator<<(std::ostream& out, const MATRIX3& m);

