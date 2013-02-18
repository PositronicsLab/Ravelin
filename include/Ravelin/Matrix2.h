/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef MATRIX2
#error This class is not to be included by the user directly. Use Matrix2d.h or Matrix2f.h instead.
#endif

/// A general 2x2 matrix
class MATRIX2
{
  public:
    MATRIX2() {}
    MATRIX2(const REAL* array);
    MATRIX2(REAL m00, REAL m01, REAL m10, REAL m11);
    MATRIX2(const MATRIX2& source) { operator=(source); }
    unsigned size() const { return 2; }
    unsigned rows() const { return 2; }
    unsigned columns() const { return 2; }
    bool is_symmetric(REAL tolerance) const;
    bool orthonormalize();
    bool is_orthonormal() const;    
    REAL det() const;
    void inverse();
    static MATRIX2 inverse(const MATRIX2& m);
    void set_rot_Z(REAL angle);
    static MATRIX2 rot_Z(REAL angle);
    static MATRIX2 transpose(const MATRIX2& m);
    void transpose();
    static bool valid_rotation(const MATRIX2& R);
    static MATRIX2 identity() { MATRIX2 m; m.set_identity(); return m; }
    static MATRIX2 zero() { MATRIX2 m; m.set_zero(); return m; }
    VECTOR2 get_row(unsigned i) const;
    VECTOR2 get_column(unsigned i) const;
    void set_identity();
    void set_zero();
    VECTOR2 mult(const VECTOR2& v) const;
    VECTOR2 transpose_mult(const VECTOR2& v) const;
    MATRIX2 mult(const MATRIX2& m) const;
    MATRIX2 transpose_mult(const MATRIX2& m) const;
    MATRIX2 transpose_mult_transpose(const MATRIX2& m) const;
    MATRIX2 mult_transpose(const MATRIX2& m) const;
    MATRIX2& operator=(const MATRIX2& source);
    MATRIX2& operator+=(const MATRIX2& m);
    MATRIX2& operator-=(const MATRIX2& m);
    MATRIX2& operator*=(const MATRIX2& m) { return *this = *this * m; }
    MATRIX2& operator*=(REAL scalar);
    MATRIX2& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    VECTOR2 operator*(const VECTOR2& v) const { return mult(v); } 
    MATRIX2 operator+(const MATRIX2& m) const { MATRIX2 n = *this; n += m; return n; }
    MATRIX2 operator-(const MATRIX2& m) const { MATRIX2 n = *this; n -= m; return n; }
    MATRIX2 operator*(const MATRIX2& m) const { return mult(m); }
    MATRIX2 operator*(REAL scalar) const { MATRIX2 m = *this; m *= scalar; return m; }
    MATRIX2 operator/(REAL scalar) const { return operator*(1.0/scalar); }
    MATRIX2 operator-() const; 
    unsigned leading_dim() const { return 2; }
    unsigned inc() const { return 1; }
    REAL xx() const { return _data[0]; }
    REAL yx() const { return _data[1]; }
    REAL xy() const { return _data[2]; }
    REAL yy() const { return _data[3]; }
    REAL& xx() { return _data[0]; }
    REAL& yx() { return _data[1]; }
    REAL& xy() { return _data[2]; }
    REAL& yy() { return _data[3]; }
    REAL& operator()(unsigned i, unsigned j);
    REAL operator()(unsigned i, unsigned j) const;
    MATRIX2& resize(unsigned rows, unsigned columns, bool preserve = false);
    ITERATOR begin();
    CONST_ITERATOR begin() const;
    ITERATOR end();
    CONST_ITERATOR end() const;
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

