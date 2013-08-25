/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef ORIGIN3
#error This class is not to be included by the user directly. Use Origin3d.h or Origin3d.h instead.
#endif

class VECTOR3;

/// A two-dimensional floating point vector used for computational geometry calculations
class ORIGIN3
{
  public:
    ORIGIN3() {}
    ORIGIN3(REAL x, REAL y, REAL z);
    ORIGIN3(const REAL* array);
    explicit ORIGIN3(const VECTOR3& v) { operator=(v); }
    REAL norm_inf() const { return std::max(std::max(std::fabs(_data[0]), std::fabs(_data[1])), std::fabs(_data[2])); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return sqr(_data[0]) + sqr(_data[1]) + sqr(_data[2]); }
    static REAL norm(const ORIGIN3& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const ORIGIN3& v) { return v.norm_sq(); }
    void set_zero() { _data[0] = _data[1] = _data[2] = (REAL) 0.0; }
    static ORIGIN3 zero() { return ORIGIN3((REAL) 0.0, (REAL) 0.0, (REAL) 0.0); }
    ORIGIN3& operator=(const ORIGIN3& o);
    ORIGIN3& operator=(const VECTOR3& v);
    ORIGIN3 operator+(const ORIGIN3& o) const;
    VECTOR3 operator+(const VECTOR3& v) const;
    ORIGIN3 operator-(const ORIGIN3& v) const;
    VECTOR3 operator-(const VECTOR3& v) const;
    ORIGIN3& operator+=(const ORIGIN3& v);
    ORIGIN3& operator-=(const ORIGIN3& v);
    ORIGIN3& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; _data[2] *= scalar; return *this; }
    ORIGIN3& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; _data[2] /= scalar; return *this; }
    ORIGIN3 operator*(REAL scalar) const { ORIGIN3 v = *this; v *= scalar; return v; }
    ORIGIN3 operator/(REAL scalar) const { ORIGIN3 v = *this; v /= scalar; return v; }
    ORIGIN3 operator-() const { return ORIGIN3(-_data[0], -_data[1], -_data[2]); }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    REAL& operator[](const unsigned i);
    const REAL& operator[](const unsigned i) const;
    REAL* data(unsigned i);
    const REAL* data(unsigned i) const;
    const REAL& x() const { return _data[0]; }
    const REAL& y() const { return _data[1]; }
    const REAL& z() const { return _data[2]; }
    REAL& x() { return _data[0]; }
    REAL& y() { return _data[1]; }
    REAL& z() { return _data[2]; }
    COLUMN_ITERATOR column_iterator_begin();
    CONST_COLUMN_ITERATOR column_iterator_begin() const;
    COLUMN_ITERATOR column_iterator_end();
    CONST_COLUMN_ITERATOR column_iterator_end() const;
    ROW_ITERATOR row_iterator_begin();
    CONST_ROW_ITERATOR row_iterator_begin() const;
    ROW_ITERATOR row_iterator_end();
    CONST_ROW_ITERATOR row_iterator_end() const;
    unsigned rows() const { return 3; }
    unsigned columns() const { return 1; }
    unsigned inc() const { return 1; }
    unsigned leading_dim() const { return 3; }
    ORIGIN3& resize(unsigned m, unsigned n, bool preserve = false);

  private:
    REAL _data[3];
    static REAL sqr(REAL x) { return x*x; }
}; // end class

inline ORIGIN3 operator*(REAL scalar, const ORIGIN3& v) { return v * scalar; }

/// Writes a ORIGIN3 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const ORIGIN3& v)
{
  out << "[" << v[0] << " " << v[1] << " " << v[2] << "] ";

  return out;
};

