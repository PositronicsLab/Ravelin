/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef ORIGIN2
#error This class is not to be included by the user directly. Use Origin2d.h or Origin2f.h instead.
#endif

class VECTOR2;

/// A two-dimensional floating point vector used for computational geometry calculations
class ORIGIN2
{
  public:
    ORIGIN2() {}
    ORIGIN2(REAL x, REAL y);
    ORIGIN2(const REAL* array);
    explicit ORIGIN2(const VECTOR2& v) { operator=(v); }
    REAL norm_inf() const { return std::max(std::fabs(_data[0]), std::fabs(_data[1])); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return sqr(_data[0]) + sqr(_data[1]); }
    static REAL norm(const ORIGIN2& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const ORIGIN2& v) { return v.norm_sq(); }
    ORIGIN2& set_zero() { _data[0] = _data[1] = 0.0; return *this; }
    ORIGIN2& set_zero(unsigned m) { assert(m==2); return set_zero(); }
    ORIGIN2& set_zero(unsigned m, unsigned n) { assert(m==2 && n==1); return set_zero(); }
    static ORIGIN2 zero() { return ORIGIN2(0.0, 0.0); }
    ORIGIN2& operator=(const VECTOR2& v);
    ORIGIN2& operator=(const ORIGIN2& v) { _data[0] = v[0]; _data[1] = v[1]; return *this; }
    ORIGIN2 operator+(const ORIGIN2& v) const { return ORIGIN2(_data[0] + v[0], _data[1] + v[1]); }
    VECTOR2 operator+(const VECTOR2& p) const;
    VECTOR2 operator-(const VECTOR2& p) const;
    ORIGIN2 operator-(const ORIGIN2& v) const { return ORIGIN2(_data[0] - v[0], _data[1] - v[1]); }
    ORIGIN2& operator+=(const ORIGIN2& v) { _data[0] += v[0]; _data[1] += v[1];  return *this; }
    ORIGIN2& operator-=(const ORIGIN2& v) { _data[0] -= v[0]; _data[1] -= v[1];  return *this; }
    ORIGIN2& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; return *this; }
    ORIGIN2& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; return *this; }
    ORIGIN2 operator*(REAL scalar) const { ORIGIN2 v = *this; v *= scalar; return v; }
    ORIGIN2 operator/(REAL scalar) const { ORIGIN2 v = *this; v /= scalar; return v; }
    ORIGIN2 operator-() const { return ORIGIN2(-_data[0], -_data[1]); }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    REAL& operator[](const unsigned i);
    const REAL& operator[](const unsigned i) const;
    REAL* data(unsigned i);
    const REAL* data(unsigned i) const;
    const REAL& x() const { return _data[0]; }
    const REAL& y() const { return _data[1]; }
    REAL& x() { return _data[0]; }
    REAL& y() { return _data[1]; }
    unsigned size() const { return 3; }
    unsigned rows() const { return 3; }
    unsigned columns() const { return 1; }
    unsigned leading_dim() const { return 3; }
    unsigned inc() const { return 1; }
    ORIGIN2& resize(unsigned m, bool preserve = false);
    ORIGIN2& resize(unsigned m, unsigned n, bool preserve = false);

  private:
    REAL _data[2];
    static REAL sqr(REAL x) { return x*x; }
}; // end class

inline ORIGIN2 operator*(REAL scalar, const ORIGIN2& v) { return v * scalar; }

/// Writes a ORIGIN2 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const ORIGIN2& v)
{
  out << '[' << v[0] << ',' << ' ' << v[1] << ']' << ' ';
  return out;
};

