/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef VECTOR3
#error This class is not to be included by the user directly. Use VECTOR3.h or Vector3f.h instead.
#endif

class MATRIX3;

/// A three-dimensional floating point vector
class VECTOR3
{
  public:
    VECTOR3() {}
    VECTOR3(REAL x, REAL y, REAL z);
    VECTOR3(const REAL* array);
    VECTOR3(const VECTOR3& source) { operator=(source); }
    REAL dot(const VECTOR3& v) const { return v[0]*_data[0] + v[1]*_data[1] + v[2]*_data[2]; }
    static REAL dot(const VECTOR3& v1, const VECTOR3& v2) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
    void normalize() { assert(norm() > std::numeric_limits<REAL>::epsilon()); operator/=(norm()); }
    void normalize_or_zero() { REAL nrm = norm(); if (nrm > std::numeric_limits<REAL>::epsilon()) operator/=(nrm); else _data[0] = _data[1] = _data[2] = (REAL) 0.0; }
    static VECTOR3 normalize(const VECTOR3& v) { VECTOR3 w = v; w.normalize(); return w; }
    unsigned size() const { return 3; }
    bool is_finite() const;
    REAL norm_inf() const { return std::max(std::fabs(_data[0]), std::max(std::fabs(_data[1]), std::fabs(_data[2]))); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return dot(*this, *this); }
    static REAL norm(const VECTOR3& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const VECTOR3& v) { return v.dot(v); }
    void set_zero() { _data[0] = _data[1] = _data[2] = 0.0; }
    void set_one() { _data[0] = _data[1] = _data[2] = 1.0; }
    static VECTOR3 zero() { return VECTOR3(0.0, 0.0, 0.0); }
    static VECTOR3 one() { return VECTOR3(1.0, 1.0, 1.0); }
    bool operator<(const VECTOR3& v) const;
    VECTOR3& operator=(const VECTOR3& v) { _data[0] = v[0]; _data[1] = v[1]; _data[2] = v[2]; return *this; }
    VECTOR3 operator+(const VECTOR3& v) const { return VECTOR3(_data[0] + v[0], _data[1] + v[1], _data[2] + v[2]); }
    VECTOR3 operator-(const VECTOR3& v) const { return VECTOR3(_data[0] - v[0], _data[1] - v[1], _data[2] - v[2]); }
    VECTOR3& operator+=(const VECTOR3& v) { _data[0] += v[0]; _data[1] += v[1]; _data[2] += v[2]; return *this; }
    VECTOR3& operator-=(const VECTOR3& v) { _data[0] -= v[0]; _data[1] -= v[1]; _data[2] -= v[2]; return *this; }
    VECTOR3& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; _data[2] *= scalar; return *this; }
    VECTOR3 operator*(REAL scalar) const { VECTOR3 v = *this; v *= scalar; return v; }
    VECTOR3 operator/(REAL scalar) const { VECTOR3 v = *this; v /= scalar; return v; }
    VECTOR3& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; _data[2] /= scalar; return *this; }
    VECTOR3 operator-() const { return VECTOR3(-_data[0], -_data[1], -_data[2]); }
    static VECTOR3 cross(const VECTOR3& v1, const VECTOR3& v2);
    static VECTOR3 determine_orthogonal_vec(const VECTOR3& v);
    static void determine_orthonormal_basis(const VECTOR3& v1, VECTOR3& v2, VECTOR3& v3);
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    unsigned rows() const { return 3; }
    unsigned columns() const { return 1; }
    VECTOR3& resize(unsigned m, unsigned n, bool preserve = false);
    unsigned inc() const { return 1; }
    unsigned leading_dim() const { return 3; }
    REAL& x() { return _data[0]; } 
    REAL x() const { return _data[0]; } 
    REAL& y() { return _data[1]; } 
    REAL y() const { return _data[1]; } 
    REAL& z() { return _data[2]; } 
    REAL z() const { return _data[2]; } 

    REAL& operator[](const unsigned i);
    REAL operator[](const unsigned i) const;
    const REAL* data(unsigned i) const;
    REAL* data(unsigned i);
    ITERATOR begin();
    CONST_ITERATOR begin() const;
    ITERATOR end();
    CONST_ITERATOR end() const;

    VECTOR3& resize(unsigned N, bool keep = true) 
    { 
      #ifndef NEXCEPT
      if (N != 3) 
        throw std::runtime_error("Can't resize a VECTOR3 to size other than 3!"); 
      #endif
      return *this; 
    }

  private:
    REAL _data[3];
}; // end class

inline VECTOR3 operator*(REAL scalar, const VECTOR3& v) { return v * scalar; }

/// Writes a VECTOR3 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const VECTOR3& v)
{
  out << "[" << v[0] << ", " << v[1] << ", " << v[2] << "] ";
  return out;
};

