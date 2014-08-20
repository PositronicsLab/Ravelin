/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef VECTOR3
#error This class is not to be included by the user directly. Use Vector3d.h or Vector3f.h instead.
#endif

class POSE3;
class MATRIX3;

/// A three-dimensional floating point vector
class VECTOR3
{
  public:
    VECTOR3(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { this->pose = pose; }
    VECTOR3(boost::shared_ptr<POSE3> pose) { this->pose = boost::const_pointer_cast<const POSE3>(pose); }
    VECTOR3(REAL x, REAL y, REAL z, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    VECTOR3(REAL x, REAL y, REAL z, boost::shared_ptr<POSE3> pose);
    VECTOR3(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    VECTOR3(const REAL* array, boost::shared_ptr<POSE3> pose);
    VECTOR3(const VECTOR3& source) { operator=(source); }
    VECTOR3(const ORIGIN3& source, boost::shared_ptr<const POSE3> pose) { this->pose = pose; operator=(source); }
    VECTOR3(const ORIGIN3& source, boost::shared_ptr<POSE3> pose) { this->pose = boost::const_pointer_cast<const POSE3>(pose); operator=(source); }
    REAL dot(const VECTOR3& v) const { return dot(*this, v); }
    static REAL dot(const VECTOR3& v1, const VECTOR3& v2);
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
    VECTOR3& set_zero() { _data[0] = _data[1] = _data[2] = 0.0; return *this; }
    VECTOR3& set_one() { _data[0] = _data[1] = _data[2] = 1.0; return *this; }
    VECTOR3& set_zero(boost::shared_ptr<const POSE3> pose) { _data[0] = _data[1] = _data[2] = 0.0; this->pose = pose; return *this; }
    VECTOR3& set_one(boost::shared_ptr<const POSE3> pose) { _data[0] = _data[1] = _data[2] = 1.0; this->pose = pose; return *this; }
    static VECTOR3 zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { return VECTOR3(0.0, 0.0, 0.0, pose); }
    static VECTOR3 one(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { return VECTOR3(1.0, 1.0, 1.0, pose); }
    bool operator<(const VECTOR3& v) const;
    VECTOR3& operator=(const ORIGIN3& o) { x() = o.x(); y() = o.y(); z() = o.z(); return *this; }
    VECTOR3& operator=(const VECTOR3& v) { pose = v.pose; x() = v.x(); y() = v.y(); z() = v.z(); return *this; }
    VECTOR3 operator+(const VECTOR3& v) const { VECTOR3 result = *this; result += v; return result; }
    VECTOR3 operator-(const VECTOR3& v) const { VECTOR3 result = *this; result -= v; return result; }
    VECTOR3 operator+(const ORIGIN3& o) const { VECTOR3 result = *this; result += o; return result; }
    VECTOR3 operator-(const ORIGIN3& o) const { VECTOR3 result = *this; result -= o; return result; }
    VECTOR3& operator+=(const VECTOR3& v);
    VECTOR3& operator-=(const VECTOR3& v);
    VECTOR3& operator+=(const ORIGIN3& o);
    VECTOR3& operator-=(const ORIGIN3& o);
    VECTOR3& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; _data[2] *= scalar; return *this; }
    VECTOR3 operator*(REAL scalar) const { VECTOR3 v = *this; v *= scalar; return v; }
    VECTOR3 operator/(REAL scalar) const { VECTOR3 v = *this; v /= scalar; return v; }
    VECTOR3& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; _data[2] /= scalar; return *this; }
    VECTOR3 operator-() const { return VECTOR3(-_data[0], -_data[1], -_data[2], pose); }
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
    const REAL& x() const { return _data[0]; } 
    REAL& y() { return _data[1]; } 
    const REAL& y() const { return _data[1]; } 
    REAL& z() { return _data[2]; } 
    const REAL& z() const { return _data[2]; } 

    REAL& operator[](const unsigned i);
    const REAL& operator[](const unsigned i) const;
    const REAL* data(unsigned i) const;
    REAL* data(unsigned i);
    COLUMN_ITERATOR begin() { return column_iterator_begin(); }
    CONST_COLUMN_ITERATOR begin() const { return column_iterator_begin(); }
    COLUMN_ITERATOR end() { return column_iterator_end(); }
    CONST_COLUMN_ITERATOR end() const { return column_iterator_end(); }
    COLUMN_ITERATOR column_iterator_begin();
    CONST_COLUMN_ITERATOR column_iterator_begin() const;
    COLUMN_ITERATOR column_iterator_end();
    CONST_COLUMN_ITERATOR column_iterator_end() const;
    ROW_ITERATOR row_iterator_begin();
    CONST_ROW_ITERATOR row_iterator_begin() const;
    ROW_ITERATOR row_iterator_end();
    CONST_ROW_ITERATOR row_iterator_end() const;

    VECTOR3& resize(unsigned N, bool keep = true) 
    { 
      #ifndef NEXCEPT
      if (N != 3) 
        throw std::runtime_error("Can't resize a VECTOR3 to size other than 3!"); 
      #endif
      return *this; 
    }

    /// The frame that this vector is defined in
    boost::shared_ptr<const POSE3> pose;

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

