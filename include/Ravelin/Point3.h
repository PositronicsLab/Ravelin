/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef POINT3
#error This class is not to be included by the user directly. Use Point3d.h or Vector3f.h instead.
#endif

class POSE3;
class ORIGIN3;
class VECTOR3;

/// A three-dimensional floating point vector
class POINT3
{
  public:
    POINT3(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { this->pose = pose; }
    POINT3(REAL x, REAL y, REAL z, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    POINT3(const ORIGIN3& o, boost::shared_ptr<const POSE3> pose) { this->pose = pose; operator=(o); }
    POINT3(const VECTOR3& p) { operator=(p); }
    POINT3(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    POINT3(const POINT3& source) { operator=(source); }
    REAL dot(const VECTOR3& v) const { return dot(*this, v); }
    REAL dot(const POINT3& p) const { return dot(*this, p); }
    static REAL dot(const POINT3& p1, const POINT3& p2);
    static REAL dot(const POINT3& p, const VECTOR3& v);
    static REAL dot(const VECTOR3& v, const POINT3& p);
    REAL norm_inf() const { return std::max(std::fabs(_data[0]), std::max(std::fabs(_data[1]), std::fabs(_data[2]))); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return sqr(_data[0]) + sqr(_data[1]) + sqr(_data[2]); }
    static REAL norm(const POINT3& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const POINT3& v) { return v.norm_sq(); }
    POINT3& set_zero() { _data[0] = _data[1] = _data[2] = 0.0; return *this; }
    POINT3& set_zero(boost::shared_ptr<const POSE3> pose) { _data[0] = _data[1] = _data[2] = 0.0; this->pose = pose; return *this; }
    static POINT3 zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { return POINT3(0.0, 0.0, 0.0, pose); }
    bool operator<(const POINT3& v) const;
    POINT3& operator=(const ORIGIN3& o);
    POINT3& operator=(const VECTOR3& v);
    POINT3& operator=(const POINT3& p);
    VECTOR3 operator+(const POINT3& v) const;
    VECTOR3 operator-(const POINT3& v) const;
    VECTOR3 operator+(const ORIGIN3& o) const;
    VECTOR3 operator-(const ORIGIN3& o) const;
    POINT3& operator+=(const ORIGIN3& o);
    POINT3& operator-=(const ORIGIN3& o);
    POINT3& operator+=(const VECTOR3& v);
    POINT3& operator-=(const VECTOR3& v);
    POINT3& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; _data[2] *= scalar; return *this; }
    POINT3 operator*(REAL scalar) const { POINT3 v = *this; v *= scalar; return v; }
    POINT3 operator/(REAL scalar) const { POINT3 v = *this; v /= scalar; return v; }
    POINT3& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; _data[2] /= scalar; return *this; }
    POINT3 operator-() const { return POINT3(-_data[0], -_data[1], -_data[2]); }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    REAL& x() { return _data[0]; } 
    const REAL& x() const { return _data[0]; } 
    REAL& y() { return _data[1]; } 
    const REAL& y() const { return _data[1]; } 
    REAL& z() { return _data[2]; } 
    const REAL& z() const { return _data[2]; } 
    ITERATOR begin();
    CONST_ITERATOR begin() const;
    ITERATOR end();
    CONST_ITERATOR end() const;

    REAL& operator[](const unsigned i);
    const REAL& operator[](const unsigned i) const;
    const REAL* data(unsigned i) const;
    REAL* data(unsigned i);

    /// The 3D pose this point is defined with respect to
    boost::shared_ptr<const POSE3> pose;

  private:
    REAL _data[3];
    static REAL sqr(REAL x) { return x*x; }
}; // end class

inline POINT3 operator*(REAL scalar, const POINT3& v) { return v * scalar; }

/// Writes a POINT3 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const POINT3& v)
{
  out << "[" << v[0] << ", " << v[1] << ", " << v[2] << "] ";
  return out;
};

