/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef POINT2
#error This class is not to be included by the user directly. Use Point2d.h or Point2f.h instead.
#endif

class POSE2;
class ORIGIN2;
class VECTOR2;

/// A point defined with respect to a two-dimensional frame 
class POINT2
{
  public:
    POINT2(boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>()) { this->pose = pose; }
    POINT2(REAL x, REAL y, boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>());
    POINT2(const REAL* array, boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>());
    POINT2(const ORIGIN2& o, boost::shared_ptr<const POSE2> pose) { this->pose = pose; operator=(o); }
    POINT2(const VECTOR2& v) { operator=(v); } 
    REAL norm_inf() const { return std::max(std::fabs(_data[0]), std::fabs(_data[1])); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return sqr(_data[0]) + sqr(_data[1]); }
    static REAL norm(const POINT2& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const POINT2& v) { return v.norm_sq(); }
    POINT2& set_zero() { _data[0] = _data[1] = 0.0; return *this; }
    POINT2& set_zero(boost::shared_ptr<const POSE2> pose) { _data[0] = _data[1] = 0.0; this->pose = pose; return *this; }
    static POINT2 zero(boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>()) { return POINT2(0.0, 0.0, pose); }
    POINT2& operator=(const ORIGIN2& o);
    POINT2& operator=(const POINT2& v);
    POINT2& operator=(const VECTOR2& v);
    POINT2 operator+(const ORIGIN2& o) const;
    POINT2 operator-(const ORIGIN2& o) const;
    VECTOR2 operator+(const POINT2& v) const;
    VECTOR2 operator-(const POINT2& v) const;
    POINT2& operator+=(const ORIGIN2& v);
    POINT2& operator-=(const ORIGIN2& v);
    POINT2& operator+=(const VECTOR2& v);
    POINT2& operator-=(const VECTOR2& v);
    POINT2& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; return *this; }
    POINT2& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; return *this; }
    POINT2 operator*(REAL scalar) const { POINT2 v = *this; v *= scalar; return v; }
    POINT2 operator/(REAL scalar) const { POINT2 v = *this; v /= scalar; return v; }
    POINT2 operator-() const { return POINT2(-_data[0], -_data[1], pose); }
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

    /// The 2D pose that this point is defined with respect to 
    boost::shared_ptr<const POSE2> pose; 

  private:
    static REAL sqr(REAL x) { return x*x; }
    REAL _data[2];
}; // end class

inline POINT2 operator*(REAL scalar, const POINT2& v) { return v * scalar; }

/// Writes a POINT2 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const POINT2& v)
{
  out << '[' << v[0] << ',' << ' ' << v[1] << ']' << ' ';
  return out;
};

