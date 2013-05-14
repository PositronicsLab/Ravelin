/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef ROT2 
#error This class is not to be included by the user directly. Use Rotd.h or Rotf.h instead.
#endif

/// Represents an orientation in 2D
class ROT2
{
  public:
    ROT2() { theta = (REAL) 0.0; }
    ROT2(REAL theta) { this->theta = theta; }
    ROT2(const ROT2& q) { theta = q.theta; }
    static ROT2 zero() { return ROT2((REAL) 0.0); } 
    static ROT2 identity() { return ROT2((REAL) 0.0); }
    ROT2& set_identity() { theta = (REAL) 0.0; return *this; }
    ROT2& set_zero() { theta = (REAL) 0.0; return *this; }
    ROT2& inverse() { theta = -theta; return *this; }
    static ROT2 invert(const ROT2& q) { return ROT2(-q.theta); }
    ROT2& operator=(const ROT2& q) { theta = q.theta; return *this; }
    ROT2 operator*(const ROT2& q) const { return ROT2(theta + q.theta); }
    ROT2 operator/(const ROT2& q) const { return ROT2(theta - q.theta); }
    ROT2 operator*(REAL scalar) const { return ROT2(theta * scalar); }
    ROT2 operator/(REAL scalar) const { return operator*(1.0/scalar); }
    ROT2& operator*=(const ROT2& q) { theta += q.theta; return *this; }
    ROT2& operator*=(REAL scalar) { theta *= scalar; return *this; }
    ROT2& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    
    ORIGIN2 operator*(const ORIGIN2& o) const
    {
      REAL cth = std::cos(theta);
      REAL sth = std::sin(theta);
      return ORIGIN2(cth*o.x() + sth*o.y(), -sth*o.x() + cth*o.y());
    }

    VECTOR2 operator*(const VECTOR2& o) const
    {
      REAL cth = std::cos(theta);
      REAL sth = std::sin(theta);
      return VECTOR2(cth*o.x() + sth*o.y(), -sth*o.x() + cth*o.y());
    }

    POINT2 operator*(const POINT2& o) const
    {
      REAL cth = std::cos(theta);
      REAL sth = std::sin(theta);
      return POINT2(cth*o.x() + sth*o.y(), -sth*o.x() + cth*o.y());
    }

    /// angle of rotation 
    REAL theta;
};

inline std::ostream& operator<<(std::ostream& out, const ROT2& q)
{
  out << "2D rotation: " << q.theta << std::endl;
  return out;
}

