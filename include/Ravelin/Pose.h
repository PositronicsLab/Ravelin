/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef POSE
#error This class is not to be included by the user directly. Use Posed.h or Posef.h instead. 
#endif

/// A rigid body pose 
class POSE
{
  public:
    POSE();
    POSE(const POSE& source) { operator=(source); }
    POSE(const AANGLE& a);
    POSE(const MATRIX3& m);
    POSE(const QUAT& q);
    POSE(const AANGLE& a, const VECTOR3& v);
    POSE(const MATRIX3& m, const VECTOR3& v);
    POSE(const QUAT& q, const VECTOR3& v);
    POSE(const VECTOR3& v);
    static POSE identity() { POSE T; T.set_identity(); return T; }
    static POSE interpolate(const POSE& m1, const POSE& m2, REAL t);
    VECTOR3 mult_point(const VECTOR3& v) const;
    VECTOR3 mult_vector(const VECTOR3& v) const;
    VECTOR3 inverse_mult_point(const VECTOR3& v) const;
    VECTOR3 inverse_mult_vector(const VECTOR3& v) const;
    POSE& set_identity();
    POSE& invert();
    POSE inverse() const { return inverse(*this); }
    static POSE inverse(const POSE& m);
    POSE& set(const AANGLE& a);
    POSE& set(const MATRIX3& m);
    POSE& set(const QUAT& q);
    POSE& set(const AANGLE& a, const VECTOR3& v);
    POSE& set(const MATRIX3& m, const VECTOR3& v);
    POSE& set(const QUAT& q, const VECTOR3& v);
    POSE& operator=(const POSE& source);
    POSE operator*(const POSE& m) const;

    /// the orientation of the pose frame
    QUAT q;

    /// the position of the pose frame
    VECTOR3 x;

    /// the pose that *this* pose is relative to (if any)
    boost::shared_ptr<POSE> rpose; 
}; // end class

std::ostream& operator<<(std::ostream& out, const POSE& m);

