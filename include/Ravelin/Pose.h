/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef POSE
#error This class is not to be included by the user directly. Use Posed.h or Posef.h instead. 
#endif

#include <boost/shared_ptr.hpp>

/// A rigid body pose 
class POSE
{
  public:
    POSE();
    POSE(const POSE& source) { operator=(source); }
    POSE(const REAL* array);
    POSE(const AANGLE& a);
    POSE(const MATRIX3& m);
    POSE(const QUAT& q);
    POSE(const AANGLE& a, const VECTOR& v);
    POSE(const MATRIX3& m, const VECTOR& v);
    POSE(const QUAT& q, const VECTOR& v);
    POSE(const VECTOR& v);
    static POSE identity() { POSE T; T.set_identity(); return T; }
    static POSE interpolate(const POSE& m1, const POSE& m2, REAL t);
    VECTOR mult_point(const VECTOR& v) const;
    VECTOR mult_vector(const VECTOR& v) const;
    VECTOR inverse_mult_point(const VECTOR& v) const;
    VECTOR transpose_mult_vector(const VECTOR& v) const;
    void set_identity();
    void invert_transform();
    POSE inverse_transform() const { return inverse_transform(*this); }
    static POSE inverse_transform(const POSE& m);
    void set(const REAL* array);
    void set(const AANGLE& a, const VECTOR&  v);
    void set(const MATRIX3& m, const VECTOR&  v);
    void set(const QUAT& q, const VECTOR&  v);
    POSE& operator=(const POSE& source);
    POSE operator*(const POSE& m) const;

    /// the orientation of the pose frame
    QUAT q;

    /// the position of the pose frame
    VECTOR x;

    /// the pose that *this* pose is relative to
    boost::shared_ptr<POSE> rpose; 
}; // end class

std::ostream& operator<<(std::ostream& out, const POSE& m);

