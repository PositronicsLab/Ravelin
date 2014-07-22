/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef POSE2
#error This class is not to be included by the user directly. Use Posed.h or Posef.h instead. 
#endif

/// A 2D rigid body pose 
class POSE2 : public boost::enable_shared_from_this<POSE2>
{
  public:
    POSE2(boost::shared_ptr<const POSE2> relative_pose = boost::shared_ptr<const POSE2>());
    POSE2(const POSE2& source) { operator=(source); }
    POSE2(const ROT2& r, const ORIGIN2& x, boost::shared_ptr<const POSE2> relative_pose = boost::shared_ptr<const POSE2>());
    POSE2(const ORIGIN2& x, boost::shared_ptr<const POSE2> relative_pose = boost::shared_ptr<const POSE2>());
    POSE2(const ROT2& r, boost::shared_ptr<const POSE2> relative_pose = boost::shared_ptr<const POSE2>());
    static POSE2 identity() { POSE2 T; T.set_identity(); return T; }
    static REAL wrap(REAL theta);
    static bool rel_equal(const POSE2& p1, const POSE2& p2, REAL tol=EPS);
    VECTOR2 transform_point(const VECTOR2& v) const;
    VECTOR2 transform_vector(const VECTOR2& v) const;
    VECTOR2 inverse_transform_point(const VECTOR2& v) const;
    VECTOR2 inverse_transform_vector(const VECTOR2& v) const;
    static VECTOR2 transform_point(boost::shared_ptr<const POSE2> target, const VECTOR2& v);
    static VECTOR2 transform_vector(boost::shared_ptr<const POSE2> target, const VECTOR2& v);
    POSE2& set_identity();
    POSE2& invert();
    POSE2 inverse() const { return invert(*this); }
    static POSE2 invert(const POSE2& m);
    POSE2& set(const ROT2& r, const ORIGIN2& v);
    POSE2& set(const ROT2& r);
    POSE2& operator=(const POSE2& source);
    POSE2 operator*(const POSE2& m) const;
    static VECTOR2 interpolate_transform_vector(const POSE2& P1, const POSE2& P2, REAL t, const ORIGIN2& o);
    static VECTOR2 interpolate_transform_point(const POSE2& P1, const POSE2& P2, REAL t, const ORIGIN2& o);

    /// the orientation of the pose frame
    ROT2 r;

    /// the position of the pose frame
    ORIGIN2 x;

    /// the pose that *this* pose is relative to (if any)
    boost::shared_ptr<const POSE2> rpose; 

  private:
    TRANSFORM2 calc_transform(boost::shared_ptr<const POSE2> p) const;
    static TRANSFORM2 calc_transform(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> target);
    static bool is_common(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> p, unsigned& i);
}; // end class

std::ostream& operator<<(std::ostream& out, const POSE2& m);

