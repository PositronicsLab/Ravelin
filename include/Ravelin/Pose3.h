/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef POSE3
#error This class is not to be included by the user directly. Use Pose3d.h or Pose3f.h instead. 
#endif

class AANGLE;
class MATRIX3;

/// A rigid body pose 
class POSE3 : public boost::enable_shared_from_this<POSE3>
{
  public:
    POSE3();
    POSE3(const POSE3& source) { operator=(source); }
    POSE3(const AANGLE& a);
    POSE3(const MATRIX3& m);
    POSE3(const QUAT& q);
    POSE3(const AANGLE& a, const ORIGIN3& v);
    POSE3(const MATRIX3& m, const ORIGIN3& v);
    POSE3(const QUAT& q, const ORIGIN3& v);
    POSE3(const ORIGIN3& v);
    static POSE3 identity() { POSE3 T; T.set_identity(); return T; }
    static POSE3 interpolate(const POSE3& m1, const POSE3& m2, REAL t);
    POINT3 transform(const POINT3& v) const;
    VECTOR3 transform(const VECTOR3& v) const;
    POINT3 inverse_transform(const POINT3& v) const;
    VECTOR3 inverse_transform(const VECTOR3& v) const;
    POINT3 transform(boost::shared_ptr<const POSE3> p, const POINT3& v) const;
    VECTOR3 transform(boost::shared_ptr<const POSE3> p, const VECTOR3& v) const;
    static POINT3 transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const POINT3& v);
    static VECTOR3 transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const VECTOR3& v);
    WRENCH transform(const WRENCH& w) const;
    WRENCH inverse_transform(const WRENCH& w) const;
    WRENCH transform(boost::shared_ptr<const POSE3> p, const WRENCH& w) const;
    static WRENCH transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const WRENCH& w);
    TWIST transform(const TWIST& t) const;
    TWIST inverse_transform(const TWIST& t) const;
    TWIST transform(boost::shared_ptr<const POSE3> p, const TWIST& t) const;
    static TWIST transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const TWIST& t);
    void set_relative_pose(boost::shared_ptr<const POSE3> p); 
    SPATIAL_RB_INERTIA transform(const SPATIAL_RB_INERTIA& j) const;
    SPATIAL_RB_INERTIA inverse_transform(const SPATIAL_RB_INERTIA& j) const;
    SPATIAL_RB_INERTIA transform(boost::shared_ptr<const POSE3> p, const SPATIAL_RB_INERTIA& j) const;
    static SPATIAL_RB_INERTIA transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const SPATIAL_RB_INERTIA& j);
    SPATIAL_AB_INERTIA transform(const SPATIAL_AB_INERTIA& j) const;
    SPATIAL_AB_INERTIA inverse_transform(const SPATIAL_AB_INERTIA& j) const;
    SPATIAL_AB_INERTIA transform(boost::shared_ptr<const POSE3> p, const SPATIAL_AB_INERTIA& j) const;
    static SPATIAL_AB_INERTIA transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const SPATIAL_AB_INERTIA& j);
    POSE3& set_identity();
    POSE3& invert();
    POSE3 inverse() const { return inverse(*this); }
    static POSE3 inverse(const POSE3& m);
    POSE3& set(const AANGLE& a);
    POSE3& set(const MATRIX3& m);
    POSE3& set(const QUAT& q);
    POSE3& set(const AANGLE& a, const ORIGIN3& v);
    POSE3& set(const MATRIX3& m, const ORIGIN3& v);
    POSE3& set(const QUAT& q, const ORIGIN3& v);
    POSE3& operator=(const POSE3& source);
    POSE3 operator*(const POSE3& m) const;

    /// the orientation of the pose frame
    QUAT q;

    /// the origin of the pose frame
    ORIGIN3 x;

    /// the pose that *this* pose is relative to (if any)
    boost::shared_ptr<const POSE3> rpose; 

  private:
    void get_r_E(VECTOR3& r, MATRIX3& E, bool inverse) const;
    std::pair<QUAT, ORIGIN3> calc_transform(boost::shared_ptr<const POSE3> p) const;
    static std::pair<QUAT, ORIGIN3> calc_transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target);
    static bool is_common(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> p, unsigned& i);
}; // end class

std::ostream& operator<<(std::ostream& out, const POSE3& m);

