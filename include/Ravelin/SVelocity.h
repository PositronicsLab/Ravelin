/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SVELOCITY 
#error This class is not to be included by the user directly. Use SVelocityf.h or SVelocityd.h instead.
#endif

class SFORCE;
class SMOMENTUM;

class SVELOCITY : public SVECTOR6 
{
  public:
    /// Constructs a velocity with zero linear and zero angular components
    SVELOCITY(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {}

    /// Constructs a velocity from a SVector6
    SVELOCITY(const SVECTOR6& v) : SVECTOR6(v.get_upper(), v.get_lower(), v.pose) {} 

    /// Constructs a velocity from six values (first three linear, next three angular) and a pose
    SVELOCITY(REAL lx, REAL ly, REAL lz, REAL ax, REAL ay, REAL az, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a velocity from six values (first three linear, next three angular) and a pose
    SVELOCITY(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[3], array[4], array[5], array[0], array[1], array[2], pose) {}

    /// Constructs a velocity from linear and angular components and a pose
    SVELOCITY(const VECTOR3& linear, const VECTOR3& angular, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(angular, linear, pose) {}

    /// Returns a zero velocity
    static SVELOCITY zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SVELOCITY t(pose); t.set_zero(); return t; }

    template <class V>
    V& transpose_to_vector(V& v) const
    {
      const unsigned SPATIAL_DIM = 6;
      v.resize(SPATIAL_DIM);
      REAL* vdata = v.data();
      const REAL* d = data();
      vdata[0] = d[3];  vdata[1] = d[4];  vdata[2] = d[5];
      vdata[3] = d[0];  vdata[4] = d[1];  vdata[5] = d[2];
      return v;
    }

    template <class V>
    static SVELOCITY from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SVELOCITY t(pose);
      REAL* tdata = t.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), tdata, 1);
      return t;
    }

    SACCEL cross(const SMOMENTUM& m) const { return SACCEL(spatial_cross(*this, m)); }
    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    SVELOCITY& operator=(const SVELOCITY& source) { SVECTOR6::operator=(source); return *this; }
    SVELOCITY& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; }
    SVELOCITY operator-() const { SVELOCITY w = *this; w.negate(); return w; }
    SVELOCITY operator+(const SVELOCITY& t) const { SVELOCITY result = *this; result += t; return result; }
    SVELOCITY operator-(const SVELOCITY& t) const { SVELOCITY result = *this; result -= t; return result; }
/*
    SVELOCITY operator*(REAL scalar) const { SVELOCITY v = *this; return v*= scalar; }
    SVELOCITY operator/(REAL scalar) const { SVELOCITY v = *this; return v/= scalar; }
    SVELOCITY& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SVELOCITY& operator*=(REAL scalar);
    SVELOCITY& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SVELOCITY& resize(unsigned rows) { assert (rows == 6); return *this; } 
    SVELOCITY& operator-=(const SVELOCITY& v);
    SVELOCITY& operator+=(const SVELOCITY& v);
    SVELOCITY operator+(const SVELOCITY& v) const;
    SVELOCITY operator-(const SVELOCITY& v) const;
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SVELOCITY& t)
{
  out << "velocity (linear= " << t.get_linear() << ", angular= " << t.get_angular() << ") frame: " << t.pose;
  return out;
}

