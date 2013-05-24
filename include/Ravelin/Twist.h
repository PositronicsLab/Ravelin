/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef TWIST 
#error This class is not to be included by the user directly. Use Twistd.h or Twistf.h instead.
#endif

class WRENCH;

class TWIST : public SVECTOR6 
{
  public:
    /// Constructs a twist with zero linear and zero angular components
    TWIST(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {}

    /// Constructs a twist from a SVector6
    TWIST(const SVECTOR6& v) : SVECTOR6(v.get_upper(), v.get_lower(), v.pose) {} 

    /// Constructs a twist from six values (first three linear, next three angular) and a pose
    TWIST(REAL lx, REAL ly, REAL lz, REAL ax, REAL ay, REAL az, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a twist from six values (first three linear, next three angular) and a pose
    TWIST(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[3], array[4], array[5], array[0], array[1], array[2], pose) {}

    /// Constructs a twist from linear and angular components and a pose
    TWIST(const VECTOR3& linear, const VECTOR3& angular, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(angular, linear, pose) {}

    /// Returns a zero twist
    static TWIST zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { TWIST t(pose); t.set_zero(); return t; }

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
    static TWIST from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      TWIST t(pose);
      REAL* tdata = t.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), tdata, 1);
      return t;
    }

    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    TWIST& operator=(const TWIST& source) { SVECTOR6::operator=(source); return *this; }
    TWIST& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; }
    TWIST operator-() const { TWIST w = *this; w.negate(); return w; }
    TWIST operator+(const TWIST& t) const { TWIST result = *this; result += t; return result; }
    TWIST operator-(const TWIST& t) const { TWIST result = *this; result -= t; return result; }
/*
    TWIST operator*(REAL scalar) const { TWIST v = *this; return v*= scalar; }
    TWIST operator/(REAL scalar) const { TWIST v = *this; return v/= scalar; }
    TWIST& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    TWIST& operator*=(REAL scalar);
    TWIST& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    TWIST& resize(unsigned rows) { assert (rows == 6); return *this; } 
    TWIST& operator-=(const TWIST& v);
    TWIST& operator+=(const TWIST& v);
    TWIST operator+(const TWIST& v) const;
    TWIST operator-(const TWIST& v) const;
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const TWIST& t)
{
  out << "Twist (linear= " << t.get_linear() << ", angular= " << t.get_angular() << ") frame: " << t.pose;
  return out;
}

