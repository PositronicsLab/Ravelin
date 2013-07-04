/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SACCEL 
#error This class is not to be included by the user directly. Use SAccelf.h or SAcceld.h instead.
#endif

class SFORCE;

class SACCEL : public SVECTOR6 
{
  public:
    /// Constructs a acceleration with zero linear and zero angular components
    SACCEL(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {}

    /// Constructs a acceleration from a SVector6
    SACCEL(const SVECTOR6& v) : SVECTOR6(v.get_upper(), v.get_lower(), v.pose) {} 

    /// Constructs a acceleration from six values (first three linear, next three angular) and a pose
    SACCEL(REAL lx, REAL ly, REAL lz, REAL ax, REAL ay, REAL az, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a acceleration from six values (first three linear, next three angular) and a pose
    SACCEL(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[3], array[4], array[5], array[0], array[1], array[2], pose) {}

    /// Constructs a acceleration from linear and angular components and a pose
    SACCEL(const VECTOR3& linear, const VECTOR3& angular, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(angular, linear, pose) {}

    /// Returns a zero acceleration
    static SACCEL zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SACCEL t(pose); t.set_zero(); return t; }

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
    static SACCEL from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SACCEL t(pose);
      REAL* tdata = t.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), tdata, 1);
      return t;
    }

    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    SACCEL& operator=(const SACCEL& source) { SVECTOR6::operator=(source); return *this; }
    SACCEL& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; }
    SACCEL operator-() const { SACCEL w = *this; w.negate(); return w; }
    SACCEL operator+(const SACCEL& t) const { SACCEL result = *this; result += t; return result; }
    SACCEL operator-(const SACCEL& t) const { SACCEL result = *this; result -= t; return result; }
/*
    SACCEL operator*(REAL scalar) const { SACCEL v = *this; return v*= scalar; }
    SACCEL operator/(REAL scalar) const { SACCEL v = *this; return v/= scalar; }
    SACCEL& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SACCEL& operator*=(REAL scalar);
    SACCEL& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SACCEL& resize(unsigned rows) { assert (rows == 6); return *this; } 
    SACCEL& operator-=(const SACCEL& v);
    SACCEL& operator+=(const SACCEL& v);
    SACCEL operator+(const SACCEL& v) const;
    SACCEL operator-(const SACCEL& v) const;
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SACCEL& t)
{
  out << "acceleration (linear= " << t.get_linear() << ", angular= " << t.get_angular() << ") frame: " << t.pose;
  return out;
}

