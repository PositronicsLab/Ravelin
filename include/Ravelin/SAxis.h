/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SAXIS 
#error This class is not to be included by the user directly. Use SAxisf.h or SAxisdd.h instead.
#endif

class SFORCE;
class SMOMENTUM;

class SAXIS : public SVECTOR6 
{
  public:
    /// Constructs a axis with zero linear and zero angular components
    SAXIS(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {}

    /// Constructs a axis from a SVector6
    explicit SAXIS(const SVECTOR6& v) : SVECTOR6(v.get_upper(), v.get_lower(), v.pose) {} 

    /// Constructs a axis from six values (first three angular, next three linear) and a pose
    SAXIS(REAL ax, REAL ay, REAL az, REAL lx, REAL ly, REAL lz, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a axis from six values (first three angular, next three linear) and a pose
    SAXIS(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a axis from linear and angular components and a pose
    SAXIS(const VECTOR3& angular, const VECTOR3& linear, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(angular, linear, pose) {}

    /// Returns a zero axis
    static SAXIS zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SAXIS t(pose); t.set_zero(); return t; }

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
    static SAXIS from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SAXIS t(pose);
      REAL* tdata = t.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), tdata, 1);
      return t;
    }

    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    SAXIS& operator=(const SAXIS& source) { SVECTOR6::operator=(source); return *this; }
    SAXIS& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; }

    /// Returns the negation of this vector
    SAXIS operator-() const
    {
      SAXIS v;
      v._data[0] = -_data[0]; 
      v._data[1] = -_data[1]; 
      v._data[2] = -_data[2]; 
      v._data[3] = -_data[3]; 
      v._data[4] = -_data[4]; 
      v._data[5] = -_data[5]; 
      v.pose = pose;

      return v;
    }

    SAXIS& operator-=(const SAXIS& v)
    {
      #ifndef NEXCEPT
      if (pose != v.pose)
        throw FrameException();
      #endif

      _data[0] -= v._data[0];
      _data[1] -= v._data[1];
      _data[2] -= v._data[2];
      _data[3] -= v._data[3];
      _data[4] -= v._data[4];
      _data[5] -= v._data[5];
     
      return *this;
    }

    SAXIS& operator+=(const SAXIS& v)
    {
      #ifndef NEXCEPT
      if (pose != v.pose)
        throw FrameException();
      #endif

      _data[0] += v._data[0];
      _data[1] += v._data[1];
      _data[2] += v._data[2];
      _data[3] += v._data[3];
      _data[4] += v._data[4];
      _data[5] += v._data[5];
     
      return *this;
    }

    REAL dot(const SMOMENTUM& m) const;
    REAL dot(const SFORCE& f) const;
    SAXIS operator+(const SAXIS& v) const { SAXIS x = *this; x += v; return x; }
    SAXIS operator-(const SAXIS& v) const { SAXIS x = *this; x -= v; return x; }
    SAXIS operator*(REAL scalar) const { SAXIS v = *this; v*= scalar; return v; }
    SAXIS operator/(REAL scalar) const { SAXIS v = *this; v/= scalar; return v; }
/*
    SAXIS& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SAXIS& operator*=(REAL scalar);
    SAXIS& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SAXIS& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SAXIS& t)
{
  out << "axis (linear= " << t.get_linear() << ", angular= " << t.get_angular() << ") frame: " << t.pose;
  return out;
}

