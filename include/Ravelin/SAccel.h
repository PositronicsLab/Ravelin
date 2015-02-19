/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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

    /// Constructs a acceleration with zero linear and zero angular components
    SACCEL(boost::shared_ptr<POSE3> pose) : SVECTOR6(pose) {}

    /// Constructs a acceleration from a SVector6
    explicit SACCEL(const SVECTOR6& v) : SVECTOR6(v.get_upper(), v.get_lower(), v.pose) {} 

    /// Constructs a acceleration from six values (first three angualr, next three linear) and a pose
    SACCEL(REAL ax, REAL ay, REAL az, REAL lx, REAL ly, REAL lz, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a acceleration from six values (first three angualr, next three linear) and a pose
    SACCEL(REAL ax, REAL ay, REAL az, REAL lx, REAL ly, REAL lz, boost::shared_ptr<POSE3> pose) : SVECTOR6(ax, ay, az, lx, ly, lz, pose) {};

    /// Constructs a acceleration from six values (first three angular, next three linear) and a pose
    SACCEL(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a acceleration from six values (first three angular, next three linear) and a pose
    SACCEL(const REAL* array, boost::shared_ptr<POSE3> pose) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a acceleration from linear and angular components and a pose
    SACCEL(const VECTOR3& angular, const VECTOR3& linear, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(angular, linear, pose) {}

    /// Constructs a acceleration from linear and angular components and a pose
    SACCEL(const VECTOR3& angular, const VECTOR3& linear, boost::shared_ptr<POSE3> pose) : SVECTOR6(angular, linear, pose) {}

    /// Returns a zero acceleration
    static SACCEL zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SACCEL t(pose); t.set_zero(); return t; }

    /// Returns a zero acceleration
    static SACCEL zero(boost::shared_ptr<POSE3> pose) { SACCEL t(pose); t.set_zero(); return t; }

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

    template <class V>
    static SACCEL from_vector(const V& v, boost::shared_ptr<POSE3> pose)
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

    /// Returns the negation of this vector
    SACCEL operator-() const
    {
      SACCEL v;
      v._data[0] = -_data[0]; 
      v._data[1] = -_data[1]; 
      v._data[2] = -_data[2]; 
      v._data[3] = -_data[3]; 
      v._data[4] = -_data[4]; 
      v._data[5] = -_data[5]; 
      v.pose = pose;

      return v;
    }

    SACCEL& operator-=(const SACCEL& v)
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

    SACCEL& operator+=(const SACCEL& v)
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

    SACCEL operator+(const SACCEL& v) const { SACCEL x = *this; x += v; return x; }
    SACCEL operator-(const SACCEL& v) const { SACCEL x = *this; x -= v; return x; }
    SACCEL operator*(REAL scalar) const { SACCEL v = *this; v*= scalar; return v;}
    SACCEL operator/(REAL scalar) const { SACCEL v = *this; v/= scalar; return v;}
/*
    SACCEL& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SACCEL& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SACCEL& t)
{
  out << "acceleration (linear= " << t.get_linear() << ", angular= " << t.get_angular() << ") frame: " << t.pose;
  return out;
}

