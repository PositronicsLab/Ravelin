/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SFORCE 
#error This class is not to be included by the user directly. Use SForced.h or SForcef.h instead.
#endif

class SAXIS;

class SFORCE : public SVECTOR6 
{
  public:
    /// Constructs a spatial force with zero force and torque components
    SFORCE(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {} 

    /// Constructs a spatial force with zero force and torque components
    SFORCE(boost::shared_ptr<POSE3> pose) : SVECTOR6(pose) {} 

    /// Constructs a spatial force from the SVector6 
    explicit SFORCE(const SVECTOR6& w) : SVECTOR6(w.get_upper(), w.get_lower(), w.pose) { }

    /// Constructs a spatial force using six values- first three force, second three torque- and a pose
    SFORCE(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(fx, fy, fz, tx, ty, tz, pose) {};

    /// Constructs a spatial force using six values- first three force, second three torque- and a pose
    SFORCE(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz, boost::shared_ptr<POSE3> pose) : SVECTOR6(fx, fy, fz, tx, ty, tz, pose) {};

    /// Constructs a spatial force using six values- first three force, second three torque0 and a pose
    SFORCE(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a spatial force using six values- first three force, second three torque0 and a pose
    SFORCE(const REAL* array, boost::shared_ptr<POSE3> pose) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a spatial force using given force and torque and pose
    SFORCE(const VECTOR3& f, const VECTOR3& t, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(f, t, pose) {}

    /// Constructs a spatial force using given force and torque and pose
    SFORCE(const VECTOR3& f, const VECTOR3& t, boost::shared_ptr<POSE3> pose) : SVECTOR6(f, t, pose) {}

    /// Constructs a zero spatial force
    static SFORCE zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SFORCE w(pose); w.set_zero(); return w; }

    /// Constructs a zero spatial force
    static SFORCE zero(boost::shared_ptr<POSE3> pose) { SFORCE w(pose); w.set_zero(); return w; }

    template <class V>
    static SFORCE from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SFORCE w(pose);
      REAL* wdata = w.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), wdata, 1);
      return w;
    }

    template <class V>
    static SFORCE from_vector(const V& v, boost::shared_ptr<POSE3> pose)
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SFORCE w(pose);
      REAL* wdata = w.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), wdata, 1);
      return w;
    }

    REAL dot(const SVELOCITY& v) const;
    REAL dot(const SAXIS& s) const;
    void set_force(const VECTOR3& f) { set_upper(f); }
    void set_torque(const VECTOR3& t) { set_lower(t); }
    VECTOR3 get_force() const { return get_upper(); }
    VECTOR3 get_torque() const { return get_lower(); }
    SFORCE& operator=(const SFORCE& source) { SVECTOR6::operator=(source); return *this; } 
    SFORCE& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; } 

    /// Returns the negation of this vector
    SFORCE operator-() const
    {
      SFORCE v;
      v._data[0] = -_data[0]; 
      v._data[1] = -_data[1]; 
      v._data[2] = -_data[2]; 
      v._data[3] = -_data[3]; 
      v._data[4] = -_data[4]; 
      v._data[5] = -_data[5]; 
      v.pose = pose;

      return v;
    }

    SFORCE& operator-=(const SFORCE& v)
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

    SFORCE& operator+=(const SFORCE& v)
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

    SFORCE operator+(const SFORCE& v) const { SFORCE x = *this; x += v; return x; }
    SFORCE operator-(const SFORCE& v) const { SFORCE x = *this; x -= v; return x; }
    SFORCE operator*(REAL scalar) const { SFORCE v = *this; v*= scalar; return v;}
    SFORCE operator/(REAL scalar) const { SFORCE v = *this; v/= scalar; return v;}
/*
    SFORCE& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SFORCE& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SFORCE& w)
{
  out << "Wrench (force = " << w.get_force() << ", torque = " << w.get_torque() << ") frame: " << w.pose;
  return out;
}

