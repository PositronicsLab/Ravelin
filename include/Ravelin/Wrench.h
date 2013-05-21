/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef WRENCH 
#error This class is not to be included by the user directly. Use Wrenchd.h or Wrenchf.h instead.
#endif

class TWIST;

class WRENCH : public SVECTOR6 
{
  public:
    /// Constructs a wrench with zero force and torque components
    WRENCH(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {} 

    /// Constructs a wrench from the SVector6 
    WRENCH(const SVECTOR6& w) : SVECTOR6(w.get_upper(), w.get_lower(), w.pose) { }

    /// Constructs a wrench using six values- first three force, second three torque- and a pose
    WRENCH(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(fx, fy, fz, tx, ty, tz, pose) {};

    /// Constructs a wrench using six values- first three force, second three torque0 and a pose
    WRENCH(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a wrench using given force and torque and pose
    WRENCH(const VECTOR3& f, const VECTOR3& t, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(f, t, pose) {}

    /// Constructs a zero wrench
    static WRENCH zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { WRENCH w(pose); w.set_zero(); return w; }

    template <class V>
    static WRENCH from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      WRENCH w(pose);
      REAL* wdata = w.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), wdata, 1);
      return w;
    }

    void set_force(const VECTOR3& f) { set_upper(f); }
    void set_torque(const VECTOR3& t) { set_lower(t); }
    VECTOR3 get_force() const { return get_upper(); }
    VECTOR3 get_torque() const { return get_lower(); }
    WRENCH& operator=(const WRENCH& source) { SVECTOR6::operator=(source); return *this; } 
    WRENCH& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; } 
    WRENCH operator-() const { WRENCH w = *this; w.negate(); return w; }
    WRENCH operator+(const WRENCH& w) const { WRENCH result = *this; result += w; return result; }
    WRENCH operator-(const WRENCH& w) const { WRENCH result = *this; result -= w; return result; }
/*
    WRENCH operator-() const;
    WRENCH operator*(REAL scalar) const { WRENCH v = *this; return v*= scalar; }
    WRENCH operator/(REAL scalar) const { WRENCH v = *this; return v/= scalar; }
    WRENCH& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    WRENCH& operator*=(REAL scalar);
    WRENCH& operator-=(const WRENCH& v);
    WRENCH& operator+=(const WRENCH& v);
    WRENCH& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    WRENCH& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const WRENCH& w)
{
  out << "Wrench (force = " << w.get_force() << ", torque = " << w.get_torque() << ") frame: " << w.pose;
  return out;
}

