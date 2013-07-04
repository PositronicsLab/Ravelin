/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SFORCE 
#error This class is not to be included by the user directly. Use SForced.h or SForcef.h instead.
#endif

class SFORCE : public SVECTOR6 
{
  public:
    /// Constructs a spatial force with zero force and torque components
    SFORCE(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {} 

    /// Constructs a spatial force from the SVector6 
    SFORCE(const SVECTOR6& w) : SVECTOR6(w.get_upper(), w.get_lower(), w.pose) { }

    /// Constructs a spatial force using six values- first three force, second three torque- and a pose
    SFORCE(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(fx, fy, fz, tx, ty, tz, pose) {};

    /// Constructs a spatial force using six values- first three force, second three torque0 and a pose
    SFORCE(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a spatial force using given force and torque and pose
    SFORCE(const VECTOR3& f, const VECTOR3& t, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(f, t, pose) {}

    /// Constructs a zero spatial force
    static SFORCE zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SFORCE w(pose); w.set_zero(); return w; }

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

    void set_force(const VECTOR3& f) { set_upper(f); }
    void set_torque(const VECTOR3& t) { set_lower(t); }
    VECTOR3 get_force() const { return get_upper(); }
    VECTOR3 get_torque() const { return get_lower(); }
    SFORCE& operator=(const SFORCE& source) { SVECTOR6::operator=(source); return *this; } 
    SFORCE& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; } 
    SFORCE operator-() const { SFORCE w = *this; w.negate(); return w; }
    SFORCE operator+(const SFORCE& w) const { SFORCE result = *this; result += w; return result; }
    SFORCE operator-(const SFORCE& w) const { SFORCE result = *this; result -= w; return result; }
/*
    SFORCE operator-() const;
    SFORCE operator*(REAL scalar) const { SFORCE v = *this; return v*= scalar; }
    SFORCE operator/(REAL scalar) const { SFORCE v = *this; return v/= scalar; }
    SFORCE& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SFORCE& operator*=(REAL scalar);
    SFORCE& operator-=(const SFORCE& v);
    SFORCE& operator+=(const SFORCE& v);
    SFORCE& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SFORCE& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SFORCE& w)
{
  out << "Wrench (force = " << w.get_force() << ", torque = " << w.get_torque() << ") frame: " << w.pose;
  return out;
}

