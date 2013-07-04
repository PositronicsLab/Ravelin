/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SMOMENTUM 
#error This class is not to be included by the user directly. Use SMomentumd.h or SMomentumf.h instead.
#endif

class SMOMENTUM : public SVECTOR6 
{
  public:
    /// Constructs a spatial momentum with zero linear and angular components
    SMOMENTUM(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(pose) {} 

    /// Constructs a spatial momentum from the SVector6 
    SMOMENTUM(const SVECTOR6& w) : SVECTOR6(w.get_upper(), w.get_lower(), w.pose) { }

    /// Constructs a spatial momentum using six values- first three linear, second three angular- and a pose
    SMOMENTUM(REAL lx, REAL ly, REAL lz, REAL ax, REAL ay, REAL az, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(lx, ly, lz, ax, ay, az, pose) {};

    /// Constructs a spatial momentum using six values- first three linear, second three angular and a pose
    SMOMENTUM(const REAL* array, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5], pose) {}

    /// Constructs a spatial momentum using given linear and angular and pose
    SMOMENTUM(const VECTOR3& l, const VECTOR3& a, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) : SVECTOR6(l, a, pose) {}

    /// Constructs a zero spatial momentum
    static SMOMENTUM zero(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>()) { SMOMENTUM w(pose); w.set_zero(); return w; }

    template <class V>
    static SMOMENTUM from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SMOMENTUM w(pose);
      REAL* wdata = w.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), wdata, 1);
      return w;
    }

    void set_linear(const VECTOR3& l) { set_upper(l); }
    void set_angular(const VECTOR3& a) { set_lower(a); }
    VECTOR3 get_linear() const { return get_upper(); }
    VECTOR3 get_angular() const { return get_lower(); }
    SMOMENTUM& operator=(const SMOMENTUM& source) { SVECTOR6::operator=(source); return *this; } 
    SMOMENTUM& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; } 
    SMOMENTUM operator-() const { SMOMENTUM w = *this; w.negate(); return w; }
    SMOMENTUM operator+(const SMOMENTUM& w) const { SMOMENTUM result = *this; result += w; return result; }
    SMOMENTUM operator-(const SMOMENTUM& w) const { SMOMENTUM result = *this; result -= w; return result; }
/*
    SMOMENTUM operator-() const;
    SMOMENTUM operator*(REAL scalar) const { SMOMENTUM v = *this; return v*= scalar; }
    SMOMENTUM operator/(REAL scalar) const { SMOMENTUM v = *this; return v/= scalar; }
    SMOMENTUM& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SMOMENTUM& operator*=(REAL scalar);
    SMOMENTUM& operator-=(const SMOMENTUM& v);
    SMOMENTUM& operator+=(const SMOMENTUM& v);
    SMOMENTUM& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SMOMENTUM& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class

inline std::ostream& operator<<(std::ostream& out, const SMOMENTUM& w)
{
  out << "Momentum (linear = " << w.get_linear() << ", angular = " << w.get_angular() << ") frame: " << w.pose;
  return out;
}

