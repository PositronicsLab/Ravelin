/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SVECTOR6
#error This class is not to be included by the user directly. Use SForced.h, SForcef.h, SMomentumd.h, SMomentumf.h, SVelocityd.h, SVelocityf.h, SAcceld.h, or SAccelf.h instead.
#endif

class POSE3;
class SFORCE;
class SVELOCITY;
class SACCEL;
class SMOMENTUM;

/// A 6-dimensional floating-point vector for use with spatial algebra
/**
 * Note that spatial algebra defines the dot product in an unusual manner: if vector x = [a; b] and 
 * vector y = [c; d] then x'y = [b'; a'][c d] = dot(b,c) + dot(a,d).
 */
class SVECTOR6
{
  public:
    SVECTOR6() {}
    SVECTOR6(boost::shared_ptr<const POSE3> pose) { this->pose = pose; }
    SVECTOR6(REAL x, REAL y, REAL z, REAL a, REAL b, REAL c);
    SVECTOR6(REAL x, REAL y, REAL z, REAL a, REAL b, REAL c, boost::shared_ptr<const POSE3> pose);
    SVECTOR6(const REAL* array);
    SVECTOR6(const REAL* array, boost::shared_ptr<const POSE3> pose);
    SVECTOR6(const VECTOR3& upper, const VECTOR3& lower);
    SVECTOR6(const VECTOR3& upper, const VECTOR3& lower, boost::shared_ptr<const POSE3> pose);
    unsigned size() const { return 6; }
    static SVECTOR6 zero() { return SVECTOR6(0,0,0,0,0,0); }
    void set_zero() { std::fill_n(_data, 6, (REAL) 0.0); }
    void set_lower(const VECTOR3& lower);
    void set_upper(const VECTOR3& upper);
    VECTOR3 get_lower() const;
    VECTOR3 get_upper() const;
    SVECTOR6& operator=(const SVECTOR6& source);
    REAL& operator[](const unsigned i) { assert(i < 6); return _data[i]; }
    const REAL& operator[](const unsigned i) const { assert(i < 6); return _data[i]; }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    SVECTOR6 operator-() const;
    SVECTOR6& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SVECTOR6& operator*=(REAL scalar);
    unsigned rows() const { return 6; }
    unsigned columns() const { return 1; }
    SVECTOR6& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SVECTOR6& resize(unsigned rows) { assert (rows == 6); return *this; } 
    ITERATOR begin();
    CONST_ITERATOR begin() const;
    ITERATOR end();
    CONST_ITERATOR end() const;
    SVECTOR6& negate() { std::transform(_data, _data+6, _data, std::negate<REAL>()); return *this; }
    unsigned inc() const { return 1; }
    unsigned leading_dim() const { return 6; }

    template <class V>
    V& to_vector(V& v) const
    {
      const unsigned SPATIAL_DIM = 6;
      v.resize(SPATIAL_DIM);
      REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, _data, 1, vdata, inc());
      return v;
    }

    template <class Vec>
    REAL dot(const Vec& v) const { return dot(*this, v); }

    template <class Vec>
    static REAL dot(const SVECTOR6& v1, const Vec& v2)
    {
      const REAL* d1 = v1.data();
      const REAL* d2 = v2.data();
      return d1[3]+d2[0] + d1[4]+d2[1] + d1[5]+d2[2]+
             d1[0]+d2[3] + d1[1]+d2[4] + d1[2]+d2[5]; 
    }

    /// The frame that this vector is defined in
    boost::shared_ptr<const POSE3> pose;

    template <class V>
    static SVECTOR6 from_vector(const V& v, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      const unsigned SPATIAL_DIM = 6;
      if (v.size() != SPATIAL_DIM)
        throw MissizeException();
      SVECTOR6 s(pose);
      REAL* sdata = s.data();
      const REAL* vdata = v.data();
      CBLAS::copy(SPATIAL_DIM, vdata, v.inc(), sdata, 1);
      return s;
    }

  protected:
    REAL _data[6];

  public:

    template <class V>
    SVECTOR6(const V& v)
    {
      const unsigned SPATIAL_DIM = 6;
      #ifndef NEXCEPT
      if (v.rows()*v.columns() != SPATIAL_DIM)
        throw MissizeException();
      #endif

      std::copy(v.begin(), v.end(), begin());
    }
}; // end class

/*
template <>
REAL SVECTOR6::dot(const SVECTOR6& v1, const SFORCE& w);

template <>
REAL SVECTOR6::dot(const SVECTOR6& v1, const SVELOCITY& t);

template <>
REAL SVECTOR6::dot(const SVECTOR6& v1, const SACCEL& t);
*/

inline std::ostream& operator<<(std::ostream& out, const SVECTOR6& v)
{
  out << "Spatial vector (upper = " << v.get_upper() << ", lower= " << v.get_lower() << ") frame: " << v.pose;
  return out;
}

