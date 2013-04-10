/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SVECTOR6
#error This class is not to be included by the user directly. Use Wrenchd.h, Wrenchf.h, Twistd.h, or Twistf.h instead.
#endif

class POSE3;

/// A 6-dimensional floating-point vector for use with spatial algebra
/**
 * Note that spatial algebra defines the dot product in an unusual manner: if vector x = [a; b] and 
 * vector y = [c; d] then x'y = [b'; a'][c d] = dot(b,c) + dot(a,d).
 */
class SVECTOR6
{
  public:
    SVECTOR6() {}
    SVECTOR6(REAL x, REAL y, REAL z, REAL a, REAL b, REAL c);
    SVECTOR6(const REAL* array);
    SVECTOR6(const VECTOR3& upper, const VECTOR3& lower);
    unsigned size() const { return 6; }
    static SVECTOR6 spatial_cross(const SVECTOR6& v1, const SVECTOR6& v2);
    static SVECTOR6 zero() { return SVECTOR6(0,0,0,0,0,0); }
    REAL dot(const SVECTOR6& v) const { return dot(*this, v); }
    static REAL dot(const SVECTOR6& v1, const SVECTOR6& v2);
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
    SVECTOR6 operator*(REAL scalar) const { SVECTOR6 v = *this; return v*= scalar; }
    SVECTOR6 operator/(REAL scalar) const { SVECTOR6 v = *this; return v/= scalar; }
    SVECTOR6& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SVECTOR6& operator*=(REAL scalar);
    SVECTOR6& operator-=(const SVECTOR6& v);
    SVECTOR6& operator+=(const SVECTOR6& v);
    SVECTOR6 operator+(const SVECTOR6& v) const;
    SVECTOR6 operator-(const SVECTOR6& v) const;
    unsigned rows() const { return 6; }
    unsigned columns() const { return 1; }
    SVECTOR6& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    SVECTOR6& resize(unsigned rows) { assert (rows == 6); return *this; } 
    ITERATOR begin();
    CONST_ITERATOR begin() const;
    ITERATOR end();
    CONST_ITERATOR end() const;

    /// The frame that this vector is defined in
    boost::shared_ptr<POSE3> pose;

  private:
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


