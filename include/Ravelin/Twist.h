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
    TWIST() : SVECTOR6() {}
 
    /// Constructs a twist from six values (first three linear, next three angular)
    TWIST(REAL lx, REAL ly, REAL lz, REAL ax, REAL ay, REAL az) : SVECTOR6(ax, ay, az, lx, ly, lz) {};

    /// Constructs a twist from six values (first three linear, next three angular)
    TWIST(const REAL* array) : SVECTOR6(array[3], array[4], array[5], array[0], array[1], array[2]) {}

    /// Constructs a twist from linear and angular components
    TWIST(const VECTOR3& linear, const VECTOR3& angular) : SVECTOR6(angular, linear) {}

    /// Returns a zero twist
    static TWIST zero() { TWIST t; t.set_zero(); return t; }

    /// Computes the dot product with a wrench
    REAL dot(const WRENCH& w) const { return SVECTOR6::dot(*this, w); }
    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    TWIST& operator=(const TWIST& source) { SVECTOR6::operator=(source); return *this; }
    TWIST& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; }
/*
    TWIST operator-() const;
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


