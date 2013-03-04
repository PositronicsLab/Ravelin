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
    WRENCH() : SVECTOR6() {} 

    /// Constructs a wrench using six values: first three force, second three torque
    WRENCH(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz) : SVECTOR6(fx, fy, fz, tx, ty, tz) {};

    /// Constructs a wrench using six values: first three force, second three torque
    WRENCH(const REAL* array) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5]) {}

    /// Constructs a wrench using given force and torque
    WRENCH(const VECTOR3& f, const VECTOR3& t) : SVECTOR6(f, t) {}

    /// Computes the dot product with a twist
    REAL dot(const TWIST& t) const { return SVECTOR6::dot(*this, t); }

    void set_force(const VECTOR3& f) { set_upper(f); }
    void set_torque(const VECTOR3& t) { set_lower(t); }
    VECTOR3 get_force() const { return get_upper(); }
    VECTOR3 get_torque() const { return get_lower(); }
    WRENCH& operator=(const WRENCH& source) { SVECTOR6::operator=(source); return *this; } 
    WRENCH& operator=(const SVECTOR6& source) { SVECTOR6::operator=(source); return *this; } 
/*
    WRENCH operator-() const;
    WRENCH operator*(REAL scalar) const { WRENCH v = *this; return v*= scalar; }
    WRENCH operator/(REAL scalar) const { WRENCH v = *this; return v/= scalar; }
    WRENCH& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    WRENCH& operator*=(REAL scalar);
    WRENCH& operator-=(const WRENCH& v);
    WRENCH& operator+=(const WRENCH& v);
    WRENCH operator+(const WRENCH& v) const;
    WRENCH operator-(const WRENCH& v) const;
    WRENCH& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    WRENCH& resize(unsigned rows) { assert (rows == 6); return *this; } 
*/
}; // end class


