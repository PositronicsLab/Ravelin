/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef WRENCH 
#error This class is not to be included by the user directly. Use Wrenchd.h or Wrenchf.h instead.
#endif

class WRENCH : public SVECTOR6 
{
  public:
    WRENCH() : SVECTOR6() {} 
    WRENCH(REAL fx, REAL fy, REAL fz, REAL tx, REAL ty, REAL tz) : SVECTOR6(fx, fy, fz, tx, ty, tz) {};
    WRENCH(const REAL* array) : SVECTOR6(array[0], array[1], array[2], array[3], array[4], array[5]) {}
    WRENCH(const VECTOR3& f, const VECTOR3& t) : SVECTOR6(f, t) {}
    REAL dot(const SVECTOR6& v) const { return dot(*this, v); }
    void set_force(const VECTOR3& f) { set_upper(f); }
    void set_torque(const VECTOR3& t) { set_lower(t); }
    VECTOR3 get_force() const { return get_upper(); }
    VECTOR3 get_torque() const { return get_lower(); }
    WRENCH& operator=(const WRENCH& source);
    WRENCH& operator=(const SVECTOR6& source);
    void transpose();
    static SVECTOR6 transpose(const WRENCH& v);
    REAL& operator[](const unsigned i) { assert(i < 6); return _data[i]; }
    REAL operator[](const unsigned i) const { assert(i < 6); return _data[i]; }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    WRENCH operator-() const;
    WRENCH operator*(REAL scalar) const { WRENCH v = *this; return v*= scalar; }
    WRENCH operator/(REAL scalar) const { WRENCH v = *this; return v/= scalar; }
    WRENCH& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    WRENCH& operator*=(REAL scalar);
    WRENCH& operator-=(const WRENCH& v);
    WRENCH& operator+=(const WRENCH& v);
    WRENCH operator+(const WRENCH& v) const;
    WRENCH operator-(const WRENCH& v) const;
    unsigned rows() const { return 6; }
    unsigned columns() const { return 1; }
    WRENCH& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    WRENCH& resize(unsigned rows) { assert (rows == 6); return *this; } 

    /// Get iterator to the start of the data
    REAL* begin() { return _data; }

    /// Get iterator to the start of the data
    const REAL* begin() const { return _data; }

    /// Get iterator to the end of the data
    REAL* end() { return _data + 6; }

    /// Get iterator to the end of the data
    const REAL* end() const { return _data + 6; }

  private:
    REAL _data[6];
}; // end class


