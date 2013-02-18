/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef TWIST 
#error This class is not to be included by the user directly. Use Twistd.h or Twistf.h instead.
#endif

class TWIST : public SVECTOR6 
{
  public:
    TWIST() : SVECTOR6() {} 
    TWIST(REAL ax, REAL ay, REAL az, REAL lx, REAL ly, REAL lz) : SVECTOR6(ax, ay, az, fx, fy, rz) {};
    TWIST(const REAL* array) : SVECTOR6(array[3], array[4], array[5], array[0], array[1], array[2]) {}
    TWIST(const VECTOR3& angular, const VECTOR3& linear) : SVECTOR6(angular, linear) {}
    REAL dot(const SVECTOR6& v) const { return dot(*this, v); }
    void set_linear(const VECTOR3& lin) { set_lower(lin); }
    void set_angular(const VECTOR3& ang) { set_upper(ang); }
    VECTOR3 get_angular() const { return get_upper(); }
    VECTOR3 get_linear() const { return get_lower(); }
    TWIST& operator=(const TWIST& source);
    TWIST& operator=(const SVECTOR6& source);
    void transpose();
    static SVECTOR6 transpose(const TWIST& v);
    REAL& operator[](const unsigned i) { assert(i < 6); return _data[i]; }
    REAL operator[](const unsigned i) const { assert(i < 6); return _data[i]; }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    TWIST operator-() const;
    TWIST operator*(REAL scalar) const { TWIST v = *this; return v*= scalar; }
    TWIST operator/(REAL scalar) const { TWIST v = *this; return v/= scalar; }
    TWIST& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    TWIST& operator*=(REAL scalar);
    TWIST& operator-=(const TWIST& v);
    TWIST& operator+=(const TWIST& v);
    TWIST operator+(const TWIST& v) const;
    TWIST operator-(const TWIST& v) const;
    unsigned rows() const { return 6; }
    unsigned columns() const { return 1; }
    TWIST& resize(unsigned rows, unsigned columns) { assert (rows == 6 && columns == 1); return *this; } 
    TWIST& resize(unsigned rows) { assert (rows == 6); return *this; } 

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


