/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef AANGLE
#error This class is not to be included by the user directly. Use AANGLE.h or AAnglef.h instead.
#endif

class QUAT;
class MATRIX3;

/// Class for representation of orientation by an angle around an axis
class AANGLE
{
  public:
    AANGLE();
    AANGLE(const AANGLE& source);
    AANGLE(REAL x, REAL y, REAL z, REAL angle);
    AANGLE(const VECTOR3& v, REAL angle);
    AANGLE(const MATRIX3& m);
    AANGLE(const QUAT& q);
    AANGLE(const MATRIX3& m, const VECTOR3& axis);
    AANGLE& set(const VECTOR3& v, REAL angle);
    AANGLE& operator=(const MATRIX3& m);
    AANGLE& operator=(const QUAT& q);
    AANGLE& operator=(const AANGLE& source);
    void set(const MATRIX3& m, const VECTOR3& axis);
    AANGLE operator*(const AANGLE& a) const;
    AANGLE& operator*=(const AANGLE& a);

    /// Does a "safe" square root - input values sufficiently close to -0 are returned as 0
    static REAL safe_sqrt(REAL x) { return (x < (REAL) 0.0 && x > -EPS) ? (REAL) 0.0 : std::sqrt(x); }

    REAL angle;
    REAL x;
    REAL y;
    REAL z;
};

std::ostream& operator<<(std::ostream& out, const AANGLE& a);


