/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SPATIAL_RB_INERTIA
#error This class is not to be included by the user directly. Use SpatialRBInertiad.h or SpatialRBInertiaf.h instead.
#endif

/// A 6x6 spatial algebra matrix used for dynamics calculations
/** 
 * The matrix is represented by:
 * | J - hx*hx*m  hx*m |
 * | -hx*m        I*m  |
 * where hx is the skew symmetric matrix determined by h and I is the identity
 * matrix.
 */
class SPATIAL_RB_INERTIA
{
  public:
    SPATIAL_RB_INERTIA();
    SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J);
    SPATIAL_RB_INERTIA(const SPATIAL_RB_INERTIA& source) { operator=(source); }
    void set_zero();
    static SPATIAL_RB_INERTIA zero() { SPATIAL_RB_INERTIA m; m.set_zero(); return m; }
    void to_matrix(REAL data[]) const;
    TWIST inverse_mult(const WRENCH& v) const;
    SPATIAL_RB_INERTIA& operator=(const SPATIAL_RB_INERTIA& source);
    SPATIAL_RB_INERTIA& operator+=(const SPATIAL_RB_INERTIA& m);
    SPATIAL_RB_INERTIA& operator-=(const SPATIAL_RB_INERTIA& m);
    SPATIAL_RB_INERTIA& operator*=(const SPATIAL_RB_INERTIA& m) { return *this = operator*(m); }
    SPATIAL_RB_INERTIA& operator*=(REAL scalar);
    SPATIAL_RB_INERTIA& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    SPATIAL_RB_INERTIA operator+(const SPATIAL_RB_INERTIA& m) const;
    SPATIAL_RB_INERTIA operator-(const SPATIAL_RB_INERTIA& m) const;
    SPATIAL_RB_INERTIA operator*(const SPATIAL_RB_INERTIA& m) const;
    SPATIAL_RB_INERTIA operator*(REAL scalar) const { SPATIAL_RB_INERTIA m = *this; m *= scalar; return m; }
    SPATIAL_RB_INERTIA operator/(REAL scalar) const { return operator*((REAL) 1.0/scalar); }
    WRENCH operator*(const TWIST& v) const;
    SPATIAL_RB_INERTIA operator-() const;

    /// The rigid body mass
    REAL m;

    /// The rigid body offset vector
    VECTOR3 h;

    /// The rigid body moment of inertia matrix
    MATRIX3 J;

    /// The pose that this inertia is defined in
    boost::shared_ptr<POSE> pose;
}; // end class

std::ostream& operator<<(std::ostream& out, const SPATIAL_RB_INERTIA& m);


