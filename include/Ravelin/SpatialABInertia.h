/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SPATIAL_AB_INERTIA
#error This class is not to be included by the user directly. Use SPATIAL_AB_INERTIAd.h or SPATIAL_AB_INERTIAf.h instead.
#endif

/// A 6x6 spatial algebra matrix typically used for dynamics calculations
class SPATIAL_AB_INERTIA 
{
  public:
    SPATIAL_AB_INERTIA();
    SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J);
    SPATIAL_AB_INERTIA(const SPATIAL_AB_INERTIA& source) { operator=(source); }
    SPATIAL_AB_INERTIA(const SPATIAL_RB_INERTIA& source) { operator=(source); }
    SPATIAL_AB_INERTIA(const MatrixN& M);
    void set_zero();
    static SPATIAL_AB_INERTIA zero() { SPATIAL_AB_INERTIA m; m.set_zero(); return m; }
    TWIST inverse_mult(const WRENCH& f) const;
    SPATIAL_AB_INERTIA& operator=(const SPATIAL_RB_INERTIA& source);
    SPATIAL_AB_INERTIA& operator=(const SPATIAL_AB_INERTIA& source);
    SPATIAL_AB_INERTIA& operator+=(const SPATIAL_AB_INERTIA& m);
    SPATIAL_AB_INERTIA& operator-=(const SPATIAL_AB_INERTIA& m);
    SPATIAL_AB_INERTIA& operator*=(const SPATIAL_AB_INERTIA& m) { return *this = operator*(m); }
    SPATIAL_AB_INERTIA& operator*=(REAL scalar);
    SPATIAL_AB_INERTIA& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    SPATIAL_AB_INERTIA operator+(const SPATIAL_RB_INERTIA& m) const;
    SPATIAL_AB_INERTIA& operator+=(const SPATIAL_RB_INERTIA& m) { return *this = *this + m; }
    SPATIAL_AB_INERTIA operator+(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator-(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator*(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator*(REAL scalar) const { SPATIAL_AB_INERTIA m = *this; m *= scalar; return m; }
    SPATIAL_AB_INERTIA operator/(REAL scalar) const { return operator*((REAL) 1.0/scalar); }
    WRENCH operator*(const TWIST& t) const;
    SPATIAL_AB_INERTIA operator-() const;
    static SPATIAL_AB_INERTIA inverse_inertia(const SPATIAL_AB_INERTIA& I);    
    static SPATIAL_AB_INERTIA mult(const MatrixN& M, const SPATIAL_AB_INERTIA& I);

    /// The upper left / lower right hand matrix 'H'
    MATRIX3 H;

    /// The lower left matrix
    MATRIX3 J;

    /// The upper right matrix
    MATRIX3 M;

    /// The pose that this inertia is defined in
    boost::shared_ptr<POSE> pose;
}; // end class

std::ostream& operator<<(std::ostream& out, const SPATIAL_AB_INERTIA& m);

