/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef QUAT
#error This class is not to be included by the user directly. Use Quatd.h or Quatf.h instead.
#endif

class VECTORN;
class MATRIX3;
class AANGLE;

/// Quaternion class used for orientation
/**
 * This class is used to represent orientation via unit quaternions.  Note, however, that the user is responsible
 * for normalizing the quaternions; normalization is not automatic.
 */
class QUAT
{
  public:
    QUAT();
    QUAT(REAL x, REAL y, REAL z, REAL w);
    QUAT(const QUAT& q);
    QUAT(const VECTORN& v);
    QUAT(const AANGLE& v);
    QUAT(const MATRIX3& v);
    bool unit() const;
    static QUAT zero();
    void conjugate();
    static QUAT conjugate(const QUAT& q);
    static QUAT identity() { return QUAT((REAL) 0.0, (REAL) 0.0, (REAL) 0.0, (REAL) 1.0); }
    void slerp(const QUAT& q, REAL alpha);
    void lerp(const QUAT& q, REAL alpha);
    static QUAT slerp(const QUAT& q1, const QUAT& q2, REAL alpha);
    static QUAT lerp(const QUAT& q1, const QUAT& q2, REAL alpha);
    QUAT& inverse();
    static QUAT invert(const QUAT& q);
    REAL norm_sq() const;
    REAL norm() const { return magnitude(); }
    REAL norm_inf() const;
    void normalize();
    static QUAT normalize(const QUAT& q); 
    QUAT& operator=(const VECTORN& v);
    QUAT& operator=(const AANGLE& a);
    QUAT& operator=(const MATRIX3& m);
    QUAT& operator=(const QUAT& q);
    void calc_generalized_force(const VECTOR3& t, REAL gt[4]) const;
    unsigned size() const { return 4; }
    REAL& operator[](unsigned i);
    const REAL& operator[](unsigned i) const;
    QUAT operator-() const;
    QUAT operator-(const QUAT& q) const;
    QUAT& operator-=(const QUAT& q);
    QUAT operator+(const QUAT& q) const;
    QUAT& operator+=(const QUAT& q);
    QUAT operator*(const QUAT& q) const;
    QUAT operator/(const QUAT& q) const;
    VECTOR3 operator*(const VECTOR3& v) const;
    QUAT operator*(REAL scalar) const;
    QUAT operator/(REAL scalar) const { return operator*(1.0/scalar); }
    QUAT& operator*=(const QUAT& q);
    QUAT& operator*=(REAL scalar);
    QUAT& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    REAL magnitude() const;
    static QUAT deriv(const QUAT& q, const VECTOR3& w);
    static QUAT dderiv(const QUAT& q, const VECTOR3& omega, const VECTOR3& alpha);
    static VECTOR3 to_omega(const QUAT& q, const QUAT& qd);
//    static VECTORN to_vect(const QUAT& q);
    VECTOR3 G_mult(REAL qw, REAL qx, REAL qy, REAL qz) const;
    QUAT G_transpose_mult(const VECTOR3& v) const;
    QUAT L_transpose_mult(const VECTOR3& v) const;
//    MATRIXN& determine_G(MATRIXN& G) const;
//    MATRIXN& determine_L(MATRIXN& L) const;
    VECTOR3 L_mult(const QUAT& q) const;
    
    /// First quaternion component
    REAL x;

    /// Second quaternion component
    REAL y;

    /// Third quaterion component
    REAL z; 

    /// Fourth quaternion component
    REAL w;

  private:
    static REAL sgn(REAL x) { return (x >= (REAL) 0.0) ? (REAL) 1.0 : (REAL) -1.0; }

    /// Computes a "safe" square root
    /**
     * If the input is negative, safe_sqrt() will return 0.
     * \note the result will be incorrect if the magnitude of the input is
     *       effectively greater than zero.
     */
    static REAL safe_sqrt(REAL x) { return (x < (REAL) 0.0) ? (REAL) 0.0 : std::sqrt(x); }
};

QUAT operator*(REAL scalar, const QUAT& q);
std::ostream& operator<<(std::ostream& out, const QUAT& q);

