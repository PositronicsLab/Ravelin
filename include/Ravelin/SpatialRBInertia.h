/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SPATIAL_RB_INERTIA
#error This class is not to be included by the user directly. Use SpatialRBInertiad.h or SpatialRBInertiaf.h instead.
#endif

class POSE3;

/// A 6x6 spatial algebra matrix used for dynamics calculations
/** 
 * The matrix is represented by:
 * | -hx*m        I*m  |
 * | J - hx*hx*m  hx*m |
 * where hx is the skew symmetric matrix determined by h (h is the vector from
 * the origin of the reference frame to the center of mass) and I is the 
 * identity matrix.
 * The inverse of this matrix is:
 * | -inv(J)*hx          inv(J)    |
 * | I/m - hx*inv(J)*hx  hx*inv(J) |
 */
class SPATIAL_RB_INERTIA
{
  public:
    SPATIAL_RB_INERTIA(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    SPATIAL_RB_INERTIA(boost::shared_ptr<POSE3> pose);
    SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J, boost::shared_ptr<POSE3> pose);
    SPATIAL_RB_INERTIA(const SPATIAL_RB_INERTIA& source) { operator=(source); }
    void set_zero();
    static SPATIAL_RB_INERTIA zero() { SPATIAL_RB_INERTIA m; m.set_zero(); return m; }
    SACCEL inverse_mult(const SFORCE& v) const;
    SVELOCITY inverse_mult(const SMOMENTUM& m) const;
    std::vector<SACCEL>& inverse_mult(const std::vector<SFORCE>& v, std::vector<SACCEL>& result) const;
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
    SFORCE operator*(const SACCEL& s) const;
    std::vector<SFORCE>& mult(const std::vector<SACCEL>& s, std::vector<SFORCE>& result) const;
    SMOMENTUM operator*(const SVELOCITY& s) const;
    SMOMENTUM mult(const SVELOCITY& s) const { return operator*(s); }
    std::vector<SMOMENTUM>& mult(const std::vector<SVELOCITY>& s, std::vector<SMOMENTUM>& result) const;
    SPATIAL_RB_INERTIA operator-() const;

    /// The rigid body mass
    REAL m;

    /// The position of the center-of-mass times the mass 
    ORIGIN3 h;

    /// The rigid body moment of inertia matrix
    MATRIX3 J;

    /// The pose that this inertia is defined in
    boost::shared_ptr<const POSE3> pose;

    /// Converts this to a matrix
    template <class Mat>
    Mat& to_matrix(Mat& M) const
    {
      const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;

      // resize the matrix
      M.resize(SPATIAL_DIM, SPATIAL_DIM);

      // precompute matrices
      ORIGIN3 hm = h*m;
      MATRIX3 hx = MATRIX3::skew_symmetric(h);
      MATRIX3 hxm = MATRIX3::skew_symmetric(hm);

      // setup the 3x3 blocks
      M.set_sub_mat(0,0, hxm, eTranspose);
      M.set_sub_mat(3,0, J-hx*hxm);
      M.set_sub_mat(0,3, MATRIX3(m, 0, 0, 0, m, 0, 0, 0, m));
      M.set_sub_mat(3,3, hxm);

      return M;
    }

    /// Converts this to a positive-definite matrix
    template <class Mat>
    Mat& to_PD_matrix(Mat& M) const
    {
      const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;

      // resize the matrix
      M.resize(SPATIAL_DIM, SPATIAL_DIM);

      // precompute matrices
      MATRIX3 hx = MATRIX3::skew_symmetric(h);
      MATRIX3 mhx = MATRIX3::skew_symmetric(h*m);

      // setup the 3x3 blocks
      M.set_sub_mat(0,3, mhx, eTranspose);
      M.set_sub_mat(3,3, J - (mhx*hx));
      M.set_sub_mat(0,0, MATRIX3(m, 0, 0, 0, m, 0, 0, 0, m));
      M.set_sub_mat(3,0, mhx);

      return M;
    }

  private:
    void mult_spatial(const SVECTOR6& t, SVECTOR6& result) const;
    void mult_spatial(const SVECTOR6& t, const MATRIX3& hxm, const MATRIX3& Jstar, SVECTOR6& result) const;
    void inverse_mult_spatial(const SVECTOR6& w, SVECTOR6& result) const;
    void inverse_mult_spatial(const SVECTOR6& w, const MATRIX3& iJ, const MATRIX3& hx, const MATRIX3& hxiJ, REAL m, SVECTOR6& result) const;

}; // end class

std::ostream& operator<<(std::ostream& out, const SPATIAL_RB_INERTIA& m);


