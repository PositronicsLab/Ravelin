/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SPATIAL_AB_INERTIA
#error This class is not to be included by the user directly. Use SPATIAL_AB_INERTIAd.h or SPATIAL_AB_INERTIAf.h instead.
#endif

class POSE3;

/// A 6x6 spatial algebra matrix typically used for dynamics calculations
/**
 * Note: this representation differs from Featherstone's new spatial AB
 * inertia representation:
 * | H' M |
 * | J  H |
 */
class SPATIAL_AB_INERTIA 
{
  public:
    SPATIAL_AB_INERTIA(boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    SPATIAL_AB_INERTIA(boost::shared_ptr<POSE3> pose);
    SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>());
    SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J, boost::shared_ptr<POSE3> pose = boost::shared_ptr<POSE3>());
    SPATIAL_AB_INERTIA(const SPATIAL_AB_INERTIA& source) { operator=(source); }
    SPATIAL_AB_INERTIA(const SPATIAL_RB_INERTIA& source) { operator=(source); }
    void set_zero();
    static SPATIAL_AB_INERTIA zero() { SPATIAL_AB_INERTIA m; m.set_zero(); return m; }
    SACCEL inverse_mult(const SFORCE& f) const;
    SVELOCITY inverse_mult(const SMOMENTUM& m) const;
    std::vector<SACCEL>& inverse_mult(const std::vector<SFORCE>& w, std::vector<SACCEL>& result) const;
    SPATIAL_AB_INERTIA& operator=(const SPATIAL_RB_INERTIA& source);
    SPATIAL_AB_INERTIA& operator=(const SPATIAL_AB_INERTIA& source);
    SPATIAL_AB_INERTIA& operator+=(const SPATIAL_AB_INERTIA& m);
    SPATIAL_AB_INERTIA& operator-=(const SPATIAL_AB_INERTIA& m);
    SPATIAL_AB_INERTIA& operator*=(const SPATIAL_AB_INERTIA& m) { return *this = operator*(m); }
    SPATIAL_AB_INERTIA& operator*=(REAL scalar);
    SPATIAL_AB_INERTIA& operator/=(REAL scalar) { return operator*=(1.0/scalar); }
    SPATIAL_RB_INERTIA to_rb_inertia() const;
    SPATIAL_AB_INERTIA operator+(const SPATIAL_RB_INERTIA& m) const;
    SPATIAL_AB_INERTIA& operator+=(const SPATIAL_RB_INERTIA& m) { return *this = *this + m; }
    SPATIAL_AB_INERTIA operator+(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator-(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator*(const SPATIAL_AB_INERTIA& m) const;
    SPATIAL_AB_INERTIA operator*(REAL scalar) const { SPATIAL_AB_INERTIA m = *this; m *= scalar; return m; }
    SPATIAL_AB_INERTIA operator/(REAL scalar) const { return operator*((REAL) 1.0/scalar); }
    SFORCE operator*(const SACCEL& s) const { return mult(s); }
    SFORCE mult(const SACCEL& s) const;
    std::vector<SFORCE>& mult(const std::vector<SACCEL>& s, std::vector<SFORCE>& result) const;
    SMOMENTUM operator*(const SVELOCITY& s) const { return mult(s); }
    SMOMENTUM mult(const SVELOCITY& s) const;
    std::vector<SMOMENTUM>& mult(const std::vector<SVELOCITY>& s, std::vector<SMOMENTUM>& result) const;
    SPATIAL_AB_INERTIA operator-() const;
    static SPATIAL_AB_INERTIA inverse_inertia(const SPATIAL_AB_INERTIA& I);    

    template <class Mat>
    static SPATIAL_AB_INERTIA from_matrix(const Mat& m, boost::shared_ptr<const POSE3> pose = boost::shared_ptr<const POSE3>())
    {
      SPATIAL_AB_INERTIA I(pose);

      #ifndef NEXCEPT
      if (m.rows() != 6 || m.columns() != 6)
        throw MissizeException();
      #endif
      m.get_sub_mat(3, 6, 3, 6, I.H);
      m.get_sub_mat(0, 3, 0, 3, I.M);   // note: using M=H' here temporarily
      I.H += MATRIX3::transpose(I.M);   //       to compute mean of the two
      I.H *= (REAL) 0.5;                //       matrices
      m.get_sub_mat(0, 3, 3, 6, I.M);
      m.get_sub_mat(3, 6, 0, 3, I.J);
      return I;
    }

    template <class Mat>
    static SPATIAL_AB_INERTIA from_matrix(const Mat& m, boost::shared_ptr<POSE3> pose = boost::shared_ptr<POSE3>())
    {
      return from_matrix(m, boost::const_pointer_cast<const POSE3>(pose));
    }

    template <class Mat>
    Mat& to_matrix(Mat& m) const
    {
      m.resize(6,6);
      m.set_sub_mat(0, 0, MATRIX3::transpose(H));
      m.set_sub_mat(0, 3, M);
      m.set_sub_mat(3, 3, H);
      m.set_sub_mat(3, 0, J);
      return m;
    }

    /// The lower right hand matrix 'H' (and transpose of upper left matrix)
    MATRIX3 H;

    /// The lower left matrix
    MATRIX3 J;

    /// The upper right matrix
    MATRIX3 M;

    /// The pose that this inertia is defined in
    boost::shared_ptr<const POSE3> pose;

  private:
    void mult_spatial(const SVECTOR6& v, SVECTOR6& result) const;
    void inverse_mult_spatial(const SFORCE& w, SVECTOR6& result) const;
    void inverse_mult_spatial(const SFORCE& w, const MATRIX3& UL, const MATRIX3& UR, const MATRIX3& LL, SVECTOR6& result) const;
}; // end class

std::ostream& operator<<(std::ostream& out, const SPATIAL_AB_INERTIA& m);

