/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using std::vector;
using boost::shared_ptr;

/// Default constructor -- constructs a zero inertia matrix
SPATIAL_AB_INERTIA::SPATIAL_AB_INERTIA(shared_ptr<const POSE3> pose)
{
  M.set_zero();
  H.set_zero();
  J.set_zero();
  this->pose = pose;
}

/// Default constructor -- constructs a zero inertia matrix
SPATIAL_AB_INERTIA::SPATIAL_AB_INERTIA(shared_ptr<POSE3> pose)
{
  M.set_zero();
  H.set_zero();
  J.set_zero();
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Constructs the spatial AB inertia from the given values 
SPATIAL_AB_INERTIA::SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J, shared_ptr<const POSE3> pose)
{
  this->M = M;
  this->H = H;
  this->J = J;
  this->pose = pose;
}

/// Constructs the spatial AB inertia from the given values 
SPATIAL_AB_INERTIA::SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J, shared_ptr<POSE3> pose)
{
  this->M = M;
  this->H = H;
  this->J = J;
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Copies a spatial AB inertia to this one 
SPATIAL_AB_INERTIA& SPATIAL_AB_INERTIA::operator=(const SPATIAL_AB_INERTIA& m)
{
  this->pose = m.pose;
  this->M = m.M;
  this->H = m.H;
  this->J = m.J;
  return *this;
}

/// Copies a spatial RB inertia to this one 
SPATIAL_AB_INERTIA& SPATIAL_AB_INERTIA::operator=(const SPATIAL_RB_INERTIA& m)
{
  // precompute some things
  MATRIX3 hx = MATRIX3::skew_symmetric(m.h);
  MATRIX3 mhx = MATRIX3::skew_symmetric(m.m * m.h);

  this->pose = m.pose;
  this->M.set_identity() *= m.m;
  this->H = mhx;
  this->J = m.J - mhx*hx;
  return *this;
}

/// Converts this matrix to a spatial RB inertia
SPATIAL_RB_INERTIA SPATIAL_AB_INERTIA::to_rb_inertia() const
{
  SPATIAL_RB_INERTIA Jx;
  Jx.pose = pose;
  Jx.m = M(0,0);
  Jx.h = MATRIX3::inverse_skew_symmetric(H/Jx.m);
  Jx.J = J + H*MATRIX3::skew_symmetric(Jx.h);
  return Jx;
}

/// Creates a zero matrix
void SPATIAL_AB_INERTIA::set_zero()
{
  M.set_zero();
  H.set_zero();
  J.set_zero();
}

/// Does spatial arithmetic
void SPATIAL_AB_INERTIA::mult_spatial(const SVECTOR6& t, SVECTOR6& result) const
{
  // get necessary components of the acceleration 
  ORIGIN3 top(t.get_upper());
  ORIGIN3 bot(t.get_lower());

  // compute top part of result
  result.pose = pose;
  result.set_upper(VECTOR3(H.transpose_mult(top) + (M * bot), pose));
  result.set_lower(VECTOR3((J * top) + (H * bot), pose));
}

/// Does inverse spatial matrix/vector multiplication
void SPATIAL_AB_INERTIA::inverse_mult_spatial(const SFORCE& w, SVECTOR6& result) const
{
  MATRIX3 nMinv = -MATRIX3::invert(M);
  MATRIX3 UR = MATRIX3::invert((H * nMinv.mult_transpose(H)) + J);
  MATRIX3 UL = UR * H * nMinv;
  MATRIX3 LL = nMinv * (H.transpose_mult(UL) - MATRIX3::identity());

  // get the components of the force 
  ORIGIN3 top(w.get_upper());
  ORIGIN3 bot(w.get_lower());
  result.pose = pose;
  result.set_upper(VECTOR3(UL*top + UR*bot, pose));
  result.set_lower(VECTOR3(LL*top + UL.transpose_mult(bot), pose)); 
}

/// Multiplies this matrix by an axis and returns the result in a momentum 
/// Does inverse spatial matrix/vector multiplication
void SPATIAL_AB_INERTIA::inverse_mult_spatial(const SFORCE& w, const MATRIX3& UL, const MATRIX3& UR, const MATRIX3& LL, SVECTOR6& result) const
{
  // get the components of the force 
  ORIGIN3 top(w.get_upper());
  ORIGIN3 bot(w.get_lower());
  result.pose = pose;
  result.set_upper(VECTOR3(UL*top + UR*bot, pose));
  result.set_lower(VECTOR3(LL*top + UL.transpose_mult(bot), pose)); 
}

/// Multiplies this matrix by an acceleration and returns the result in a force 
SFORCE SPATIAL_AB_INERTIA::mult(const SACCEL& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  SFORCE result;
  mult_spatial(t, result);
  return result;
}

/// Multiplies this matrix by a velocity and returns the result in a momentum 
SMOMENTUM SPATIAL_AB_INERTIA::mult(const SVELOCITY& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  SMOMENTUM result;
  mult_spatial(t, result);
  return result;
}

/// Multiplies this matrix by a vector of accelerations and returns the result in a vector of forces 
vector<SFORCE>& SPATIAL_AB_INERTIA::mult(const vector<SACCEL>& t, vector<SFORCE>& result) const
{
  result.resize(t.size());

  // get necessary components of the acceleration 
  for (unsigned i=0; i< t.size(); i++)
  { 
    #ifndef NEXCEPT
    if (pose != t[i].pose)
      throw FrameException();
    #endif

    mult_spatial(t[i], result[i]);
  }

  return result;
}

/// Multiplies this inertia by a vector of velocities and returns the result in a vector of momenta 
vector<SMOMENTUM>& SPATIAL_AB_INERTIA::mult(const vector<SVELOCITY>& t, vector<SMOMENTUM>& result) const
{
  result.resize(t.size());

  // get necessary components of the velocity 
  for (unsigned i=0; i< t.size(); i++)
  { 
    #ifndef NEXCEPT
    if (pose != t[i].pose)
      throw FrameException();
    #endif

    mult_spatial(t[i], result[i]);
  }

  return result;
}

/// Multiplies this matrix by a scalar in place
SPATIAL_AB_INERTIA& SPATIAL_AB_INERTIA::operator*=(REAL scalar)
{
  M *= scalar;
  H *= scalar;
  J *= scalar;
  return *this;
}

/// Returns the negation of this matrix
SPATIAL_AB_INERTIA SPATIAL_AB_INERTIA::operator-() const
{
  SPATIAL_AB_INERTIA result;
  result.M = -this->M;
  result.H = -this->H;
  result.J = -this->J;
  result.pose = this->pose;
  return result;
}

/// Adds a spatial articulated body inertia and a spatial rigid body inertia 
SPATIAL_AB_INERTIA SPATIAL_AB_INERTIA::operator+(const SPATIAL_RB_INERTIA& m) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  // precompute some things
  MATRIX3 hxm = MATRIX3::skew_symmetric(m.h*m.m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(m.h)*hxm;

  // do some preliminary calculations
  SPATIAL_AB_INERTIA result(pose);
  result.M = M;
  result.H = H + hxm;
  result.J = J + m.J - hxhxm;

  // update M with mass
  result.M(X,X) += m.m;
  result.M(Y,Y) += m.m;
  result.M(Z,Z) += m.m;

  return result;
}

/// Adds two spatial matrices
SPATIAL_AB_INERTIA SPATIAL_AB_INERTIA::operator+(const SPATIAL_AB_INERTIA& m) const
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  SPATIAL_AB_INERTIA result(pose);
  result.M = this->M + m.M;
  result.H = this->H + m.H;
  result.J = this->J + m.J;
  return result;
}

/// Subtracts two spatial matrices
SPATIAL_AB_INERTIA SPATIAL_AB_INERTIA::operator-(const SPATIAL_AB_INERTIA& m) const
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  SPATIAL_AB_INERTIA result(pose);
  result.M = this->M - m.M;
  result.H = this->H - m.H;
  result.J = this->J - m.J;
  return result;
}

/// Adds m to this in place
SPATIAL_AB_INERTIA& SPATIAL_AB_INERTIA::operator+=(const SPATIAL_AB_INERTIA& m)
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  this->M += m.M;
  this->H += m.H;
  this->J += m.J;
  return *this;
}

/// Subtracts m from this in place
SPATIAL_AB_INERTIA& SPATIAL_AB_INERTIA::operator-=(const SPATIAL_AB_INERTIA& m)
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  this->M -= m.M;
  this->H -= m.H;
  this->J -= m.J;
  return *this;
}

/// Multiplies the inverse of this spatial AB inertia by a force to yield an acceleration 
SACCEL SPATIAL_AB_INERTIA::inverse_mult(const SFORCE& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  SACCEL result;
  inverse_mult_spatial(w, result);
  return result;
}

/// Multiplies the inverse of this spatial RB inertia by a momentum to yield a velocity 
SVELOCITY SPATIAL_AB_INERTIA::inverse_mult(const SMOMENTUM& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  MATRIX3 nMinv = -MATRIX3::invert(M);
  MATRIX3 UR = MATRIX3::invert((H * nMinv.mult_transpose(H)) + J);
  MATRIX3 UL = UR * H * nMinv;
  MATRIX3 LL = nMinv * (H.transpose_mult(UL) - MATRIX3::identity());

  // get the components of the momentum
  ORIGIN3 top(w.get_upper());
  ORIGIN3 bot(w.get_lower());

  // result is set in the same order as the force based version
  // (theory indicates this should be the case)
  SVELOCITY result(pose);
  result.set_upper(VECTOR3(UL*top + UR*bot, pose));
  result.set_lower(VECTOR3(LL*top + UL.transpose_mult(bot), pose)); 
  return result;
}

/// Multiplies the inverse of this spatial AB inertia by a force to yield an accel 
vector<SACCEL>& SPATIAL_AB_INERTIA::inverse_mult(const std::vector<SFORCE>& w, vector<SACCEL>& result) const
{
  result.resize(w.size());
  if (result.empty())
    return result;

  // do precomputation
  MATRIX3 nMinv = -MATRIX3::invert(M);
  MATRIX3 UR = MATRIX3::invert((H * nMinv.mult_transpose(H)) + J);
  MATRIX3 UL = UR * H * nMinv;
  MATRIX3 LL = nMinv * (H.transpose_mult(UL) - MATRIX3::identity());

  // loop
  for (unsigned i=0; i< w.size(); i++)
  {
    #ifndef NEXCEPT
    if (pose != w[i].pose)
      throw FrameException();
    #endif

    inverse_mult_spatial(w[i], UL, UR, LL, result[i]); 
  }

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const SPATIAL_AB_INERTIA& m) 
{
  out << "spatial AB H:" << std::endl << m.H;
  out << "spatial AB M:" << std::endl << m.M;
  out << "spatial AB J:" << std::endl << m.J;
  out << "pose: " << m.pose << std::endl;
   
  return out;
}

