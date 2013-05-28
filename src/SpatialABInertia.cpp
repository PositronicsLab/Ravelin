/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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

/// Constructs the spatial AB inertia from the given values 
SPATIAL_AB_INERTIA::SPATIAL_AB_INERTIA(const MATRIX3& M, const MATRIX3& H, const MATRIX3& J, shared_ptr<const POSE3> pose)
{
  this->M = M;
  this->H = H;
  this->J = J;
  this->pose = pose;
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
  this->pose = m.pose;
  this->M.set_identity() *= m.m;
  this->H = MATRIX3::skew_symmetric(m.h);
  this->J = m.J;
  return *this;
}

/// Creates a zero matrix
void SPATIAL_AB_INERTIA::set_zero()
{
  M.set_zero();
  H.set_zero();
  J.set_zero();
}

/// Multiplies this matrix by a twist and returns the result in a wrench 
WRENCH SPATIAL_AB_INERTIA::mult(const TWIST& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  // get necessary components of the twist 
  ORIGIN3 top(t.get_angular());
  ORIGIN3 bot(t.get_linear());

  // compute top part of result
  VECTOR3 rtop(H.transpose_mult(top) + (M * bot), pose);
  VECTOR3 rbot((J * top) + (H * bot), pose);

  return WRENCH(rtop, rbot, pose); 
}

/// Multiplies this matrix by a vector of twists and returns the result in a vector of wrenches 
vector<WRENCH>& SPATIAL_AB_INERTIA::mult(const vector<TWIST>& t, vector<WRENCH>& result) const
{
  result.resize(t.size());

  // get necessary components of the twist 
  for (unsigned i=0; i< t.size(); i++)
  { 
    #ifndef NEXCEPT
    if (pose != t[i].pose)
      throw FrameException();
    #endif

    ORIGIN3 top(t[i].get_angular());
    ORIGIN3 bot(t[i].get_linear());
    VECTOR3 wtop(H.transpose_mult(top) + (M * bot), pose);
    VECTOR3 wbot((J * top) + (H * bot), pose);
    result[i] = WRENCH(wtop, wbot, pose); 
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

  // do some preliminary calculations
  SPATIAL_AB_INERTIA result(pose);
  result.M = M;
  result.H = H;
  result.J = m.J + J;

  // update M with mass
  result.M(X,X) += m.m;
  result.M(Y,Y) += m.m;
  result.M(Z,Z) += m.m;

  // update H
  MATRIX3 hx = MATRIX3::skew_symmetric(m.h);
  result.H += hx;

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

/// Multiplies the inverse of this spatial AB inertia by a wrench to yield a twist
TWIST SPATIAL_AB_INERTIA::inverse_mult(const WRENCH& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  MATRIX3 nMinv = -MATRIX3::invert(M);
  MATRIX3 UR = MATRIX3::invert((H * nMinv.mult_transpose(H)) + J);
  MATRIX3 UL = UR * H * nMinv;
  MATRIX3 LL = nMinv * (H.transpose_mult(UL) - MATRIX3::identity());

  // get the components of the wrench 
  ORIGIN3 top(w.get_force());
  ORIGIN3 bot(w.get_torque());
  VECTOR3 ttop(UL*top + UR*bot, pose);
  VECTOR3 tbot(LL*top + UL.transpose_mult(bot), pose); 

  // do the arithmetic
  return TWIST(ttop, tbot, pose);
}

/// Multiplies the inverse of this spatial AB inertia by a wrench to yield a twist
vector<TWIST>& SPATIAL_AB_INERTIA::inverse_mult(const std::vector<WRENCH>& w, vector<TWIST>& result) const
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

    // get the components of the wrench 
    ORIGIN3 top(w[i].get_force());
    ORIGIN3 bot(w[i].get_torque());
    VECTOR3 ttop(UL*top + UR*bot, pose);
    VECTOR3 tbot(LL*top + UL.transpose_mult(bot), pose);

    // do the arithmetic
    result[i] = TWIST(ttop, tbot, pose);
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

