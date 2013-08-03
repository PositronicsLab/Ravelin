/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

using boost::shared_ptr;

/// Default constructor -- constructs a zero inertia matrix
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(shared_ptr<const POSE3> pose)
{
  m = (REAL) 0.0;
  h.set_zero();
  J.set_zero();
  this->pose = pose;
}

/// Constructs the 6x6 spatial matrix from the given values 
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J, shared_ptr<const POSE3> pose)
{
  this->m = m;
  this->h = h;
  this->J = J;
  this->pose = pose;
}

/// Copies a spatial matrix to this one 
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator=(const SPATIAL_RB_INERTIA& m)
{
  this->m = m.m;
  this->h = m.h;
  this->J = m.J;
  return *this;
}

/// Creates a zero matrix
void SPATIAL_RB_INERTIA::set_zero()
{
  m = (REAL) 0.0;
  h.set_zero();
  J.set_zero();
}

/// Multiplies the inverse of this spatial matrix by a force 
SACCEL SPATIAL_RB_INERTIA::inverse_mult(const SFORCE& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  // compute the inverse of the inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute hx * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // get the components of the force 
  ORIGIN3 top(w.get_force());
  ORIGIN3 bot(w.get_torque());

  // do the arithmetic
  VECTOR3 ttop(hxiJ.transpose_mult(top) + iJ*bot, pose); 
  VECTOR3 tbot(hxiJ*hx.transpose_mult(top) + top*m + hxiJ*bot, pose); 

  return SACCEL(ttop, tbot, pose);
}

/// Multiplies the inverse of this spatial matrix by a momentum 
SVELOCITY SPATIAL_RB_INERTIA::inverse_mult(const SMOMENTUM& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  // compute the inverse of the inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute hx * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // get the components of the force 
  ORIGIN3 top(w.get_linear());
  ORIGIN3 bot(w.get_angular());

  // do the arithmetic
  VECTOR3 ttop(hxiJ.transpose_mult(top) + iJ*bot, pose); 
  VECTOR3 tbot(hxiJ*hx.transpose_mult(top) + top*m + hxiJ*bot, pose); 

  return SVELOCITY(ttop, tbot, pose);
}

/// Multiplies the inverse of this spatial matrix by a force 
std::vector<SACCEL>& SPATIAL_RB_INERTIA::inverse_mult(const std::vector<SFORCE>& w, std::vector<SACCEL>& result) const
{
  result.resize(w.size());
  if (result.empty())
    return result;

  // compute the inverse of the inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute hx * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // get the components of the force 
  for (unsigned i=0; i< w.size(); i++)
  {
    #ifndef NEXCEPT
    if (pose != w[i].pose)
      throw FrameException();
    #endif

    // get the components of the force 
    ORIGIN3 top(w[i].get_force());
    ORIGIN3 bot(w[i].get_torque());

    // do the arithmetic
    VECTOR3 ttop(hxiJ.transpose_mult(top) + iJ*bot, pose); 
    VECTOR3 tbot(hxiJ*hx.transpose_mult(top) + top*m + hxiJ*bot, pose); 
    result[i] = SACCEL(ttop, tbot, pose);
  }

  return result;
}

/// Multiplies the inverse of this spatial matrix by a momentum 
std::vector<SVELOCITY>& SPATIAL_RB_INERTIA::inverse_mult(const std::vector<SMOMENTUM>& w, std::vector<SVELOCITY>& result) const
{
  result.resize(w.size());
  if (result.empty())
    return result;

  // compute the inverse of the inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute hx * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // get the components of the force 
  for (unsigned i=0; i< w.size(); i++)
  {
    #ifndef NEXCEPT
    if (pose != w[i].pose)
      throw FrameException();
    #endif

    // get the components of the momentum 
    ORIGIN3 top(w[i].get_linear());
    ORIGIN3 bot(w[i].get_angular());

    // do the arithmetic
    VECTOR3 ttop(hxiJ.transpose_mult(top) + iJ*bot, pose); 
    VECTOR3 tbot(hxiJ*hx.transpose_mult(top) + top*m + hxiJ*bot, pose); 
    result[i] = SVELOCITY(ttop, tbot, pose);
  }

  return result;
}

/// Multiplies this inertia by an acceleration and returns a force 
SFORCE SPATIAL_RB_INERTIA::operator*(const SACCEL& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  // get necessary components of t
  ORIGIN3 ttop(t.get_angular());
  ORIGIN3 tbot(t.get_linear());

  // do some precomputation
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(h) * hxm; 

  // compute result
  VECTOR3 rtop((tbot * m) - (hxm * ttop), pose);
  VECTOR3 rbot(((J - hxhxm) * ttop) + (hxm * tbot), pose);
  return SFORCE(rtop, rbot, pose); 
}

/// Multiplies this inertia by a velocity and returns a momentum 
SMOMENTUM SPATIAL_RB_INERTIA::operator*(const SVELOCITY& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  // get necessary components of t
  ORIGIN3 ttop(t.get_angular());
  ORIGIN3 tbot(t.get_linear());

  // do some precomputation
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(h) * hxm; 

  // compute result
  VECTOR3 rtop((tbot * m) - (hxm * ttop), pose);
  VECTOR3 rbot(((J - hxhxm) * ttop) + (hxm * tbot), pose);
  return SMOMENTUM(rtop, rbot, pose); 
}

/// Multiplies this matrix by a scalar in place
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator*=(REAL scalar)
{
  m *= scalar;
  h *= scalar;
  J *= scalar;
  return *this;
}

/// Returns the negation of this matrix
SPATIAL_RB_INERTIA SPATIAL_RB_INERTIA::operator-() const
{
  SPATIAL_RB_INERTIA result;
  result.m = -this->m;
  result.h = -this->h;
  result.J = -this->J;
  result.pose = pose;
  return result;
}

/// Adds two spatial matrices
SPATIAL_RB_INERTIA SPATIAL_RB_INERTIA::operator+(const SPATIAL_RB_INERTIA& m) const
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  SPATIAL_RB_INERTIA result;
  result.m = this->m + m.m;
  result.h = this->h + m.h;
  result.J = this->J + m.J;
  return result;
}

/// Subtracts two spatial matrices
SPATIAL_RB_INERTIA SPATIAL_RB_INERTIA::operator-(const SPATIAL_RB_INERTIA& m) const
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  SPATIAL_RB_INERTIA result;
  result.m = this->m - m.m;
  result.h = this->h - m.h;
  result.J = this->J - m.J;
  result.pose = pose;
  return result;
}

/// Adds m to this in place
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator+=(const SPATIAL_RB_INERTIA& m)
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  this->m += m.m;
  this->h += m.h;
  this->J += m.J;
  return *this;
}

/// Subtracts m from this in place
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator-=(const SPATIAL_RB_INERTIA& m)
{
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  this->m -= m.m;
  this->h -= m.h;
  this->J -= m.J;
  return *this;
}

/// Multiplies this inertia by a vector of accelerations and returns a vector of forces
std::vector<SFORCE>& SPATIAL_RB_INERTIA::mult(const std::vector<SACCEL>& t, std::vector<SFORCE>& result) const
{
  // get number of accels 
  const unsigned N = t.size(); 

  // resize the result
  result.resize(N);

  // look for empty result
  if (N == 0)
    return result;

  // do some precomputation
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(h) * hxm; 

  // carry out multiplication one column at a time
  for (unsigned i=0; i< N; i++)
  {
    const SACCEL& accel = t[i];

    #ifndef NEXCEPT
    if (pose != accel.pose)
      throw FrameException();
    #endif

    // compute result
    ORIGIN3 ttop(accel.get_angular());
    ORIGIN3 tbot(accel.get_linear());
    VECTOR3 rtop((tbot * m) - (hxm * ttop), pose);
    VECTOR3 rbot(((J - hxhxm) * ttop) + (hxm * tbot), pose);
    result[i] = SFORCE(rtop, rbot, pose);
  } 

  return result;
}

/// Multiplies this inertia by a vector of accelerations and returns a vector of forces
std::vector<SMOMENTUM>& SPATIAL_RB_INERTIA::mult(const std::vector<SVELOCITY>& t, std::vector<SMOMENTUM>& result) const
{
  // get number of accels 
  const unsigned N = t.size(); 

  // resize the result
  result.resize(N);

  // look for empty result
  if (N == 0)
    return result;

  // do some precomputation
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(h) * hxm; 

  // carry out multiplication one column at a time
  for (unsigned i=0; i< N; i++)
  {
    const SVELOCITY& vel = t[i];

    #ifndef NEXCEPT
    if (pose != vel.pose)
      throw FrameException();
    #endif

    // compute result
    ORIGIN3 ttop(vel.get_angular());
    ORIGIN3 tbot(vel.get_linear());
    VECTOR3 rtop((tbot * m) - (hxm * ttop), pose);
    VECTOR3 rbot(((J - hxhxm) * ttop) + (hxm * tbot), pose);
    result[i] = SMOMENTUM(rtop, rbot, pose);
  } 

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const SPATIAL_RB_INERTIA& m) 
{
  out << "spatial rigid body mass=" << m.m << " h = " << m.h << " J = ";
  out << std::endl << m.J;
  out << " pose: " << m.pose << std::endl;
   
  return out;
}

