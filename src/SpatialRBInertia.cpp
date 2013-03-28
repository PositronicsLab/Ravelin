/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Default constructor -- constructs a zero inertia matrix
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA()
{
  m = (REAL) 0.0;
  h.set_zero();
  J.set_zero();
}

/// Constructs the 6x6 spatial matrix from the given values 
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J)
{
  this->m = m;
  this->h = h;
  this->J = J;
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

/// Multiplies the inverse of this spatial matrix by a wrench 
TWIST SPATIAL_RB_INERTIA::inverse_mult(const WRENCH& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // compute the components of the matrix
  MATRIX3 UR = MATRIX3::inverse((hx * hx * inv_m) + J);
  MATRIX3 UL = UR * hx * -inv_m;
  MATRIX3 LL = ((hx * UL) - MATRIX3::identity()) * inv_m;

  // get the components of the wrench 
  VECTOR3 top = w.get_force();
  VECTOR3 bot = w.get_torque();

  // do the arithmetic
  TWIST t(UL*top + UR*bot, LL*top + UL.transpose_mult(bot));
  t.pose = pose;
  return t;
}

/// Multiplies the inverse of this spatial matrix by a wrench 
std::vector<TWIST>& SPATIAL_RB_INERTIA::inverse_mult(const std::vector<WRENCH>& w, std::vector<TWIST>& result) const
{
  result.resize(w.size());
  if (result.empty())
    return result;

  // compute the skew symmetric version of h
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute inverse mass
  REAL inv_m = (REAL) 1.0/m;

  // compute the components of the matrix
  MATRIX3 UR = MATRIX3::inverse((hx * hx * inv_m) + J);
  MATRIX3 UL = UR * hx * -inv_m;
  MATRIX3 LL = ((hx * UL) - MATRIX3::identity()) * inv_m;

  // get the components of the wrench 
  for (unsigned i=0; i< w.size(); i++)
  {
    #ifndef NEXCEPT
    if (pose != w[i].pose)
      throw FrameException();
    #endif

    VECTOR3 top = w[i].get_force();
    VECTOR3 bot = w[i].get_torque();

    // do the arithmetic
    result[i] = TWIST(UL*top + UR*bot, LL*top + UL.transpose_mult(bot));
    result[i].pose = pose;
  }

  return result;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
WRENCH SPATIAL_RB_INERTIA::operator*(const TWIST& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  // get necessary components of t
  VECTOR3 ttop = t.get_angular();
  VECTOR3 tbot = t.get_linear();

  // do some precomputation
  MATRIX3 hX = MATRIX3::skew_symmetric(h);

  // compute top part of result
  VECTOR3 rtop = (tbot * m) - (hX * ttop);
  VECTOR3 rbot = (J * ttop) + (hX * tbot);

  WRENCH w(rtop, rbot); 
  w.pose = pose;
  return w;
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

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
std::vector<WRENCH>& SPATIAL_RB_INERTIA::mult(const std::vector<TWIST>& t, std::vector<WRENCH>& result) const
{
  // get number of twists 
  const unsigned N = t.size(); 

  // resize the result
  result.resize(N);

  // look for empty result
  if (N == 0)
    return result;

  // compute the skew symmetric matrix corresponding to h
  MATRIX3 hX = MATRIX3::skew_symmetric(this->h);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< N; i++)
  {
    const TWIST& twist = t[i];

    #ifndef NEXCEPT
    if (pose != twist.pose)
      throw FrameException();
    #endif

    VECTOR3 top = twist.get_angular();
    VECTOR3 bot = twist.get_linear();
    WRENCH w((bot * this->m)-(hX * top), (this->J * top)+(hX * bot));
    result[i] = w;
    result[i].pose = pose;
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
