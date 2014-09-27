/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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

/// Default constructor -- constructs a zero inertia matrix
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(shared_ptr<POSE3> pose)
{
  m = (REAL) 0.0;
  h.set_zero();
  J.set_zero();
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Constructs the 6x6 spatial matrix from the given values 
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J, shared_ptr<const POSE3> pose)
{
  this->m = m;
  this->h = h;
  this->J = J;
  this->pose = pose;
}

/// Constructs the 6x6 spatial matrix from the given values 
SPATIAL_RB_INERTIA::SPATIAL_RB_INERTIA(REAL m, const VECTOR3& h, const MATRIX3& J, shared_ptr<POSE3> pose)
{
  this->m = m;
  this->h = h;
  this->J = J;
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Copies a spatial matrix to this one 
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator=(const SPATIAL_RB_INERTIA& m)
{
  this->pose = m.pose;
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

/// Multiplies a spatial vector
void SPATIAL_RB_INERTIA::mult_spatial(const SVECTOR6& t, SVECTOR6& result) const
{
  // get necessary components of t
  ORIGIN3 ttop(t.get_upper());
  ORIGIN3 tbot(t.get_lower());

  // do some precomputation
  MATRIX3 hx = MATRIX3::skew_symmetric(h);
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);

  // compute result
  VECTOR3 rtop(tbot*m + hxm.transpose_mult(ttop), pose);
  VECTOR3 rbot((J-hx*hxm)*ttop + hxm*tbot, pose);
  result.pose = pose;
  result.set_upper(rtop);
  result.set_lower(rbot);  
}

/// Multiplies a spatial vector
void SPATIAL_RB_INERTIA::mult_spatial(const SVECTOR6& t, const MATRIX3& hxm, const MATRIX3& Jmod, SVECTOR6& result) const
{
  // get necessary components of t
  ORIGIN3 ttop(t.get_upper());
  ORIGIN3 tbot(t.get_lower());

  // compute result
  VECTOR3 rtop(tbot*m + hxm.transpose_mult(ttop), pose);
  VECTOR3 rbot(Jmod*ttop + hxm*tbot, pose);
  result.pose = pose;
  result.set_upper(rtop);
  result.set_lower(rbot);  
}

/// Multiplies the inverse of this inertia by a spatial vector
void SPATIAL_RB_INERTIA::inverse_mult_spatial(const SVECTOR6& w, SVECTOR6& result) const
{
  // precompute some things
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute the inverse of the c.o.m. inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute h * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // setup the result
  inverse_mult_spatial(w, iJ, hx, hxiJ, 1.0/m, result); 
}

void SPATIAL_RB_INERTIA::inverse_mult_spatial(const SVECTOR6& w, const MATRIX3& iJ, const MATRIX3& hx, const MATRIX3& hxiJ, REAL inv_m, SVECTOR6& result) const
{
  // get the components of the force 
  ORIGIN3 top(w.get_upper());
  ORIGIN3 bot(w.get_lower());

  // do the arithmetic
  VECTOR3 ttop(hxiJ.transpose_mult(top) + iJ*bot, pose); 
  VECTOR3 tbot(top*inv_m + hx*(hxiJ.transpose_mult(top)) + hxiJ.mult(bot), pose);
  result.pose = pose;
  result.set_upper(ttop);
  result.set_lower(tbot);
}

/// Multiplies the inverse of this spatial matrix by a force 
SACCEL SPATIAL_RB_INERTIA::inverse_mult(const SFORCE& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  SACCEL a;
  inverse_mult_spatial(w, a);
  return a;
}

/// Multiplies the inverse of this spatial matrix by a momentum 
SVELOCITY SPATIAL_RB_INERTIA::inverse_mult(const SMOMENTUM& w) const
{
  #ifndef NEXCEPT
  if (pose != w.pose)
    throw FrameException();
  #endif

  SVELOCITY v;
  inverse_mult_spatial(w, v);

  return v;
}

/// Multiplies the inverse of this spatial matrix by a force 
std::vector<SACCEL>& SPATIAL_RB_INERTIA::inverse_mult(const std::vector<SFORCE>& w, std::vector<SACCEL>& result) const
{
  result.resize(w.size());
  if (result.empty())
    return result;

  // precompute some things
  MATRIX3 hx = MATRIX3::skew_symmetric(h);

  // compute the inverse of the c.o.m. inertia matrix
  MATRIX3 iJ = MATRIX3::invert(J);

  // compute h * inv(J)
  MATRIX3 hxiJ = hx * iJ;

  // compute 1/m
  REAL inv_m = ((REAL) 1.0)/m;

  // get the components of the force 
  for (unsigned i=0; i< w.size(); i++)
  {
    #ifndef NEXCEPT
    if (pose != w[i].pose)
      throw FrameException();
    #endif

    // do the spatial arithmetic
    inverse_mult_spatial(w[i], iJ, hx, hxiJ, inv_m, result[i]);
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

  // do the spatial arithmetic
  SFORCE f;
  mult_spatial(t, f);
  return f;
}

/// Multiplies this inertia by a velocity and returns a momentum 
SMOMENTUM SPATIAL_RB_INERTIA::operator*(const SVELOCITY& t) const
{
  #ifndef NEXCEPT
  if (pose != t.pose)
    throw FrameException();
  #endif

  // do the spatial arithmetic
  SMOMENTUM f;
  mult_spatial(t, f);
  return f;
}

//// Multiplies this matrix by a scalar in place
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
SPATIAL_RB_INERTIA SPATIAL_RB_INERTIA::operator+(const SPATIAL_RB_INERTIA& Jx) const
{
  #ifndef NEXCEPT
  if (pose != Jx.pose)
    throw FrameException();
  #endif

  SPATIAL_RB_INERTIA result;
  result.m = m + Jx.m;
  if (result.m > 0.0)
    result.h = (h*m + Jx.h*Jx.m)/(m + Jx.m);
  else
    result.h.set_zero();
  MATRIX3 hrskew = MATRIX3::skew_symmetric(result.h);
  MATRIX3 mrhrskew = MATRIX3::skew_symmetric(result.h * result.m);
  MATRIX3 hskew = MATRIX3::skew_symmetric(h);
  MATRIX3 mhskew = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxskew = MATRIX3::skew_symmetric(Jx.h);
  MATRIX3 mxhxskew = MATRIX3::skew_symmetric(Jx.h*Jx.m);
  result.J = J + Jx.J - (mhskew*hskew + hxskew*mxhxskew - hrskew*mrhrskew);
  result.pose = pose;
  return result;
}

/// Subtracts two spatial matrices
SPATIAL_RB_INERTIA SPATIAL_RB_INERTIA::operator-(const SPATIAL_RB_INERTIA& Jx) const
{
  #ifndef NEXCEPT
  if (pose != Jx.pose)
    throw FrameException();
  #endif

  SPATIAL_RB_INERTIA result;
  result.m = m - Jx.m;
  if (result.m > 0.0)
    result.h = (h*m - Jx.h*Jx.m)/(m - Jx.m);
  else
    result.h.set_zero();
  MATRIX3 hrskew = MATRIX3::skew_symmetric(result.h);
  MATRIX3 mrhrskew = MATRIX3::skew_symmetric(result.h * result.m);
  MATRIX3 hskew = MATRIX3::skew_symmetric(h);
  MATRIX3 mhskew = MATRIX3::skew_symmetric(h*m);
  MATRIX3 hxskew = MATRIX3::skew_symmetric(Jx.h);
  MATRIX3 mxhxskew = MATRIX3::skew_symmetric(Jx.h*Jx.m);
  result.J = J - Jx.J - mhskew*hskew + hxskew*mxhxskew + hrskew*mrhrskew;
  result.pose = pose;
  return result;
}

/// Adds m to this in place
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator+=(const SPATIAL_RB_INERTIA& Jx)
{
  #ifndef NEXCEPT
  if (pose != Jx.pose)
    throw FrameException();
  #endif

  *this = *this + Jx;
  return *this;
}

/// Subtracts m from this in place
SPATIAL_RB_INERTIA& SPATIAL_RB_INERTIA::operator-=(const SPATIAL_RB_INERTIA& Jx)
{
  #ifndef NEXCEPT
  if (pose != Jx.pose)
    throw FrameException();
  #endif

  *this = *this - Jx;
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
  MATRIX3 hx = MATRIX3::skew_symmetric(h);
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 Jstar = J - (hx*hxm);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< N; i++)
  {
    #ifndef NEXCEPT
    if (pose != t[i].pose)
      throw FrameException();
    #endif

    // compute result
    mult_spatial(t[i], hx, Jstar, result[i]);
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
  MATRIX3 hx = MATRIX3::skew_symmetric(h);
  MATRIX3 hxm = MATRIX3::skew_symmetric(h*m);
  MATRIX3 Jstar = J - (hx*hxm);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< N; i++)
  {
    #ifndef NEXCEPT
    if (pose != t[i].pose)
      throw FrameException();
    #endif

    // compute result
    mult_spatial(t[i], hxm, Jstar, result[i]);
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

