/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Constants.h>
#include <Moby/SMatrix6N.h>
#include <Moby/SpatialABInertiaf.h>

using namespace Moby;

/// Default constructor -- constructs a zero inertia matrix
SpatialABInertiaf::SpatialABInertiaf()
{
  M.set_zero();
  H.set_zero();
  J.set_zero();
}

/// Constructs a spatial AB inertia from the given MatrixNf object
SpatialABInertiaf::SpatialABInertiaf(const MatrixNf& m)
{
  m.get_sub_mat(0, 3, 3, 6, M);
  m.get_sub_mat(3, 6, 3, 6, H);
  m.get_sub_mat(3, 6, 0, 3, J);
}

/// Constructs the spatial AB inertia from the given values 
SpatialABInertiaf::SpatialABInertiaf(const Matrix3f& M, const Matrix3f& H, const Matrix3f& J)
{
  this->M = M;
  this->H = H;
  this->J = J;
}

/// Copies a spatial AB inertia to this one 
SpatialABInertiaf& SpatialABInertiaf::operator=(const SpatialABInertiaf& m)
{
  this->M = m.M;
  this->H = m.H;
  this->J = m.J;
  return *this;
}

/// Copies a spatial RB inertia to this one 
SpatialABInertiaf& SpatialABInertiaf::operator=(const SpatialRBInertia& m)
{
  this->M.set_identity() *=m.m;
  this->H = Matrix3f::skew_symmetric(m.h);
  this->J = m.J;
  return *this;
}

/// Creates a zero matrix
void SpatialABInertiaf::set_zero()
{
  M = H = J.set_zero();
}

/// Multiplies this matrix by a spatial twist and returns the result in a new vector
Wrenchf SpatialABInertiaf::operator*(const Twistf& t) const
{
  // get necessary components of v
  Vector3f ang = t.get_angular();
  Vector3f lin = t.get_linear();

  // do some precomputation
  Matrix3f HT = Matrix3f::transpose(H);

  // compute top part of result
  Vector3f force = HT * ang+ (M * lin);
  Vector3f torque = (J * ang) + (H * lin);

  return Wrenchf(force, torque); 
}

/// Multiplies this matrix by a scalar in place
SpatialABInertiaf& SpatialABInertiaf::operator*=(Real scalar)
{
  M *= scalar;
  H *= scalar;
  J *= scalar;
  return *this;
}

/// Returns the negation of this matrix
SpatialABInertiaf SpatialABInertiaf::operator-() const
{
  SpatialABInertiaf result;
  result.M = -this->M;
  result.H = -this->H;
  result.J = -this->J;
  return result;
}

/// Adds a spatial articulated body inertia and a spatial rigid body inertia 
SpatialABInertiaf SpatialABInertiaf::operator+(const SpatialRBInertia& m) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // do some preliminary calculations
  SpatialABInertiaf result;
  result.M = M;
  result.H = H;
  result.J = m.J + J;

  // update M with mass
  result.M(X,X) += m.m;
  result.M(Y,Y) += m.m;
  result.M(Z,Z) += m.m;

  // update H
  Matrix3f hx = Matrix3f::skew_symmetric(m.h);
  result.H += hx;

  return result;
}

/// Adds two spatial matrices
SpatialABInertiaf SpatialABInertiaf::operator+(const SpatialABInertiaf& m) const
{
  SpatialABInertiaf result;
  result.M = this->M + m.M;
  result.H = this->H + m.H;
  result.J = this->J + m.J;
  return result;
}

/// Subtracts two spatial matrices
SpatialABInertiaf SpatialABInertiaf::operator-(const SpatialABInertiaf& m) const
{
  SpatialABInertiaf result;
  result.M = this->M - m.M;
  result.H = this->H - m.H;
  result.J = this->J - m.J;
  return result;
}

/// Adds m to this in place
SpatialABInertiaf& SpatialABInertiaf::operator+=(const SpatialABInertiaf& m)
{
  this->M += m.M;
  this->H += m.H;
  this->J += m.J;
  return *this;
}

/// Subtracts m from this in place
SpatialABInertiaf& SpatialABInertiaf::operator-=(const SpatialABInertiaf& m)
{
  this->M -= m.M;
  this->H -= m.H;
  this->J -= m.J;
  return *this;
}

/// Multiplies the inverse of this spatial AB inertia by a wrench 
Twistf SpatialABInertiaf::inverse_mult(const Wrenchf& w) const
{
  Matrix3f nMinv = -Matrix3f::inverse(M);
  Matrix3f HT = Matrix3f::transpose(H);
  Matrix3f UR = Matrix3f::inverse((H * nMinv * HT) + J);
  Matrix3f UL = UR * H * nMinv;
  Matrix3f LR = Matrix3f::transpose(UL);
  Matrix3f LL = nMinv * ((HT * UL) - IDENTITY_3x3);

  // get the components of w 
  Vector3f force = w.get_force();
  Vector3f torque = w.get_torque();

  // do the arithmetic
  return Twistf(UL*force + UR*torque, LL*force + LR*torque);
}

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
SMatrix6N& SpatialABInertiaf::mult(const SMatrix6N& m, SMatrix6N& result) const
{
  const unsigned SPATIAL_DIM = 6;

  // get number of columns of m
  const unsigned NCOLS = m.columns();

  // resize the result
  result.resize(SPATIAL_DIM, NCOLS);

  // look for empty result
  if (NCOLS == 0)
  {
    result.set_zero();
    return result;
  }

  // compute the transpose of H
  Matrix3f HT = Matrix3f::transpose(this->H);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< NCOLS; i++)
  {
    SVector6 v = m.get_column(i);
    Vector3f vtop = v.get_upper();
    Vector3f vbot = v.get_lower();
    v = SVector6((this->M * vbot)+(HT * vtop), (this->J * vtop)+(H * vbot));
    result.set_column(i, v);
  } 

  return result;
}

/// Multiplies a 6x6 matrix by a Spatial AB inertia matrix
SpatialABInertiaf SpatialABInertiaf::mult(const MatrixNf& m, const SpatialABInertiaf& I)
{
  assert(m.rows() == 6 && m.columns() == 6);

  // get the components of m
  Matrix3f UL, UR, LL, LR;
  m.get_sub_mat(0,3,0,3,UL);
  m.get_sub_mat(0,3,3,6,UR);
  m.get_sub_mat(3,6,0,3,LL);
  m.get_sub_mat(3,6,3,6,LR);
 
  // multiply by components of I
  SpatialABInertiaf result;
  result.M = (UL*I.M) + (UR*I.H);
  result.J = (LL*Matrix3f::transpose(I.H)) + (LR*I.J);
  result.H = (LL*I.M) + (LR*I.H);
  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const SpatialABInertiaf& m) 
{
  out << "spatial AB H:" << std::endl << m.H;
  out << "spatial AB M:" << std::endl << m.M;
  out << "spatial AB J:" << std::endl << m.J;
   
  return out;
}

