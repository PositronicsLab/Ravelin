/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/Opsd.h>

using namespace Ravelin;

Matrix3d to_Matrix3d(const Quatd& q)
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // verify that the quaternion is normalized
  assert(std::fabs(q.magnitude()) - 1 < EPS_FLOAT);

  // setup repeated products
  const double xx = q.x*q.x;
  const double xy = q.x*q.y;
  const double xz = q.x*q.z;
  const double xw = q.x*q.w;
  const double yy = q.y*q.y;
  const double yz = q.y*q.z;
  const double yw = q.y*q.w;
  const double zz = q.z*q.z;
  const double zw = q.z*q.w; 
  const double ww = q.w*q.w;

  Matrix3d m;
  m(X,X) = 2*(xx + ww) - 1;
  m(X,Y) = 2*(xy - zw);
  m(X,Z) = 2*(xz + yw);
  m(Y,X) = 2*(xy + zw);
  m(Y,Y) = 2*(yy + ww) - 1;
  m(Y,Z) = 2*(yz - xw);
  m(Z,X) = 2*(xz - yw);
  m(Z,Y) = 2*(yz + xw);
  m(Z,Z) = 2*(zz + ww) - 1;
  return m;
}

Matrix3d to_Matrix3(const AAngled& a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  double x = a.x;
  double y = a.y;
  double z = a.z;
  double ca = cos(a.angle);
  double sa = sin(a.angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const double SOMEWHAT_EPS_FLOAT = 1e-2;
  double va = (std::fabs(a.angle) > SOMEWHAT_EPS_FLOAT) ? 1 - ca : (sa*sa)/(1+ca);

  // setup the matrix
  Matrix3d m;
  m(X,X) = x*x*va + ca;
  m(X,Y) = x*y*va - z*sa;
  m(X,Z) = x*z*va + y*sa;
  m(Y,X) = x*y*va + z*sa;
  m(Y,Y) = y*y*va + ca;
  m(Y,Z) = y*z*va - x*sa;
  m(Z,X) = x*z*va - y*sa;
  m(Z,Y) = y*z*va + x*sa;
  m(Z,Z) = z*z*va + ca;
  return m;
}

/// Converts a rotation matrix to a unit Quaternion
Quatd to_Quat(const Matrix3d& m)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // core computation
  Quatd q;
  q.w = std::sqrt(std::max(0.0, 1.0 + m(X,X) + m(Y,Y) + m(Z,Z))) * 0.5;
  q.x = std::sqrt(std::max(0.0, 1.0 + m(X,X) - m(Y,Y) - m(Z,Z))) * 0.5;
  q.y = std::sqrt(std::max(0.0, 1.0 - m(X,X) + m(Y,Y) - m(Z,Z))) * 0.5;
  q.z = std::sqrt(std::max(0.0, 1.0 - m(X,X) - m(Y,Y) + m(Z,Z))) * 0.5;

  // sign computation
  if (m(Z,Y) - m(Y,Z) < 0.0)
    q.x = -q.x;
  if (m(X,Z) - m(Z,X) < 0.0)
    q.y = -q.y;
  if (m(Y,X) - m(X,Y) < 0.0)
    q.z = -q.z;

  #ifndef NEXCEPT
  if (!q.unit())
    std::cerr << "Quatd::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}

/// Sets quaternion to that represented by an axis-angle representation
Quatd to_Quat(const AAngled& a)
{
  const double half = a.angle*0.5;
  double sina = std::sin(half);
  Quatd q;
  q.x = a.x * sina;
  q.y = a.y * sina;
  q.z = a.z * sina;
  q.w = std::cos(half);

  #ifndef NEXCEPT
  if (!q.unit())
    std::cerr << "Quatd::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}

/// Sets axis angle from unit quaternion 
/**
 * \todo test this method
 */
AAngled to_AAngle(const Quatd& q)
{
  Quatd qn = Quatd::normalize(q);
  AAngled a;
  a.x = qn.x;
  a.y = qn.y;
  a.z = qn.z;
  a.angle = std::acos(qn.w) * 2.0;

  // verify that [x y z] normalized
  double nrm = AAngled::safe_sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  if (std::fabs(nrm) < EPS_FLOAT)
  {
    // axis is zero; set it to [1 0 0] arbitrarily
    a.x = 1.0;
    a.y = a.z = 0.0;
  }
  else if (std::fabs(nrm - 1.0) > EPS_FLOAT)
  {
    a.x /= nrm;
    a.y /= nrm;
    a.z /= nrm;
  }

  return a;
}

/// Sets the object from a rotation matrix
AAngled to_AAngle(const Matrix3d& m)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // get the arc-cosine of the angle (clip it as necessary)
  double acosangle = (m(X,X) + m(Y,Y) + m(Z,Z) - 1)*0.5; 
  if (acosangle > 1.0)
    acosangle = 1.0;
  else if (acosangle < -1.0)
    acosangle = -1.0;

  // compute angle of rotation
  AAngled a;
  a.angle = std::acos(acosangle);
  
  // if angle is 0, then axis is arbitrary
  if (a.angle < std::numeric_limits<double>::epsilon())
  {
    a.x = 1.0;
    a.y = a.z = 0.0;
    return a;
  }
  
  // if angle is pi then must solve for rx, ry, rz
  if (std::fabs(a.angle-M_PI) < std::numeric_limits<double>::epsilon())
  {
    a.x = AAngled::safe_sqrt((m(X,X)+1)*0.5);
    a.y = AAngled::safe_sqrt((m(Y,Y)+1)*0.5);
    a.z = AAngled::safe_sqrt((m(Z,Z)+1)*0.5);
    assert(!std::isnan(a.x));
    assert(!std::isnan(a.y));
    assert(!std::isnan(a.z));
    return a;
  }
  
  // standard case
  double constant = 1.0/(2.0*std::sin(a.angle));
  a.x = constant * (m(Z,Y) - m(Y,Z));
  a.y = constant * (m(X,Z) - m(Z,X));
  a.z = constant * (m(Y,X) - m(X,Y));
  
  // normalize the axis (generally not necessary, but safe...)
  double len = AAngled::safe_sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  assert(len != 0.0);
  if (std::fabs(len - 1.0) > EPS_FLOAT)
  {
    double ilen = 1.0/len;
    a.x *= ilen;
    a.y *= ilen;
    a.z *= ilen;
  }

  assert(!std::isnan(a.angle));
  assert(!std::isnan(a.x));
  assert(!std::isnan(a.y));
  assert(!std::isnan(a.z));

  return a;
}

