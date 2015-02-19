/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Ravelin/Opsf.h>

using namespace Ravelin;

Matrix3f to_Matrix3f(const Quatf& q)
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // verify that the quaternion is normalized
  assert(std::fabs(q.magnitude()) - 1 < EPS_FLOAT);

  // setup repeated products
  const float xx = q.x*q.x;
  const float xy = q.x*q.y;
  const float xz = q.x*q.z;
  const float xw = q.x*q.w;
  const float yy = q.y*q.y;
  const float yz = q.y*q.z;
  const float yw = q.y*q.w;
  const float zz = q.z*q.z;
  const float zw = q.z*q.w; 
  const float ww = q.w*q.w;

  Matrix3f m;
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

Matrix3f to_Matrix3(const AAnglef& a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  float x = a.x;
  float y = a.y;
  float z = a.z;
  float ca = cos(a.angle);
  float sa = sin(a.angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const float SOMEWHAT_EPS_FLOAT = 1e-2;
  float va = (std::fabs(a.angle) > SOMEWHAT_EPS_FLOAT) ? 1 - ca : (sa*sa)/(1+ca);

  // setup the matrix
  Matrix3f m;
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
Quatf to_Quat(const Matrix3f& m)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // core computation
  Quatf q;
  q.w = std::sqrt(std::max(0.0f, 1.0f + m(X,X) + m(Y,Y) + m(Z,Z))) * 0.5f;
  q.x = std::sqrt(std::max(0.0f, 1.0f + m(X,X) - m(Y,Y) - m(Z,Z))) * 0.5f;
  q.y = std::sqrt(std::max(0.0f, 1.0f - m(X,X) + m(Y,Y) - m(Z,Z))) * 0.5f;
  q.z = std::sqrt(std::max(0.0f, 1.0f - m(X,X) - m(Y,Y) + m(Z,Z))) * 0.5f;

  // sign computation
  if (m(Z,Y) - m(Y,Z) < 0.0f)
    q.x = -q.x;
  if (m(X,Z) - m(Z,X) < 0.0f)
    q.y = -q.y;
  if (m(Y,X) - m(X,Y) < 0.0f)
    q.z = -q.z;

  #ifndef NEXCEPT
  if (!q.unit())
    std::cerr << "Quatf::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}

/// Sets quaternion to that represented by an axis-angle representation
Quatf to_Quat(const AAnglef& a)
{
  const float half = a.angle*0.5f;
  float sina = std::sin(half);
  Quatf q;
  q.x = a.x * sina;
  q.y = a.y * sina;
  q.z = a.z * sina;
  q.w = std::cos(half);

  #ifndef NEXCEPT
  if (!q.unit())
    std::cerr << "Quatf::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}

/// Sets axis angle from unit quaternion 
/**
 * \todo test this method
 */
AAnglef to_AAngle(const Quatf& q)
{
  Quatf qn = Quatf::normalize(q);
  AAnglef a;
  a.x = qn.x;
  a.y = qn.y;
  a.z = qn.z;
  a.angle = std::acos(qn.w) * 2.0f;

  // verify that [x y z] normalized
  float nrm = AAnglef::safe_sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  if (std::fabs(nrm) < EPS_FLOAT)
  {
    // axis is zero; set it to [1 0 0] arbitrarily
    a.x = 1.0f;
    a.y = a.z = 0.0f;
  }
  else if (std::fabs(nrm - 1.0f) > EPS_FLOAT)
  {
    a.x /= nrm;
    a.y /= nrm;
    a.z /= nrm;
  }

  return a;
}

/// Sets the object from a rotation matrix
AAnglef to_AAngle(const Matrix3f& m)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // get the arc-cosine of the angle (clip it as necessary)
  float acosangle = (m(X,X) + m(Y,Y) + m(Z,Z) - 1)*0.5f; 
  if (acosangle > 1.0f)
    acosangle = 1.0f;
  else if (acosangle < -1.0f)
    acosangle = -1.0f;

  // compute angle of rotation
  AAnglef a;
  a.angle = std::acos(acosangle);
  
  // if angle is 0, then axis is arbitrary
  if (a.angle < std::numeric_limits<float>::epsilon())
  {
    a.x = 1.0f;
    a.y = a.z = 0.0f;
    return a;
  }
  
  // if angle is pi then must solve for rx, ry, rz
  if (std::fabs(a.angle-M_PI) < std::numeric_limits<float>::epsilon())
  {
    a.x = AAnglef::safe_sqrt((m(X,X)+1)*0.5f);
    a.y = AAnglef::safe_sqrt((m(Y,Y)+1)*0.5f);
    a.z = AAnglef::safe_sqrt((m(Z,Z)+1)*0.5f);
    assert(!std::isnan(a.x));
    assert(!std::isnan(a.y));
    assert(!std::isnan(a.z));
    return a;
  }
  
  // standard case
  float constant = 1.0f/(2.0f*std::sin(a.angle));
  a.x = constant * (m(Z,Y) - m(Y,Z));
  a.y = constant * (m(X,Z) - m(Z,X));
  a.z = constant * (m(Y,X) - m(X,Y));
  
  // normalize the axis (generally not necessary, but safe...)
  float len = AAnglef::safe_sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  assert(len != 0.0f);
  if (std::fabs(len - 1.0f) > EPS_FLOAT)
  {
    float ilen = 1.0f/len;
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

