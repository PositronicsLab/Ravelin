/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Computes the dot product between a velocity and a momentum
REAL SVELOCITY::dot(const SMOMENTUM& v2) const
{
  // verify that both vectors are defined in the same frame
  #ifndef NEXCEPT
  if (pose != v2.pose)
    throw FrameException();
  #endif

  const REAL* d1 = data();
  const REAL* d2 = v2.data();
  return d1[3]*d2[0] + d1[4]*d2[1] + d1[5]*d2[2]+
         d1[0]*d2[3] + d1[1]*d2[4] + d1[2]*d2[5];
}

/// Computes the dot product between a velocity and a force 
REAL SVELOCITY::dot(const SFORCE& v2) const
{
  // verify that both vectors are defined in the same frame
  #ifndef NEXCEPT
  if (pose != v2.pose)
    throw FrameException();
  #endif

  const REAL* d1 = data();
  const REAL* d2 = v2.data();
  return d1[3]*d2[0] + d1[4]*d2[1] + d1[5]*d2[2]+
         d1[0]*d2[3] + d1[1]*d2[4] + d1[2]*d2[5];
}

/// Returns the spatial cross product between two velocity vectors
SVELOCITY SVELOCITY::cross(const SVELOCITY& v) const
{
  // verify that both vectors are defined in the same frame
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  VECTOR3 ax = get_angular();
  VECTOR3 bx = get_linear();

  // multiply
  VECTOR3 vtop = v.get_angular();
  VECTOR3 vbot = v.get_linear();
  VECTOR3 top = VECTOR3::cross(ax, vtop);
  VECTOR3 bot = VECTOR3::cross(bx, vtop) + VECTOR3::cross(ax, vbot);
  return SVELOCITY(top, bot, pose);
}

/// Returns the spatial cross product between a velocity and an acceleration
SFORCE SVELOCITY::cross(const SMOMENTUM& m) const
{
  // verify that both vectors are defined in the same frame
  #ifndef NEXCEPT
  if (pose != m.pose)
    throw FrameException();
  #endif

  VECTOR3 ax = get_angular();
  VECTOR3 bx = get_linear();

  // multiply
  VECTOR3 vtop = m.get_angular();
  VECTOR3 vbot = m.get_linear();
  VECTOR3 bot = VECTOR3::cross(ax, vtop);
  VECTOR3 top = VECTOR3::cross(bx, vtop) + VECTOR3::cross(ax, vbot);
  return SFORCE(top, bot, pose);
}

