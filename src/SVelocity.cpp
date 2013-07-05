/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

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
  VECTOR3 vtop = m.get_linear();
  VECTOR3 vbot = m.get_angular();
  VECTOR3 top = VECTOR3::cross(ax, vtop);
  VECTOR3 bot = VECTOR3::cross(bx, vtop) + VECTOR3::cross(ax, vbot);
  return SFORCE(top, bot, pose);
}

