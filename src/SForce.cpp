/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

REAL SFORCE::dot(const SVELOCITY& v2) const
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

