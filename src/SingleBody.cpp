/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;

/// Gets the "super" body
shared_ptr<DYNAMIC_BODY> SINGLE_BODY::get_super_body() const
{
  shared_ptr<ARTICULATED_BODY> ab = get_articulated_body();
  if (ab)
    return boost::dynamic_pointer_cast<DYNAMIC_BODY>(ab);
  else
  {
    shared_ptr<DYNAMIC_BODY> base = boost::const_pointer_cast<DYNAMIC_BODY>(shared_from_this());
    return boost::dynamic_pointer_cast<DYNAMIC_BODY>(base);
  }
}

REAL SINGLE_BODY::calc_point_vel(const VECTOR3& point, const VECTOR3& dir)
{
  VECTOR3 pv = calc_point_vel(point);
  return POSE3::transform_vector(dir.pose, pv).dot(dir);
}

