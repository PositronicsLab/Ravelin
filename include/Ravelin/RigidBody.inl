/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

template <class OutputIterator>
OutputIterator RIGIDBODY::get_parent_links(OutputIterator begin) const
{
  BOOST_FOREACH(boost::shared_ptr<JOINT> j, _inner_joints)
    *begin++ = get_parent_link(j);

  return begin;
}

template <class OutputIterator>
OutputIterator RIGIDBODY::get_child_links(OutputIterator begin) const
{
  BOOST_FOREACH(boost::shared_ptr<JOINT> j, _outer_joints)
    *begin++ = get_child_link(j);

  return begin;
}

/// Sets generalized coordinates using a templated vector
template <class V>
void RIGIDBODY::get_generalized_coordinates_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gc) 
{
  const unsigned N_SPATIAL = 6, N_EULER = 7;
  const boost::shared_ptr<const POSE3> GLOBAL;

  // special case: disabled body
  if (!_enabled)
  {
    gc.resize(0);
    return;
  }

  // resize vector
  switch (gctype)
  {
    case DYNAMICBODY::eEuler:   gc.resize(N_EULER); break;
    case DYNAMICBODY::eSpatial: gc.resize(N_SPATIAL); break;
  }

  // get current inertial pose 
  POSE3 P = *_F;
  P.update_relative_pose(GLOBAL);

  // get linear components
  gc[0] = P.x[0];
  gc[1] = P.x[1];
  gc[2] = P.x[2];

  // get angular components 
  if (gctype == DYNAMICBODY::eSpatial)
    P.q.to_rpy(gc[3], gc[4], gc[5]);
  else
  {
    // return the generalized position using Euler parameters
    assert(gctype == DYNAMICBODY::eEuler);
    gc[3] = P.q.x;
    gc[4] = P.q.y;
    gc[5] = P.q.z;
    gc[6] = P.q.w;
  }
}

/// Sets the generalized coordinates of this rigid body (does not call articulated body)
template <class V>
void RIGIDBODY::set_generalized_coordinates_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gc)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // do easiest case first 
  if (gctype == DYNAMICBODY::eSpatial)
  {
    // this isn't correct
    assert(false);
/*
    // note: generalized coordinates in eSpatial are always set with regard to 
    // the global frame
    Origin3d x(gc[0], gc[1], gc[2]);
    QUAT q = QUAT::rpy(gc[3], gc[4], gc[5]); 

    // convert the pose to the correct relative frame
    POSE3 P(q, x);


    P.update_relative_pose(_F->rpose);

    // set the transform
    set_pose(P);
*/
  }
  else
  {
    assert(gctype == DYNAMICBODY::eEuler);

    // get the position
    ORIGIN3 x(gc[0], gc[1], gc[2]);

    // get the unit quaternion
    QUAT q;
    q.x = gc[3];
    q.y = gc[4];
    q.z = gc[5];
    q.w = gc[6];

    // normalize the unit quaternion, just in case
    q.normalize();

    POSE3 P(q, x);
    set_pose(P);
/*
    // get the transform from the link pose to the inertial pose
    Transform3d lTm = POSE3::calc_relative_pose(_jF, _F);

    // set the pose 
    *_jF = P;
    *_F = lTm.apply_transform();
    _jF->update_relative_pose(_F);
*/
    

    // invalidate the pose vectors 
    invalidate_pose_vectors();    
  }
}

/// Sets the generalized velocity of this rigid body (does not call articulated body version)
template <class V>
void RIGIDBODY::set_generalized_velocity_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gv)
{
  const boost::shared_ptr<const POSE3> GLOBAL;

  // special case: disabled body
  if (!_enabled)
    return;

  // get the velocity
  SVELOCITY xd;
  xd.pose = _F2;
 
  // set the linear velocity first
  xd.set_linear(VECTOR3(gv[0], gv[1], gv[2]));
 
  // simplest case: spatial coordinates
  if (gctype == DYNAMICBODY::eSpatial)
    xd.set_angular(VECTOR3(gv[3], gv[4], gv[5]));
  else
  {
    assert(gctype == DYNAMICBODY::eEuler);

    // get the quaternion derivatives
    QUAT qd;
    qd.x = gv[3] * (REAL) 2.0;
    qd.y = gv[4] * (REAL) 2.0;
    qd.z = gv[5] * (REAL) 2.0;
    qd.w = gv[6] * (REAL) 2.0;

    // setup the pose
    POSE3 F = *_F;
    F.update_relative_pose(GLOBAL);

    // setup the angular component
    xd.set_angular(F.q.G_mult(qd.x, qd.y, qd.z, qd.w));
  }

  // set the velocity
  set_velocity(xd);
}

/// Gets the generalized velocity of this rigid body (does not call articulated body version)
template <class V>
void RIGIDBODY::get_generalized_velocity_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gv) 
{
  const unsigned N_SPATIAL = 6, N_EULER = 7;
  const boost::shared_ptr<const POSE3> GLOBAL;

  // special case: disabled body
  if (!_enabled)
  {
    gv.resize(0);
    return;
  }

  // resize the generalized velocity vector
  switch (gctype)
  {
    case DYNAMICBODY::eEuler:   gv.resize(N_EULER); break;
    case DYNAMICBODY::eSpatial: gv.resize(N_SPATIAL); break;
  }

  // get/set linear components of velocity
  VECTOR3 lv = _xdcom.get_linear();
  gv[0] = lv[0];
  gv[1] = lv[1];
  gv[2] = lv[2];

  // determine the proper generalized coordinate type
  if (gctype == DYNAMICBODY::eSpatial)
  {
    // get/set angular components of velocity
    VECTOR3 av = _xdcom.get_angular();
    gv[3] = av[0];
    gv[4] = av[1];
    gv[5] = av[2];
  }
  else
  {
    assert(gctype == DYNAMICBODY::eEuler);

    // going to need Euler coordinate derivatives
    POSE3 F = *_F;
    F.update_relative_pose(GLOBAL);
    QUAT qd = F.q.G_transpose_mult(_xdcom.get_angular()) * (REAL) 0.5;

    // setup the angular components 
    gv[3] = qd.x;
    gv[4] = qd.y;
    gv[5] = qd.z;
    gv[6] = qd.w; 
  }
}

/// Gets the generalized acceleration of this body (does not call articulated body version)
template <class V>
void RIGIDBODY::get_generalized_acceleration_generic(V& ga)
{
  const unsigned N_SPATIAL = 6;

  // special case: body is disabled
  if (!_enabled)
  {
    ga.resize(0);
    return;
  } 

  // setup the linear components
  ga.resize(N_SPATIAL);

  // get linear and angular components
  VECTOR3 la = _xddcom.get_linear();
  VECTOR3 aa = _xddcom.get_angular();

  // set linear components
  ga[0] = la[0];
  ga[1] = la[1];
  ga[2] = la[2];
  ga[3] = aa[0];
  ga[4] = aa[1];
  ga[5] = aa[2];
}

/// Sets the generalized acceleration of this rigid body (does not call articulated body version)
template <class V>
void RIGIDBODY::set_generalized_acceleration_generic(V& ga)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // get the acceleration
  SACCEL xdd;
  xdd.pose = _F2;
 
  // set the linear acceleration first
  xdd.set_linear(VECTOR3(ga[0], ga[1], ga[2]));
  xdd.set_angular(VECTOR3(ga[3], ga[4], ga[5]));

  // set the acceleration
  set_accel(xdd);
}


