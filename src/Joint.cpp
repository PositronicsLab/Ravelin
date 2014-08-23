/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
JOINT::JOINT()
{
  // initialize _q_tare
  _q_tare.resize(0);

  // initialize the two frames
  _F = shared_ptr<POSE3>(new POSE3);
  _Fb = shared_ptr<POSE3>(new POSE3);
  _Fprime = shared_ptr<POSE3>(new POSE3);
  _Fprime->rpose = _F;
}

/// Determines q tare
void JOINT::determine_q_tare()
{
  // determine q tare
  determine_q(_q_tare);
}

/// Evaluates the time derivative of the constraint
void JOINT::evaluate_constraints_dot(REAL C[6])
{
  REAL Cx[6];

/*
  // get the inboard and outboard links
  RigidBodyPtr in = get_inboard_link();
  RigidBodyPtr out = get_outboard_link();

  // get the linear angular velocities
  const SVELOCITY& inv = in->get_velocity();
  const SVELOCITY& outv = out->get_velocity();
  VECTOR3 lvi = inv.get_linear();
  VECTOR3 lvo = outv.get_linear();
  VECTOR3 avi = inv.get_angular();
  VECTOR3 avo = outv.get_angular();

  // compute
  const unsigned NEQ = num_constraint_eqns();
  for (unsigned i=0; i< NEQ; i++)
  {
    // TODO: fix this to do frame calculations
    calc_constraint_jacobian(DynamicBody::eSpatial, in, i, Cx);
    VECTOR3 lv(Cx[0], Cx[1], Cx[2]);
    VECTOR3 av(Cx[3], Cx[4], Cx[5]);
    C[i] = lv.dot(lvi) + av.dot(avi);
    calc_constraint_jacobian(DynamicBody::eSpatial, out, i, Cx);
    lv = VECTOR3(Cx[0], Cx[1], Cx[2]);
    av = VECTOR3(Cx[3], Cx[4], Cx[5]);
    C[i] += -lv.dot(lvo) - av.dot(avo);
  }
*/
}

/// Abstract method to update the local spatial axes
/**
 * Only applicable for reduced-coordinate articulated bodies
 */
void JOINT::update_spatial_axes()
{
  // setup the spatial axis vector and frame
  _s.resize(num_dof());
  for (unsigned i=0; i< _s.size(); i++)
    _s[i].pose = get_pose();
}

/// Sets the number of degrees-of-freedom for this joint
/**
 * \note resets all joint values (q) to zero
 */
void JOINT::init_data()
{
  const unsigned NDOF = num_dof();
  const unsigned NEQ = num_constraint_eqns();

  q.set_zero(NDOF);
  qd.set_zero(NDOF);
  _q_tare.set_zero(NDOF);
  _s.resize(NDOF);
}

/// Sets the inboard pose on the joint
void JOINT::set_inboard_pose(shared_ptr<const POSE3> pose, bool update_joint_pose) 
{
  if (update_joint_pose)
  {
    _F->update_relative_pose(pose);
    shared_ptr<const POSE3> old_rpose = _Fb->rpose;
    *_Fb = *_F;
    _Fb->update_relative_pose(old_rpose);    
  }
  else
    _F->rpose = pose;

  // update spatial axes if both poses are set
  if (_F->rpose && _Fb->rpose)
    update_spatial_axes();
}

/// Sets the outboard pose on the joint
void JOINT::set_outboard_pose(shared_ptr<POSE3> pose, bool update_joint_pose) 
{
  if (update_joint_pose)
  {
    // update Fb
    _Fb->update_relative_pose(pose);
    shared_ptr<const POSE3> old_rpose = _F->rpose;
    *_F = *_Fb;

    // now update F
    _F->update_relative_pose(old_rpose);    
  }
  else
    _Fb->rpose = pose;

  // set the pose to be relative to Fprime
  pose->update_relative_pose(_Fprime);

  // update spatial axes if both poses are set
  if (_F->rpose && _Fb->rpose)
    update_spatial_axes();
}

/// Sets the location of this joint
void JOINT::set_location(const VECTOR3& point) 
{
  // verify inboard and outboard poses are set
  if (!_F->rpose)
    throw std::runtime_error("JOINT::set_location() called and inboard pose not set");
  if (!_Fb->rpose)
    throw std::runtime_error("JOINT::set_location() called and outboard pose not set");

  // convert p to the inboard and outboard links' frames
  VECTOR3 pi = POSE3::transform_point(_F->rpose, point);
  VECTOR3 po = POSE3::transform_point(_Fb->rpose, point);

  // set _F's and Fb's origins
  _F->x = ORIGIN3(pi);
  _Fb->x = ORIGIN3(po);
}

/// Gets the location of this joint
/**
 * \param use_outboard if <b>true</b> then the joint position is calculated 
 *        using the outboard link rather than inboard link; the position will
 *        not be identical if the joint constraint is violated (therefore,
 *        this method will behave identically for reduced-coordinate 
 *        articulated bodies)
 */
VECTOR3 JOINT::get_location(bool use_outboard) const
{
  // compute the global position
  if (!use_outboard)
  {
    // joint is defined with respect to inboard frame
    VECTOR3 p(_F);
    p.set_zero();
    return p;
  }
  else
  {
    VECTOR3 p(_Fb);
    p.set_zero();
    return p;
  }
}

/// Gets the spatial axes for this joint
/**
 * Spatial axes describe the motion of the joint. Note that for rftype = eLink,
 * spatial axes are given in outboard link's frame. 
 */
const vector<SVELOCITY>& JOINT::get_spatial_axes()
{
  return _s;
}


