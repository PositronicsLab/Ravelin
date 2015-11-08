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

  // initialize the constraint type to unknown
  _constraint_type = eUnknown;
}

/// Resets the force on the joint
void JOINT::reset_force()
{
  force.set_zero(num_dof());
}

/// Adds a force to the joint
void JOINT::add_force(const VECTORN& f)
{
  force += f;
}

/// Sets inboard link
void JOINT::set_inboard_link(shared_ptr<RIGIDBODY> inboard, bool update_pose)
{
  _inboard_link = inboard;
  if (!inboard)
    return;

  // setup F's pose relative to the inboard
  set_inboard_pose(inboard->get_pose(), update_pose);
}

/// Sets the pointer to the outboard link for this joint
/**
 * \note also points the outboard link to this joint
 */
void JOINT::set_outboard_link(shared_ptr<RIGIDBODY> outboard, bool update_pose)
{
  _outboard_link = outboard;
  if (!outboard)
    return;

  // get the outboard pose
  if (outboard->_F->rpose)
    throw std::runtime_error("Joint::set_inboard_link() - relative pose on inboard link already set");

  // setup Fb's pose relative to the outboard 
  set_outboard_pose(outboard->_F, update_pose);
}

/// Determines q tare
void JOINT::determine_q_tare()
{
  // determine q tare
  determine_q(_q_tare);
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
  qdd.set_zero(NDOF);
  force.set_zero(NDOF);
  lambda.set_zero(NEQ);
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

  // update spatial axes if both poses are set
  if (_F->rpose && _Fb->rpose)
    update_spatial_axes();
}

/// Sets the location of this joint
void JOINT::set_location(const VECTOR3& point, shared_ptr<RIGIDBODY> inboard, shared_ptr<RIGIDBODY> outboard) 
{
  assert(inboard && outboard);

  // convert p to the inboard and outboard links' frames
  VECTOR3 pi = POSE3::transform_point(inboard->get_pose(), point);
  VECTOR3 po = POSE3::transform_point(outboard->get_pose(), point);

  // set _F's and Fb's origins
  _F->x = ORIGIN3(pi);
  _Fb->x = ORIGIN3(po);

  // invalidate all outboard pose vectors
  outboard->invalidate_pose_vectors();

  // set inboard and outboard links
  set_inboard_link(inboard, false);
  set_outboard_link(outboard, false);
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

/// Transforms a Jacobian (if necessary)
/**
 * If the inboard/outboard link is part of an articulated body, transforms
 * the Jacobian to the coordinates of that body.
 * \return true if any transformation is done (will be found in output)
 */
bool JOINT::transform_jacobian(MATRIXN& J, bool use_inboard, MATRIXN& output)
{
  const shared_ptr<const POSE3> GLOBAL_3D;
  MATRIXN Jm;

  // get the appropriate link
  shared_ptr<RIGIDBODY> rb = (use_inboard) ? get_inboard_link() : get_outboard_link();

  // see whether it is part of an articulated body
  shared_ptr<ARTICULATED_BODY> ab = rb->get_articulated_body();
  if (ab)
  {
    shared_ptr<RC_ARTICULATED_BODY> rcab = dynamic_pointer_cast<RC_ARTICULATED_BODY>(ab);
    if (rcab)
    {
      rcab->calc_jacobian(GLOBAL_3D, rb, Jm);
      J.mult(Jm, output);
      return true;
    }
  }

  // no transformation necessary
  return false;
}

/// Computes the constraint jacobian with respect to a body *numerically*
void JOINT::calc_constraint_jacobian_numeric(bool inboard, MATRIXN& Cq)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;
  MATRIXN tmp;
  REAL C[SPATIAL_DIM];

  // get the inboard and outboard links
  shared_ptr<RIGIDBODY> board = (inboard) ? get_inboard_link() : get_outboard_link();

  // if the body is disabled, quit now
  if (!board->is_enabled())
  {
    Cq.resize(num_constraint_eqns(), 0);
    return;
  }

  // get the super body and the number of generalized coordinates 
  shared_ptr<DYNAMIC_BODY> super = board->get_super_body();
  const unsigned NGC = super->num_generalized_coordinates(DYNAMIC_BODY::eSpatial);

  // initialize Cq
  Cq.set_zero(num_constraint_eqns(), NGC);
  
  // store the current generalized velocity
  VECTORN vsave;
  super->get_generalized_velocity(DYNAMIC_BODY::eSpatial, vsave);

  // we want to find J such that J*v = \dot{\phi}, where \dot{\phi} is the
  // constraint velocities 
  // J*v = \dot{\phi}
  // J = inv(V)*\dot{\Phi}
  // if V is the identity matrix, no inversion is necessary
  VECTORN v;
  for (unsigned i=0; i< NGC; i++)
  {
    // set the generalized velocity to all zeros and a single one
    v.set_zero(NGC);
    v[i] = 1.0;
    super->set_generalized_velocity(DYNAMIC_BODY::eSpatial, v);     

    // evaluate the time derivative of the constraint equation and set
    // the appropriate column of the matrix 
    evaluate_constraints_dot(C);
    for (unsigned j=0; j< num_constraint_eqns(); j++)
      Cq(j,i) = C[j];
  }

  // restore the generalized velocity
  super->set_generalized_velocity(DYNAMIC_BODY::eSpatial, vsave);     
}

/// Evaluates the time derivative of the constraint equations
void JOINT::evaluate_constraints_dot(REAL C_dot[])
{
  const unsigned SPATIAL_D = 6;
  REAL C_before[SPATIAL_D];
  const REAL DT = 1e-5;
  vector<VECTORN> q;
  vector<shared_ptr<DYNAMIC_BODY> > supers;

  // get the inboard and outboard links
  shared_ptr<RIGIDBODY> inboard = get_inboard_link();
  shared_ptr<RIGIDBODY> outboard = get_outboard_link();

  // get the inboard and outboard super bodies
  shared_ptr<DYNAMIC_BODY> inboard_sb = inboard->get_super_body();
  shared_ptr<DYNAMIC_BODY> outboard_sb = outboard->get_super_body();

  // make super bodies unique, if necessary
  supers.push_back(inboard_sb);
  supers.push_back(outboard_sb);
  std::sort(supers.begin(), supers.end());
  supers.erase(std::unique(supers.begin(), supers.end()), supers.end());

  // save configurations of super bodies
  for (unsigned i=0; i< supers.size(); i++)
  {
    q.push_back(VECTORN());
    supers[i]->get_generalized_coordinates_euler(q.back());
  }

  // evaluate constraint equations
  evaluate_constraints(C_before);
  
  // integrate configurations of super bodies forward
  for (unsigned i=0; i< supers.size(); i++)
  {
    VECTORN qd;
    supers[i]->get_generalized_velocity(DYNAMIC_BODY::eEuler, qd);
    qd *= DT;
    qd += q[i];
    supers[i]->set_generalized_coordinates_euler(qd);
  }

  // re-evaluate constraint equations
  evaluate_constraints(C_dot);

  // compute time derivative
  for (unsigned i=0; i< num_constraint_eqns(); i++)
  {
    C_dot[i] -= C_before[i];
    C_dot[i] /= DT;
  }

  // restore configurations of super bodies
  for (unsigned i=0; i< supers.size(); i++)
    supers[i]->set_generalized_coordinates_euler(q[i]);
}

