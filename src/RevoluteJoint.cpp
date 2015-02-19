/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

/// Initializes the joint
/**
 * The axis of rotation is set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
REVOLUTEJOINT::REVOLUTEJOINT() : JOINT()
{
  // init joint data
  init_data();

  // init the joint axes
  _u.set_zero(_F);
  _ui.set_zero();
  _uj.set_zero();
  _v2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}

/// Sets the axis of rotation for this joint
void REVOLUTEJOINT::set_axis(const VECTOR3& axis) 
{ 
  // normalize the joint axis, in case the caller didn't 
  VECTOR3 naxis = VECTOR3::normalize(axis);

  // transform the axis as necessary
  _u = POSE3::transform_vector(_F, naxis);

  // update the spatial axes
  update_spatial_axes(); 

  // setup ui and uj
  VECTOR3::determine_orthonormal_basis(_u, _ui, _uj);
 
  // make sure that the two poss are set 
  shared_ptr<const POSE3> b1 = get_inboard_pose();
  shared_ptr<const POSE3> b2 = get_outboard_pose();
  if (!b1 || !b2)
    throw std::runtime_error("Attempt to set joint axis without setting inboard and outboard links first!");
 
  // compute joint axis in outer link frame 
  _v2 = POSE3::transform_vector(b2, naxis); 
}        

/// Updates the spatial axis for this joint
void REVOLUTEJOINT::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;
  const VECTOR3 ZEROS_3((REAL) 0.0, (REAL) 0.0, (REAL) 0.0, get_pose());

  // call parent method
  JOINT::update_spatial_axes();

  // update the spatial axis in joint coordinates
  _s[0].set_angular(_u);
  _s[0].set_linear(ZEROS_3);
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void REVOLUTEJOINT::determine_q(VECTORN& q)
{
  shared_ptr<const POSE3> GLOBAL;

  // get the outboard link pointer
  shared_ptr<const POSE3> Fo = get_outboard_pose();
  
  // verify that the outboard link is set
  if (!Fo)
    throw std::runtime_error("determine_q() called on NULL outboard link!");

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - (REAL) 1.0) > EPS)
    throw UndefinedAxisException();

  // get the poses of the joint and outboard link
  shared_ptr<const POSE3> Fj = get_pose();

  // compute transforms
  TRANSFORM3 wTo = POSE3::calc_relative_pose(Fo, GLOBAL); 
  TRANSFORM3 jTw = POSE3::calc_relative_pose(GLOBAL, Fj);
  TRANSFORM3 jTo = jTw * wTo;

  // determine the joint transformation
  MATRIX3 R = jTo.q;
  AANGLE aa(R, _u);

  // set q 
  q.resize(num_dof());
  q[DOF_1] = aa.angle;
}

/// Gets the pose for this joint
shared_ptr<const POSE3> REVOLUTEJOINT::get_induced_pose()
{
  assert(std::fabs(_u.norm() - (REAL) 1.0) < EPS);

  // note that translation is set to zero in the constructors
  _Fprime->q = AANGLE(_u, this->q[DOF_1]+this->_q_tare[DOF_1]);

  return _Fprime;
}

/// Gets the derivative for the spatial axes for this joint
const std::vector<SVELOCITY>& REVOLUTEJOINT::get_spatial_axes_dot()
{
  return _s_dot;
}

/// Evaluates the constraint equations
void REVOLUTEJOINT::evaluate_constraints(REAL C[])
{
  shared_ptr<const POSE3> GLOBAL;
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two poses 
  shared_ptr<const POSE3> b1 = get_inboard_pose();
  shared_ptr<const POSE3> b2 = get_outboard_pose();

  // This code was developed using [Shabana, 2003], p. 435-436; variable names
  // have been altered however

  // determine v1, v1i, v1j, and v2 (all in global coordinates)
  VECTOR3 v1i, v1j;
  VECTOR3 v1 = POSE3::transform_vector(GLOBAL, _u);
  VECTOR3 v2 = POSE3::transform_vector(GLOBAL, _v2);
  VECTOR3::determine_orthonormal_basis(v1, v1i, v1j);

  // determine the global positions of the attachment points and subtract them
  VECTOR3 r1 = get_location(false);
  VECTOR3 r2 = get_location(true);
  VECTOR3 r12 = r1 - r2; 

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = v1i.dot(v2);
  C[4] = v1j.dot(v2); 
}


