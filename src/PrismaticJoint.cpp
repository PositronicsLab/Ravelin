/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Ravelin/UndefinedAxisException.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

/// Initializes the joint
/**
 * The axis of rotation is set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
PRISMATICJOINT::PRISMATICJOINT() : JOINT()
{
  // init the joint data
  init_data();

  // set the axes and associated vectors to zeros initially
  _u.set_zero();
  _v2.set_zero();
  _ui.set_zero();
  _uj.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.clear();
}

/// Sets the axis of translation for this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void PRISMATICJOINT::set_axis(const VECTOR3& axis) 
{
  // check that axis is ok 
  if (std::fabs(axis.norm() - (REAL) 1.0) > EPS)
    throw UndefinedAxisException(); 
 
  // normalize the axis, in case caller did not 
  VECTOR3 naxis = VECTOR3::normalize(axis); 

  // transform axis to joint frame
  _u = POSE3::transform_vector(get_pose(), naxis);

  // set the joint axis in the inner link frame
  update_spatial_axes(); 
/*
  // set the joint axis in the outer link frame and setup associated
  // vectors needed for maximal coordinate articulated bodies
  _v2 = inner->get_transform_vector().mult_vector(naxis);
  _v2 = outer->get_transform_vector().transpose_mult_vector(_v2);
  VECTOR3::determine_orthonormal_basis(_u, _ui, _uj);
*/
}        

/// Updates the spatial axis for this joint
void PRISMATICJOINT::update_spatial_axes()
{
  const VECTOR3 ZEROS_3((REAL) 0.0, (REAL) 0.0, (REAL) 0.0, get_pose());

  // call parent method
  JOINT::update_spatial_axes();

  // if the axis is not normal, return
  if (std::fabs(_u.norm_sq() - (REAL) 1.0) > EPS)
    return;

  // update the spatial axis in link coordinates
  _s[0].set_linear(_u);
  _s[0].set_angular(ZEROS_3);
}

/// Determines (and sets) the value of Q from the axis and the inboard link and outboard link transforms
void PRISMATICJOINT::determine_q(VECTORN& q)
{
  shared_ptr<const POSE3> GLOBAL; 
  shared_ptr<const POSE3> Fo = get_outboard_pose();

  // verify that the outboard pose is set
  if (!Fo)
    throw std::runtime_error("determine_q() called on when outboard pose is NULL!");

  // if axis is not defined, can't use this method
  if (std::fabs(_u.norm() - (REAL) 1.0) > EPS)
    throw UndefinedAxisException();

  // get the poses of the joint and outboard link
  shared_ptr<const POSE3> Fj = get_pose();

  // compute transforms
  TRANSFORM3 wTo = POSE3::calc_relative_pose(Fo, GLOBAL); 
  TRANSFORM3 jTw = POSE3::calc_relative_pose(GLOBAL, Fj);
  TRANSFORM3 jTo = jTw * wTo;

  // get the vector of translation
  VECTOR3 x(jTo.x, _u.pose);
  q.resize(num_dof());
  q[DOF_1] = x.norm();

  // see whether to reverse q
  if (x.dot(_u) < (REAL) 0.0)
    q[DOF_1] = -q[DOF_1];
}

/// Gets the (local) transform for this joint
shared_ptr<const POSE3> PRISMATICJOINT::get_induced_pose()
{
  _Fprime->x = ORIGIN3(_u * (this->q[DOF_1] + this->_q_tare[DOF_1]));
  return _Fprime;
}

/// Gets the derivative fo the spatial axes for this joint
const vector<SVELOCITY>& PRISMATICJOINT::get_spatial_axes_dot()
{
  return _s_dot;
}

/// Evaluates the constraint equations
void PRISMATICJOINT::evaluate_constraints(REAL C[])
{
  shared_ptr<const POSE3> GLOBAL;

  // get the two links
  shared_ptr<const POSE3> Fi = get_inboard_pose();
  shared_ptr<const POSE3> Fo = get_outboard_pose();

  // This code was developed using [Shabana, 2003], p. 437; some variable names
  // have been altered

  // get axis (with respect to inboard) 
  VECTOR3 v1 = get_axis();

  // get axis (with respect to outboard)
  VECTOR3 v2 = POSE3::transform_vector(GLOBAL, _v2);

  // determine v1i, v1j
  VECTOR3 v1i, v1j;
  VECTOR3::determine_orthonormal_basis(get_axis(), v1i, v1j);

  // determine h1 and h2
  VECTOR3 h1 = Fi->transform_vector(GLOBAL, _ui);
  VECTOR3 h2 = Fo->transform_vector(GLOBAL, _uj);

  // determine the global positions of the attachment points and subtract them
  const VECTOR3 p1 = POSE3::transform_point(GLOBAL, get_location(false)); 
  const VECTOR3 p2 = POSE3::transform_point(GLOBAL, get_location(true)); 
  VECTOR3 r12 = p1 - p2; 

  // evaluate the constraint equations
  C[0] = v1i.dot(v2);
  C[1] = v1j.dot(v2);
  C[2] = v1i.dot(r12);
  C[3] = v1j.dot(r12);
  C[4] = h1.dot(h2);
}

