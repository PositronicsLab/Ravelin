/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Ravelin/UndefinedAxisException.h>

using boost::shared_ptr;
using std::vector;
using namespace Ravelin;

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
UNIVERSALJOINT::UNIVERSALJOINT() : JOINT()
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1].set_zero();
  _u[eAxis2].set_zero();
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _h2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;
}

/// Gets the axis for this joint
VECTOR3 UNIVERSALJOINT::get_axis(Axis a) const
{
  const unsigned X = 0;

  // axis one is already set 
  if (a == eAxis1)
    return _u[DOF_1];
  else
    return _u[DOF_2];
}

/// Sets an axis of this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void UNIVERSALJOINT::set_axis(const VECTOR3& axis, Axis a) 
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // normalize the axis in case the user did not
  VECTOR3 naxis = VECTOR3::normalize(axis); 

  // set the axis
  _u[a] = POSE3::transform_vector(get_pose(), naxis); 

  // update the spatial axes
  update_spatial_axes(); 

  // set the second axis in outer link coordinates, if we just changed it
  if (a == eAxis2)
  {
    shared_ptr<const POSE3> Fi = get_inboard_pose();
    shared_ptr<const POSE3> Fo = get_outboard_pose();
    _h2 = POSE3::transform_vector(Fo, naxis);
  }
}        

/// Updates the spatial axis for this joint
void UNIVERSALJOINT::update_spatial_axes()
{
  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // call the parent method
  JOINT::update_spatial_axes();

  // update the spatial axes in joint pose 
  _s[0].set_angular(_u[DOF_1]);
  _s[0].set_linear(ZEROS_3);

  // update the spatial axes in link coordinates
  _s_dot[0].set_angular(ZEROS_3);
  _s_dot[0].set_linear(ZEROS_3);
}

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& UNIVERSALJOINT::get_spatial_axes()
{
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes() called with NULL outboard link");

  // get the transformed second axis
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 u1(R * ORIGIN3(_u[DOF_2]), get_pose());

  // get the inner link pose
  boost::shared_ptr<const POSE3> pose = get_inboard_pose(); 

  // get the second spatial axis
  VECTOR3 u2(R * ORIGIN3(_u[1]), get_pose());

  // update the spatial axes in link coordinates
  _s[0].set_angular(_u[0]);
  _s[0].set_linear(ZEROS_3);
  _s[1].set_angular(u2);

  // update the second spatial axes
  _s[1].set_angular(u1);
  _s[1].set_linear(ZEROS_3);

  // use the JOINT function to do the rest
  return JOINT::get_spatial_axes();
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& UNIVERSALJOINT::get_spatial_axes_dot()
{
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes_dot() called with NULL outboard link");

  // get the transformed second axis
  const double S1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double C1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  REAL qd1 = qd[DOF_1];
  VECTOR3 u2(0.0, -S1*qd1, C1*qd1, get_pose());

  // get the second spatial axis time derivative. This is:
  // R * \dot{u2} + \dot{R} * u
  // note that \dot{u2} is zero
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 omega = _u[DOF_1] * qd1;
  ORIGIN3 u1 = ORIGIN3::cross(ORIGIN3(omega), R * ORIGIN3(_u[DOF_2]));
  _s_dot[1].set_angular(VECTOR3(u1, get_pose()));
  _s_dot[1].set_linear(ZEROS_3);

/*
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;
  const VECTORN& qd = this->qd;

  // compute some needed quantities
  const REAL c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const REAL s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);

  // form the time derivative of the spatial axis for the second DOF; note that spatial
  // axis for first DOF is constant, so time-derivative is zero 
  MATRIX3 Rx = MATRIX3::zero();
  REAL qd1 = qd[DOF_1];
  ORIGIN3 v(0.0, -s1*qd1, c1*qd1); 
  Rx.set_column(Y, v);
  ORIGIN3 cross = ORIGIN3::cross(ORIGIN3(1.0, 0.0, 0.0), v); 
  Rx.set_column(Z, cross);

  // compute u
  ORIGIN3 ux = (_R * Rx * MATRIX3::transpose(_R)).get_column(Y);
  VECTOR3 u(ux, get_pose());

  // get the time derivative of second spatial axis
  REAL qd1 = qd[DOF_1];
  VECTOR3 u2(0.0, -s1*qd1, c1*qd1, get_pose());

  // get the second spatial axis time derivative. This is:
  // R * \dot{u2} + \dot{R} * u
  // note that \dot{u2} is zero
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 omega = _u[DOF_1] * qd1;
  VECTOR3 dotRu = VECTOR3::cross(omega, VECTOR3(R * ORIGIN3(_u[DOF_2]), get_pose()));

  // update the spatial axis in link coordinates; note that axis 1 is always
  // set to zero (init'd in constructor)
  _s_dot[1].set_upper(dotRu);
  _s_dot[1].set_lower(ZEROS_3);
*/

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void UNIVERSALJOINT::determine_q(VECTORN& q)
{
  std::cerr << "SPHERICALJOINT::determine_q() warning- determine_q(.) is not currently functional" << std::endl;
}

/// Gets the (local) transform for this joint
MATRIX3 UNIVERSALJOINT::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q and _q_tare
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // compute some needed quantities
  const REAL Q1 = q[DOF_1]+q_tare[DOF_1];
  const REAL Q2 = q[DOF_2]+q_tare[DOF_2];

  // get the first rotation matrix
  ORIGIN3 u1(_u[DOF_1]);
  MATRIX3 R = AANGLE(u1, Q1);

  // get the second axis 
  ORIGIN3 u2 = R * ORIGIN3(_u[DOF_2]);

  // transform using orientation matrix
  return R * AANGLE(u2, Q2);
}

/// Gets the transform induced by this joint
shared_ptr<const POSE3> UNIVERSALJOINT::get_induced_pose()
{
  // get the rotation
  _Fprime->q = get_rotation();
  return _Fprime;
}

/// Computes the constraint jacobian with respect to a body
void UNIVERSALJOINT::calc_constraint_jacobian(bool inboard, SHAREDMATRIXN& Cq)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;

  // get the two links
  shared_ptr<const POSE3> Pi = get_inboard_pose();
  shared_ptr<const POSE3> Po = get_outboard_pose();

  // get the attachment frames in the global frame
  TRANSFORM3 wPi = POSE3::calc_relative_pose(_F, GLOBAL);
  TRANSFORM3 wPo = POSE3::calc_relative_pose(_Fb, GLOBAL);

  // get the vector from the inboard and outboard poses to the joint pose 
  VECTOR3 joint_pos(0.0, 0.0, 0.0, _F);
  VECTOR3 ui = POSE3::transform_vector(Pi, joint_pos);
  VECTOR3 uo = POSE3::transform_vector(Po, joint_pos);

  // setup the constraint equations (from Shabana, p. 432)
  if (inboard)
  {
    // get the vector from the inboard pose to the joint pose 
    ORIGIN3 u = ORIGIN3(POSE3::transform_point(Pi, joint_pos));

    // get the information necessary to compute the constraint equations
    QUAT& R = wPi.q;
    ORIGIN3 Ru = R*u;

    // get positional components of Cq
    SHAREDMATRIXN Cq_trans = Cq.block(0, 3, 0, 3);
    Cq_trans.set_identity();
    SHAREDMATRIXN Cq_rot = Cq.block(0, 3, 3, 6);
    Cq_rot = MATRIX3::skew_symmetric(-Ru);

    // Jacobian of dot(Ri*_u[DOF_1], Ro*_h2) w.r.t. inboard is
    // dot(w x Ri*_u[DOF_1], Ro*h2) = (h2*Ro)' * (-Ri*u) x w 
    // transpose((h2 * Ro)' * skew(-Ri*u)) = skew(Ri * u) * Ro * h2 
    ORIGIN3 h2w = wPo.q * ORIGIN3(_h2);
    ORIGIN3 result = MATRIX3::skew_symmetric(R*ORIGIN3(_u[DOF_1])) * h2w;
    SHAREDVECTORN last_row = Cq.row(3);
    last_row.segment(0, 3).set_zero();
    last_row.segment(4, 6) = result; 
  }
  else
  {
    // get the vector from the outboard pose to the joint pose 
    ORIGIN3 u = ORIGIN3(POSE3::transform_point(Po, joint_pos));

    // get the information necessary to compute the constraint equations
    MATRIX3 R = wPo.q;
    ORIGIN3 Ru = R*u;

    // get positional components of Cq
    SHAREDMATRIXN Cq_trans = Cq.block(0, 3, 0, 3);
    Cq_trans.set_identity();
    Cq_trans *= -1.0;
    SHAREDMATRIXN Cq_rot = Cq.block(0, 3, 3, 6);
    Cq_rot = MATRIX3::skew_symmetric(Ru);

    // Jacobian of dot(Ri*_u[DOF_1], Ro*_h2) w.r.t. outboard is
    // dot(Ri*_u[DOF_1], Ro*h2) = (Ri*_h1)' * (-Ro*u) x w 
    // transpose((Ri * _h1)' * skew(-Ro*u)) = skew(Ro * u) * Ri * h1 
    ORIGIN3 h1w = R * ORIGIN3(_u[DOF_1]);
    ORIGIN3 result = MATRIX3::skew_symmetric(R*ORIGIN3(_h2)) * h1w;
    SHAREDVECTORN last_row = Cq.row(3);
    last_row.segment(0, 3).set_zero();
    last_row.segment(4, 6) = result; 
  }
}

/// Evaluates the constraint equations
void UNIVERSALJOINT::evaluate_constraints(REAL C[])
{
  shared_ptr<const POSE3> GLOBAL;
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  shared_ptr<const POSE3> inner = get_inboard_pose();
  shared_ptr<const POSE3> outer = get_outboard_pose();

  // This code was developed using [Shabana, 2003], p. 438; variable names
  // have been altered however

  // determine h1 and h2 in global coordinates
  VECTOR3 h1 = inner->transform_vector(GLOBAL, _u[DOF_1]);
  VECTOR3 h2 = outer->transform_vector(GLOBAL, _h2);

  // determine the global positions of the attachment points and subtract them
  VECTOR3 r1 = get_location(false);
  VECTOR3 r2 = get_location(true);
  VECTOR3 r12 = r1 - r2;

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = h1.dot(h2);
}


