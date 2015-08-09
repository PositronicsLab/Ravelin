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

  // setup v1i and v1j 
  VECTOR3::determine_orthonormal_basis(_u, _v1i, _v1j);

  // set _ui
  VECTOR3 outboard_origin(0.0, 0.0, 0.0, _Fb);
  _ui = POSE3::transform_point(_F, outboard_origin);

  // set _uj
  _uj = POSE3::transform_vector(_Fb, _v1i); 

  // set the joint axis in the inner link frame
  update_spatial_axes(); 
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

/// Computes the constraint jacobian with respect to a body
void PRISMATICJOINT::calc_constraint_jacobian(bool inboard, SHAREDMATRIXN& Cq)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;

  // get the two links
  shared_ptr<const POSE3> Pi = get_inboard_pose();
  shared_ptr<const POSE3> Po = get_outboard_pose();

  // get the attachment frames in the global frame
  TRANSFORM3 wPi = POSE3::calc_relative_pose(_F, GLOBAL);
  TRANSFORM3 wPo = POSE3::calc_relative_pose(_Fb, GLOBAL);

  // get inboard and outboard rotations
  const QUAT& Ri = wPi.q;
  const QUAT& Ro = wPo.q;

  // get the vector from the inboard and outboard poses to the joint pose 
  VECTOR3 joint_pos(0.0, 0.0, 0.0, _F);
  ORIGIN3 ti(POSE3::transform_vector(Pi, joint_pos));
  ORIGIN3 to(POSE3::transform_vector(Po, joint_pos));

  // compute r12
  VECTOR3 joint_pos2(0.0, 0.0, 0.0, _Fb);
  ORIGIN3 r12 = ORIGIN3(POSE3::transform_point(GLOBAL, joint_pos)) -
                ORIGIN3(POSE3::transform_point(GLOBAL, joint_pos2));

  // setup the constraint equations (from Shabana, p. 432)
  if (inboard)
  {
    // get the vector from the inboard pose to the joint pose 
    ORIGIN3 u = ORIGIN3(POSE3::transform_point(Pi, joint_pos));

    // get the information necessary to compute the constraint equations
    MATRIX3 R = wPi.q;
    ORIGIN3 Ru = R*u;

    // Jacobian of dot(Ri* _v1i, Ro*_v2) w.r.t. inboard is
    // dot(w x Ri*_v1i, Ro*h2) = (_v2*Ro)' * (-Ri*_v1i) x w 
    // transpose((_v2 * Ro)' * skew(-Ri*_v1i)) = skew(Ri * _v1i) * Ro * _v2 
    ORIGIN3 v2w = wPo.q * ORIGIN3(_v2);
    ORIGIN3 result = MATRIX3::skew_symmetric(R*ORIGIN3(_v1i)) * v2w;
    SHAREDVECTORN first_row = Cq.row(0);
    first_row.segment(0, 3).set_zero();
    first_row.segment(4, 6) = result; 

    // Jacobian of dot(Ri* _v1j, Ro*_v2) w.r.t. inboard is
    // dot(w x Ri*_v1j, Ro*_v2) = (_v2*Ro)' * (-Ri*_v1j) x w 
    // transpose((_v2 * Ro)' * skew(-Ri*_v1j)) = skew(Ri * _v1j) * Ro * _v2 
    result = MATRIX3::skew_symmetric(R*ORIGIN3(_v1j)) * v2w;
    SHAREDVECTORN second_row = Cq.row(1);
    second_row.segment(0, 3).set_zero();
    second_row.segment(4, 6) = result; 

    // Jacobian of dot(Ri * _v1i, xi + Ri*ti - xo - Ro*to) w.r.t. inboard
    // dot(wi x Ri * _v1i, r12) + dot(Ri * _v1i, \dot{xi} + wi x Ri*ti)
    // = (skew(Ri * _v1i) * r12 + skew(Ri*t1)*Ri*_v1i) * wi + 
    //      dot(Ri * _v1i)*dot{xi}
    SHAREDVECTORN third_row = Cq.row(2);
    third_row.segment(0,3) = Ri * ORIGIN3(_v1i);
    third_row.segment(4,6) = MATRIX3::skew_symmetric(Ri * ORIGIN3(_v1i)) * r12 +
                             MATRIX3::skew_symmetric(Ri * ti)*(Ri*ORIGIN3(_v1i)); 

    // Jacobian of dot(Ri * _v1j, xi + Ri*ti - xo - Ro*to) w.r.t. inboard
    // dot(wi x Ri * _v1j, r12) + dot(Ri * _v1j, \dot{xi} + wi x Ri*ti)
    // = (skew(Ri * _v1j) * r12 + skew(Ri*t1)*Ri*_v1j) * wi + 
    //      dot(Ri * _v1j)*dot{xi}
    SHAREDVECTORN fourth_row = Cq.row(3);
    fourth_row.segment(0,3) = Ri * ORIGIN3(_v1j);
    fourth_row.segment(4,6) = MATRIX3::skew_symmetric(Ri * ORIGIN3(_v1j)) * r12 +
                              MATRIX3::skew_symmetric(Ri * ti)*(Ri*ORIGIN3(_v1j)); 

    // Jacobian of dot(Ri* _ui, Ro*_uj) w.r.t. inboard is
    // dot(w x Ri*_ui, Ro*_uj) = (_uj*Ro)' * (-Ri*_ui) x w 
    // transpose((_uj * Ro)' * skew(-Ri*_ui)) = skew(Ri * _ui) * Ro * _uj 
    ORIGIN3 ujw = wPo.q * ORIGIN3(_uj);
    result = MATRIX3::skew_symmetric(R*ORIGIN3(_ui)) * ujw;
    SHAREDVECTORN last_row = Cq.row(4);
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

    // Jacobian of dot(Ri*_v1i, Ro*_v2) w.r.t. outboard is
    // dot(Ri*_v1i, Ro*_v2) = (Ri*_v1i)' * (-Ro*_v2) x w 
    // transpose((Ri * _v1i)' * skew(-Ro*_v2)) = skew(Ro * _v2) * Ri * _v1i 
    ORIGIN3 v1w = R * ORIGIN3(_v1i);
    ORIGIN3 result = MATRIX3::skew_symmetric(R*ORIGIN3(_v2)) * v1w;
    SHAREDVECTORN first_row = Cq.row(0);
    first_row.segment(0, 3).set_zero();
    first_row.segment(4, 6) = result; 

    // Jacobian of dot(Ri*_v1j, Ro*_v2) w.r.t. outboard is
    // dot(Ri*_v1j, Ro*_v2) = (Ri*_v1j)' * (-Ro*_v2) x w 
    // transpose((Ri * _v1j)' * skew(-Ro*_v2)) = skew(Ro * _v2) * Ri * _v1j 
    v1w = R * ORIGIN3(_v1j);
    result = MATRIX3::skew_symmetric(R*ORIGIN3(_v2)) * v1w;
    SHAREDVECTORN second_row = Cq.row(1);
    second_row.segment(0, 3).set_zero();
    second_row.segment(4, 6) = result; 

    // Jacobian of dot(Ri * _v1i, xi + Ri*ti - xo - Ro*to) w.r.t. outboard
    // dot(Ri * _v1i, -\dot{xo} - wo x Ro*to)
    // = (skew(-Ro * to)*Ri*_v1i*wo + (Ri*_v1i)*\dot{xo}
    SHAREDVECTORN third_row = Cq.row(2);
    third_row.segment(0, 3) = Ri*ORIGIN3(_v1i);
    third_row.segment(4, 6) = MATRIX3::skew_symmetric(-(Ro * to))*(Ri*ORIGIN3(_v1i));

    // Jacobian of dot(Ri * _v1j, xi + Ri*ti - xo - Ro*to) w.r.t. outboard
    // dot(Ri * _v1j, -\dot{xo} - wo x Ro*to)
    // = (skew(-Ro * to)*Ri*_v1j*wo + (Ri*_v1j)*\dot{xo}
    SHAREDVECTORN fourth_row = Cq.row(3);
    fourth_row.segment(0, 3) = Ri*ORIGIN3(_v1j);
    fourth_row.segment(4, 6) = MATRIX3::skew_symmetric(-(Ro * to))*(Ri*ORIGIN3(_v1j));

    // Jacobian of dot(Ri*_ui, Ro*_uj) w.r.t. outboard is
    // dot(Ri*_ui, Ro*_uj) = (Ri*_ui)' * (-Ro*_uj) x w 
    // transpose((Ri * _ui)' * skew(-Ro*_uj)) = skew(Ro * _uj) * Ri * _ui 
    ORIGIN3 uiw = R * ORIGIN3(_ui);
    result = MATRIX3::skew_symmetric(R*ORIGIN3(_uj)) * uiw;
    SHAREDVECTORN last_row = Cq.row(4);
    last_row.segment(0, 3).set_zero();
    last_row.segment(4, 6) = result; 
  }
}


