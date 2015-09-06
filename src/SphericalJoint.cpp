/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
SPHERICALJOINT::SPHERICALJOINT() : JOINT()
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1].set_zero();
  _u[eAxis2].set_zero();
  _u[eAxis3].set_zero();
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _u[eAxis3].pose = _F;

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;

  // assign the spherical joint tolerance
  SINGULAR_TOL = (REAL) 1e-2;
}

/// Determines whether two values are relatively equal
bool SPHERICALJOINT::rel_equal(REAL x, REAL y)
{
  return (std::fabs(x - y) <= EPS * std::max(std::fabs(x), std::max(std::fabs(y), (REAL) 1.0)));
}

/// Gets the global axis for this joint
/**
 * The global axis for this joint takes the orientation of the inboard link into account; thus, if the orientation
 * of the inboard link changes, then the global axis changes.
 * \sa getAxisLocal()
 * \sa setAxisLocal()
 */
VECTOR3 SPHERICALJOINT::get_axis(Axis a) const
{
  const unsigned X = 0;

  // axis one is easy 
  if (a == eAxis1)
  {
    return _u[DOF_1];
  }

  // for both axes 2 and 3 we need cos and sin of q(1)
  const REAL c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const REAL s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);

  // axis two is obtained by rotating around x axis
  if (a == eAxis2)
  {
    VECTOR3 u(0,c1,s1);
    u.pose = _u[DOF_1].pose;
    return u;
  }
    
  // axis 3, requires the rotation matrix induced by the axis-angle
  // representation for axis 2..  much simpler to just use the rotation matrix from the
  // universal joint induced transform
  const REAL c2 = std::cos(q[DOF_2]+_q_tare[DOF_2]);
  const REAL s2 = std::sin(q[DOF_2]+_q_tare[DOF_2]);
  assert (a == eAxis3);
  VECTOR3 u(s2, -c2*s1, c1*c2);
  u.pose = _u[DOF_1].pose;
  return u;
}

/// Sets the local axis of translation for this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void SPHERICALJOINT::set_axis(const VECTOR3& axis, Axis a) 
{ 
  // normalize the axis in case the caller did not
  VECTOR3 naxis = VECTOR3::normalize(axis);
  _u[a] = POSE3::transform_vector(get_pose(), naxis); 
  update_spatial_axes(); 
}        

/// Updates the spatial axis for this joint
void SPHERICALJOINT::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // call the parent method 
  JOINT::update_spatial_axes();

  // assign axes, if possible
  if (!assign_axes())
    return;

  // set first spatial axis 
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // update the spatial axes in joint pose 
  _s[DOF_1].set_angular(_u[DOF_1]);
  _s[DOF_1].set_linear(ZEROS_3);

  // update the spatial axes in link coordinates
  _s_dot[DOF_1].set_angular(ZEROS_3);
  _s_dot[DOF_1].set_linear(ZEROS_3);

  // verify that all are unit-length and they are orthogonal
  assert(std::fabs(_u[eAxis1].norm() - (REAL) 1.0) < EPS);
  assert(std::fabs(_u[eAxis2].norm() - (REAL) 1.0) < EPS);
  assert(std::fabs(_u[eAxis3].norm() - (REAL) 1.0) < EPS);
  assert(std::fabs(VECTOR3::dot(_u[eAxis1], _u[eAxis2])) < EPS);
  assert(std::fabs(VECTOR3::dot(_u[eAxis1], _u[eAxis3])) < EPS);
  assert(std::fabs(VECTOR3::dot(_u[eAxis2], _u[eAxis3])) < EPS);
}

/// Attempts to compute unassigned axes
bool SPHERICALJOINT::assign_axes()
{
  // compute the rotation matrix necessary to transform the first axis to the
  // z-axis; can only compute this matrix if all three axes are set
  if (_u[eAxis1].norm() < EPS)
  {
    // check whether axis2 and/or 3 are set
    if (_u[eAxis2].norm() < EPS)
    {
      // axis2 is not set; if axis3 is not set, we can do nothing...
      if (_u[eAxis3].norm() < EPS)
        return false;

      // otherwise, set axes1 and 2 from 3
      _u[eAxis3].normalize();
      VECTOR3::determine_orthonormal_basis(_u[eAxis3], _u[eAxis1], _u[eAxis2]);
    }
    else
    {
      // check for axis2 set and axes 1 and 3 not set
      if (_u[eAxis3].norm() < EPS)
      {
        // axis3 is not set; set axes1/3 from 2
        _u[eAxis2].normalize();
        VECTOR3::determine_orthonormal_basis(_u[eAxis2], _u[eAxis3], _u[eAxis1]);
      }
      else
      {
        // axis2 and axis3 are set; set axis1
        _u[eAxis2].normalize();
        _u[eAxis3].normalize();
        _u[eAxis1] = VECTOR3::cross(_u[eAxis2], _u[eAxis3]);      
      }
    }
  }
  else
  {
    // check whether axis2 and/or 3 are set
    if (_u[eAxis2].norm() < EPS)
    {
      // axis2 is not set; if axis3 is not set, set them arbitrarily
      if (_u[eAxis3].norm() < EPS)
      {
        _u[eAxis1].normalize();
        VECTOR3::determine_orthonormal_basis(_u[eAxis1], _u[eAxis2], _u[eAxis3]);
      }
      else
      {
        // axis3 is set; just need to set axis2
        _u[eAxis1].normalize();
        _u[eAxis3].normalize();
        _u[eAxis2] = VECTOR3::cross(_u[eAxis3], _u[eAxis1]);
      }
    }
    else
    {
      // axis1 and 2 are set; see whether we need to set axis3
      if (_u[eAxis3].norm() < EPS)
      {
        _u[eAxis1].normalize();        
        _u[eAxis2].normalize();        
        _u[eAxis3] = VECTOR3::cross(_u[eAxis1], _u[eAxis2]);
      }
    }
  }

  return true;
} 

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& SPHERICALJOINT::get_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;
  const VECTOR3 ZEROS_3((REAL) 0.0, (REAL) 0.0, (REAL) 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("SPHERICALJOINT::get_spatial_axes() called with NULL inboard pose");
  if (!outboard)
    throw std::runtime_error("SPHERICALJOINT::get_spatial_axes() called with NULL outboard pose");

  // get current values of q
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // get the transformed second axis
  MATRIX3 R1 = AANGLE(_u[eAxis1], q[DOF_1] + q_tare[DOF_1]);
  VECTOR3 u2(R1 * ORIGIN3(_u[eAxis2]), get_pose());

  // get the transformed third axis
  MATRIX3 R2 = AANGLE(_u[eAxis2], q[DOF_2] + q_tare[DOF_2]);
  VECTOR3 u3(R1 * (R2 * ORIGIN3(_u[eAxis3])), get_pose());

  // update the spatial axis in link coordinates
  _s[DOF_2].set_angular(u2);
  _s[DOF_2].set_linear(ZEROS_3);
  _s[DOF_3].set_angular(u3);
  _s[DOF_3].set_linear(ZEROS_3);

  // use the JOINT function to do the rest
  return JOINT::get_spatial_axes();
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& SPHERICALJOINT::get_spatial_axes_dot()
{
  const VECTOR3 ZEROS_3((REAL) 0.0, (REAL) 0.0, (REAL) 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("SPHERICALJOINT::get_spatial_axes_dot() called with NULL inboard pose");
  if (!outboard)
    throw std::runtime_error("SPHERICALJOINT::get_spatial_axes_dot() called with NULL outboard pose");

  // get q, _q_tare, and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;
  const VECTORN& qd = this->qd;

  // get the two transformed axes
  REAL c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  REAL c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  REAL s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  REAL s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  REAL qd1 = qd[DOF_1];
  REAL qd2 = qd[DOF_2];

  // get the transformed second axis time derivative. This is:
  // R1 * \dot{u2} + \dot{R} * u2
  // note that \dot{u2} is zero
  MATRIX3 R1 = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 omega1 = _u[DOF_1] * qd1;
  ORIGIN3 u2 = ORIGIN3::cross(ORIGIN3(omega1), R1 * ORIGIN3(_u[DOF_2]));
  _s_dot[DOF_2].set_angular(VECTOR3(u2, get_pose()));
  _s_dot[DOF_2].set_linear(ZEROS_3);

  // get the transformed third axis time derivative. This is:
  // R1 * R2 * \dot{u3} + \dot{R1} * R2 * u3 + R1 * \dot{R2} * u3
  // note that \dot{u3} is zero
  // compute R2 and omega2 
  MATRIX3 R2 = AANGLE(_u[eAxis2], q[DOF_2] + q_tare[DOF_2]);
  VECTOR3 omega2 = _u[DOF_2] * qd2;
  ORIGIN3 u3 = ORIGIN3::cross(ORIGIN3(omega1), R1 * R2 * ORIGIN3(_u[DOF_3])) +
               R1 * (ORIGIN3::cross(ORIGIN3(omega2), R2 * ORIGIN3(_u[DOF_3]))); 
  _s_dot[DOF_3].set_angular(VECTOR3(u3, get_pose()));
  _s_dot[DOF_3].set_linear(ZEROS_3);

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void SPHERICALJOINT::determine_q(VECTORN& q)
{
  std::cerr << "SPHERICALJOINT::determine_q() warning- determine_q(.) is not currently functional" << std::endl;
}

/// Gets the (local) rotation induced by this joint
MATRIX3 SPHERICALJOINT::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q, _q_tare
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // get the transformed second axis
  MATRIX3 R1 = AANGLE(_u[eAxis1], q[DOF_1] + q_tare[DOF_1]);
  VECTOR3 u2(R1 * ORIGIN3(_u[eAxis2]), get_pose());

  // get the transformed third axis
  MATRIX3 R2 = AANGLE(_u[eAxis2], q[DOF_2] + q_tare[DOF_2]);
  VECTOR3 u3(R1 * (R2 * ORIGIN3(_u[eAxis3])), get_pose());

  // return the rotation
  return R1 * R2 * AANGLE(_u[eAxis3], q[DOF_3] + q_tare[DOF_3]);
}

/// Gets the (local) transform for this joint
shared_ptr<const POSE3> SPHERICALJOINT::get_induced_pose()
{
  // note that translation is zero by default 
  _Fprime->q = get_rotation();
  return _Fprime;
}

/// Evaluates the constraint equations
void SPHERICALJOINT::evaluate_constraints(REAL C[])
{
  const unsigned X = 0, Y = 1, Z = 2;
  const shared_ptr<const POSE3> GLOBAL;

  // This code was developed using [Shabana, 2003], pp. 430-431; variable
  // names have been altered, however

  // get the attachment frames in the global frame
  TRANSFORM3 wFi = POSE3::calc_relative_pose(_F, GLOBAL);
  TRANSFORM3 wFo = POSE3::calc_relative_pose(_Fb, GLOBAL);

  // determine the global positions of the attachment points and subtract them
  ORIGIN3 r12 = wFi.x - wFo.x; 

  // copy values
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
}

/// Computes the constraint jacobian with respect to a body
void SPHERICALJOINT::calc_constraint_jacobian(bool inboard, MATRIXN& Cq)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;
  MATRIXN tmp;

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
    MATRIX3 R = wPi.q;
    ORIGIN3 Ru = R*u;

    // get positional components of Cq
    SHAREDMATRIXN Cq_trans = Cq.block(0, 3, 0, 3);
    Cq_trans.set_identity();
    SHAREDMATRIXN Cq_rot = Cq.block(0, 3, 3, 6);
    Cq_rot = MATRIX3::skew_symmetric(-Ru);
  }
  else
  {
    // get the vector from the outboard pose to the joint pose 
    ORIGIN3 u = ORIGIN3(POSE3::transform_point(Po, joint_pos));

    // get the information necessary to compute the constraint equations
    MATRIX3 R = wPi.q;
    ORIGIN3 Ru = R*u;

    // get positional components of Cq
    SHAREDMATRIXN Cq_trans = Cq.block(0, 3, 0, 3);
    Cq_trans.set_identity();
    Cq_trans *= -1.0;
    SHAREDMATRIXN Cq_rot = Cq.block(0, 3, 3, 6);
    Cq_rot = MATRIX3::skew_symmetric(Ru);
  }

  // transform the Jacobian as necessary
  if (transform_jacobian(Cq, inboard, tmp))
    Cq = tmp;
}

/// Computes the time derivative of the constraint jacobian with respect to a body
/**
 * TODO: implement this
 */
void SPHERICALJOINT::calc_constraint_jacobian_dot(bool inboard, MATRIXN& Cq)
{
  throw std::runtime_error("Implementation required");
}

