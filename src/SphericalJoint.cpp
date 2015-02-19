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
    return _u[0];
  }

  // for both axes 2 and 3 we need cos and sin of q(1)
  const REAL c1 = std::cos(q[DOF_1]+_q_tare[DOF_1]);
  const REAL s1 = std::sin(q[DOF_1]+_q_tare[DOF_1]);

  // axis two is obtained by rotating around x axis
  if (a == eAxis2)
  {
    VECTOR3 u(0,c1,s1);
    u.pose = _u[0].pose;
    return u;
  }
    
  // axis 3, requires the rotation matrix induced by the axis-angle
  // representation for axis 2..  much simpler to just use the rotation matrix from the
  // universal joint induced transform
  const REAL c2 = std::cos(q[DOF_2]+_q_tare[DOF_2]);
  const REAL s2 = std::sin(q[DOF_2]+_q_tare[DOF_2]);
  assert (a == eAxis3);
  VECTOR3 u(s2, -c2*s1, c1*c2);
  u.pose = _u[0].pose;
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

  // get the set of spatial axes
  REAL c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  REAL c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  REAL s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  REAL s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);

  // form untransformed spatial axes -- this are the vectors describing each axis, after
  // rotation by preceding axis/axes; note that first axis always points toward 1,0,0
  VECTOR3 uu2(0, c1, s1);
  VECTOR3 uu3(s2, -c2*s1, c1*c2);

  // transform the spatial axes into the joint frame 
  VECTOR3 u1 = _u[0];
  VECTOR3 u2 =  uu2;
  VECTOR3 u3 =  uu3;

  // setup relative poses for all three
  u1.pose = _F;
  u2.pose = _F;
  u3.pose = _F;

  // update the spatial axis in link coordinates
  _s[0].set_angular(u1);
  _s[0].set_linear(ZEROS_3);
  _s[1].set_angular(u2);
  _s[1].set_linear(ZEROS_3);
  _s[2].set_angular(u3);
  _s[2].set_linear(ZEROS_3);

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

  // form the time derivatives of the non-constant spatial axes (untransformed) 
  VECTOR3 uu2(0, -s1*qd1, c1*qd1);
  VECTOR3 uu3(c2*qd2, -c2*c1*qd1 + s2*s1*qd2, -c1*s2*qd2 - s1*c2*qd1);

  // transform the axes into outer link coordinates
  VECTOR3 u2 = uu2; 
  VECTOR3 u3 = uu3; 

  // update the spatial axis in joint coordinates; note that third column of spatial axis
  // derivative set to zero in constructor and is never modified
  _s_dot[1].set_angular(u2);
  _s_dot[1].set_lower(ZEROS_3);
  _s_dot[2].set_angular(u3);
  _s_dot[2].set_linear(ZEROS_3);

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void SPHERICALJOINT::determine_q(VECTORN& q)
{
  shared_ptr<const POSE3> GLOBAL;
  const unsigned X = 0, Y = 1, Z = 2;

  // get the outboard link
  shared_ptr<const POSE3> Fo = get_outboard_pose();

  // verify that the outboard pose is set
  if (!Fo)
    throw std::runtime_error("determine_q() called with NULL outboard pose!");

  // if any of the axes are not defined, can't use this method
  if (std::fabs(_u[0].norm_sq() - (REAL) 1.0) > EPS ||
      std::fabs(_u[1].norm_sq() - (REAL) 1.0) > EPS ||
      std::fabs(_u[2].norm_sq() - (REAL) 1.0) > EPS)
    throw UndefinedAxisException();

  // set proper size for q
  q.resize(num_dof());

  // get the poses of the joint and outboard link
  shared_ptr<const POSE3> Fj = get_pose();

  // compute transforms
  TRANSFORM3 wTo = POSE3::calc_relative_pose(Fo, GLOBAL); 
  TRANSFORM3 jTw = POSE3::calc_relative_pose(GLOBAL, Fj);
  TRANSFORM3 jTo = jTw * wTo;

  // determine the joint transformation
  MATRIX3 R = jTo.q;

  // determine cos and sin values for q1, q2,  and q3
  REAL s2 = R(X,Z);
  REAL c2 = std::cos(std::asin(s2));
  REAL s1, c1, s3, c3;
  if (std::fabs(c2) > EPS)
  {
    s1 = -R(Y,Z)/c2;
    c1 = R(Z,Z)/c2;
    s3 = -R(X,Y)/c2;
    c3 = R(X,X)/c2;
    assert(!std::isnan(s1));
    assert(!std::isnan(c1));
    assert(!std::isnan(s3));
    assert(!std::isnan(c3));
  }
  else
  {
    // singular, we can pick any value for s1, c1, s3, c3 as long as the
    // following conditions are satisfied
    // c1*s3 + s1*c3*s2 = R(Y,X)
    // c1*c3 - s1*s3*s2 = R(Y,Y)
    // s1*s3 - c1*c3*s2 = R(Z,X)
    // s1*c3 + c1*s3*s2 = R(Z,Y)
    // so, we'll set q1 to zero (arbitrarily) and obtain
    s1 = 0;
    c1 = 1;
    s3 = R(Y,X);
    c3 = R(Y,Y);
  }

  // now determine q; only q2 can be determined without ambiguity
  if (std::fabs(s1) < EPS)
    q[DOF_2] = std::atan2(R(X,Z), R(Z,Z)/c1);
  else
    q[DOF_2] = std::atan2(R(X,Z), -R(Y,Z)/s1);
  assert(!std::isnan(q[DOF_2]));

  // if cos(q2) is not singular, proceed easily from here..
  if (std::fabs(c2) > EPS)
  {
    q[DOF_1] = std::atan2(-R(Y,Z)/c2, R(Z,Z)/c2);
    q[DOF_3] = std::atan2(-R(X,Y)/c2, R(X,X)/c2);
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
  else
  {
    if (std::fabs(c1) > EPS)
      q[DOF_3] = std::atan2((R(Y,X) - s1*s2*c3)/c1, (R(Y,Y) + s1*s2*s3)/c1);
    else
      q[DOF_3] = std::atan2((R(Z,X) + c1*s2*c3)/s1, (R(Z,Y) - c1*s2*s3)/s1);
    if (std::fabs(c3) > EPS)
      q[DOF_1] = std::atan2((R(Y,X) - c1*s3)/(s2*c3), (-R(Y,X) + s1*s3)/(s2*c3));
    else
      q[DOF_1] = std::atan2((-R(Y,Y) + c1*c3)/(s2*s3), (R(Z,Y) - s1*c3)/(s2*s3));
    assert(!std::isnan(q[DOF_1]));
    assert(!std::isnan(q[DOF_3]));
  }
}

/// Gets the (local) rotation induced by this joint
MATRIX3 SPHERICALJOINT::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q, _q_tare
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // compute some needed quantities
  const REAL c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const REAL s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const REAL c2 = std::cos(q[DOF_2]+q_tare[DOF_2]);
  const REAL s2 = std::sin(q[DOF_2]+q_tare[DOF_2]);
  const REAL c3 = std::cos(q[DOF_3]+q_tare[DOF_3]);
  const REAL s3 = std::sin(q[DOF_3]+q_tare[DOF_3]);

  // determine rotation
  // this is just the rotation matrix induced by using Tait-Bryan angles
  MATRIX3 R;
  R(X,X) = c2*c3;              R(X,Y) = -c2*s3;              R(X,Z) = s2;
  R(Y,X) = s1*s2*c3 + c1*s3;   R(Y,Y) = -s1*s2*s3 + c1*c3;   R(Y,Z) = -c2*s1;
  R(Z,X) = -c1*s2*c3 + s1*s3;  R(Z,Y) = c1*s2*s3 + s1*c3;    R(Z,Z) = c2*c1;

  return R;
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
/*
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // This code was developed using [Shabana, 2003], pp. 430-431; variable
  // names have been altered, however

  // determine the global positions of the attachment points and subtract them
  VECTOR3 r1 = get_position_global(false);
  VECTOR3 r2 = get_position_global(true);
  VECTOR3 r12 = r1 - r2; 

  // copy values
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
*/
}

/*
/// Computes the constraint jacobian with respect to a body
void SPHERICALJOINT::calc_constraint_jacobian(RigidBodyPtr body, unsigned index, REAL Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (REAL) 0.0;
    return;
  }

  // setup the constraint equations (from Shabana, p. 432)
  if (body == inner)
  {
    // get the information necessary to compute the constraint equations
    const QUAT& q = inner->get_orientation();
    const VECTOR3& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
    const REAL qx = q.x;
    const REAL qy = q.y;
    const REAL qz = q.z;
    const REAL qw = q.w;
    const REAL px = p[X];
    const REAL py = p[Y];
    const REAL pz = p[Z];

    switch (index)
    {
      case 0:
        Cq[0] = 1.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*px*qw + 2*pz*qy - 2*py*qz;
        Cq[4] = 4*qx*px + 2*qy*py + 2*qz*pz;
        Cq[5] = 2*pz*qw + 2*py*qx;
        Cq[6] = 2*pz*qx - 2*py*qw;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 1.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*py*qw - 2*pz*qx + 2*px*qz;
        Cq[4] = 2*qy*px - 2*qw*pz;
        Cq[5] = 2*px*qx + 4*py*qy + 2*pz*qz;
        Cq[6] = 2*px*qw + 2*pz*qy;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 1.0;
        Cq[3] = 4*pz*qw + 2*py*qx - 2*px*qy;
        Cq[4] = 2*qz*px + 2*qw*py;
        Cq[5] = 2*py*qz - 2*px*qw;
        Cq[6] = 4*pz*qz + 2*py*qy + 2*px*qx;
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // get the information necessary to compute the constraint equations
    const QUAT& q = outer->get_orientation();
    const VECTOR3& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
    const REAL qx = q.x;
    const REAL qy = q.y;
    const REAL qz = q.z;
    const REAL qw = q.w;
    const REAL px = -p[X];
    const REAL py = -p[Y];
    const REAL pz = -p[Z];

    switch (index)
    {
      case 0:
        Cq[0] = -1.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -(4*px*qw + 2*pz*qy - 2*py*qz);
        Cq[4] = -(4*qx*px + 2*qy*py + 2*qz*pz);
        Cq[5] = -(2*pz*qw + 2*py*qx);
        Cq[6] = -(2*pz*qx - 2*py*qw);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = -1.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*py*qw - 2*pz*qx + 2*px*qz);
        Cq[4] = -(2*qy*px - 2*qw*pz);
        Cq[5] = -(2*px*qx + 4*py*qy + 2*pz*qz);
        Cq[6] = -(2*px*qw + 2*pz*qy);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = -1.0;
        Cq[3] = -(4*pz*qw + 2*py*qx - 2*px*qy);
        Cq[4] = -(2*qz*px + 2*qw*py);
        Cq[5] = -(2*py*qz - 2*px*qw);
        Cq[6] = -(4*pz*qz + 2*py*qy + 2*px*qx);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Computes the time derivative of the constraint jacobian with respect to a body
void SPHERICALJOINT::calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, REAL Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  RigidBodyPtr inner = get_inboard_link();
  RigidBodyPtr outer = get_outboard_link();

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (REAL) 0.0;
    return;
  }

  // setup the constraint equations (from Shabana, p. 432)
  if (body == inner)
  {
    // get the information necessary to compute the constraint equations
    const QUAT& q = inner->get_orientation();
    const VECTOR3& p = inner->get_outer_joint_data(outer).com_to_joint_vec;
    const REAL px = p[X];
    const REAL py = p[Y];
    const REAL pz = p[Z];
    QUAT qd = QUAT::deriv(q, inner->get_avel());
    const REAL dqw = qd.w;
    const REAL dqx = qd.x;
    const REAL dqy = qd.y;
    const REAL dqz = qd.z;

    switch (index)
    {
      case 0:
        Cq[0] = 0.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*px*dqw + 2*pz*dqy - 2*py*dqz;
        Cq[4] = 4*dqx*px + 2*dqy*py + 2*dqz*pz;
        Cq[5] = 2*pz*dqw + 2*py*dqx;
        Cq[6] = 2*pz*dqx - 2*py*dqw;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*py*dqw - 2*pz*dqx + 2*px*dqz;
        Cq[4] = 2*dqy*px - 2*dqw*pz;
        Cq[5] = 2*px*dqx + 4*py*dqy + 2*pz*dqz;
        Cq[6] = 2*px*dqw + 2*pz*dqy;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = 4*pz*dqw + 2*py*dqx - 2*px*dqy;
        Cq[4] = 2*dqz*px + 2*dqw*py;
        Cq[5] = 2*py*dqz - 2*px*dqw;
        Cq[6] = 4*pz*dqz + 2*py*dqy + 2*px*dqx;
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    // get the information necessary to compute the constraint equations
    const QUAT& q = outer->get_orientation();
    const VECTOR3& p = body->get_inner_joint_data(inner).joint_to_com_vec_of;
    const REAL px = -p[X];
    const REAL py = -p[Y];
    const REAL pz = -p[Z];
    QUAT qd = QUAT::deriv(q, outer->get_avel());
    const REAL dqw = qd.w;
    const REAL dqx = qd.x;
    const REAL dqy = qd.y;
    const REAL dqz = qd.z;

    switch (index)
    {
      case 0:
        Cq[0] = 0.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -(4*px*dqw + 2*pz*dqy - 2*py*dqz);
        Cq[4] = -(4*dqx*px + 2*dqy*py + 2*dqz*pz);
        Cq[5] = -(2*pz*dqw + 2*py*dqx);
        Cq[6] = -(2*pz*dqx - 2*py*dqw);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = 0.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*py*dqw - 2*pz*dqx + 2*px*dqz);
        Cq[4] = -(2*dqy*px - 2*dqw*pz);
        Cq[5] = -(2*px*dqx + 4*py*dqy + 2*pz*dqz);
        Cq[6] = -(2*px*dqw + 2*pz*dqy);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = -(4*pz*dqw + 2*py*dqx - 2*px*dqy);
        Cq[4] = -(2*dqz*px + 2*dqw*py);
        Cq[5] = -(2*py*dqz - 2*px*dqw);
        Cq[6] = -(4*pz*dqz + 2*py*dqy + 2*px*dqx);
        break;

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}
*/

