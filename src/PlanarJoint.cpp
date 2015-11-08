/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
PLANARJOINT::PLANARJOINT() : JOINT()
{
  const unsigned SPATIAL_DIM = 6;

  // set the offset to zero by default
  _offset = 0.0;

  // init joint data
  init_data();

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
  {
    _s_dot[i].set_zero();
    _s_dot[i].pose = _F;
  }
}

/// Sets the normal vector to the plane
void PLANARJOINT::set_normal(const VECTOR3& normal)
{
  _normal = normal;
  update_spatial_axes();
}

/// Updates the offset value to the plane
void PLANARJOINT::update_offset()
{
  // setup the constraint vector
  REAL C[3];

  // evaluate the constraints
  evaluate_constraints(C);

  // set the offset
  _offset = C[0];
}

/// Updates the spatial axis for this joint
void PLANARJOINT::update_spatial_axes()
{
  const unsigned X = 0, Y = 1, Z = 2;
  shared_ptr<const POSE3> GLOBAL_3D;

  // call the parent method 
  JOINT::update_spatial_axes();

  // set first spatial axis 
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get two vectors orthogonal to the normal vector in the outboard link's
  // frame
  VECTOR3::determine_orthonormal_basis(_normal, _vi, _vj);
  _vi = POSE3::transform_vector(_Fb, _vi);
  _vj = POSE3::transform_vector(_Fb, _vj);

  // get the normal and tangent directions in the joint pose
  TRANSFORM3 jT0 = POSE3::calc_relative_pose(GLOBAL_3D, get_pose());
  VECTOR3 normal = jT0.transform_vector(_normal);
  VECTOR3 tan1 = jT0.transform_vector(_tan1);
  VECTOR3 tan2 = jT0.transform_vector(_tan2);

  // update the spatial axes in joint pose 
  _s[DOF_1].set_angular(ZEROS_3);
  _s[DOF_1].set_linear(tan1);
  _s[DOF_2].set_angular(ZEROS_3);
  _s[DOF_2].set_linear(tan2);
  _s[DOF_3].set_angular(normal);
  _s[DOF_3].set_linear(ZEROS_3);

  // update the offset
  update_offset();
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void PLANARJOINT::determine_q(VECTORN& q)
{
  std::cerr << "PLANARJOINT::determine_q() warning- determine_q(.) is not currently functional" << std::endl;
}

/// Gets the (local) transform for this joint
shared_ptr<const POSE3> PLANARJOINT::get_induced_pose()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q, _q_tare
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // setup the translation
  VECTOR3 xlat = _tan1*(q[DOF_1] + q_tare[DOF_1]) +
                 _tan2*(q[DOF_2] + q_tare[DOF_2]) + _normal*_offset;

  // setup the rotation 
  MATRIX3 R = AANGLE(_normal, q[DOF_3] + q_tare[DOF_3]);

  // setup the pose
  _Fprime->x = ORIGIN3(xlat);
  _Fprime->q = R;

  return _Fprime;
}

/// Evaluates the constraint equations
void PLANARJOINT::evaluate_constraints(REAL C[])
{
  const unsigned X = 0, Y = 1, Z = 2;
  const shared_ptr<const POSE3> GLOBAL_3D;

  // This code was developed using [Shabana, 2003], pp. 430-431; variable
  // names have been altered, however

  // get the attachment frames in the global frame
  TRANSFORM3 wPo = POSE3::calc_relative_pose(_Fb, GLOBAL_3D);

  // equation 1: translational deviation of outboard from plane
  C[0] = wPo.x.dot(ORIGIN3(_normal)) - _offset;

  // get two orthogonal vectors in the global frame
  // frame
  VECTOR3 vi0 = POSE3::transform_vector(GLOBAL_3D, _vi);
  VECTOR3 vj0 = POSE3::transform_vector(GLOBAL_3D, _vj);

  // equation 2: pitch of outboard
  C[1] = _normal.dot(vi0);

  // equation 3: roll of outboard
  C[2] = _normal.dot(vj0);
}

/// Computes the constraint jacobian with respect to a body
void PLANARJOINT::calc_constraint_jacobian(bool inboard, MATRIXN& Cq)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;
  MATRIXN tmp;

  // resize the matrix
  Cq.resize(num_constraint_eqns(), SPATIAL_DIM);

  // inboard will always have no Jacobian (it should be fixed)
  assert(!get_inboard_link()->is_enabled());

  // get the attachment frame in the global frame
  TRANSFORM3 wPo = POSE3::calc_relative_pose(_Fb, GLOBAL);

  // setup the constraint equations (from Shabana, p. 432)
  if (inboard)
  {
    Cq.set_zero();
    return;
  }
  else
  {
    // get the information necessary to compute the constraint equations
    MATRIX3 R = wPo.q;

    // time derivative of wPo.x.dot(_normal) - _offset:
    // identity * _normal
    SHAREDVECTORN first_row = Cq.row(0);
    first_row.segment(0,3) = _normal;
    first_row.segment(3,6).set_zero();  

    // Jacobian of dot(normal, R*_vi) w.r.t. outboard is
    // normal' * (-R*_vi) x w = normal' * skew(-R*_vi)
    ORIGIN3 result = MATRIX3::skew_symmetric(R * ORIGIN3(_vi)) * ORIGIN3(_normal);
    SHAREDVECTORN second_to_last_row = Cq.row(1);
    second_to_last_row.segment(0,3).set_zero();
    second_to_last_row.segment(3,6) = result;

    // Jacobian of dot(normal, R*_vj) w.r.t. outboard is
    // normal' * (-R*_vj) x w = normal' * skew(-R*_vj)
    result = MATRIX3::skew_symmetric(R * ORIGIN3(_vj)) * ORIGIN3(_normal);
    SHAREDVECTORN last_row = Cq.row(2);
    last_row.segment(0,3).set_zero();
    last_row.segment(3,6) = result;
  }

  // transform the Jacobian as necessary
  if (transform_jacobian(Cq, inboard, tmp))
    Cq = tmp;
}

/// Computes the time derivative of the constraint jacobian with respect to a body
/**
 * TODO: implement this
 */
void PLANARJOINT::calc_constraint_jacobian_dot(bool inboard, MATRIXN& Cq)
{
  throw std::runtime_error("Implementation required");
}

