/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

/// Initializes the joint
/**
 * The inboard and outboard links are set to NULL.
 */
FIXEDJOINT::FIXEDJOINT() : JOINT()
{
  // init joint data
  init_data();

  // setup the spatial axis time derivative to zero
  _s_dot.clear();

  // initialize frames
  _F1 = shared_ptr<POSE3>(new POSE3);
  _F2 = shared_ptr<POSE3>(new POSE3);
}

/// Sets spatial axes to zero 
void FIXEDJOINT::update_spatial_axes()
{
  // call parent method
  JOINT::update_spatial_axes();
}

/// Setup the transform from one joint to another
void FIXEDJOINT::setup_joint()
{
  const unsigned X = 0, Y = 1, Z = 2;
  shared_ptr<const POSE3> GLOBAL;

  // get the two poses 
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard || !outboard)
    return;

  // get the rotation matrices
  MATRIX3 Ri = inboard->q;
  MATRIX3 Ro = outboard->q;

  // TODO: fix this
  // compute the relative transform
/*
  MATRIX3 Rrel = MATRIX3::transpose(Ro) * Ri;
  _T->q = Rrel;
  _T->x.set_zero();
  _T->source = ;

  // compute the vector from the inner link to the outer link in inner link
  // frame
  VECTOR3 ox(outboard->x, To);
//  _ui = inboard->inverse_transform(ox);
*/
  // compute the constant orientation term
  _rconst[X] = VECTOR3(Ri.get_row(X), GLOBAL).dot(VECTOR3(Ro.get_row(X), GLOBAL));
  _rconst[Y] = VECTOR3(Ri.get_row(Y), GLOBAL).dot(VECTOR3(Ro.get_row(Y), GLOBAL));
  _rconst[Z] = VECTOR3(Ri.get_row(Z), GLOBAL).dot(VECTOR3(Ro.get_row(Z), GLOBAL));
}

/// Sets the inboard link
void FIXEDJOINT::set_inboard_pose(shared_ptr<const POSE3> inboard_pose, bool update_joint_pose)
{
  // call parent method since it does all of the work
  JOINT::set_inboard_pose(inboard_pose, update_joint_pose);
  setup_joint();
}

/// Sets the outboard link
void FIXEDJOINT::set_outboard_pose(shared_ptr<POSE3> outboard_pose, bool update_joint_pose)
{
  // call parent method since it does all of the work
  JOINT::set_outboard_pose(outboard_pose, update_joint_pose);
  setup_joint();
}

/// Gets the (local) transform for this joint (constant)
shared_ptr<const POSE3> FIXEDJOINT::get_induced_pose()
{
  // get the link transform (always identity)
  return _Fprime;
}

/// Gets the derivative for the spatial axes for this joint
const vector<SVELOCITY>& FIXEDJOINT::get_spatial_axes_dot()
{
  return _s_dot;
}

/// Computes the constraint Jacobian
void FIXEDJOINT::calc_constraint_jacobian(bool inboard, SHAREDMATRIXN& Cq)
{
  // set the index
  Cq.set_identity();
  if (!inboard)
    Cq *= -1.0;
}

/// Computes the constraint Jacobian
void FIXEDJOINT::calc_constraint_jacobian_dot(bool inboard, SHAREDMATRIXN& Cq)
{
  // set the index
  Cq.set_zero();
}

/// Evaluates the constraint equations
void FIXEDJOINT::evaluate_constraints(REAL C[])
{
  const unsigned X = 0, Y = 1, Z = 2;
  const shared_ptr<const POSE3> GLOBAL;

  // get ui in the global frame
  shared_ptr<const POSE3> P1 = get_inboard_pose();
  VECTOR3 ui0 = P1->transform_vector(_ui);

  // get the poses for the two links
  POSE3 wP1 = *P1;
  POSE3 wP2 = *get_outboard_pose();

  // get the poses relative to the global frame
  wP1.update_relative_pose(GLOBAL);
  wP2.update_relative_pose(GLOBAL);

  // get the transforms and orientations for the two links
  MATRIX3 R1 = wP1.q;
  MATRIX3 R2 = wP2.q;

  // evaluate the relative position
  ORIGIN3 rpos = ORIGIN3(ui0) + wP1.x - wP2.x; 

  const REAL XX = VECTOR3(R1.get_row(X), GLOBAL).dot(VECTOR3(R2.get_row(X), GLOBAL));
  const REAL YY = VECTOR3(R1.get_row(Y), GLOBAL).dot(VECTOR3(R2.get_row(Y), GLOBAL));
  const REAL ZZ = VECTOR3(R1.get_row(Z), GLOBAL).dot(VECTOR3(R2.get_row(Z), GLOBAL));

  // setup C
  C[0] = rpos[0];
  C[1] = rpos[1];
  C[2] = rpos[2];
  C[3] = XX - _rconst[X];
  C[4] = YY - _rconst[Y];
  C[5] = ZZ - _rconst[Z];
}


