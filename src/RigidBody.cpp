/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::cerr;
using std::endl;
using std::map;
using std::list;
using std::queue;

/// Default constructor
/**
 * Constructs a rigid body with zero mass, zero inertia tensor, and center
 * of mass at [0,0,0] with position at [0,0,0], identity orientation, and zero
 * linear and angular velocity.  Body is enabled by default.
 */
RIGIDBODY::RIGIDBODY()
{
  const unsigned SPATIAL_DIM = 6;

  // setup reference pose
  _F = shared_ptr<POSE3>(new POSE3(POSE3::identity()));
  _Ji.pose = _xdi.pose = _xddi.pose = _forcei.pose = _F;

  // setup inertial pose
  _jF = shared_ptr<POSE3>(new POSE3(POSE3::identity()));
  _jF->rpose = _F;
  _Jm.pose = _xdm.pose = _xddm.pose = _forcem.pose = _jF;

  // setup c.o.m. frame link pose
  _F2 = shared_ptr<POSE3>(new POSE3);
  _Jcom.pose = _xdcom.pose = _xddcom.pose = _forcecom.pose = _F2;

  // invalidate everything
  _forcei_valid = false;
  _forcej_valid = false;
  _forcecom_valid = false;
  _force0_valid = false;
  _xdi_valid = false;
  _xdj_valid = false;
  _xdcom_valid = false;
  _xd0_valid = false;
  _xddi_valid = false;
  _xddj_valid = false;
  _xddcom_valid = false;
  _xdd0_valid = false;
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // use link c.o.m. frame by default
  _rftype = eLinkCOM;

  // set everything else
  _enabled = true;
  _link_idx = std::numeric_limits<unsigned>::max();

  // setup the default limit bound expansion
  compliance = eRigid;
}

/// Gets the frame in which kinematics and dynamics computations occur
shared_ptr<const POSE3> RIGIDBODY::get_computation_frame() const
{
  switch (_rftype)
  {
    case eLink:
      return _F;

    case eLinkInertia:
      return _jF;

    case eLinkCOM:
      return _F2;

    case eGlobal:
      return shared_ptr<const POSE3>();

    case eJoint:
      return (_abody.expired() || is_base()) ? _F : get_inner_joint_explicit()->get_pose();

    default:
      assert(false);
  }

  return shared_ptr<const POSE3>();
}

/// Rotates the rigid body
void RIGIDBODY::rotate(const QUAT& q)
{
  const shared_ptr<const POSE3> GLOBAL;

  // save the current relative pose
  shared_ptr<const POSE3> Frel = _F->rpose;

  // update the rotation
  _F->update_relative_pose(GLOBAL);
  _F->q *= q;
  _F->update_relative_pose(Frel);

  // update the mixed pose
  update_mixed_pose();

  // invalidate vector quantities
  _forcei_valid = _forcecom_valid = _force0_valid = false;
  _xdi_valid = _xdcom_valid = _xd0_valid = false;
  _xddi_valid = _xddcom_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<shared_ptr<RIGIDBODY> > outer;
  BOOST_FOREACH(shared_ptr<JOINT> j, _outer_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<shared_ptr<RIGIDBODY> >::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<shared_ptr<RIGIDBODY> >::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// Gets the time derivative of the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose
MATRIXN& RIGIDBODY::calc_jacobian_dot(shared_ptr<const POSE3> source_pose, shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMIC_BODY> body, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;

  if (body != shared_from_this())
    throw std::runtime_error("RigidBody::calc_jacobian_dot() called with wrong body!");

  if (!is_enabled())
  {
    J.set_zero(SPATIAL_DIM, 0);
    return J;
  }

  // construct the spatial transform
  POSE3::dot_spatial_transform_to_matrix2(source_pose, target_pose, J);

  return J;
}

/// Gets the time derivative of the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose
MATRIXN& RIGIDBODY::calc_jacobian(shared_ptr<const POSE3> source_pose, shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMIC_BODY> body, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;
  const shared_ptr<const POSE3> GLOBAL;

  if (body != shared_from_this())
    throw std::runtime_error("RigidBody::calc_jacobian() called with wrong body!");

  // if the body is disabled, do not compute a Jacobian
  if (!is_enabled())
  {
    J.set_zero(6, 0);
    return J;
  }

  // construct the spatial transform
  POSE3::spatial_transform_to_matrix2(source_pose, target_pose, J);

  FILE_LOG(LOG_DYNAMICS) << "RigidBody::calc_jacobian() entered" << std::endl;
  if (LOGGING(LOG_DYNAMICS))
  {
    POSE3 P;
    if (target_pose)
    {
      P.rpose = target_pose;
      P.update_relative_pose(GLOBAL);
      FILE_LOG(LOG_DYNAMICS) << "  pose: " << P << std::endl;
    }
    else
      FILE_LOG(LOG_DYNAMICS) << "  pose: " << GLOBAL << std::endl;
  } 

  return J;
}

/// Translates the rigid body
void RIGIDBODY::translate(const ORIGIN3& x)
{
  const shared_ptr<const POSE3> GLOBAL;

  // save the current relative pose
  shared_ptr<const POSE3> Frel = _F->rpose;

  // update the translation
  _F->update_relative_pose(GLOBAL);
  _F->x += x;
  _F->update_relative_pose(Frel);

  // update the mixed pose
  update_mixed_pose();

  // invalidate vector quantities
  _forcei_valid = _forcecom_valid = _force0_valid = false;
  _xdi_valid = _xdcom_valid = _xd0_valid = false;
  _xddi_valid = _xddcom_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<shared_ptr<RIGIDBODY> > outer;
  BOOST_FOREACH(shared_ptr<JOINT> j, _outer_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<shared_ptr<RIGIDBODY> >::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<shared_ptr<RIGIDBODY> >::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// (Re)sets the computation frame
void RIGIDBODY::set_computation_frame_type(ReferenceFrameType rftype)
{
  // correct rftype if necessary
  if (_abody.expired() && rftype == eJoint)
    rftype = eLink;

  // store the new reference frame type
  _rftype = rftype;
}

/// Computes the forward dynamics for this body
void RIGIDBODY::calc_fwd_dyn()
{
  // if the body is free, just compute linear and angular acceleration via
  // Newton's and Euler's laws
  if (_abody.expired())
  {
    // don't do anything if the body is enabled
    if (!is_enabled())
      return;

    // make sure that the inertia is reasonable
    const SPATIAL_RB_INERTIA& J = get_inertia();
    #ifndef NDEBUG
    if (J.m <= 0.0 || J.J.norm_inf() <= 0.0)
      throw std::runtime_error("Tried to calculate forward dynamics on body with zero mass/inertia");
    #endif

    // otherwise, calculate forward dynamics
    SFORCE f = sum_forces() - calc_euler_torques();
    SACCEL xdd = J.inverse_mult(f);

    // get current velocity
    SVELOCITY xd = get_velocity();

    FILE_LOG(LOG_DYNAMICS) << "Dynamics: " << POSE3::transform(_F2, xdd, _xdcom) << std::endl;
    // set the acceleration
    switch (_rftype)
    {
      case eGlobal:
        _xdd0 = xdd;
        _xddm = POSE3::transform(_jF, xdd, xd);
        _xddi_valid = _xddj_valid = _xddcom_valid = false;
        break;

      case eLinkCOM:
        _xddcom = xdd;
        _xddcom_valid = true;
        _xddm = POSE3::transform(_jF, xdd, xd);
        _xddi_valid = _xddj_valid = _xdd0_valid = false;
        break;

      case eLink:
        _xddi = xdd;
        _xddm = POSE3::transform(_jF, xdd, xd);
        _xddi_valid = true;
        _xddj_valid = _xddcom_valid = _xdd0_valid = false;
        break;

      case eLinkInertia:
        _xddm = xdd;
        _xddcom_valid = _xddi_valid = _xddj_valid = _xdd0_valid = false;
        break;

      case eJoint:
        _xddj = xdd;
        _xddcom = POSE3::transform(_jF, xdd, xd);
        _xddj_valid = true;
        _xddi_valid = _xddcom_valid = _xdd0_valid = false;
        break;

      default:
        assert(false);
    }
  }
  else
  {
    // otherwise, need to call forward dynamics on the articulated body
    shared_ptr<ARTICULATED_BODY> abody(_abody);

    // calculate forward dynamics on it
    abody->calc_fwd_dyn();
  }
}

/// Sets the body to enabled / disabled.
/**
 * If the body is disabled, the linear and angular velocity are set to zero,
 * and the body will not be updated if it is attempted to integrate its
 * equations of motion.
 */
void RIGIDBODY::set_enabled(bool flag)
{
  // mark as enabled / disabled
  _enabled = flag;

  // if disabled, then zero the velocities and accelerations
  if (!_enabled)
  {
    _xdcom.set_zero();
    _xddcom.set_zero();
    _xd0.set_zero();
    _xdd0.set_zero();
    _xdi.set_zero();
    _xddi.set_zero();
    _xdj.set_zero();
    _xddj.set_zero();
    _xdm.set_zero();
    _xddm.set_zero();
    _xdi_valid = _xdj_valid = _xdcom_valid = _xd0_valid = true;
    _xddi_valid = _xddj_valid = _xddcom_valid = _xdd0_valid = true;
  }
}

/// Sets the velocity of this body
void RIGIDBODY::set_velocity(const SVELOCITY& xd)
{
  const shared_ptr<const POSE3> GLOBAL;

  // set the velocity
  _xdm = POSE3::transform(dynamic_pointer_cast<const POSE3>(_jF), xd);

  // invalidate the remaining velocities
  _xdi_valid = _xdj_valid = _xdcom_valid = _xd0_valid = false;

  // see whether we can re-validate a velocity
  if (xd.pose == _F)
  {
    _xdi_valid = true;
    _xdi = xd;
  }
  else if (xd.pose == _F2)
  {
    _xdcom_valid = true;
    _xdcom = xd;
  }
  else if (!is_base() && xd.pose == get_inner_joint_explicit()->get_pose())
  {
    _xdj_valid = true;
    _xdj = xd;
  }
  else if (xd.pose == GLOBAL)
  {
    _xd0_valid = true;
    _xd0 = xd;
  }
}

/// Sets the acceleration of this body
void RIGIDBODY::set_accel(const SACCEL& xdd)
{
  const shared_ptr<const POSE3> GLOBAL;

  // get the velocity
  SVELOCITY xd = POSE3::transform(xdd.pose, _xdm);

  // set the acceleration
  _xddm = POSE3::transform(_jF, xdd, xd);

  // invalidate the remaining accelerations
  _xddi_valid = _xddj_valid = _xddcom_valid = _xdd0_valid = false;

  // see whether we can re-validate an acceleration
  if (xdd.pose == _F)
  {
    _xddi_valid = true;
    _xddi = xdd;
  }
  else if (xdd.pose == _F2)
  {
    _xddcom_valid = true;
    _xddcom = xdd;
  }
  else if (!is_base() && xdd.pose == get_inner_joint_explicit()->get_pose())
  {
    _xddj_valid = true;
    _xddj = xdd;
  }
  else if (xdd.pose == GLOBAL)
  {
    _xdd0_valid = true;
    _xdd0 = xdd;
  }
}

/// Sets the rigid body inertia for this body
void RIGIDBODY::set_inertia(const SPATIAL_RB_INERTIA& inertia)
{
  const shared_ptr<const POSE3> GLOBAL;

  // set the inertia
  _Jm = POSE3::transform(_jF, inertia);

  // invalidate the remaining inertias
  _Ji_valid = _Jj_valid = _J0_valid = _Jcom_valid = false;

  // see whether we can re-validate an acceleration
  if (inertia.pose == _F)
  {
    _Ji_valid = true;
    _Ji = inertia;
  }
  else if (inertia.pose == _F2)
  {
    _Jcom_valid = true;
    _Jcom = inertia;
  }
  else if (!is_base() && inertia.pose == get_inner_joint_explicit()->get_pose())
  {
    _Jj_valid = true;
    _Jj = inertia;
  }
  else if (inertia.pose == GLOBAL)
  {
    _J0_valid = true;
    _J0 = inertia;
  }
}

/// Sets the inertial pose for this rigid body
/**
 * Inertial pose should be defined relative to the rigid body pose
 */
void RIGIDBODY::set_inertial_pose(const POSE3& P)
{
  // verify that pose is set relative to body pose
  if (P.rpose != _F)
    throw std::runtime_error("RigidBody::set_inertial_pose() - inertial pose not defined relative to body pose");

  // update P to refer to _jF's pose
  POSE3 Q = P;
  Q.update_relative_pose(_jF->rpose);
  *_jF = Q;

  // update the mixed pose
  update_mixed_pose();

  // invalidate vectors using inertial frame
  _xdcom_valid = _xddcom_valid = _forcecom_valid = false;
}

/// Gets the current sum of forces on this body
const SFORCE& RIGIDBODY::sum_forces()
{
  const shared_ptr<const POSE3> GLOBAL;

  switch (_rftype)
  {
    case eGlobal:
      if (!_force0_valid)
        _force0 = POSE3::transform(GLOBAL, _forcem);
      _force0_valid = true;
      return _force0;

    case eLink:
      if (!_forcei_valid)
        _forcei = POSE3::transform(_F, _forcem);
      _forcei_valid = true;
      return _forcei;

    case eLinkCOM:
      if (!_forcecom_valid)
        _forcecom = POSE3::transform(_F2, _forcem);
      _forcecom_valid = true;
      return _forcecom;

    case eLinkInertia:
      return _forcem;

    case eJoint:
      if (!_forcej_valid)
        _forcej = POSE3::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _forcem);
      _forcej_valid = true;
      return _forcej;

    default:
      assert(false);
  }
}

/// Gets the current body velocity
const SVELOCITY& RIGIDBODY::get_velocity()
{
  const shared_ptr<const POSE3> GLOBAL;

  switch (_rftype)
  {
    case eLinkInertia:
      return _xdm;

    case eLink:
      if (!_xdi_valid)
        _xdi = POSE3::transform(const_pointer_cast<const POSE3>(_F), _xdcom);
      _xdi_valid = true;
      return _xdi;

    case eLinkCOM:
      if (!_xdcom_valid)
        _xdcom = POSE3::transform(const_pointer_cast<const POSE3>(_F2), _xdm);
      _xdcom_valid = true;
      return _xdcom;

    case eJoint:
      if (!_xdj_valid)
        _xdj = POSE3::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _xdm);
      _xdj_valid = true;
      return _xdj;

    case eGlobal:
      if (!_xd0_valid)
        _xd0 = POSE3::transform(GLOBAL, _xdm);
      _xd0_valid = true;
      return _xd0;

    default:
      assert(false);
  }
}

/// Gets the body inertia
const SPATIAL_RB_INERTIA& RIGIDBODY::get_inertia()
{
  const shared_ptr<const POSE3> GLOBAL;

  switch (_rftype)
  {
    case eGlobal:
      if (!_J0_valid)
        _J0 = POSE3::transform(GLOBAL, _Jm);
      _J0_valid = true;
      return _J0;

    case eLink:
      if (!_Ji_valid)
        _Ji = POSE3::transform(_F, _Jm);
      _Ji_valid = true;
      return _Ji;

    case eLinkCOM:
      if (!_Jcom_valid)
        _Jcom = POSE3::transform(_F2, _Jm);
      _Jcom_valid = true;
      return _Jcom;

    case eLinkInertia:
      return _Jm;

    case eJoint:
      if (!_Jj_valid)
        _Jj = POSE3::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _Jm);
      _Jj_valid = true;
      return _Jj;

    default:
      assert(false);
  }
}

/// Gets the current body acceleration
const SACCEL& RIGIDBODY::get_accel()
{
  const shared_ptr<const POSE3> GLOBAL;

  // do simplified case where body is disabled
  if (!is_enabled())
  {
    _xddi_valid = _xddcom_valid = _xddj_valid = _xdd0_valid = true;
  }

  switch (_rftype)
  {
    case eLinkInertia:
      return _xddm;

    case eLink:
      if (!_xddi_valid)
        _xddi = POSE3::transform(_F, _xddm, _xdm);
      _xddi_valid = true;
      return _xddi;

    case eLinkCOM:
      if (!_xddcom_valid)
        _xddcom = POSE3::transform(_F2, _xddm, _xdm);
      _xddcom_valid = true;
      return _xddcom;

    case eJoint:
      if (!_xddj_valid)
        _xddj = POSE3::transform((is_base()) ? _F : get_inner_joint_explicit()->get_pose(), _xddm, _xdm);
      _xddj_valid = true;
      return _xddj;

    case eGlobal:
      if (!_xdd0_valid)
        _xdd0 = POSE3::transform(GLOBAL, _xddm, _xdm);
      _xdd0_valid = true;
      return _xdd0;

    default:
      assert(false);
  }
}

/// Resets the force accumulators on this body
void RIGIDBODY::reset_accumulators()
{
  // clear forces
  _force0.set_zero();
  _forcei.set_zero();
  _forcem.set_zero();
  _forcej.set_zero();
  _forcecom.set_zero();

  // validate all forces
  _forcei_valid = true;
  _forcecom_valid = true;
  _forcej_valid = true;
  _force0_valid = true;
}

/// Computes the torques (w x Jw) that come from the Euler component of the Newton-Euler equations
SFORCE RIGIDBODY::calc_euler_torques()
{
  FILE_LOG(LOG_DYNAMICS) << "Calculating Euler torques for " << body_id << std::endl;
  VECTOR3 omega = get_velocity().get_angular();
  MATRIX3 J = get_inertia().J;
  ORIGIN3 w(omega);

  SFORCE f(omega.pose);
  f.set_zero();
  f.set_torque(VECTOR3(ORIGIN3::cross(w, J * w), omega.pose));
  return f;

/*
  VECTOR3 xd = get_velocity();
  return xd.cross(get_inertia() * xd);
*/
}

/// Updates the center-of-mass center / global aligned frame
/**
 * \note this function is called by set_pose() and 
 *       RCArticulatedBody::update_link_poses(.)
 */
void RIGIDBODY::update_mixed_pose()
{
  const shared_ptr<const POSE3> GLOBAL;

  // update the mixed pose
  _F2->set_identity();
  _F2->rpose = _jF;
  _F2->update_relative_pose(GLOBAL);
  _F2->q.set_identity();

  // invalidate pose vectors
  invalidate_pose_vectors();
}

/// Sets the current 3D pose for this rigid body
/**
 * Also updates the transforms for associated visualization and collision data.
 */
void RIGIDBODY::set_pose(const POSE3& p)
{
  // verify that the two poses are relative to the same pose
  if (p.rpose != _F->rpose)
    throw std::runtime_error("RigidBody::set_pose() - relative pose is not correct");

  // update the pose
  *_F = p;

  // update the mixed pose
  update_mixed_pose();
}

/// Invalidates pose quantities
void RIGIDBODY::invalidate_pose_vectors()
{
  // invalidate vector quantities
  _forcei_valid = _forcecom_valid = _force0_valid = false;
  _xdi_valid = _xdcom_valid = _xd0_valid = false;
  _xddi_valid = _xddcom_valid = _xdd0_valid = false;

  // invalidate inertias
  _Ji_valid = false;
  _Jj_valid = false;
  _J0_valid = false;
  _Jcom_valid = false;

  // invalidate every outer rigid body
  vector<shared_ptr<RIGIDBODY> > outer;
  BOOST_FOREACH(shared_ptr<JOINT> j, _outer_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
      outer.push_back(j->get_outboard_link());
  }
  vector<shared_ptr<RIGIDBODY> >::const_iterator end = std::unique(outer.begin(), outer.end());
  for (vector<shared_ptr<RIGIDBODY> >::const_iterator i = outer.begin(); i != end; i++)
    (*i)->invalidate_pose_vectors();
}

/// Gets the desired child link
shared_ptr<RIGIDBODY> RIGIDBODY::get_child_link(shared_ptr<JOINT> j) const
{
  assert(_outer_joints.find(j) != _outer_joints.end());
  return j->get_outboard_link();
}

/// Sets the force on the body (this function is private b/c I can't imagine where it should be called by the user)
void RIGIDBODY::set_force(const SFORCE& w)
{
  const shared_ptr<const POSE3> GLOBAL;

  // do not add forces to disabled bodies
  if (!_enabled)
    return;

  // update the force
  _forcem = POSE3::transform(_jF, w);

  // see whether we update a force
  if (w.pose == _F)
  {
    _forcei_valid = true;
    _forcei = w;

    // invalidate the remaining forces
    _forcej_valid = _forcecom_valid = _force0_valid = false;
  }
  else if (w.pose == _F2)
  {
    _forcecom_valid = true;
    _forcecom = w;

    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _force0_valid = false;
  }
  else if (!is_base() && w.pose == get_inner_joint_explicit()->get_pose())
  {
    _forcej_valid = true;
    _forcej = w;

    // invalidate the remaining forces
    _forcei_valid = _forcecom_valid = _force0_valid = false;
  }
  else if (w.pose == GLOBAL)
  {
    _force0_valid = true;
    _force0 = w;

    // invalidate the remaining forces
    _forcei_valid = _forcecom_valid = _forcej_valid = false;
  }
  else
    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _forcecom_valid = _force0_valid = false;
}

/// Adds a force to the body
void RIGIDBODY::add_force(const SFORCE& w)
{
  const shared_ptr<const POSE3> GLOBAL;

  // do not add forces to disabled bodies
  if (!_enabled)
    return;

  // update the force
  _forcem += POSE3::transform(_jF, w);

  // see whether we update a force
  if (w.pose == _F)
  {
    if (_forcei_valid)
      _forcei += w;

    // invalidate the remaining forces
    _forcej_valid = _forcecom_valid = _force0_valid = false;
  }
  else if (w.pose == _F2)
  {
    if (_forcecom_valid)
      _forcecom += w;

    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _force0_valid = false;
  }
  else if (!is_base() && w.pose == get_inner_joint_explicit()->get_pose())
  {
    if (_forcej_valid)
      _forcej += w;

    // invalidate the remaining forces
    _forcei_valid = _forcecom_valid = _force0_valid  = false;
  }
  else if (w.pose == GLOBAL)
  {
    if (_force0_valid)
      _force0 += w;

    // invalidate the remaining forces
    _forcei_valid = _forcecom_valid = _forcej_valid  = false;
  }
  else
    // invalidate the remaining forces
    _forcei_valid = _forcej_valid = _forcecom_valid = _force0_valid  = false;
}

/// Calculates the velocity of a point on this rigid body in the body frame
VECTOR3 RIGIDBODY::calc_point_vel(const VECTOR3& point) const
{
  // if the body is disabled, point velocity is zero
  if (!_enabled)
    return VECTOR3::zero(_F);

  // convert point to a vector in the body frame
  VECTOR3 r = POSE3::transform_point(_F, point);

  // get the velocity in the body frame
  SVELOCITY xd = POSE3::transform(const_pointer_cast<const POSE3>(_F), _xd0);

  // compute the point velocity - in the body frame
  VECTOR3 pv = xd.get_linear() + VECTOR3::cross(xd.get_angular(), r);
  pv.pose = _F;
  return pv;
}

/// Adds an inner joint for this link
/**
 * \param parent the outer link of the parent
 * \param j the joint connecting parent and this
 */
void RIGIDBODY::add_inner_joint(shared_ptr<JOINT> j)
{
  _inner_joints.insert(j);

  // update the spatial axes
  j->update_spatial_axes();

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(shared_ptr<ARTICULATED_BODY>(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> abody1 = j->get_articulated_body();
    shared_ptr<ARTICULATED_BODY> abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Adds an outer joint for this link
/**
 * \param j the joint connecting this and child
 * \note replaces the outer joint if it is already attached to this link
 */
void RIGIDBODY::add_outer_joint(shared_ptr<JOINT> j)
{
  // add the outer joint
  _outer_joints.insert(j);

  // update the spatial axes
  j->update_spatial_axes();

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(shared_ptr<ARTICULATED_BODY>(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> abody1 = j->get_articulated_body();
    shared_ptr<ARTICULATED_BODY> abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Determines whether the given link is a child link of this
bool RIGIDBODY::is_child_link(shared_ptr<const RIGIDBODY> query) const
{
  BOOST_FOREACH(shared_ptr<JOINT> j, _outer_joints)
    if (shared_ptr<RIGIDBODY>(j->get_outboard_link()) == query)
      return true;

  return false;
}

/// Determines whether the given link is a descendant of this
/**
 * \note returns <b>true</b> if query == this
 */
bool RIGIDBODY::is_descendant_link(shared_ptr<const RIGIDBODY> query) const
{
  queue<shared_ptr<const RIGIDBODY> > q;

  // check for query == this
  if (query == shared_from_this())
    return true;

  // add all children to the queue
  BOOST_FOREACH(shared_ptr<JOINT> j, _outer_joints)
    q.push(shared_ptr<const RIGIDBODY>(j->get_outboard_link()));

  // continue processing children until no more children are able to be processed
  while (!q.empty())
  {
    shared_ptr<const RIGIDBODY> link = q.front();
    q.pop();
    if (link == query)
      return true;
    BOOST_FOREACH(shared_ptr<JOINT> j, link->_outer_joints)
      q.push(shared_ptr<const RIGIDBODY>(j->get_outboard_link()));
  }

  return false;
}

/// Removes the specified outer joint from this link
void RIGIDBODY::remove_outer_joint(shared_ptr<JOINT> joint)
{
  _outer_joints.erase(joint);
}

/// Removes the specified outer joint from this link
/**
 * Returns true if the link was found.
 */
void RIGIDBODY::remove_inner_joint(shared_ptr<JOINT> joint)
{
  _inner_joints.erase(joint);
}

/// Applies a impulse to this link
/**
 * \param w the impulse as a force
 */
void RIGIDBODY::apply_impulse(const SMOMENTUM& w)
{
  const shared_ptr<const POSE3> GLOBAL;

  // if this is not an articulated body, just update linear and angular
  // momenta and velocites
  if (_abody.expired())
  {
    if (!_enabled)
      return;

    // get velocity update
    SMOMENTUM wx = POSE3::transform(get_computation_frame(), w);
    SVELOCITY dxd = get_inertia().inverse_mult(wx);

    // update linear and angular velocities
    _xdm += POSE3::transform(const_pointer_cast<const POSE3>(_jF), dxd);

    // see whether we can update any velocities
    if (dxd.pose == _F)
    {
      if (_xdi_valid)
        _xdi += dxd;
      _xdj_valid = _xdcom_valid = _xd0_valid = false;
    }
    else if (dxd.pose == _F2)
    {
      if (_xdcom_valid)
        _xdcom += dxd;
      _xdj_valid = _xdi_valid = _xd0_valid = false;
    }
    else if (!is_base() && dxd.pose == get_inner_joint_explicit()->get_pose())
    {
      if (_xdj_valid)
        _xdj += dxd;
      _xdcom_valid = _xdi_valid = _xd0_valid = false;
    }
    else if (dxd.pose == GLOBAL)
    {
      if (_xd0_valid)
        _xd0 += dxd;
      _xdcom_valid = _xdi_valid = _xdj_valid = false;
    }
    else
      _xdcom_valid = _xdi_valid = _xdj_valid = _xd0_valid = false;


    // reset the force and torque accumulators for this body
    reset_accumulators();
  }
  else
  {
    // get the articulated body
    shared_ptr<ARTICULATED_BODY> abody(_abody);

    // apply the impulse to the articulated body
    abody->apply_impulse(w, get_this());
  }
}

/// Gets the generalized inertia of this rigid body
unsigned RIGIDBODY::num_generalized_coordinates(GeneralizedCoordinateType gctype) const
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->num_generalized_coordinates(gctype);
  }
  else
    return num_generalized_coordinates_single(gctype);
}

/// Sets the generalized forces on the rigid body
void RIGIDBODY::set_generalized_forces(const Ravelin::SHAREDVECTORN& gf)
{
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->add_generalized_force(gf);
    return;
  }

  // if we're still here, this is only an individual body
  assert(gf.size() == num_generalized_coordinates(DYNAMIC_BODY::eSpatial));
  SFORCE w;

  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  // set the pose for w
  w.pose = _F2;

  // get the force and torque
  w.set_force(VECTOR3(gf[0], gf[1], gf[2]));
  w.set_torque(VECTOR3(gf[3], gf[4], gf[5]));

  // add the force to the sum of forces
  set_force(w);
}

/// Adds a generalized force to this rigid body
void RIGIDBODY::add_generalized_force(const SHAREDVECTORN& gf)
{
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->add_generalized_force(gf);
    return;
  }

  // if we're still here, this is only an individual body
  assert(gf.size() == num_generalized_coordinates(DYNAMIC_BODY::eSpatial));
  SFORCE w;

  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  // set the pose for w
  w.pose = _F2;

  // get the force and torque
  w.set_force(VECTOR3(gf[0], gf[1], gf[2]));
  w.set_torque(VECTOR3(gf[3], gf[4], gf[5]));

  // add the force to the sum of forces
  add_force(w);
}

/// Applies a generalized impulse to this rigid body
void RIGIDBODY::apply_generalized_impulse(const SHAREDVECTORN& gj)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->apply_generalized_impulse(gj);
    return;
  }
  else
    apply_generalized_impulse_single(gj);
}

/// Applies a generalized impulse to this rigid body
void RIGIDBODY::apply_generalized_impulse_single(const SHAREDVECTORN& gj)
{
  SMOMENTUM w;

  // don't do anything if this body is disabled
  if (!_enabled)
    return;

  // simple error check...
  assert(gj.size() == num_generalized_coordinates(DYNAMIC_BODY::eSpatial));

  // clear the force accumulators (and validate them all)
  reset_accumulators();

  // get the impulses
  w.set_linear(VECTOR3(gj[0], gj[1], gj[2]));
  w.set_angular(VECTOR3(gj[3], gj[4], gj[5]));
  w.pose = _F2;

  // get the impulse in the inertial frame
  SMOMENTUM wj = POSE3::transform(_jF, w);

  // get the current velocity in the inertial frame
  SVELOCITY v = _xdm;

  // update the velocity
  v += _Jm.inverse_mult(wj);

  set_velocity(v);
}

/// Solves using the generalized inertia matrix
SHAREDMATRIXN& RIGIDBODY::transpose_solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->transpose_solve_generalized_inertia(B, X);
  }
  else
  {
    // get proper generalized inertia matrix
    const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
    MATRIXN M;
    M.resize(NGC, NGC);
    SHAREDMATRIXN Mshared = M.block(0, NGC, 0, NGC);
    get_generalized_inertia_inverse(Mshared);
    M.mult_transpose(B, X);
    return X;
  }
}

/// Solves using the generalized inertia matrix
SHAREDMATRIXN& RIGIDBODY::solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->solve_generalized_inertia(B, X);
  }
  else
    return solve_generalized_inertia_single(B, X);
}

/// Solves using the generalized inertia matrix (does not call articulated body version)
SHAREDMATRIXN& RIGIDBODY::transpose_solve_generalized_inertia_single(const SHAREDMATRIXN& B, SHAREDMATRIXN& X)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  MATRIXN M;
  M.resize(NGC, NGC);
  SHAREDMATRIXN Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult_transpose(B, X);

  return X;
}

/// Solves using the generalized inertia matrix (does not call articulated body version)
SHAREDMATRIXN& RIGIDBODY::solve_generalized_inertia_single(const SHAREDMATRIXN& B, SHAREDMATRIXN& X)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  MATRIXN M;
  M.resize(NGC, NGC);
  SHAREDMATRIXN Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult(B, X);

  return X;
}

/// Solves using the generalized inertia matrix
SHAREDVECTORN& RIGIDBODY::solve_generalized_inertia(const SHAREDVECTORN& b, SHAREDVECTORN& x)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->solve_generalized_inertia(b, x);
  }
  else
    return solve_generalized_inertia_single(b, x);
}

/// Solves using the generalized inertia matrix
SHAREDVECTORN& RIGIDBODY::solve_generalized_inertia_single(const SHAREDVECTORN& b, SHAREDVECTORN& x)
{
  // get proper generalized inertia matrix
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  MATRIXN M;
  M.resize(NGC, NGC);
  SHAREDMATRIXN Mshared = M.block(0, NGC, 0, NGC);
  get_generalized_inertia_inverse(Mshared);
  M.mult(b, x);

  return x;
}

/// Gets the generalized position of this rigid body
SHAREDVECTORN& RIGIDBODY::get_generalized_coordinates(GeneralizedCoordinateType gctype, SHAREDVECTORN& gc)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->get_generalized_coordinates(gctype, gc);
  }
  else
    get_generalized_coordinates_generic(gctype, gc);

  return gc;
}

/// Sets the generalized coordinates of this rigid body
void RIGIDBODY::set_generalized_coordinates(GeneralizedCoordinateType gctype, const SHAREDVECTORN& gc)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->set_generalized_coordinates(gctype, gc);
  }
  else
    set_generalized_coordinates_generic(gctype, gc);
}

/// Sets the generalized velocity of this rigid body
void RIGIDBODY::set_generalized_velocity(GeneralizedCoordinateType gctype, const SHAREDVECTORN& gv)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->set_generalized_velocity(gctype, gv);
  }
  else
    set_generalized_velocity_generic(gctype, gv);
}

/// Sets the generalized acceleration of this rigid body
void RIGIDBODY::set_generalized_acceleration(const SHAREDVECTORN& ga)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->set_generalized_acceleration(ga);
  }
  else
    set_generalized_acceleration_generic(ga);
}

/// Gets the generalized velocity of this rigid body
SHAREDVECTORN& RIGIDBODY::get_generalized_velocity(GeneralizedCoordinateType gctype, SHAREDVECTORN& gv)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->get_generalized_velocity(gctype, gv);
  }
  else
    get_generalized_velocity_generic(gctype, gv);

  return gv;
}

/// Gets the generalized acceleration of this body
SHAREDVECTORN& RIGIDBODY::get_generalized_acceleration(SHAREDVECTORN& ga)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    ab->get_generalized_acceleration(ga);
  }
  else
    get_generalized_acceleration_generic(ga);

  return ga;
}

/// Gets the generalized inertia of this rigid body
SHAREDMATRIXN& RIGIDBODY::get_generalized_inertia(SHAREDMATRIXN& M)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->get_generalized_inertia(M);
  }
  else
    return get_generalized_inertia_single(M);
}

/// Gets the generalized inertia of this rigid body (does not call articulated body version)
SHAREDMATRIXN& RIGIDBODY::get_generalized_inertia_single(SHAREDMATRIXN& M)
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;

  // special case: disabled body
  if (!_enabled)
    return M.resize(0,0);

  // get the inertia
  SPATIAL_RB_INERTIA J = POSE3::transform(_F2, _Jm);

  // precompute some things
  MATRIX3 hxm = MATRIX3::skew_symmetric(J.h * J.m);
  MATRIX3 hxhxm = MATRIX3::skew_symmetric(J.h) * hxm;

  // arrange the matrix the way we want it: mass upper left, inertia lower right
  M.resize(SPATIAL_DIM, SPATIAL_DIM);
  M.set_sub_mat(0, 0, MATRIX3(J.m, 0, 0, 0, J.m, 0, 0, 0, J.m));
  M.set_sub_mat(3, 0, hxm);
  M.set_sub_mat(0, 3, hxm, Ravelin::eTranspose);
  M.set_sub_mat(3, 3, J.J - hxhxm);

  return M;
}

/// Gets the generalized inertia of this rigid body (does not call articulated body version)
SHAREDMATRIXN& RIGIDBODY::get_generalized_inertia_inverse(SHAREDMATRIXN& M) const
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 6;
  static LINALG _LA;

  // don't invert inertia for disabled bodies
  if (!_enabled)
  {
    M.resize(0,0);
    return M;
  }

  // get the inertia
  SPATIAL_RB_INERTIA J = POSE3::transform(_F2, _Jm);

  // setup the matrix
  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
  MATRIX3 hxm = MATRIX3::skew_symmetric(J.h*J.m);
  M.resize(6,6);
  M.set_sub_mat(0,3, hxm, eTranspose);
  M.set_sub_mat(3,3, J.J - hx*hxm);
  M.set_sub_mat(0,0, MATRIX3(J.m, 0, 0, 0, J.m, 0, 0, 0, J.m));
  M.set_sub_mat(3,0, hxm);

  // invert the matrix
  _LA.invert(M);

  return M;
}

/// Gets the generalized inertia of this rigid body
SHAREDVECTORN& RIGIDBODY::get_generalized_forces(SHAREDVECTORN& gf)
{
  // if this body part of an articulated body, call that function instead
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->get_generalized_forces(gf);
  }
  else
    return get_generalized_forces_single(gf);
}

/// Gets the generalized external forces (does not call articulated body version)
SHAREDVECTORN& RIGIDBODY::get_generalized_forces_single(SHAREDVECTORN& gf)
{
  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // resize the generalized forces vector
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  gf.resize(NGC);

  // get force and torque
  VECTOR3 f = _forcecom.get_force();
  VECTOR3 t = _forcecom.get_torque();

  // setup the linear components of f
  gf[0] = f[0];
  gf[1] = f[1];
  gf[2] = f[2];
  gf[3] = t[0];
  gf[4] = t[1];
  gf[5] = t[2];

  return gf;
}

/// Converts a force to a generalized force
SHAREDVECTORN& RIGIDBODY::convert_to_generalized_force(shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf)
{
  // if this belongs to an articulated body, call the articulated body method
  if (!_abody.expired())
  {
    shared_ptr<ARTICULATED_BODY> ab(_abody);
    return ab->convert_to_generalized_force(body, w, gf);
  }
  else
    return convert_to_generalized_force_single(body, w, gf);
}

/// Converts a force to a generalized force (does not call articulated body version)
SHAREDVECTORN& RIGIDBODY::convert_to_generalized_force_single(shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf)
{
  // verify that body == this
  assert(body.get() == this);

  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // transform w to computation frame
  SFORCE wt = POSE3::transform(_F2, w);

  // get linear and angular components of wt
  VECTOR3 f = wt.get_force();
  VECTOR3 t = wt.get_torque();

  // resize gf
  gf.resize(num_generalized_coordinates(DYNAMIC_BODY::eSpatial));

  // setup the linear components
  gf[0] = f[0];
  gf[1] = f[1];
  gf[2] = f[2];
  gf[3] = t[0];
  gf[4] = t[1];
  gf[5] = t[2];

  return gf;
}

/// Calculates the kinetic energy of the body with respect to the inertial frame
REAL RIGIDBODY::calc_kinetic_energy(shared_ptr<const POSE3> P)
{
  if (!_enabled)
    return (REAL) 0.0;

/*
  const SVELOCITY& xd = get_velocity();
  const SPATIAL_RB_INERTIA& J = get_inertia();
*/
  SVELOCITY xd = POSE3::transform(P, get_velocity());
  SPATIAL_RB_INERTIA J = POSE3::transform(P, get_inertia());

  VECTOR3 v = xd.get_linear();
  VECTOR3 w = xd.get_angular();
  VECTOR3 wx = VECTOR3(J.J*ORIGIN3(w), w.pose);
  return (v.norm_sq()*J.m + w.dot(wx))*(REAL) 0.5;

// return xd.dot(J * xd) * 0.5;
}

/// Gets the number of generalized coordinates
unsigned RIGIDBODY::num_generalized_coordinates_single(DYNAMIC_BODY::GeneralizedCoordinateType gctype) const
{
  const unsigned NGC_EULER = 7, NGC_SPATIAL = 6;

  // no generalized coordinates if this body is disabled
  if (!_enabled)
    return 0;

  // return the proper number of coordinates
  switch (gctype)
  {
    case DYNAMIC_BODY::eEuler:
      return NGC_EULER;

    case DYNAMIC_BODY::eSpatial:
      return NGC_SPATIAL;

    default:
      assert(false);
  }

  // make compiler happy
  assert(false);
  return 0;
}

/// Gets the first parent link of this link; returns NULL if there is no parent
shared_ptr<RIGIDBODY> RIGIDBODY::get_parent_link() const
{
  // special case (no parent!)
  if (_inner_joints.empty())
    return shared_ptr<RIGIDBODY>();

  // iterate through joints, looking for explicit constraint
  BOOST_FOREACH(shared_ptr<JOINT> j, _inner_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
      return shared_ptr<RIGIDBODY>(j->get_inboard_link());
  }

  // should not still be here
  assert(false);
  return shared_ptr<RIGIDBODY>();
}

/// Gets the explicit inner joint of this link; returns NULL if there is no explicit inner joint
/**
 * Throws an exception if this link has multiple explicit inner joints
 */
shared_ptr<JOINT> RIGIDBODY::get_inner_joint_explicit() const
{
  shared_ptr<JOINT> ij;
  BOOST_FOREACH(shared_ptr<JOINT> j, _inner_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
    {
      if (ij)
        throw std::runtime_error("Multiple explicit joints detected for a single link!");
      else
        ij = j;
    }
  }

  return ij;
}

/// Determines whether this link is a "ground" (fixed link)
bool RIGIDBODY::is_ground() const
{
  // clear easy cases
  if (!_enabled)
    return true;

  // can't be a ground if not disabled and not part of an articulated body
  if (_abody.expired())
    return false;

  // now, case will differ depending on what type of articulated body this is
  shared_ptr<ARTICULATED_BODY> ab(_abody);
  shared_ptr<RC_ARTICULATED_BODY> rcab = dynamic_pointer_cast<RC_ARTICULATED_BODY>(ab);
  if (rcab)
  {
    // check whether inner explicit joints are present (if none are present,
    // this is a base link)
    bool is_base = true;
    BOOST_FOREACH(shared_ptr<JOINT> j, _inner_joints)
    {
      if (j->get_constraint_type() == JOINT::eExplicit)
      {
        is_base = false;
        break;
      }
    }

    // if this link is a base and the base is fixed, it is a ground
    if (is_base && !rcab->is_floating_base())
      return true;
  }

  // still here? can't be a ground link
  return false;
}

/// Determines whether this link is the base
bool RIGIDBODY::is_base() const
{
  // clear easy cases
  if (_abody.expired())
    return true;

  // check whether no explicit joints are present
  BOOST_FOREACH(shared_ptr<JOINT> j, _inner_joints)
  {
    if (j->get_constraint_type() == JOINT::eExplicit)
      return false;
  }

  // no explicit joints... it's the base
  return true;
}


