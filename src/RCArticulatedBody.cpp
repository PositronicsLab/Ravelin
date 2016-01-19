/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::queue;
using std::list;
using std::map;
using std::string;

/// Default constructor
/**
 * Constructs a reduced-coordinate articulated body with no joints and no links.
 */
RC_ARTICULATED_BODY::RC_ARTICULATED_BODY()
{
  _floating_base = false;
  _n_joint_DOF_explicit = 0;

  // create the linear algebra object
  _LA = shared_ptr<LINALG>(new LINALG);
  _fsab._LA = _LA;
  _crb._LA = _LA;

  // set default algorithm to CRB and computation frame to link c.o.m. 
  algorithm_type = eCRB;
  set_computation_frame_type(eLinkCOM);

  // invalidate position quanitites
  _position_invalidated = true;
}

/// Validates position variables
void RC_ARTICULATED_BODY::validate_position_variables()
{
  _position_invalidated = false;
}

/// Gets the frame used for generalized coordinate calculations
shared_ptr<const POSE3> RC_ARTICULATED_BODY::get_gc_pose() const
{
  if (_links.empty())
    throw std::runtime_error("Cannot return gc pose when articulated body has no links");

  return _links.front()->get_gc_pose();
}

/// Sets the computation frame type
void RC_ARTICULATED_BODY::set_computation_frame_type(ReferenceFrameType rftype)
{
  // set the reference frame
  _rftype = rftype;

  // invalidate
  _position_invalidated = true;

  // set the reference frame type for all links
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->set_computation_frame_type(rftype);
}

/// Determines whether all of the children of a link have been processed
bool RC_ARTICULATED_BODY::all_children_processed(shared_ptr<RIGIDBODY> link) const
{
  const std::set<boost::shared_ptr<JOINT> >& joints = link->get_outer_joints();
  BOOST_FOREACH(boost::shared_ptr<JOINT> j, joints)
  {
    shared_ptr<RIGIDBODY> child = j->get_outboard_link();
    if (!_processed[child->get_index()])
      return false;
  }

  return true;
}

/// Gets the number of generalized coordinates for this body
unsigned RC_ARTICULATED_BODY::num_generalized_coordinates(DYNAMIC_BODY::GeneralizedCoordinateType gctype) const
{
  // look for trivial case
  if (_links.empty())
    return 0;

  if (!_floating_base)
    return _n_joint_DOF_explicit;
  else
    return _n_joint_DOF_explicit + _links.front()->num_generalized_coordinates_single(gctype);
}

/// Updates inverse generalized inertia matrix, as necessary
void RC_ARTICULATED_BODY::update_factorized_generalized_inertia()
{
  // see whether we need to update
  if (!_position_invalidated)
    return;

  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body = dynamic_pointer_cast<RC_ARTICULATED_BODY>(shared_from_this());

  // do precalculation on the body
  if (algorithm_type == eFeatherstone)
    _fsab.calc_spatial_inertias(body);
  else
    _crb.precalc(body);

  // indicate factorized inertia is valid
  validate_position_variables(); 
}

/// Solves using a generalized inertia matrix
SHAREDVECTORN& RC_ARTICULATED_BODY::solve_generalized_inertia(const SHAREDVECTORN& v, SHAREDVECTORN& result)
{
  // store the body's computation reference frame type
  ReferenceFrameType rftype = get_computation_frame_type();

  // set the reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(eLinkCOM);

  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // make x/b one vector
    result = v;

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // make x/b one vector
    result = v;

    // solve once
    _crb.M_solve_noprecalc(result);
  }

  // revert the link reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(rftype); 

  return result;
}

/// Solves the transpose using a generalized inertia matrix
SHAREDMATRIXN& RC_ARTICULATED_BODY::transpose_solve_generalized_inertia(const SHAREDMATRIXN& m, SHAREDMATRIXN& result)
{
  // store the body's computation reference frame type
  ReferenceFrameType rftype = get_computation_frame_type();

  // set the reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(eLinkCOM);

  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    MATRIXN::transpose(m, result);

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    MATRIXN::transpose(m, result);

    // solve
    _crb.M_solve_noprecalc(result);
  }

  // revert the link reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(rftype); 

  return result;
}

/// Solves using a generalized inertia matrix
SHAREDMATRIXN& RC_ARTICULATED_BODY::solve_generalized_inertia(const SHAREDMATRIXN& m, SHAREDMATRIXN& result)
{
  // store the body's computation reference frame type
  ReferenceFrameType rftype = get_computation_frame_type();

  // set the reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(eLinkCOM);

  if (algorithm_type == eFeatherstone)
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    result = m;

    // solve
    _fsab.solve_generalized_inertia_noprecalc(result);
  }
  else
  {
    // update the inverse / factorized inertia (if necessary)
    update_factorized_generalized_inertia();

    // setup the result
    result = m;

    // solve
    _crb.M_solve_noprecalc(result);
  }

  // revert the link reference frame type
  if (rftype != eLinkCOM && is_floating_base())
    set_computation_frame_type(rftype); 

  return result;
}

/// Applies a generalized impulse to the articulated body
void RC_ARTICULATED_BODY::apply_generalized_impulse(const SHAREDVECTORN& gj)
{
  if (algorithm_type == eFeatherstone)
    _fsab.apply_generalized_impulse(gj);
  else
  {
    assert(algorithm_type == eCRB);

    // setup work variables
    static VECTORN gv, gv_delta;

    // get the current generalized velocity
    get_generalized_velocity(DYNAMIC_BODY::eSpatial, gv);

    // we'll solve for the change in generalized velocity
    DYNAMIC_BODY::solve_generalized_inertia(gj, gv_delta);

    // apply the change in generalized velocity
    gv += gv_delta;
    set_generalized_velocity(DYNAMIC_BODY::eSpatial, gv);
  }

  // reset the force and torque accumulators
  reset_accumulators();
}

/// Sets the generalized forces for the articulated body
void RC_ARTICULATED_BODY::set_generalized_forces(const SHAREDVECTORN& gf)
{
  unsigned index = 0;
  SFORCE f0;

  if (_floating_base)
  {
    // get the base
    shared_ptr<RIGIDBODY> base = _links.front();

    // first, get the force on the base link
    gf.get_sub_vec(num_joint_dof_explicit(), gf.size(), f0);

    // add the force to the base
    SFORCE fx = POSE3::transform(base->get_gc_pose(), f0);
    base->set_force(fx);
  }

  // add to joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    CONST_SHAREDVECTORN f = gf.segment(idx, idx+_ejoints[i]->num_dof());
    _ejoints[i]->force = f;
  }
}

/// Adds a generalized force to the articulated body
void RC_ARTICULATED_BODY::add_generalized_force(const SHAREDVECTORN& gf)
{
  unsigned index = 0;

  if (_floating_base)
  {
    // get the base
    shared_ptr<RIGIDBODY> base = _links.front();

    // first, get the force on the base link
    SFORCE f0;
    gf.get_sub_vec(num_joint_dof_explicit(), gf.size(), f0);

    // add the force to the base
    f0.pose = base->get_gc_pose();
    base->add_force(f0);
  }

  // add to joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    CONST_SHAREDVECTORN f = gf.segment(idx, idx+_ejoints[i]->num_dof());
    _ejoints[i]->force += f;
  }
}

/// Determines whether the link is effectively a leaf link
bool RC_ARTICULATED_BODY::treat_link_as_leaf(shared_ptr<RIGIDBODY> link) const
{
  // if all children have lower link indices, link treatable as end-effector
  const std::set<boost::shared_ptr<JOINT> >& joints = link->get_outer_joints();
  BOOST_FOREACH(boost::shared_ptr<JOINT> j, joints)
  {
    shared_ptr<RIGIDBODY> child = j->get_outboard_link();
    if (child->get_index() > link->get_index())
      return false;
  }

  return true;
}

/// Sets whether the base of this body is "floating" (or fixed)
void RC_ARTICULATED_BODY::set_floating_base(bool flag)
{
  _floating_base = flag;
  compile();
}

/// Compiles this body (updates the link transforms and velocities)
void RC_ARTICULATED_BODY::compile()
{
  // call parent method first
  ARTICULATED_BODY::compile();

  // verify all links are enabled
  if (!is_floating_base())
  {
    for (unsigned i=1; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("Only first link can be disabled in a reduced coordinate body with fixed-base!");
  }
  else
  {
    for (unsigned i=0; i< _links.size(); i++)
      if (!_links[i]->is_enabled())
        throw std::runtime_error("No links can be disabled in a reduced coordinate body with floating-base!");
  }

  // update processed vector size
  _processed.resize(_links.size());

  // setup explicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx = 0; i< _ejoints.size(); i++)
  {
    _ejoints[i]->set_coord_index(cidx);
    _ejoints[i]->set_constraint_index(ridx);
    cidx += _ejoints[i]->num_dof();
    ridx += _ejoints[i]->num_constraint_eqns();
  }

  // setup implicit joint generalized coordinate and constraint indices
  for (unsigned i=0, cidx = 0, ridx=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->set_coord_index(cidx);
    _ijoints[i]->set_constraint_index(ridx);
    cidx += _ijoints[i]->num_dof();
    ridx += _ijoints[i]->num_constraint_eqns();
  }

  // store all joint values and reset to zero; this is necessary because the
  // links are expected to be initialized at the joint zero positions
  vector<VECTORN> q_save(_ejoints.size());
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    q_save[i] = _ejoints[i]->q;
    _ejoints[i]->q.set_zero();
  }

  // update relative poses for explicit joints only
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    shared_ptr<RIGIDBODY> outboard = _ejoints[i]->get_outboard_link();
    shared_ptr<POSE3> pose = boost::const_pointer_cast<POSE3>(outboard->get_pose());
    pose->update_relative_pose(_ejoints[i]->get_induced_pose());
  }

  // restore all joint values
  for (unsigned i=0; i< _ejoints.size(); i++)
    _ejoints[i]->q = q_save[i];

  // point both algorithms to this body
  _crb.set_body(get_this());
  _fsab.set_body(get_this());

  // update link transforms and velocities
  update_link_poses();
  update_link_velocities();
}

/// Sets the vector of links and joints
void RC_ARTICULATED_BODY::set_links_and_joints(const vector<shared_ptr<RIGIDBODY> >& links, const vector<boost::shared_ptr<JOINT> >& joints)
{
  // setup the processed vector
  _processed.resize(links.size());

  // setup pointers
  for (unsigned i=0; i< joints.size(); i++)
  {
    // set pointers for inner and outer joints
    shared_ptr<RIGIDBODY> inboard = joints[i]->get_inboard_link();
    shared_ptr<RIGIDBODY> outboard = joints[i]->get_outboard_link();
    outboard->_inner_joints.insert(joints[i]);
    inboard->_outer_joints.insert(joints[i]);

    // set pointer for articulated body
    joints[i]->set_articulated_body(dynamic_pointer_cast<RC_ARTICULATED_BODY>(shared_from_this()));

    // set frames for outboard
    outboard->_xdj.pose = joints[i]->get_pose();
    outboard->_xddj.pose = joints[i]->get_pose();
    outboard->_Jj.pose = joints[i]->get_pose();
    outboard->_forcej.pose = joints[i]->get_pose();
  }

  // clear the vectors of joints
  _ejoints.clear();
  _ijoints.clear();
  
  // start processed at the base link
  map<shared_ptr<RIGIDBODY>, bool> processed;
  BOOST_FOREACH(boost::shared_ptr<RIGIDBODY> link, links)
  {
    // if the link has already been processed, no need to process it again
    if (processed[link])
      continue;

    // get all outer joints for this link
    BOOST_FOREACH(boost::shared_ptr<JOINT> joint, link->get_outer_joints())
    {
      // see whether the child has already been processed
      shared_ptr<RIGIDBODY> child(joint->get_outboard_link());
      if (processed[child])
        _ijoints.push_back(joint);
      else
      {
        _ejoints.push_back(joint);
      }
    }

    // indicate that the link has been processed
    processed[link] = true;
  }
  
  // recalculate the explicit joint degrees-of-freedom of this body
  _n_joint_DOF_explicit = 0;
  for (unsigned i=0; i< _ejoints.size(); i++)
    _n_joint_DOF_explicit += _ejoints[i]->num_dof();

  // mark joints as the correct type
  for (unsigned i=0; i< _ejoints.size(); i++)
    _ejoints[i]->set_constraint_type(JOINT::eExplicit);
  for (unsigned i=0; i< _ijoints.size(); i++)
    _ijoints[i]->set_constraint_type(JOINT::eImplicit);

  // set the computation frame type on all links
  for (unsigned i=0; i< links.size(); i++)
    links[i]->set_computation_frame_type(_rftype);

  // find the base link and setup a processed map
  shared_ptr<RIGIDBODY> base;
  for (unsigned i=0; i< joints.size(); i++)
  {
    shared_ptr<RIGIDBODY> inboard = joints[i]->get_inboard_link();
    shared_ptr<RIGIDBODY> outboard = joints[i]->get_outboard_link();
    shared_ptr<RIGIDBODY> parent = inboard->get_parent_link();
    if (!parent)
    {
      if (base && inboard != base)
        throw std::runtime_error("Multiple base links detected!");
      base = inboard;
    }
  }

  // if there is no clearly defined base link, there are no links *and* there
  // are joints, we can't do anything
  if (!base)
  {
    if (joints.empty() && links.size() == 1)
      base = links.front();
    else
    {
      // look for a disabled link
      for (unsigned i=0; i< links.size(); i++)
      {
        if (!links[i]->is_enabled())
        {
          if (base)
            throw std::runtime_error("Could not find unique base link!");
          else
            base = links[i];
        }
      }
    }
    if (!base)
      throw std::runtime_error("Could not find base link!");
  }

  // check to see whether user's numbering scheme is acceptable
  for (unsigned i=1; i< links.size(); i++)
  {
    // look for an unknown constraint
    BOOST_FOREACH(boost::shared_ptr<JOINT> joint, links[i]->get_inner_joints())
      if (joint->get_constraint_type() == JOINT::eUnknown)
        throw std::runtime_error("Unknown constraint type found!");

    // no unknown constraint; look for an explicit constraint
    if (!links[i]->get_inner_joint_explicit())
      throw std::runtime_error("Nonzero link does not have an inner explicit joint!");
  }

  // look whether it's a floating base
  _floating_base = base->is_enabled();

  // copy the vector
  _links = links;

  // setup the link in the map 
  for (unsigned i=0; i< _links.size(); i++)
  {
    _links[i]->set_index(i);
    _links[i]->set_articulated_body(get_this());
  }

  // set vector of joints
  _joints = joints;

  // iterate over each joint
  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->set_index(i);

  // compile the body
  compile();
}

/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as
 * well as the base momentum; therefore, the derivative of the state vector is
 * composed of the joint-space accelerations and base forces (and torques).
 */
/// Gets the derivative of the velocity state vector for this articulated body
/**
 * The state vector consists of the joint-space velocities of the robot as
 * well as the base momentum; therefore, the derivative of the state vector is
 * composed of the joint-space accelerations and base forces (and torques).
 */
SHAREDVECTORN& RC_ARTICULATED_BODY::get_generalized_acceleration(SHAREDVECTORN& ga)
{
  get_generalized_acceleration_generic(ga);
  return ga;
}

/// Updates the transforms of the links based on the current joint positions
/**
 * \note this doesn't actually calculate other than the joint positions; all
 *       links are defined with respect to the joints, which are defined
 *       with respect to their inner link
 */
void RC_ARTICULATED_BODY::update_link_poses()
{
  // indicate factorized inertia matrix is no longer valid
  _position_invalidated = true;

  // update all joint poses
  for (unsigned i=0; i< _joints.size(); i++){
    _joints[i]->get_induced_pose();
  }

  // update the center-of-mass centered / global aligned frame in all links
  for (unsigned i=0; i< _links.size(); i++){
    _links[i]->update_mixed_pose();
  }

  // print all link poses and joint poses
  if (LOGGING(LOG_DYNAMICS))
  {
    const shared_ptr<const POSE3> GLOBAL;
    for (unsigned i=0; i< _links.size(); i++)
    {
      TRANSFORM3 Tx = POSE3::calc_relative_pose(_links[i]->get_pose(), GLOBAL);
      FILE_LOG(LOG_DYNAMICS) << "  link " << _links[i]->body_id << " pose (relative to global frame): " << Tx.x << " " << AANGLE(Tx.q) << std::endl;
    }
    for (unsigned i=0; i< _joints.size(); i++)
    {
      TRANSFORM3 Tx = POSE3::calc_relative_pose(_joints[i]->get_pose(), GLOBAL);
      FILE_LOG(LOG_DYNAMICS) << "  joint " << _joints[i]->joint_id << " pose (relative to global frame): " << Tx.x << " " << AANGLE(Tx.q) << std::endl;
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "RC_ARTICULATED_BODY::update_link_poses() exited" << std::endl;
}

/// Updates the link velocities
void RC_ARTICULATED_BODY::update_link_velocities()
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  vector<SVELOCITY> sprime;

  // look for easy exit
  if (_links.empty() || _joints.empty())
    return;

  // get the base link
  shared_ptr<RIGIDBODY> base = _links.front();

  FILE_LOG(LOG_DYNAMICS) << "RC_ARTICULATED_BODY::update_link_velocities() entered" << std::endl;
  if (LOGGING(LOG_DYNAMICS))
  {
    for (unsigned i=0; i< _links.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << " link " << i << " " << _links[i]->body_id << std::endl;
    for (unsigned i=0; i< _joints.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << " joint " << i << " " <<  _joints[i]->joint_id << std::endl;
  }

  // add all children of the base to the link queue
  list<shared_ptr<RIGIDBODY> > child_links;
  base->get_child_links(std::back_inserter(child_links));
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
    link_queue.push(rb);

  // reset processed vector
  for (unsigned i=0; i< _links.size(); i++)
    _processed[i] = false;
  _processed[base->get_index()] = true;

  // propagate link velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> outboard = link_queue.front();
    shared_ptr<const POSE3> opose = outboard->get_computation_frame();
    link_queue.pop();
    unsigned i = outboard->get_index();

    // and push all children of the link onto the queue
    list<shared_ptr<RIGIDBODY> > child_links;
    outboard->get_child_links(std::back_inserter(child_links));
    BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
      if (!_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the inboard link
    boost::shared_ptr<JOINT> joint(outboard->get_inner_joint_explicit());
    shared_ptr<RIGIDBODY> inboard(joint->get_inboard_link());
    shared_ptr<const POSE3> ipose = inboard->get_computation_frame();
    unsigned h = inboard->get_index();

    // set this link's velocity to the parent's link velocity
    outboard->set_velocity(POSE3::transform(opose, inboard->get_velocity()));

    // get the (transformed) link spatial axis
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    POSE3::transform(opose, s, sprime);

    // determine the link velocity due to the parent velocity + joint velocity
    if (sprime.empty())
      outboard->set_velocity(outboard->get_velocity());
    else
      outboard->set_velocity(outboard->get_velocity() + SPARITH::mult(sprime, joint->qd));

    // indicate that the link has been processed
    _processed[i] = true;

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << outboard->body_id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent velocity: " << inboard->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- link velocity : " << outboard->get_velocity()  << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RC_ARTICULATED_BODY::update_link_velocities() exited" << std::endl;
}

/// Calculates the column of a Jacobian matrix for a floating base with respect to a given point
/**
 * \param point the point in 3D space which the Jacobian is calculated against
 * \return a pointer to a 6x6 matrix; top three dimensions will be linear
 *         velocity and bottom three dimensions will be angular velocity
 */
/*
MATRIXN& RC_ARTICULATED_BODY::calc_jacobian_floating_base(const VECTOR3& point, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;
  static vector<SVELOCITY> sbase, sbase_prime;
  static shared_ptr<POSE3> P(new POSE3);

  // get the base link and the base pose
  shared_ptr<RIGIDBODY> base = get_base_link();
  shared_ptr<const POSE3> baseP = base->get_gc_pose();

  // setup the twists at the base link - first three vectors are linear motion
  // and second three are angular
  sbase.clear();
  sbase.resize(SPATIAL_DIM, SVELOCITY::zero(baseP));
  sbase[3][0] = 1.0;  sbase[4][1] = 1.0;  sbase[5][2] = 1.0;
  sbase[0][3] = 1.0;  sbase[1][4] = 1.0;  sbase[2][5] = 1.0;

  // convert the poses to the point frame
  P->x = ORIGIN3(point);
  POSE3::transform(P, sbase, sbase_prime);

  // init the base Jacobian
  J.resize(SPATIAL_DIM, SPATIAL_DIM);
  for (unsigned i=0; i< SPATIAL_DIM; i++)
  {
    SHAREDVECTORN Ji = J.column(i);
    sbase_prime[i].transpose_to_vector(Ji);
  }

  return J;
}

/// Calculates the Jacobian for the current robot configuration at a given point and with respect to a given link
MATRIXN& RC_ARTICULATED_BODY::calc_jacobian(const VECTOR3& p, shared_ptr<RIGIDBODY> link, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;
  static MATRIXN Jsub;

  // resize the Jacobian
  J.set_zero(SPATIAL_DIM, num_generalized_coordinates(DYNAMIC_BODY::eSpatial));

  // get the base link
  shared_ptr<RIGIDBODY> base = get_base_link();

  if (is_floating_base())
  {
shared_ptr<POSE3> frame(new POSE3);
frame->rpose = p.pose;
frame->q.set_identity();
frame->x = ORIGIN3(p);

    // construct the spatial transform
    POSE3::spatial_transform_to_matrix2(base->get_gc_pose(), frame, Jsub);

    // setup the floating base
    J.set_sub_mat(0, num_joint_dof_explicit(),Jsub);
  }

  // calculate all relevant columns
  while (link != base)
  {
    boost::shared_ptr<JOINT> joint = link->get_inner_joint_explicit();
    calc_jacobian_column(joint, p, Jsub);
    J.set_sub_mat(0, joint->get_coord_index(), Jsub);
    link = link->get_parent_link();
  }

  return J;
}
*/

/// Calculates the column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \param base_pose the transform to use for the base
 * \param q a mapping from joints to joint positions to use; joints without
 *        positions specified will be set to their current values
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 * \note this method works by changing the joint values, recomputing all link
 *       transforms, computing the Jacobian, restoring the old joint values and
 *       link transforms; thus, it will be slower than a special purpose method
 */
/*
MATRIXN& RC_ARTICULATED_BODY::calc_jacobian(const VECTOR3& point, const POSE3& base_pose, const map<boost::shared_ptr<JOINT>, VECTORN>& q, shared_ptr<RIGIDBODY> link, MATRIXN& Jc)
{
  // store current joint values
  map<boost::shared_ptr<JOINT>, VECTORN> currentQ;
  for (unsigned i=0; i< _ejoints.size(); i++)
    currentQ[_ejoints[i]] = _ejoints[i]->q;

  // overwrite current joint values
  for (map<boost::shared_ptr<JOINT>, VECTORN>::const_iterator i = q.begin(); i != q.end(); i++)
    i->first->q = i->second;

  // save the current base pose and set it to that desired
  shared_ptr<RIGIDBODY> base = get_base_link();
  POSE3 saved_base_pose = *base->get_pose();
  base->set_pose(base_pose);

  // update link transforms
  update_link_poses();

  // compute and store the Jacobian
  calc_jacobian(point, link, Jc);

  // restore joint values
  for (map<boost::shared_ptr<JOINT>, VECTORN>::const_iterator i = currentQ.begin(); i != currentQ.end(); i++)
    i->first->q = i->second;

  // restore the base pose and transforms
  base->set_pose(saved_base_pose);
  update_link_poses();

  return Jc;
}
*/

/// Calculates column(s) of a Jacobian matrix
/*
 * \param joint the joint with which the Jacobian will be calculated
 * \param point the reference point in 3D space used to calculate the Jacobian
 * \return a 6xN matrix, where N is the number of DOF of the joint; the top
 *         three dimensions will be the contribution to linear velocity, and
 *         the bottom three dimensions will be the contribution to angular
 *         velocity
 */
MATRIXN& RC_ARTICULATED_BODY::calc_jacobian_column(boost::shared_ptr<JOINT> joint, const VECTOR3& point, MATRIXN& Jc)
{
  const unsigned SPATIAL_DIM = 6;
  static shared_ptr<POSE3> target(new POSE3);
  static vector<SVELOCITY> sprime;

  // NOTE: spatial algebra provides us with a simple means to compute the
  // Jacobian of a joint with respect to a point.  The spatial axis of the
  // joint transforms joint velocity to spatial velocity. The spatial
  // transform is then used to transform the spatial velocity into a
  // convenient frame: the orientation of that frame is aligned with the
  // global frame, and the origin of the frame is aligned with the desired
  // point.  Applying the spatial transform from the spatial axis (either in
  // link or global frame) to the new frame will give us the desired vector.

  // get the spatial axis of the joint
  const vector<SVELOCITY>& s = joint->get_spatial_axes();
  Jc.resize(SPATIAL_DIM, s.size());

  // compute a spatial transformation using the point as the target frame
  target->rpose = point.pose;
  target->x = ORIGIN3(point);
  target->q.set_identity();
  POSE3::transform(target, s, sprime);

  // calculate the Jacobian column
  for (unsigned i=0; i< sprime.size(); i++)
  {
    VECTOR3 top = sprime[i].get_upper();
    VECTOR3 bot = sprime[i].get_lower();
    Jc(0,i) = bot[0];
    Jc(1,i) = bot[1];
    Jc(2,i) = bot[2];
    Jc(3,i) = top[0];
    Jc(4,i) = top[1];
    Jc(5,i) = top[2];
  }

  return Jc;
}

/// Resets the force and torque accumulators for all links and joints in the rigid body
void RC_ARTICULATED_BODY::reset_accumulators()
{
  for (unsigned i=0; i< _links.size(); i++)
    _links[i]->reset_accumulators();

  for (unsigned i=0; i< _joints.size(); i++)
    _joints[i]->reset_force();
}

/// The signum function
REAL RC_ARTICULATED_BODY::sgn(REAL x)
{
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<REAL>::epsilon());
  if (x < -NEAR_ZERO)
    return (REAL) -1.0;
  else if (x > NEAR_ZERO)
    return (REAL) 1.0;
  else
    return (REAL) 0.0;
}

/// Computes the forward dynamics
/**
 * Given the joint positions and velocities, joint forces, and external
 * forces on the links, determines the joint and link accelerations as well as
 * floating base accelerations (if applicable).  The joint
 * accelerations are stored in the individual joints and the link accelerations
 * are stored in the individual links.
 * \note only computes the forward dynamics if the state-derivative is no longer valid
 */
void RC_ARTICULATED_BODY::calc_fwd_dyn()
{
  static VECTORN ff;

  FILE_LOG(LOG_DYNAMICS) << "RC_ARTICULATED_BODY::calc_fwd_dyn() entered" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "  computing forward dynamics in ";
  if (get_computation_frame_type() == eGlobal)
    FILE_LOG(LOG_DYNAMICS) << "global ";
  else if (get_computation_frame_type() == eLink)
    FILE_LOG(LOG_DYNAMICS) << "link ";
  else if (get_computation_frame_type() == eLinkCOM)
    FILE_LOG(LOG_DYNAMICS) << "link c.o.m. ";
  else
    FILE_LOG(LOG_DYNAMICS) << "joint ";
  FILE_LOG(LOG_DYNAMICS) << "coordinate system" << std::endl;

  // use the proper dynamics algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      if (!_position_invalidated)
        _fsab.calc_fwd_dyn_special();
      else
        _fsab.calc_fwd_dyn();
      break;

    case eCRB:
      if (!_position_invalidated)
        _crb.calc_fwd_dyn_special();
      else
        _crb.calc_fwd_dyn();
      break;

    default:
      assert(false);
  }

  FILE_LOG(LOG_DYNAMICS) << "RC_ARTICULATED_BODY::calc_fwd_dyn() exited" << std::endl;
}

/// Computes the forward dynamics with loops
/*
void RC_ARTICULATED_BODY::calc_fwd_dyn_loops()
{
  REAL Cx[6];
  static MATRIXN Jx_iM_JxT, U, V;
  static MATRIXN Jx_dot, iM_JxT;
  static VECTORN v, fext, C, alpha_x, beta_x, Dx_v, Jx_v, Jx_dot_v, workv;
  static VECTORN iM_fext, a, S;

  // get the generalized velocity, generalized forces, and inverse generalized
  // inertia matrix
  get_generalized_velocity(eSpatial, v);
  get_generalized_forces(fext);

  // determine how many implicit constraint equations
  unsigned N_EXPLICIT_CONSTRAINT_EQNS = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    N_EXPLICIT_CONSTRAINT_EQNS = _ijoints[i]->num_constraint_eqns();

  // evaluate implicit constraints
  C.resize(N_EXPLICIT_CONSTRAINT_EQNS);
  for (unsigned i=0, r=0, s=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->evaluate_constraints(Cx);
    for (unsigned j=0; j< _ijoints[i]->num_constraint_eqns(); j++)
      C[r++] = Cx[j];
  }

  // get the implicit constraint Jacobian and its time derivative
  determine_implicit_constraint_jacobian(_Jx);
  determine_implicit_constraint_jacobian_dot(Jx_dot);

  // compute Jx * v and Jx_dot * v
  _Jx.mult(v, Jx_v);
  Jx_dot.mult(v, Jx_dot_v);

  // get movement Jacobian for implicit constraints and compute velocities
  determine_implicit_constraint_movement_jacobian(_Dx);
  _Dx.mult(v, Dx_v);
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    Dx_v.get_sub_vec(k, k+_ijoints[i]->num_dof(), _ijoints[i]->qd);
    k += _ijoints[i]->num_dof();
  }

  // add in implicit actuator forces
  beta_x.resize(_Dx.rows());
  for (unsigned i=0, k=0; i< _ijoints.size(); i++)
  {
    _ijoints[i]->get_scaled_force(workv);
    for (unsigned j=0; j< _ijoints[i]->num_dof(); j++)
      beta_x[k++] = workv[j];
  }
  _Dx.transpose_mult(beta_x, workv);
  fext += workv;

  // compute the constraint forces
  DYNAMIC_BODY::solve_generalized_inertia(fext, iM_fext);
  _Jx.mult(iM_fext, alpha_x) += Jx_dot_v;
  _Jx.mult(v, workv) *= ((REAL) 2.0 * b_alpha);
  alpha_x += workv;
  C *= (b_beta*b_beta);
  alpha_x += C;
  DYNAMIC_BODY::transpose_solve_generalized_inertia(_Jx, iM_JxT);
  _Jx.mult(iM_JxT, Jx_iM_JxT);
  _LA->svd(Jx_iM_JxT, U, S, V);
  _LA->solve_LS_fast(U, S, V, alpha_x);

  // compute generalized acceleration
  fext -= _Jx.transpose_mult(alpha_x, workv);
  DYNAMIC_BODY::solve_generalized_inertia(fext, a);
  set_generalized_acceleration(a);
}
*/

/// Determines the constraint Jacobian for implicit constraints
void RC_ARTICULATED_BODY::determine_implicit_constraint_jacobian(MATRIXN& J)
{
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  const unsigned NSPATIAL = 6;
  MATRIXN tmp2;

  // determine the number of implicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NEQ += _ijoints[i]->num_constraint_eqns();

  // resize J
  J.set_zero(NEQ, NGC);

  // get the base link
  shared_ptr<RIGIDBODY> base = get_base_link();

  // setup a temporary matrix
  MATRIXN tmp(NSPATIAL, NGC);

  // compute the constraint Jacobian for all implicit constraints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the constraint index for this joint
    unsigned j = _ijoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    shared_ptr<RIGIDBODY> rbi = _ijoints[i]->get_inboard_link();
    shared_ptr<RIGIDBODY> rbo = _ijoints[i]->get_outboard_link();

    // get the number of constraint equations for this joint 
    const unsigned THIS_EQ = _joints[i]->num_constraint_eqns();

    // resize the temporary matrix
    tmp.resize(THIS_EQ, NGC);

    // setup the appropriate block of the matrices
    SHAREDMATRIXN Cqi = J.block(j, j+THIS_EQ, 0, NGC);
    SHAREDMATRIXN Cqo = tmp.block(0, tmp.rows(), 0, NGC);

    // get constraint equations for inner and outer links
    _ijoints[i]->calc_constraint_jacobian(true, tmp2); Cqi = tmp2;
    _ijoints[i]->calc_constraint_jacobian(false, tmp2); Cqo = tmp2;

    // add Cqo to Cqi
    Cqi += Cqo; 
  }
}

/// Determines the time derivative of the constraint Jacobian for implicit constraints
/**
 * Because this matrix is the product of two matrices, each of which is a
 * function of time, we have to add the results together.
 */
/// Determines the constraint Jacobian for implicit constraints
void RC_ARTICULATED_BODY::determine_implicit_constraint_jacobian_dot(MATRIXN& J)
{
  const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  const unsigned NSPATIAL = 6;
  MATRIXN tmp2;

  // determine the number of implicit constraint equations
  unsigned NEQ = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    NEQ += _ijoints[i]->num_constraint_eqns();

  // resize J
  J.set_zero(NEQ, NGC);

  // get the base link
  shared_ptr<RIGIDBODY> base = get_base_link();

  // setup a temporary matrix
  MATRIXN tmp(NSPATIAL, NGC);

  // compute the constraint Jacobian for all implicit constraints
  for (unsigned i=0; i< _ijoints.size(); i++)
  {
    // get the constraint index for this joint
    unsigned j = _ijoints[i]->get_constraint_index();

    // get the rigid bodies of this joint
    shared_ptr<RIGIDBODY> rbi = _ijoints[i]->get_inboard_link();
    shared_ptr<RIGIDBODY> rbo = _ijoints[i]->get_outboard_link();

    // get the number of constraint equations for this joint 
    const unsigned THIS_EQ = _joints[i]->num_constraint_eqns();

    // resize the temporary matrix
    tmp.resize(THIS_EQ, NGC);

    // setup the appropriate block of the matrices
    SHAREDMATRIXN Cqi = J.block(j, j+THIS_EQ, 0, NGC);
    SHAREDMATRIXN Cqo = tmp.block(0, tmp.rows(), 0, NGC);

    // get constraint equations for inner and outer links
    _ijoints[i]->calc_constraint_jacobian_dot(true, tmp2); Cqi = tmp2;
    _ijoints[i]->calc_constraint_jacobian_dot(false, tmp2); Cqo = tmp2;

    // add Cqo to Cqi
    Cqi += Cqo; 
  }
}

/// Sets the generalized acceleration for this body
void RC_ARTICULATED_BODY::set_generalized_acceleration(const SHAREDVECTORN& a)
{
  static VECTORN base_a;

  if (_floating_base)
  {
    a.get_sub_vec(num_joint_dof_explicit(), a.size(), base_a);
    shared_ptr<RIGIDBODY> base = _links.front();
    SACCEL xdd;
    xdd.set_linear(VECTOR3(base_a[0], base_a[1], base_a[2]));
    xdd.set_angular(VECTOR3(base_a[3], base_a[4], base_a[5]));
    base->set_accel(xdd);
  }

  // set joint accelerations
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    a.get_sub_vec(idx, idx+_ejoints[i]->num_dof(), _ejoints[i]->qdd);
  }
}

/// Applies an impulsive force at the specified point and propagates it up the articulated body chain
/**
 * \param w the impulse
 * \param link the link that the impulse is applied to
 * \param fsab_algo if non-null, already computed quanities from dynamics
 *         algorithm will be used; otherwise, forward dynamics algorithm will
 *         be called using FSAB algorithm
 */
void RC_ARTICULATED_BODY::apply_impulse(const SMOMENTUM& w, shared_ptr<RIGIDBODY> link)
{
  // compute the forward dynamics, given the algorithm
  switch (algorithm_type)
  {
    case eFeatherstone:
      _fsab.apply_impulse(w, link);
      break;

    case eCRB:
      _crb.apply_impulse(w, link);
      break;

    default:
      assert(false);
  }
}

/// Gets the generalized coordinates of this body
SHAREDVECTORN& RC_ARTICULATED_BODY::get_generalized_coordinates_euler(SHAREDVECTORN& gc)
{
  get_generalized_coordinates_euler_generic(gc);
  return gc;
}

/// Sets the generalized position of this body
void RC_ARTICULATED_BODY::set_generalized_coordinates_euler(const SHAREDVECTORN& gc)
{
  set_generalized_coordinates_euler_generic(gc);
}

/// Sets the generalized velocity of this body
void RC_ARTICULATED_BODY::set_generalized_velocity(DYNAMIC_BODY::GeneralizedCoordinateType gctype, const SHAREDVECTORN& gv)
{
  set_generalized_velocity_generic(gctype, gv);
}

/// Gets the generalized velocity of this body
SHAREDVECTORN& RC_ARTICULATED_BODY::get_generalized_velocity(DYNAMIC_BODY::GeneralizedCoordinateType gctype, SHAREDVECTORN& gv)
{
  get_generalized_velocity_generic(gctype, gv);
  return gv;
}

/// Gets the generalized inertia of this body
SHAREDMATRIXN& RC_ARTICULATED_BODY::get_generalized_inertia(SHAREDMATRIXN& M)
{
  // calculate the generalized inertia matrix
  _crb.calc_generalized_inertia(M);

  return M;
}

/// Gets the number of degrees of freedom for implicit constraints
unsigned RC_ARTICULATED_BODY::num_joint_dof_implicit() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _ijoints.size(); i++)
    ndof += _ijoints[i]->num_dof();
  return ndof;
}

/// Gets the generalized forces on this body
/**
 * \note does not add forces from implicit joint constraints!
 */
SHAREDVECTORN& RC_ARTICULATED_BODY::get_generalized_forces(SHAREDVECTORN& f)
{
  const unsigned SPATIAL_DIM = 6, X = 0, Y = 1, Z = 2, A = 3, B = 4, C = 5;

  // compute the generalized forces
  SFORCE f0;
  VECTORN CmQ;
  _crb.calc_generalized_forces(f0, CmQ);

  // determine the vector of joint forces
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    unsigned idx = _ejoints[i]->get_coord_index();
    for (unsigned m=0; m< _ejoints[i]->num_dof(); m++, idx++)
      f[idx] = _ejoints[i]->force[m] - CmQ[idx];
  }

  // setup joint space part of f
  if (_floating_base)
  {
    // determine external and inertial forces on base
    shared_ptr<RIGIDBODY> base = _links.front();
    unsigned idx = CmQ.size();
    SFORCE fb = POSE3::transform(base->get_gc_pose(), f0);
    f[idx++] = -fb[X];  f[idx++] = -fb[Y];  f[idx++] = -fb[Z];
    f[idx++] = -fb[A];  f[idx++] = -fb[B];  f[idx++] = -fb[C];
  }

  return f;
}

/// Converts a force to a generalized force
SHAREDVECTORN& RC_ARTICULATED_BODY::convert_to_generalized_force(shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf)
{
  const unsigned SPATIAL_DIM = 6;
  static vector<SVELOCITY> J;
  static vector<SVELOCITY> sprime;

  // get the body as a rigid body
  shared_ptr<RIGIDBODY> link = dynamic_pointer_cast<RIGIDBODY>(body);
  assert(link);

  // get the gc frame
  shared_ptr<const POSE3> P = _links.front()->get_gc_pose();

  // get w in the mixed frame
  SFORCE wP = POSE3::transform(P, w);

  // clear the Jacobian
  J.resize(num_joint_dof_explicit());
  for (unsigned i=0; i< J.size(); i++)
    J[i] = SVELOCITY::zero(P);

  // compute the Jacobian in w's frame
  for (unsigned i=0; i< _ejoints.size(); i++)
  {
    // if link is not a descendent of this joint's outboard, keep looping
    shared_ptr<RIGIDBODY> outboard = _ejoints[i]->get_outboard_link();
    if (!outboard->is_descendant_link(link))
      continue;

    // transform the Jacobian
    const vector<SVELOCITY>& s = _ejoints[i]->get_spatial_axes();
    POSE3::transform(P, s, sprime);
    for (unsigned j=0, k=_ejoints[i]->get_coord_index(); j < s.size(); j++, k++)
      J[k] = sprime[j];
  }

  // resize gf
  gf.resize(num_generalized_coordinates(DYNAMIC_BODY::eSpatial));

  // get the torque on the joints
  SHAREDVECTORN jf = gf.segment(0, J.size());
  SPARITH::transpose_mult(J, wP, jf);

  // determine the generalized force on the base, if the base is floating
  if (_floating_base)
  {
    shared_ptr<RIGIDBODY> base = _links.front();
    SHAREDVECTORN gfbase = gf.segment(J.size(), gf.size());
    wP.to_vector(gfbase);
  }

  // return the generalized force vector
  return gf;
}

/// Determines whether a joint supports a link
bool RC_ARTICULATED_BODY::supports(boost::shared_ptr<JOINT> joint, shared_ptr<RIGIDBODY> link)
{
  // save the original link
  shared_ptr<RIGIDBODY> l = link;

  do
  {
    boost::shared_ptr<JOINT> j = l->get_inner_joint_explicit();
    // check for l is base
    if (!j)
      return false;
    if (j == joint)
      return true;

    // proceed up the chain
    l = j->get_inboard_link();
  }
  while (link != l);

  return false;
}



