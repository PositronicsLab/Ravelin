/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using std::set;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::list;
using std::vector;
using std::map;
using std::string;
using std::queue;

ARTICULATEDBODY::ARTICULATEDBODY()
{
}

/// "Compiles" the articulated body
void ARTICULATEDBODY::compile()
{
}

/// Gets the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose 
MATRIXN& ARTICULATEDBODY::calc_jacobian_dot(boost::shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMICBODY> body, MATRIXN& J)
{
  const shared_ptr<const POSE3> GLOBAL;

  // get the generalized coordinate pose
  if (is_floating_base())
  {
    boost::shared_ptr<const POSE3> gc_pose = get_base_link()->get_pose();
    return calc_jacobian_dot(gc_pose, target_pose, body, J);
  }
  else
    return calc_jacobian_dot(GLOBAL, target_pose, body, J);
}

/// Gets the time derivative of the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose 
MATRIXN& ARTICULATEDBODY::calc_jacobian_dot(shared_ptr<const POSE3> source_pose, shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMICBODY> body, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NEXP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NEXP_DOF + SPATIAL_DIM : NEXP_DOF;

  // setup the Jacobian
  J.set_zero(SPATIAL_DIM, NDOF); 

  // get the current link
  shared_ptr<RIGIDBODY> link = dynamic_pointer_cast<RIGIDBODY>(body);
  
  // get the base link
  shared_ptr<RIGIDBODY> base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    shared_ptr<JOINT> joint = link->get_inner_joint_explicit();

    // get the parent link
    shared_ptr<RIGIDBODY> parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVELOCITY>& s = joint->get_spatial_axes_dot();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SHAREDVECTORN v = J.column(CIDX+i);
      POSE3::transform(target_pose, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  // NOTE: we do not even check for a floating base, because the 
  // time-derivative of its Jacobian will always be zero

  return J;
}

/// Gets the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose 
MATRIXN& ARTICULATEDBODY::calc_jacobian(boost::shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMICBODY> body, MATRIXN& J)
{
  const shared_ptr<const POSE3> GLOBAL;

  // get the generalized coordinate pose
  if (is_floating_base())
  {
    boost::shared_ptr<const POSE3> gc_pose = get_base_link()->get_pose();
    return calc_jacobian(gc_pose, target_pose, body, J);
  }
  else
    return calc_jacobian(GLOBAL, target_pose, body, J);
}

/// Gets the Jacobian that converts velocities from this body in the source pose to velocities of the particular link in the target pose 
MATRIXN& ARTICULATEDBODY::calc_jacobian(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, shared_ptr<DYNAMICBODY> body, MATRIXN& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NEXP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NEXP_DOF + SPATIAL_DIM : NEXP_DOF;

  // setup the Jacobian
  J.set_zero(SPATIAL_DIM, NDOF); 

  // get the current link
  shared_ptr<RIGIDBODY> link = dynamic_pointer_cast<RIGIDBODY>(body);
  
  // get the base link
  shared_ptr<RIGIDBODY> base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    shared_ptr<JOINT> joint = link->get_inner_joint_explicit();

    // get the parent link
    shared_ptr<RIGIDBODY> parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVELOCITY>& s = joint->get_spatial_axes();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SHAREDVECTORN v = J.column(CIDX+i);
      POSE3::transform(target_pose, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  // if base is floating, setup Jacobian columns at the end
  if (is_floating_base())
  {
    SHAREDMATRIXN Jbase = J.block(0, SPATIAL_DIM, NEXP_DOF, NEXP_DOF+SPATIAL_DIM);
    POSE3::spatial_transform_to_matrix2(source_pose, target_pose, Jbase);
  }

  return J;
}

/// Determines the loop indices corresponding to each joint and the vector of links for each joint
void ARTICULATEDBODY::find_loops(vector<unsigned>& loop_indices, vector<vector<unsigned> >& loop_links) const
{
  vector<shared_ptr<JOINT> > loop_joints, implicit_joints;
  queue<shared_ptr<RIGIDBODY> > q;

  // clear vectors
  loop_indices.resize(_joints.size());
  implicit_joints.clear();

  // get all implicit joints
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == JOINT::eImplicit)
      implicit_joints.push_back(_joints[i]);

  // set all loop indices to INF (indicates no loop) initially
  for (unsigned i=0; i< _joints.size(); i++)
    loop_indices[i] = std::numeric_limits<unsigned>::max();

  // look for early exit
  if (implicit_joints.empty())
    return;

  // two cases: 1) body uses *only* implicit joints and 2) body uses 
  // explicit and implicit joints
  if (_joints.size() == implicit_joints.size())
  {
    // we're going to reset implicit_joints to hold only the joints that
    // complete loops
    implicit_joints.clear();
    for (unsigned i=0; i< _joints.size(); i++)
    {
      shared_ptr<RIGIDBODY> inboard = _joints[i]->get_inboard_link();
      shared_ptr<RIGIDBODY> outboard = _joints[i]->get_outboard_link();

      // check for obvious loop closure
      if (inboard->get_index() > outboard->get_index())
      {
        implicit_joints.push_back(_joints[i]);
        continue;
      }

      // check for not-so-obvious loop closure
      if (!outboard->is_enabled())
      {
        // outboard is fixed; look back to see whether one of the predecessor
        // links is fixed as well
        q.push(inboard);
        while (!q.empty())
        {
          shared_ptr<RIGIDBODY> link = q.front();
          q.pop();
          if (!link->is_enabled())
          {
            implicit_joints.push_back(_joints[i]);
            break;
          }
          const set<shared_ptr<JOINT> >& ij = link->get_inner_joints();
          BOOST_FOREACH(shared_ptr<JOINT> j, ij)
            q.push(shared_ptr<RIGIDBODY>(j->get_inboard_link()));
         }
       }
     }
  }

  // reset loop links
  loop_links.clear();
  loop_links.resize(implicit_joints.size());

  // for every kinematic loop
  for (unsigned k=0; k< implicit_joints.size(); k++)
  {
    // get the implicit joint
    shared_ptr<JOINT> ejoint = implicit_joints[k];
    shared_ptr<RIGIDBODY> outboard = ejoint->get_outboard_link();
    bool ground_outboard = outboard->is_ground();

    // determine all joints and links in the loop by iterating backward until
    // we get back to the first link in the loop
    loop_joints.clear();
    loop_joints.push_back(ejoint);
    loop_links[k].push_back(outboard->get_index());
    shared_ptr<RIGIDBODY> inboard = ejoint->get_inboard_link();
    while (true)
    {
      shared_ptr<JOINT> jx = inboard->get_inner_joint_explicit();
      loop_joints.push_back(jx);
      loop_links[k].push_back(inboard->get_index());
      inboard = jx->get_inboard_link();
 
      // check for loop termination
      if (inboard == outboard)
        break;
      if (ground_outboard && inboard->is_ground())
        break;
    }

    // reverse the vector of links so that it is (almost) sorted (all but
    // last link)
    std::reverse(loop_links[k].begin(), loop_links[k].end());
    #ifndef NDEBUG
    for (unsigned i=1; i< loop_links[k].size()-1; i++)
      assert(loop_links[k][i] > loop_links[k][i-1]);
    #endif

    // setup loop indices for each joint in the loop
    for (unsigned i=0; i< loop_joints.size(); i++)
      loop_indices[loop_joints[i]->get_index()] = k;
  }
}

/// Sets the vectors of links and joints
void ARTICULATEDBODY::set_links_and_joints(const vector<shared_ptr<RIGIDBODY> >& links, const vector<shared_ptr<JOINT> >& joints)
{
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
  {
    _joints[i]->set_index(i);
    _joints[i]->set_articulated_body(get_this());
  }

  compile();
}

/// Gets the number of explicit joint constraint equations
unsigned ARTICULATEDBODY::num_constraint_eqns_explicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == JOINT::eExplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of implicit joint constraint equations
unsigned ARTICULATEDBODY::num_constraint_eqns_implicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == JOINT::eImplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of joint degrees of freedom permitted by both implicit and explicit joint constraints
unsigned ARTICULATEDBODY::num_joint_dof() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    ndof += _joints[i]->num_dof();
  return ndof;
}

/// Finds the joint with the given name
/**
 * \return NULL if the joint wasshared_ptr<void> not found
 */
shared_ptr<JOINT> ARTICULATEDBODY::find_joint(const string& jointname) const
{
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->id == jointname)
      return _joints[i];
      
  return shared_ptr<JOINT>();
}

/// Gets the adjacent links
void ARTICULATEDBODY::get_adjacent_links(list<sorted_pair<shared_ptr<RIGIDBODY> > >& links) const
{
  for (unsigned i=0; i< _joints.size(); i++)
  {
    shared_ptr<RIGIDBODY> ib = _joints[i]->get_inboard_link();
    shared_ptr<RIGIDBODY> ob = _joints[i]->get_outboard_link();
    if (ib && ob)
      links.push_back(make_sorted_pair(ib, ob));
  }
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ARTICULATEDBODY::translate(const ORIGIN3& x)
{
  // apply transform to all links
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, _links)
    rb->translate(x);
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ARTICULATEDBODY::rotate(const QUAT& q)
{
  // apply transform to all links
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, _links)
    rb->rotate(q);
}

/// Calculates the combined kinetic energy of all links in this body with respect to the base frame
REAL ARTICULATEDBODY::calc_kinetic_energy() 
{
  // TODO: fix this to do computation in the base frame
  REAL KE = 0;
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, _links)
    KE += rb->calc_kinetic_energy();

  return KE;
}

/// Resets force and torque accumulators on the body
/**
 * Force and torque accumulators on all links are reset.
 */
void ARTICULATEDBODY::reset_accumulators()
{
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, _links)
    rb->reset_accumulators();
}

/// Finds the link with the given name
shared_ptr<RIGIDBODY> ARTICULATEDBODY::find_link(const string& linkid) const
{
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, _links)
    if (rb->id == linkid)
      return rb;

  return shared_ptr<RIGIDBODY>();
}


