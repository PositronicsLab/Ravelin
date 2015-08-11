/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;
using std::vector;
using std::map;
using std::endl;
using std::queue;
using std::list;

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics
/**
 * Computed joint actuator forces are stored in inv_dyn_data.
 */
map<shared_ptr<JOINT>, VECTORN> RNE_ALGORITHM::calc_inv_dyn(shared_ptr<RC_ARTICULATED_BODY> body, const map<shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data)
{
  if (!body->is_floating_base())
    return calc_inv_dyn_fixed_base(body, inv_dyn_data);
  else
    return calc_inv_dyn_floating_base(body, inv_dyn_data);
}

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a fixed base
/**
 * Computed joint actuator forces are stored in inv_dyn_data.
 */
map<shared_ptr<JOINT>, VECTORN> RNE_ALGORITHM::calc_inv_dyn_fixed_base(shared_ptr<RC_ARTICULATED_BODY> body, const map<shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  map<shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>::const_iterator idd_iter;
  vector<SVELOCITY> sprime;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_fixed_base() entered" << endl;

  // get the reference frame for computation
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // ** STEP 0: compute isolated inertias 

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // ** STEP 1: compute velocities and accelerations

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // setup the vector of link accelerations
  vector<SACCEL> a(links.size(), SACCEL::zero());

  // set frames for accelerations 
  for (unsigned i=0; i< links.size(); i++)
    a[i].pose = links[i]->get_computation_frame(); 
  
  // add all child links of the base to the processing queue
  list<shared_ptr<RIGIDBODY> > child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
    link_queue.push(rb);
  
  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // push all children of the link onto the queue
    child_links.clear();
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
      link_queue.push(rb);

    // get the link's parent
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    unsigned h = parent->get_index();

    // get the joint for this link
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());

    // get the spatial link velocity
    const SVELOCITY& v = link->get_velocity(); 

    // get the current joint velocity
    const VECTORN& qd = joint->qd;

    // **** compute acceleration

    // get the desired joint acceleration
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VECTORN& qdd_des = idd_iter->second.qdd;  

    // get the spatial axes and time derivative
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    const vector<SVELOCITY>& sdot = joint->get_spatial_axes_dot();

    // put s into the proper frame (that of v/a) (if necessary)
    POSE3::transform(v.pose, s, sprime);  

    // add this link's contribution
    SVELOCITY sqd = SPARITH::mult(sprime, qd);
    a[i] += SACCEL(v.cross(sqd));
    a[i] += SACCEL(SPARITH::mult(s, qdd_des));

    // put s into the proper frame (that of v/a) (if necessary)
    POSE3::transform(v.pose, sdot, sprime);  
    a[i] += SACCEL(SPARITH::mult(sprime, qd));

    // now add parent's contribution
    a[i] += SPARITH::transform_accel(a[i].pose, a[h]);

    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->id << endl;
//    FILE_LOG(LOG_DYNAMICS) << "  spatial axes: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << v << endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << a[i] << endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<SFORCE> f(links.size(), SFORCE::zero());

  // all forces should be in the same frame as the individual links
  for (unsigned i=0; i< links.size(); i++)
    f[i].pose = links[i]->get_computation_frame();

   // add all leaf links to the queue
  for (unsigned i=1; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // if this link has already been processed, do not process it again
    if (processed[i])
      continue;

    // if the link is the base, continue the loop
    if (link->is_base())
      continue;
    
    // link is not the base; add the parent to the queue for processing
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    link_queue.push(parent);
    unsigned h = parent->get_index();
   
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << f[i] << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (link->get_inertia() * a[i]) << endl;

    // get the spatial velocity for this link
    const SVELOCITY& v = link->get_velocity();

    // add I*a to the link force
    f[i] += link->get_inertia() * a[i];

    // update the determined force to this link w/Coriolis + centrifugal terms
    f[i] += v.cross(link->get_inertia() * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << f[i] << endl;    

    // determine external forces in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const SFORCE& fx = idd_iter->second.wext;  

    // add in external forces
    f[i] -= fx; 

    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << f[i] << endl;

    // indicate that this link has been processed
    processed[i] = true;

    // update the parent force, if the parent is not the base
    if (parent->is_base())
      continue;
    else 
      f[h] += POSE3::transform(f[h].pose, f[i]); 
  }
  
  // ** STEP 3: compute joint forces

  // setup a map from joints to actuator forces
  map<shared_ptr<JOINT>, VECTORN> actuator_forces;

  // compute actuator forces
  for (unsigned j=1; j< links.size(); j++)
  {
    shared_ptr<RIGIDBODY> link = links[j];
    const unsigned i = link->get_index();
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    VECTORN& Q = actuator_forces[joint]; 
    SFORCE w = POSE3::transform(joint->get_pose(), f[i]); 
    SPARITH::transpose_mult(s, w, Q);
  
    FILE_LOG(LOG_DYNAMICS) << "joint " << joint->id << " inner joint force: " << actuator_forces[joint] << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_fixed_base() exited" << endl;

  return actuator_forces;
}

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a fixed base
/**
 * \pre Uses joint accelerations computed by forward dynamics, so the 
 *      appropriate forward dynamics function must be run first.
 */
/*
void RNE_ALGORITHM::calc_constraint_forces(shared_ptr<RC_ARTICULATED_BODY> body)
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  vector<SPATIAL_RB_INERTIA> Iiso;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_constraint_forces() entered" << endl;

  // get the reference frame for computation
  ReferenceFrameType rftype = body->computation_frame_type;

  // ** STEP 0: compute isolated inertias 

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // get the isolated inertiae
  Iiso.resize(links.size());
  for (unsigned i=1; i< links.size(); i++)
  {
    unsigned idx = links[i]->get_index();
    Iiso[idx] = links[i]->get_spatial_iso_inertia(rftype); 
  }

   // ** STEP 1: compute velocities and accelerations

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // setup the vector of link accelerations
  vector<SACCEL> accels(links.size(), SACCEL::zero());
  
  // add all child links of the base to the processing queue
  list<shared_ptr<RIGIDBODY> > child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
    link_queue.push(rb);

  // ** STEP 1: compute link forces -- backward recursion
  vector<bool> processed(links.size(), false);
  vector<SFORCE> forces(links.size(), SForce::zero());

  // add all leaf links to the queue
  for (unsigned i=1; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue and push all children of the link onto the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned idx = link->get_index();

    // if this link has already been processed, do not process it again
    if (processed[idx])
      continue;

    // if the link is the base, continue the loop
    if (link->is_base())
      continue;
    
    // link is not the base; add the parent to the queue for processing
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    link_queue.push(parent);
    unsigned pidx = parent->get_index();

    // get the forces for this link and this link's parent
    SFORCE& fi = forces[idx];
    SFORCE& fim1 = forces[pidx];

    // get this link's acceleration
    SACCEL& a = link->get_accel(rftype);
    
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << fi << endl;    
    FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (Iiso[idx] * a) << endl;

    // get the spatial velocity for this link
    const SVELOCITY& v = link->velocity();

    // add I*a to the link force
    fi += Iiso[idx] * a;

    // update the determined force to this link w/Coriolis + centrifugal terms
    fi += v.cross(Iiso[idx] * v);

    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a): " << fi << endl;    

    // determine external forces in link frame
    const SFORCE& wext = link->sum_force(); 
    shared_ptr<const POSE3>& T = link->get_pose();
    SFORCE fx(T.transpose_mult_vector(fext), T.transpose_mult_vector(text));

    // subtract external forces in the appropriate frame
    if (rftype == eGlobal)
    {
      SpatialTransform X_0_i = link->get_spatial_transform_link_to_global();
      fi -= X_0_i.transform(fx);
    }
    else
      fi -= fx;

    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << fi << endl;

    // indicate that this link has been processed
    processed[idx] = true;

    // update the parent force, if the parent is not the base
    if (parent->is_base())
      continue;
    else
    { 
      if (rftype == eGlobal)
        fim1 += fi;
      else
        fim1 += link->get_spatial_transform_backward().transform(fi);
    }
  }
  
  // ** STEP 2: compute constraint forces

  // compute actuator forces
  for (unsigned i=1; i< links.size(); i++)
  {
    shared_ptr<RIGIDBODY> link = links[i];
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    joint->get_spatial_constraints(rftype, s);
    transpose_mult(s, forces[link->get_index()], joint->lambda);
  
    FILE_LOG(LOG_DYNAMICS) << "joint " << joint->id << " constraint force: " << joint->lambda << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_constraint_forces() exited" << endl;
}
*/

/// Executes the Recursive Newton-Euler algorithm for inverse dynamics for a floating base
/**
 * \param inv_dyn_data a mapping from links to the external forces (and
 *        torques) applied to them and to the desired inner joint
 *        accelerations; note that all links in the robot should be included
 *        in this map (even the base link, although inner joint acceleration
 *        is not applicable in that case and will be ignored for it)
 * \return a mapping from joints to actuator forces
 */
map<shared_ptr<JOINT>, VECTORN> RNE_ALGORITHM::calc_inv_dyn_floating_base(shared_ptr<RC_ARTICULATED_BODY> body, const map<shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data) const
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  map<shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>::const_iterator idd_iter;
  vector<SPATIAL_RB_INERTIA> I;
  vector<SVELOCITY> v;
  vector<SACCEL> a;
  vector<SFORCE> Z;
  vector<SVELOCITY> sprime;

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() entered" << endl;

  // get the reference frame type
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // ** STEP 0: compute isolated inertias 

  // ** STEP 1: compute velocities and relative accelerations

  // set all desired velocities and accelerations (latter relative to the base)
  // to zero initially
  v.resize(links.size());
  a.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    v[i] = SVELOCITY::zero();
    a[i] = SACCEL::zero();
    v[i].pose = a[i].pose = links[i]->get_computation_frame();
  }  

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();
  
  // set velocity for the base
  v.front() = base->get_velocity();

  // add all child links of the base to the processing queue
  list<shared_ptr<RIGIDBODY> > child_links;
  base->get_child_links(std::back_inserter(child_links)); 
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
    link_queue.push(rb);
    
  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();
    
    // add all child links of the link to the processing queue
    child_links.clear();
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
      link_queue.push(rb);
    
    // get the parent link and inner joint
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    
    // get the index of this link and its parent
    const unsigned i = link->get_index();
    const unsigned h = parent->get_index();

    // get spatial axes and derivatives
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    const vector<SVELOCITY>& sdot = joint->get_spatial_axes_dot();

    // put s into the proper frame (that of v/a) 
    POSE3::transform(link->get_computation_frame(), s, sprime);

    // compute s * qdot
    SVELOCITY sqd = SPARITH::mult(sprime, joint->qd);
    
    // get the desired acceleration for the current link
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const VECTORN& qdd_des = idd_iter->second.qdd;

    // compute parent contributions to velocity and relative acceleration
    v[i] = POSE3::transform(v[i].pose, v[h]);

    // compute velocity and relative acceleration
    v[i] += sqd;
    a[i] += SACCEL(SPARITH::mult(sprime, qdd_des) + v[i].cross(sqd));

    // compute time derivative of spatial axes contributions (if any)
    POSE3::transform(link->get_computation_frame(), sdot, sprime);
    a[i] += SACCEL(SPARITH::mult(sprime, joint->qd));

//    FILE_LOG(LOG_DYNAMICS) << "  s: " << s << endl;
    FILE_LOG(LOG_DYNAMICS) << "  velocity for link " << links[i]->id << ": " << v[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  relative accel for link " << links[i]->id << ": " << a[i] << endl;
  }
  
  // ** STEP 2: compute composite inertias and Z.A. forces

  // resize vectors of spatial inertias and Z.A. vectors
  I.resize(links.size());
  Z.resize(links.size());

  // zero out I and Z and set poses
  for (unsigned i=0; i< links.size(); i++)
  {
    I[i].set_zero();
    Z[i].set_zero();
    I[i].pose = Z[i].pose = links[i]->get_computation_frame();
  }

  // set all spatial isolated inertias and Z.A. forces
  for (unsigned j=0; j< links.size(); j++)
  {
    // get the i'th link
    shared_ptr<RIGIDBODY> link = links[j];
    const unsigned i = link->get_index();

    // add the spatial isolated inertia for this link to the composite inertia
    I[i] += link->get_inertia();

    // setup forces due to (relative) acceleration on link
    Z[i] = link->get_inertia() * a[i];

    // add coriolis and centrifugal forces on link
    Z[i] += v[i].cross(link->get_inertia() * v[i]);

    // determine external forces on the link in link frame
    idd_iter = inv_dyn_data.find(link);
    assert(idd_iter != inv_dyn_data.end());
    const SFORCE& wx = idd_iter->second.wext;

    // subtract external force from Z.A. vector
    Z[i] -= wx;

    FILE_LOG(LOG_DYNAMICS) << " -- processing link " << link->id << endl;
    FILE_LOG(LOG_DYNAMICS) << "   external force: " << wx << endl;
    FILE_LOG(LOG_DYNAMICS) << "   ZA vector: " << Z[i] << endl;
  }
  
  // *** compute composite inertias and zero acceleration vectors

  // setup vector that indicates when links have been processed
  vector<bool> processed(links.size(), false);

  // put all leaf links into the queue
  for (unsigned i=0; i< links.size(); i++)
    if (links[i]->num_child_links() == 0)
      link_queue.push(links[i]);

  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the index for this link
    const unsigned i = link->get_index();
    
    // see whether this link has already been processed
    if (processed[i])
      continue;
    
    // process the parent link, if possible
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (parent)
    {
      // put the parent on the queue
      link_queue.push(parent);
    
      // get the parent index
      const unsigned h = parent->get_index();
    
      // add this inertia and Z.A. force to its parent
      I[h] += POSE3::transform(I[h].pose, I[i]);
      Z[h] += POSE3::transform(Z[h].pose, Z[i]);

      // indicate that the link has been processed
      processed[i] = true;
    }
  }

  // ** STEP 3: compute base acceleration
  a.front() = I.front().inverse_mult(-Z.front());
  
  FILE_LOG(LOG_DYNAMICS) << "  Composite inertia for the base: " << endl << I.front();
  FILE_LOG(LOG_DYNAMICS) << "  ZA vector for the base: " << Z.front() << endl;
  FILE_LOG(LOG_DYNAMICS) << "  Determined base acceleration: " << a.front() << endl;

  // ** STEP 4: compute joint forces
  
  // setup the map of actuator forces
  map<shared_ptr<JOINT>, VECTORN> actuator_forces;

  // compute the forces
  for (unsigned j=1; j< links.size(); j++)
  {
    const unsigned i = links[j]->get_index();
    shared_ptr<JOINT> joint(links[j]->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    VECTORN& Q = actuator_forces[joint];
    SFORCE w = I[i] * a.front() + Z[i];
    SPARITH::transpose_mult(s, POSE3::transform(joint->get_pose(), w), Q);

    FILE_LOG(LOG_DYNAMICS) << "  processing link: " << links[j]->id << endl;
//    FILE_LOG(LOG_DYNAMICS) << "    spatial axis: " << endl << s;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << endl << I[i];
    FILE_LOG(LOG_DYNAMICS) << "    Z: " << endl << Z[i];
    FILE_LOG(LOG_DYNAMICS) << "    actuator force: " << actuator_forces[joint] << endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "RNEAlgorithm::calc_inv_dyn_floating_base() exited" << endl;

  return actuator_forces;
}

