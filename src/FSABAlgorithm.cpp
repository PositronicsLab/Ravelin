/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;
using std::set;
using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;

FSAB_ALGORITHM::FSAB_ALGORITHM()
{
}

/// Solves the equation MX = B, where M is the generalized inertia matrix
/**
 * \param Y the matrix B on entry, the matrix X on return
 * \pre spatial inertias already computed for the body's current configuration 
 */
void FSAB_ALGORITHM::solve_generalized_inertia_noprecalc(SHAREDMATRIXN& Y)
{
  VECTORN gv;

  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // store the current generalized velocity
  body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, gv);

  // get the number of generalized coords
  const unsigned NGC = body->num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
  if (Y.rows() != NGC)
    throw MissizeException();

  // get the number of columns 
  const unsigned NC = Y.columns(); 

  // resize inv(M)
  _workv.set_zero(NGC);

  // solve using the matrix
  for (unsigned i=0; i< NC; i++)
  {
    // setup the generalized impulse
    Y.get_column(i, _workv);

    // get the current generalized velocity
    body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv2);

    // apply the generalized impulse
    apply_generalized_impulse(_workv);

    // get the new velocity out
    body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv);
    _workv -= _workv2;

    // set the appropriate column of Y 
    Y.set_column(i, _workv);
  } 

  // restore the current generalized velocity
  body->set_generalized_velocity(DYNAMIC_BODY::eSpatial, gv);
}

/// Solves the equation Mx = b, where M is the generalized inertia matrix
/**
 * \param v the vector b on entry, the vector x on return
 * \pre spatial inertias already computed for the body's current configuration 
 */
void FSAB_ALGORITHM::solve_generalized_inertia_noprecalc(SHAREDVECTORN& v)
{
  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // store the current generalized velocity
  body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv);
  if (v.rows() != _workv.rows())
    throw MissizeException();

  // apply the generalized impulse
  _workv2 = v;
  apply_generalized_impulse(_workv2);

  // get the new velocity out
  body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv2);
  _workv2 -= _workv;

  // restore the current generalized velocity
  body->set_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv);

  // store the change in velocity
  v = _workv2;
}

/// Calculates the inverse generalized inertia matrix
void FSAB_ALGORITHM::calc_inverse_generalized_inertia_noprecalc(MATRIXN& iM)
{
  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // store the current generalized velocity
  body->get_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv2);

  // get the number of generalized coords
  const unsigned NGC = body->num_generalized_coordinates(DYNAMIC_BODY::eSpatial);

  // resize inv(M)
  iM.resize(NGC, NGC);
  _workv.set_zero(NGC);

  // compute the inverse inertia
  for (unsigned i=0; i< NGC; i++)
  {
    // setup the generalized impulse
    _workv.set_zero(NGC);
    _workv[i] = (REAL) 1.0;

    // apply the generalized impulse
    apply_generalized_impulse(i, _workv);

    // set the appropriate column of inv(M)
    iM.set_column(i, _workv);
  } 

  // make inverse(M) symmetric
  for (unsigned i=0; i< NGC; i++)
  {
    CONST_ROW_ITERATOR source = iM.row(i).row_iterator_begin();
    COLUMN_ITERATOR target = iM.column(i).column_iterator_begin();
    for (unsigned j=0; j< i; j++)
      *target++ = *source++;
  }

  // restore the current generalized velocity
  body->set_generalized_velocity(DYNAMIC_BODY::eSpatial, _workv2);
 
  FILE_LOG(LOG_DYNAMICS) << "inverse M: " << std::endl << iM;
}

/// Applies a generalized impulse using the algorithm of Drumwright
void FSAB_ALGORITHM::apply_generalized_impulse(unsigned index, VECTORN& vgj)
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  const unsigned SPATIAL_DIM = 6;
  const unsigned BASE_IDX = 0;
  vector<SVELOCITY> sprime;

  // get the body as well as sets of links and joints
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<boost::shared_ptr<JOINT> >& joints = body->get_explicit_joints();

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = (body->is_floating_base()) ? SPATIAL_DIM : 0;

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << vgj << endl;

  // clear values for vectors
  _Y.resize(links.size());
  _processed.resize(links.size());

  // reset Y and processed
  for (unsigned j=0; j< links.size(); j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
    _processed[i] = false;
  }

  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  for (unsigned i=0; i< links.size(); i++)
    if (i > BASE_IDX && links[i]->num_child_links() == 0)
      link_queue.push(links[i]);

  // backward recursion
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();
    unsigned i = link->get_index();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (_processed[i])
      continue;
  
    // indicate that this link has been processed
    _processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (parent->get_index() != BASE_IDX)
    {
      link_queue.push(parent);
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent << " to queue for processing" << endl;
    }
    unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();

    // get I, Y, and mu
    const SPATIAL_AB_INERTIA& I = _I[i];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && h == BASE_IDX)
      continue;

    // don't update parent Y using joint DOF if there is no joint DOF
    if (joint->num_dof() == 0)
    {
      _Y[h] += POSE3::transform(_Y[h].pose, _Y[i]);
      continue;
    } 

    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    vgj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[i]);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu[i] << endl;
    
    // compute the qm subexpression
    POSE3::transform(_Y[i].pose, s, sprime);
    _mu[i] -= SPARITH::transpose_mult(sprime, _Y[i], _workv2);

    // update parent impulsive force 
    const vector<SMOMENTUM>& Is = _Is[i];
    solve_sIs(i, _mu[i], _sIsmu);
    SPARITH::mult(Is, _sIsmu, _workv2);
    SMOMENTUM Y = _Y[i] + SMOMENTUM::from_vector(_workv2, _Y[i].pose);
    _Y[h] += POSE3::transform(_Y[h].pose, Y);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    parent Y is now: " << endl << _Y[h] << endl;
  }

  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // update base components 
    shared_ptr<RIGIDBODY> base = links.front();

    // get the start of the generalized coordinates for the base 
    const unsigned S = body->_n_joint_DOF_explicit;

    // generalized impulse on base will be in mixed frame, by Moby convention 
    SMOMENTUM w(vgj[S+0], vgj[S+1], vgj[S+2], vgj[S+3], vgj[S+4], vgj[S+5], base->get_gc_pose());
    _Y.front() -= POSE3::transform(_Y.front().pose, w);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
    FILE_LOG(LOG_DYNAMICS) << endl;
    for (unsigned i=st_idx; i< links.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << "Articulated zero-velocity delta vector for " << links[i] << ": " << _Y[links[i]->get_index()] << endl;  
    FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_velocity_deltas() ended" << endl;
  }

  // ************************************************************
  // determine the changes in joint and link velocities
  // ************************************************************

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();
  const unsigned NUM_LINKS = links.size();

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  change in velocity: " << _dv.front() << endl;

    // figure out where the base starts
    const unsigned BASE_START = body->_n_joint_DOF_explicit;

    // determine change in generalized velocities -- NOTE: we have to reverse
    // linear and angular components
    vgj[BASE_START+0] = _dv.front()[3+0];
    vgj[BASE_START+1] = _dv.front()[3+1];
    vgj[BASE_START+2] = _dv.front()[3+2];
    vgj[BASE_START+3] = _dv.front()[0+0];
    vgj[BASE_START+4] = _dv.front()[0+1];
    vgj[BASE_START+5] = _dv.front()[0+2];
  }
  else 
    _dv.front().set_zero();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the parent link
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());

    // place all children on the link queue
    push_children(link, link_queue);

    // check the starting index of this joint
    const unsigned CSTART = joint->get_coord_index();
    if (CSTART > index)
      continue; 
    
    // get spatial axes of the inner joint for link i
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent << endl;

    // do the first update of the link velocity`
    _dv[i] = POSE3::transform(_dv[i].pose, _dv[h]);

    // check whether we can do a simple update of link velocity
    if (joint->num_dof() == 0)
      continue;
    
    // transform s to Y's frame
    POSE3::transform(_Y[i].pose, s, sprime); 

    // compute dq = inv(s'*I*s) * (mu - compute s'*I*dv[h] + s'*Y)
    SMOMENTUM Y = _Y[i] + _I[i]*POSE3::transform(_dv[i].pose, _dv[h]); 
    SPARITH::transpose_mult(sprime, Y, _workv2);
    _workv2.negate();
    _workv2 += _mu[i]; 
    solve_sIs(i, _workv2, _qd_delta);

    // update the link velocity   
    _dv[i] += SPARITH::mult(sprime, _qd_delta);

    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- I: " << endl << _I[i];
    FILE_LOG(LOG_DYNAMICS) << "    -- Qi: " << _mu[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- Qi - (s'I*dv + s'Y): " << endl << _workv2;
    FILE_LOG(LOG_DYNAMICS) << "    -- [Qi - (s'I*dv + s'Y)]/s'Is = qd_delta: " << _qd_delta << endl; 
    FILE_LOG(LOG_DYNAMICS) << "    -- dv[parent]: " << _dv[h] << endl;    
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // update vgj
    for (unsigned k=0; k< joint->num_dof(); k++)
      vgj[CSTART+k] = _qd_delta[k];
  }

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_generalized_impulse() exited" << endl;
}

/// Applies a generalized impulse using the algorithm of Drumwright
void FSAB_ALGORITHM::apply_generalized_impulse(const VECTORN& gj)
{
  static VECTORN tmp, tmp2;
  static vector<VECTORN> mu;
  queue<shared_ptr<RIGIDBODY> > link_queue;
  vector<SVELOCITY> sprime;

  // determine the number of generalized coordinates for the base
  const unsigned N_BASE_GC = 6;

  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_generalized_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << "gj: " << gj << endl;

  // get the sets of links and joints
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<boost::shared_ptr<JOINT> >& joints = body->get_explicit_joints();

  // init mu and Y 
  _mu.resize(links.size());
  _Y.resize(links.size());
  for (unsigned j=0; j< links.size(); j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
  }

  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
 
  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index(); 
 
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;
 
    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue;
 
    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
      FILE_LOG(LOG_DYNAMICS) << "added link " << parent << " to queue for processing" << endl;
    }
    unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();

    // get I
    const SPATIAL_AB_INERTIA& I = _I[i];

    // don't update parent Y for direct descendants of non-floating bases
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _mu[i]);
    FILE_LOG(LOG_DYNAMICS) << " gj subvector for this link: " << _mu[i] << endl;
    
    // compute the qm subexpression
    POSE3::transform(_Y[i].pose, s, sprime);
    _mu[i] -= SPARITH::transpose_mult(sprime, _Y[i], tmp2);

    // get Is
    const vector<SMOMENTUM>& Is = _Is[i];

    // prepare to update parent Y
    solve_sIs(i, _mu[i], _sIsmu);
    SMOMENTUM uY = _Y[i] + SMOMENTUM::from_vector(SPARITH::mult(Is, _sIsmu, tmp), _Y[i].pose);

    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << _mu[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Y: " << _Y[i] << endl;

    // update parent force
    _Y[h] += POSE3::transform(_Y[h].pose, _Y[i]); 
  }

  // if we're dealing with a floating base
  if (body->is_floating_base())
  {
    // update base components 
    shared_ptr<RIGIDBODY> base = links.front();

    // get the start of the generalized coordinates for the base
    const unsigned S = body->_n_joint_DOF_explicit;

    // momentum is in mixed frame by Moby convention 
    SMOMENTUM basew(gj[S+0], gj[S+1], gj[S+2], gj[S+3], gj[S+4], gj[S+5], base->get_gc_pose());
    _Y.front() -= POSE3::transform(_Y.front().pose, basew);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
    FILE_LOG(LOG_DYNAMICS) << endl;
    for (unsigned i=st_idx; i< links.size(); i++)
      FILE_LOG(LOG_DYNAMICS) << "Articulated zero-velocity delta vector for " << links[i] << ": " << _Y[links[i]->get_index()] << endl;  
    FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_velocity_deltas() ended" << endl;
  }

  // ************************************************************
  // determine the new joint and link velocities
  // ************************************************************

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();
  const unsigned NUM_LINKS = links.size();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base velocity: " << base->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  change in base velocity: " << _dv.front() << endl;

    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << _Y.front() << endl;

    // update the base velocity
    base->set_velocity(base->get_velocity() + _dv.front());

    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << base->get_velocity() << endl;
  }
  else 
    _dv.front().set_zero();
  
  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the parent link
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    
    // get spatial axes of the inner joint for link i
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    
    // get articulated body inertia for link i
    const SPATIAL_AB_INERTIA& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent << endl;
    
    // determine appropriate components of gj
    const unsigned CSTART = joint->get_coord_index(); 
    gj.get_sub_vec(CSTART,CSTART+joint->num_dof(), _Qi);

    // determine the joint and link velocity updates
    POSE3::transform(_Y[i].pose, s, sprime);
    SMOMENTUM w1 = _I[i] * POSE3::transform(_Y[i].pose, _dv[h]);
    SPARITH::transpose_mult(sprime, w1 + _Y[i], tmp2).negate();
    tmp2 += _Qi;
    solve_sIs(i, tmp2, _qd_delta);
    _dv[i] = POSE3::transform(_Y[i].pose, _dv[h]);
    if (joint->num_dof() > 0)
      _dv[i] += SPARITH::mult(sprime, _qd_delta);
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    link->set_velocity(link->get_velocity() + _dv[i]);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_generalized_impulse() exited" << endl;
}

/// Computes the combined spatial coriolis / centrifugal forces vectors for the body
void FSAB_ALGORITHM::calc_spatial_coriolis_vectors(shared_ptr<RC_ARTICULATED_BODY> body)
{
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_coriolis_vectors() entered" << endl;
  vector<SVELOCITY> sprime;

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // clear spatial values for all links
  _c.resize(links.size());

  // process all links except the base
  for (unsigned i=1; i< links.size(); i++)
  {
    // get the link
    shared_ptr<RIGIDBODY> link = links[i];
    unsigned idx = link->get_index();

    // get the link's joint
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());

    // get the spatial axis and transform it
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    POSE3::transform(link->get_computation_frame(), s, sprime);

    // compute the coriolis vector
    if (sprime.empty())
      _c[idx].set_zero(link->get_computation_frame());
    else
    {
      SVELOCITY sqd = SPARITH::mult(sprime, joint->qd);
      _c[idx] = link->get_velocity().cross(sqd);

      FILE_LOG(LOG_DYNAMICS) << "processing link: " << link << endl;
      FILE_LOG(LOG_DYNAMICS) << "v: " << link->get_velocity() << endl;
      FILE_LOG(LOG_DYNAMICS) << "s * qdot: " << sqd << endl;
      FILE_LOG(LOG_DYNAMICS) << "c: " << _c[idx] << endl;
    }
  }
}

/// Computes articulated body zero acceleration forces used for computing forward dynamics
void FSAB_ALGORITHM::calc_spatial_zero_accelerations(shared_ptr<RC_ARTICULATED_BODY> body)
{
  VECTORN tmp, workv;
  MATRIXN workM;
  queue<shared_ptr<RIGIDBODY> > link_queue;
  vector<SVELOCITY> sprime;

  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() entered" << endl;

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // clear spatial values for all links
  _Z.resize(links.size());
  _mu.resize(links.size());

  // ************************************************
  // form the isolated z.a. vectors
  // ************************************************

  // process all links (including the base, if it's floating)
  unsigned start_idx = (body->is_floating_base()) ? 0 : 1;

  // process all links except the base
  for (unsigned j=start_idx; j< links.size(); j++)
  {
    // get the link
    shared_ptr<RIGIDBODY> link = links[j];
    const unsigned i = link->get_index();

    // calculate inertial and Coriolis forces 
    const SVELOCITY& v = link->get_velocity();

    // set 6-dimensional spatial isolated zero-acceleration vector of link  
    _Z[i] = v.cross(link->get_inertia() * v) - link->sum_forces();
    
    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso ZA: " << endl << _Z[i] << endl;
  }
 
  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false; 

  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index(); 
 
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;
  
    // verify that all children have already been processed
    if (!body->all_children_processed(link))
      continue;

    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (!parent->is_base())
      link_pqueue.push(parent->get_index());
    const unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    POSE3::transform(link->get_computation_frame(), s, sprime);

    // get I, c, and Z
    const SPATIAL_AB_INERTIA& I = _I[i];
    const SACCEL& c = _c[i];
    const SFORCE& Z = _Z[i];
    
    // compute the qm subexpression
    _mu[i] = joint->force;
    _mu[i] -= SPARITH::transpose_mult(sprime, Z + I*c, workv);

    // get Is
    const vector<SMOMENTUM>& Is = _Is[i];

    // get the qm subexpression
    const VECTORN& mu = _mu[i];
  
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    parent link: " << parent << endl;
    if (LOGGING(LOG_DYNAMICS) && !sprime.empty()) 
    FILE_LOG(LOG_DYNAMICS) << "    s': " << sprime.front() << endl;  
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;
    if (!Is.empty())
      FILE_LOG(LOG_DYNAMICS) << "    Is: " << Is.front() << endl;
    FILE_LOG(LOG_DYNAMICS) << "    c: " << c << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm subexp: " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    recursive Z: " << Z << endl;

    // don't update Z for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // compute a couple of necessary matrices
    solve_sIs(i, mu, _sIsmu);
    SFORCE uZ = Z + (I*c) + SFORCE::from_vector(SPARITH::mult(Is, _sIsmu, workv), Z.pose);

    // update the parent zero acceleration
    _Z[h] += POSE3::transform(_Z[h].pose, uZ);
  }

  FILE_LOG(LOG_DYNAMICS) << endl;
  unsigned st_idx = (body->is_floating_base()) ? 0 : 1;
  for (unsigned i=st_idx; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << "Articulated zero-acceleration vector for " << links[i] << ": " << _Z[links[i]->get_index()] << endl;  
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() ended" << endl;
}

/// Computes articulated body inertias used for computing forward dynamics
void FSAB_ALGORITHM::calc_spatial_inertias(shared_ptr<RC_ARTICULATED_BODY> body)
{
  FILE_LOG(LOG_DYNAMICS) << "calc_spatial_zero_accelerations() entered" << endl;
  vector<SVELOCITY> sprime;
  MATRIXN tmp, tmp2, tmp3;

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // clear spatial values for all links
  _rank_deficient.resize(links.size());
  _I.resize(links.size());
  _Is.resize(links.size());
  _sIs.resize(links.size());
  _usIs.resize(links.size());
  _ssIs.resize(links.size());
  _vsIs.resize(links.size());

  // ************************************************
  // form the isolated spatial inertia vectors
  // ************************************************

  // process all links (including the base, if it's floating)
  unsigned start_idx = (body->is_floating_base()) ? 0 : 1;

  // process all links
  for (unsigned j=start_idx; j< links.size(); j++)
  {
    // get the link
    shared_ptr<RIGIDBODY> link = links[j];
    const unsigned i = link->get_index();

    // set the articulated body inertia for this link to be its isolated
    // spatial inertia (this will be updated in the phase below)
   _I[i] = link->get_inertia();

    // check for degenerate inertia
    #ifndef NDEBUG
    if (link->is_base() && body->is_floating_base() && 
        (link->get_inertia().m <= 0.0 || link->get_inertia().J.norm_inf() <= 0.0))
      throw std::runtime_error("Attempted to compute dynamics given degenerate inertia for a floating base body");
    #endif

    FILE_LOG(LOG_DYNAMICS) << "  processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    Link spatial iso inertia: " << endl << _I[i];
  }
 
  // indicate that no links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
 
  // doing a recursion backward from the end-effectors; add all leaf links to the link_queue
  std::priority_queue<unsigned> link_pqueue;
  for (unsigned i=0; i< links.size(); i++)
    if (!links[i]->is_base() && links[i]->num_child_links() == 0)
      link_pqueue.push(i);

  // backward recursion
  while (!link_pqueue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = links[link_pqueue.top()];
    link_pqueue.pop();
    const unsigned i = link->get_index();
  
    // see whether this link has already been processed (because two different children can have the same parent)
    if (body->_processed[i])
      continue;

    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue; 
  
    // indicate that this link has been processed
    body->_processed[i] = true;

    // push the parent of the link onto the queue, *unless* the parent is the base
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (!parent->is_base())
    {
      link_pqueue.push(parent->get_index());
FILE_LOG(LOG_DYNAMICS) << "added link " << parent << " to queue for processing" << endl;
    }
    const unsigned h = parent->get_index();

    // get the inner joint and the spatial axis
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    POSE3::transform(_I[i].pose, s, sprime);

    // get I
    const SPATIAL_AB_INERTIA& I = _I[i];
    
    // compute Is
    SPARITH::mult(I, sprime, _Is[i]);

    // compute sIs
    SPARITH::transpose_mult(sprime, _Is[i], _sIs[i]);

    // get whether s is rank deficient
    _rank_deficient[i] = joint->is_singular_config();

    // if the joint is not rank deficient, compute a Cholesky factorization 
    // of sIs
    if (_sIs[i].rows() == 1)
      _sIs[i].data()[0] = 1.0/_sIs[i].data()[0];
    else
    { 
      if (!_rank_deficient[i])
        _LA->factor_chol(_sIs[i]);
      else
        _LA->svd(_sIs[i], _usIs[i], _ssIs[i], _vsIs[i]);
    }

    // get Is
    const vector<SMOMENTUM>& Is = _Is[i];
 
    FILE_LOG(LOG_DYNAMICS) << "  *** Backward recursion processing link " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    I: " << I << endl;

    // don't update I for direct descendants of the base if the base is 
    // not floating
    if (!body->is_floating_base() && parent->is_base())
      continue;
 
    // compute a couple of necessary matrices
    transpose_solve_sIs(i, sprime, _sIss);
    SPARITH::mult(Is, _sIss, tmp);
    I.to_matrix(tmp2);
    MATRIXN::mult(tmp, tmp2, tmp3);
    SPATIAL_AB_INERTIA uI = I - SPATIAL_AB_INERTIA::from_matrix(tmp3, I.pose);

    // output the updates
    if (LOGGING(LOG_DYNAMICS) && _Is[i].size() > 0)
      FILE_LOG(LOG_DYNAMICS) << "  Is: " << _Is[i][0] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  s/(s'Is): " << _sIss << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  Is*s/(s'Is): " << std::endl << tmp;
    FILE_LOG(LOG_DYNAMICS) << "  Is*s/(s'Is)*I: " << std::endl << tmp3;
    FILE_LOG(LOG_DYNAMICS) << "  inertial update: " << uI << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  transformed I: " << POSE3::transform(_I[h].pose, uI) << std::endl;

    // update the parent inertia
    _I[h] += POSE3::transform(_I[h].pose, uI);
  }
}

/// Computes joint and spatial link accelerations 
void FSAB_ALGORITHM::calc_spatial_accelerations(shared_ptr<RC_ARTICULATED_BODY> body)
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  VECTORN result;
  vector<SVELOCITY> sprime, sdotprime;

  // get the links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // verify that base is index 0
  assert(base->get_index() == 0);

  // set spatial acceleration of base
  if (!body->is_floating_base())
  {
    base->set_accel(SACCEL::zero(base->get_computation_frame()));
    FILE_LOG(LOG_DYNAMICS) << "  base acceleration: (zero)" << endl;
  }
  else
  {
    SACCEL a0 = _I.front().inverse_mult(-_Z.front());
    base->set_accel(a0);
    FILE_LOG(LOG_DYNAMICS) << "  articulated base inertia: " << POSE3::transform(base->get_mixed_pose(), _I.front()) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  negated base Z: " << POSE3::transform(base->get_mixed_pose(), -_Z.front()) << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base acceleration: " << SPARITH::transform_accel(base->get_mixed_pose(), a0) << endl;
  }
  
  // *****************************************************************
  // compute joint accelerations
  // *****************************************************************

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // process all links for forward recursion
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();
    const unsigned i = link->get_index();

    // push all children of the link onto the queue
    push_children(link, link_queue);
      
    // get the parent link and the inner joint
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());

    // get the parent index
    const unsigned h = parent->get_index();
    
    // compute transformed parent link acceleration
    SACCEL ah = SPARITH::transform_accel(link->get_computation_frame(), parent->get_accel()); 

    // get the spatial axis and its derivative
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    const vector<SVELOCITY>& sdot = joint->get_spatial_axes_dot();

    // transform spatial axes
    POSE3::transform(link->get_computation_frame(), s, sprime);
    POSE3::transform(link->get_computation_frame(), sdot, sdotprime);

    // get the Is and qm subexpressions
    const VECTORN& mu = _mu[i];    
    const SACCEL& c = _c[i];

    // compute joint i acceleration
    SFORCE w = _I[i] * ah;
    SPARITH::transpose_mult(sprime, w, result);
    result.negate();
    result += mu;
    solve_sIs(i, result, joint->qdd);
    
    // compute link i spatial acceleration
    SACCEL ai = ah + c;
    if (!sprime.empty())
      ai += SACCEL(SPARITH::mult(sprime, joint->qdd));
    if (!sdotprime.empty())
      ai += SACCEL(SPARITH::mult(sdotprime, joint->qd));
    link->set_accel(ai);

    FILE_LOG(LOG_DYNAMICS) << endl << endl << "  *** Forward recursion processing link " << link << endl;  
    FILE_LOG(LOG_DYNAMICS) << "    a[h]: " << ah << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qm(subexp): " << mu << endl;
    FILE_LOG(LOG_DYNAMICS) << "    qdd: " << joint->qdd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    spatial acceleration: " << link->get_accel() << endl;
  }
  
  FILE_LOG(LOG_DYNAMICS) << "    joint accel: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << boost::shared_ptr<JOINT>(links[i]->get_inner_joint_explicit())->qdd << " ";
  FILE_LOG(LOG_DYNAMICS) << "joint q: ";
  for (unsigned i=1; i< links.size(); i++)
    FILE_LOG(LOG_DYNAMICS) << boost::shared_ptr<JOINT>(links[i]->get_inner_joint_explicit())->q << " ";
  FILE_LOG(LOG_DYNAMICS) << endl;
  FILE_LOG(LOG_DYNAMICS) << "calc_joint_accels() exited" << endl;

  FILE_LOG(LOG_DYNAMICS) << endl << "Outputting link data" << endl;
  for (unsigned i=0; i< links.size(); i++)
  {
    FILE_LOG(LOG_DYNAMICS) << "*** Link: " << links[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  v: " << links[i]->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  a: " << links[i]->get_accel() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  c: " << _c[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  Z: " << _Z[links[i]->get_index()] << endl;
    FILE_LOG(LOG_DYNAMICS) << "  I: " << endl << _I[links[i]->get_index()];
  }
}

/// Computes the joint accelerations (forward dynamics) for an articulated body
/**
 * Featherstone Algorithm taken from Mirtich's thesis (p. 113).  Mirtich's 
 * numbering is a little funny, so I decrement his joint indices by one, while 
 * leaving his link indices intact.
 */
void FSAB_ALGORITHM::calc_fwd_dyn()
{
  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() entered" << endl;

  // get the body and the reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  if (!body->_ijoints.empty())
    throw std::runtime_error("FSAB_ALGORITHM cannot process bodies with kinematic loops!");

  // get the links and joints for the body
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<boost::shared_ptr<JOINT> >& joints = body->get_explicit_joints();

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();
 
  // compute spatial coriolis vectors
  calc_spatial_coriolis_vectors(body);

  // compute spatial articulated body inertias
  calc_spatial_inertias(body);

  // compute spatial ZAs
  calc_spatial_zero_accelerations(body);
    
  // compute spatial accelerations
  calc_spatial_accelerations(body);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() exited" << endl;
}

/// Computes the joint accelerations (forward dynamics) when position variables have not been changed 
/**
 * Featherstone Algorithm taken from Mirtich's thesis (p. 113).  Mirtich's 
 * numbering is a little funny, so I decrement his joint indices by one, while 
 * leaving his link indices intact.
 */
void FSAB_ALGORITHM::calc_fwd_dyn_special()
{
  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() entered" << endl;

  // get the body and the reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  if (!body->_ijoints.empty())
    throw std::runtime_error("FSAB_ALGORITHM cannot process bodies with kinematic loops!");

  // get the links and joints for the body
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<boost::shared_ptr<JOINT> >& joints = body->get_explicit_joints();

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();
 
  // compute spatial coriolis vectors
  calc_spatial_coriolis_vectors(body);

  // compute spatial ZAs
  calc_spatial_zero_accelerations(body);
    
  // compute spatial accelerations
  calc_spatial_accelerations(body);

  FILE_LOG(LOG_DYNAMICS) << "FSABAlgorith::calc_fwd_dyn() exited" << endl;
}

/// Signum function with NEAR_ZERO cutoff
REAL FSAB_ALGORITHM::sgn(REAL x)
{
  const REAL NEAR_ZERO = std::sqrt(std::numeric_limits<REAL>::epsilon());
  if (x > NEAR_ZERO)
    return (REAL) 1.0;
  else if (x < -NEAR_ZERO)
    return (REAL) -1.0;
  else
    return (REAL) 0.0;
}

/// Pushes all children of the given link onto the given queue
void FSAB_ALGORITHM::push_children(shared_ptr<RIGIDBODY> link, queue<shared_ptr<RIGIDBODY> >& q)
{
  const set<boost::shared_ptr<JOINT> >& ojs = link->get_outer_joints();
  BOOST_FOREACH(boost::shared_ptr<JOINT> j, ojs)
  {
    shared_ptr<RIGIDBODY> child(j->get_outboard_link());
    q.push(child);
  }
}

/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
/**
 * \pre spatial inertias already computed for the body's current configuration 
 */
void FSAB_ALGORITHM::apply_impulse(const SMOMENTUM& w, shared_ptr<RIGIDBODY> link)
{
  static MATRIXN tmp;
  static VECTORN tmp2, workv;
  vector<SVELOCITY> sprime;

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_impulse() entered" << endl;
  FILE_LOG(LOG_DYNAMICS) << " -- applying impulse " << w << endl;

  // get the computation reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);   
  if (!body->_ijoints.empty())
    throw std::runtime_error("FSAB_ALGORITHM cannot process bodies with kinematic loops!");

  // initialize spatial zero velocity deltas to zeros
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const unsigned NUM_LINKS = links.size();
  _Y.resize(NUM_LINKS);
  for (unsigned j=0; j< NUM_LINKS; j++)
  {
    const unsigned i = links[j]->get_index();
    _Y[i].set_zero();
    _Y[i].pose = links[j]->get_computation_frame();
  } 

  // **********************************************************************
  // NOTE: uses articulated body inertias and spatial axes already computed 
  // **********************************************************************


  // get the index for this link
  unsigned i = link->get_index();
  
  // transform the impulse
  _Y[i] = -POSE3::transform(link->get_computation_frame(), w);
  
  FILE_LOG(LOG_DYNAMICS) << "  -- impulse applied to link " << link << " = " << w << endl;
  FILE_LOG(LOG_DYNAMICS) << "  -- recursing backward" << endl;

  // do not recurse backward if this is the base
  if (!link->is_base())
  {
    // recurse backward
    while (!link->is_base())
    {
      // get the parent link
      shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    
      // get the parent index
      const unsigned h = parent->get_index();

      // get spatial axes of the inner joint for link i
      boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
      const vector<SVELOCITY>& s = joint->get_spatial_axes();
      POSE3::transform(link->get_computation_frame(), s, sprime);
    
      // get Is for link i
      const vector<SMOMENTUM>& Is = _Is[i];
    
      // compute Is * inv(sIs) * s'
      transpose_solve_sIs(i, sprime, _sIss);
      SPARITH::mult(Is, _sIss, tmp); 
      tmp.mult(_Y[i], workv);

      // compute impulse for h in i's frame (or global frame)
      SMOMENTUM Yi = _Y[i] - SMOMENTUM::from_vector(workv,  _Y[i].pose);

      // transform the spatial impulse, if necessary
      _Y[h] = POSE3::transform(_Y[h].pose, Yi);
   
      FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- this transformed impulse is: " << _Y[i] << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- parent is link: " << h << endl;
      FILE_LOG(LOG_DYNAMICS) << "    -- transformed spatial impulse for parent: " << _Y[h] << endl; 
    
      // update the link to be the parent
      link = parent;
      i = h;
    }
  }

  // ************************************************************
  // determine the new joint and link velocities
  // ************************************************************
  
  // setup a queue for recursing forward through the links
  queue<shared_ptr<RIGIDBODY> > link_queue;

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // add all children of the base to the link queue
  push_children(base, link_queue);

  // setup a vector of link velocity updates
  _dv.resize(NUM_LINKS);
  
  // if floating base, apply spatial impulse
  if (body->is_floating_base())
  {
    // determine the change in velocity
    _dv.front() = _I.front().inverse_mult(-_Y.front());

    FILE_LOG(LOG_DYNAMICS) << "base is floating..." << endl;
    FILE_LOG(LOG_DYNAMICS) << "  base transform: " << endl << base->get_pose();
    FILE_LOG(LOG_DYNAMICS) << "  current base velocity: " << base->get_velocity() << endl;
    FILE_LOG(LOG_DYNAMICS) << "  impulse on the base: " << _Y.front() << endl;

    // update the base velocity
    base->set_velocity(base->get_velocity() + _dv.front());
    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << base->get_velocity() << endl;
  }
  else 
    _dv.front().set_zero();
  
  // update link and joint velocities
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the parent link
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    
    // get the index of the link and its parent
    unsigned i = link->get_index(); 
    unsigned h = parent->get_index(); 
    
    // get the inboard joint
    boost::shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    
    // get spatial axes of the inner joint for link i
    const vector<SVELOCITY>& s = joint->get_spatial_axes();
    POSE3::transform(link->get_computation_frame(), s, sprime);
    
    // get articulated body inertia for link i
    const SPATIAL_AB_INERTIA& I = _I[i];    
    
    FILE_LOG(LOG_DYNAMICS) << "  -- processing link: " << link << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- parent is link " << parent << endl;
    
    // determine the joint and link velocity updates
    SVELOCITY dvh = POSE3::transform(_dv[i].pose, _dv[h]);
    SMOMENTUM f = (I * dvh) + _Y[i];
    SPARITH::transpose_mult(sprime, (I * dvh) + _Y[i], tmp2);
    solve_sIs(i, tmp2, _qd_delta).negate();
    _dv[i] = dvh + SPARITH::mult(s, _qd_delta);

    FILE_LOG(LOG_DYNAMICS) << "    -- I * dv[parent]: " << _dv[h] << endl;
    
    // update the joint velocity
    joint->qd += _qd_delta;

    // update the link velocity
    link->set_velocity(link->get_velocity() + _dv[i]);
 
    FILE_LOG(LOG_DYNAMICS) << "    -- cumulative transformed impulse on this link: " << _Y[i] << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta qd: " << _qd_delta << "  qd: " << joint->qd << endl;
    FILE_LOG(LOG_DYNAMICS) << "    -- delta v: " << _dv[i] << endl;
   
    // place all children on the link queue
    push_children(link, link_queue);
  }

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "FSAB_ALGORITHM::apply_impulse() exited" << endl;
}

/// Solves a system for sIs*x = m' using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MATRIXN& FSAB_ALGORITHM::transpose_solve_sIs(unsigned i, const vector<SVELOCITY>& m, MATRIXN& result) const
{
  // transpose m
  SPARITH::transpose_to_matrix(m, result);   

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);
 
  return result;
}

/// Solves a system for sIs*x = m using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
MATRIXN& FSAB_ALGORITHM::solve_sIs(unsigned i, const MATRIXN& m, MATRIXN& result) const
{
  result = m;

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);
    
  return result;
}

/// Solves a system for sIs*x = v using a factorization (if sIs is nonsingular) or the pseudo-inverse of sIs otherwise
VECTORN& FSAB_ALGORITHM::solve_sIs(unsigned i, const VECTORN& v, VECTORN& result) const
{
  result = v;

  // look for simplest case
  if (_sIs[i].rows() == 1)
  {
    result *= _sIs[i].data()[0];
    return result;
  } 

  // determine whether we are dealing with a rank-deficient sIs
  if (_rank_deficient[i])
    _LA->solve_LS_fast(_usIs[i], _ssIs[i], _vsIs[i], result);
  else
    _LA->solve_chol_fast(_sIs[i], result);

  return result;
}

