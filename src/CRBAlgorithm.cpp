/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using std::list;
using std::map;
using std::vector;
using std::queue;
using std::endl;
using boost::shared_array;
using boost::shared_ptr;

CRB_ALGORITHM::CRB_ALGORITHM()
{
}

/// Computes the parent array for sparse Cholesky factorization
void CRB_ALGORITHM::setup_parent_array()
{
  // get the number of generalized coordinates
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  const unsigned N = body->num_generalized_coordinates(DYNAMIC_BODY::eSpatial);

  // get explicit joints
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();

  // determine parent array (lambda)
  _lambda.resize(N);

  // set all values of lambda to inf initially
  for (unsigned i=0; i< N; i++)
    _lambda[i] = std::numeric_limits<unsigned>::max();
  
  // loop over all explicit joints
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get the index of this joint
    unsigned idx = ejoints[i]->get_coord_index();

    // get the parent joint and its index
    shared_ptr<RIGIDBODY> inboard = ejoints[i]->get_inboard_link();
    shared_ptr<JOINT> parent = inboard->get_inner_joint_explicit();
    if (!parent)
      continue;
    unsigned pidx = parent->get_coord_index();

    // now set the elements of lambda
    for (unsigned j=0; j< ejoints[i]->num_dof(); j++)
    {
      _lambda[idx] = pidx;
      pidx = idx;
      idx++;
    }
  }
}

/// Factorizes (Cholesky) the generalized inertia matrix, exploiting sparsity
bool CRB_ALGORITHM::factorize_cholesky(MATRIXN& M)
{
  // check whether the parent array has been setup
  if (_lambda.size() == 0)
    setup_parent_array(); 

  // get the number of degrees of freedom
  const unsigned n = M.rows();

  // loop
  for (unsigned kk=n; kk> 0; kk--)
  {
    unsigned k = kk - 1;

    // get a column iterator
    COLUMN_ITERATOR kiter = M.column_iterator_begin();

    // set it to row k
    kiter += n*k;

    if (kiter[k] < (double) 0.0)
    {
      assert(false);
      return false;
    }
    kiter[k] = std::sqrt(kiter[k]);
    unsigned i = _lambda[k];
    while (i != std::numeric_limits<unsigned>::max())
    {
      kiter[i] /= kiter[k];
      i = _lambda[i];
    }
    i = _lambda[k];
    while (i != std::numeric_limits<unsigned>::max())
    {
      unsigned j=i;

      // get a column iterator
      COLUMN_ITERATOR iiter = M.column_iterator_begin()+(i*n);

      // iterate
      while (j != std::numeric_limits<unsigned>::max())
      {
        iiter[j] -= kiter[i]*kiter[j];
        j = _lambda[j];
      }
      i = _lambda[i];
    }
  }

  return true;
}

// Transforms (as necessary) and multiplies
void CRB_ALGORITHM::transform_and_mult(shared_ptr<const POSE3> target, const SPATIAL_RB_INERTIA& I, const vector<SVELOCITY>& s, vector<SMOMENTUM>& Is)
{
  // transform s
  POSE3::transform(I.pose, s, _sprime);

  // do the multiplication
  SPARITH::mult(I, _sprime, _Isprime);

  // do the transformation
  POSE3::transform(target, _Isprime, Is);
}

/// Gets the frame in which kinematics and dynamics computations occur
/**
 * Only valid for bodies with floating bases.
 */
shared_ptr<const POSE3> CRB_ALGORITHM::get_computation_frame(shared_ptr<RC_ARTICULATED_BODY> body)
{
  assert(body->is_floating_base());

  // get the base 
  shared_ptr<RIGIDBODY> base = body->get_base_link();

  switch (body->get_computation_frame_type())
  {
    case eLink:
    case eJoint:
      return base->get_pose();

    case eLinkInertia:
      return base->get_inertial_pose();

    case eLinkCOM:
      return base->get_gc_pose();

    case eGlobal:
      return shared_ptr<const POSE3>();

    default:
      assert(false);
  }

  return shared_ptr<const POSE3>();
}

/// Calculates the generalized inertia of this body
/**
 * Specialized function for use with the CRB algorithm
 */
void CRB_ALGORITHM::calc_generalized_inertia(shared_ptr<RC_ARTICULATED_BODY> body)
{
  const unsigned SPATIAL_DIM = 6;

  // get the appropriate M
  MATRIXN& M = this->_M;

  // get the set of links
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // get explicit joints
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();

  // compute the joint space inertia
  calc_joint_space_inertia(body, _H, _Ic);

  // get the number of base degrees-of-freedom
  const unsigned N_BASE_DOF = (body->is_floating_base()) ? 6 : 0;

  // resize M
  M.resize(N_BASE_DOF + body->num_joint_dof_explicit(), N_BASE_DOF + body->num_joint_dof_explicit());

  // set appropriate part of H
  M.set_sub_mat(0, 0, _H);

  // see whether we are done
  if (!body->is_floating_base())
    return;

  // ************************************************************************
  // floating base
  // ************************************************************************

  // setup the indices for the base
  const unsigned BASE_START = body->num_joint_dof_explicit();

  // get components of M
  SHAREDMATRIXN Ic0 = M.block(BASE_START, M.rows(), BASE_START, M.columns());
  SHAREDMATRIXN KS = M.block(0, BASE_START, BASE_START, M.columns()); 
  SHAREDMATRIXN K = M.block(BASE_START, M.rows(), 0, BASE_START); 

  // get composite inertia in desired frame
  shared_ptr<const POSE3> P = get_computation_frame(body);
  POSE3::transform(P, _Ic.front()).to_PD_matrix(Ic0);

  FILE_LOG(LOG_DYNAMICS) << "Ic0: " << std::endl << Ic0;

  // compute K
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get the spatial axes for the joint
    shared_ptr<JOINT> joint = ejoints[i];
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes();
    if (joint->num_dof() == 0)
      continue;

    // get the index for this joint
    unsigned jidx = joint->get_coord_index();

    // get the outboard link and link index
    shared_ptr<RIGIDBODY> outboard = joint->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // transform and multiply
    transform_and_mult(P, _Ic[oidx], s, _Is);

    // compute the requisite columns of K
    SHAREDMATRIXN Kb = K.block(0, SPATIAL_DIM, jidx, jidx+joint->num_dof()); 
    SHAREDMATRIXN KSb = KS.block(jidx, jidx+joint->num_dof(), 0, SPATIAL_DIM); 
    SPARITH::to_matrix(_Is, Kb); 
    MATRIXN::transpose(Kb, KSb); 
  }

  FILE_LOG(LOG_DYNAMICS) << "[H K'; K Ic0] (permuted): " << std::endl << M;
}

/// Computes *just* the joint space inertia matrix
void CRB_ALGORITHM::calc_joint_space_inertia(shared_ptr<RC_ARTICULATED_BODY> body, MATRIXN& H, vector<SPATIAL_RB_INERTIA>& Ic)
{
  queue<shared_ptr<RIGIDBODY> > link_queue;
  const unsigned SPATIAL_DIM = 6;

  // get the reference frame
  ReferenceFrameType rftype = body->get_computation_frame_type();

  // get the sets of links and joints
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();
  const vector<shared_ptr<JOINT> >& joints = body->get_joints();

  // set the composite inertias to the isolated inertias initially 
  Ic.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    Ic[links[i]->get_index()] = links[i]->get_inertia();

    // check for degenerate inertia
    #ifndef NDEBUG
    SPATIAL_RB_INERTIA& J = Ic[links[i]->get_index()];
    if (links[i]->is_base() && body->is_floating_base() && 
        (J.m <= 0.0 || J.J.norm_inf() <= 0.0))
      throw std::runtime_error("Attempted to compute dynamics given degenerate inertia for a floating base body");
    #endif

    if (LOGGING(LOG_DYNAMICS))
    {
      MATRIXN X;
      Ic[links[i]->get_index()].to_matrix(X);
      FILE_LOG(LOG_DYNAMICS) << "isolated inertia for link " << links[i]->body_id << ": " << endl << X;
    }
  }

  // ************************************************************************
  // first, determine the supports for the joints and the number of joint DOF
  // ************************************************************************
 
  // create and initialize the supports array
  _supports.resize(joints.size());
  for (unsigned i=0; i< joints.size(); i++)
  {
    _supports[i].resize(links.size());
    for (unsigned j=0; j< links.size(); j++)
      _supports[i][j] = false;
  }

  // add all leaf links to the link queue
  for (unsigned i=1; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);

  // process until the queue is empty
  while (!link_queue.empty())
  {
    // get the element out of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the explicit inner joint for this link
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    unsigned jidx = joint->get_index();
    assert(joint);

    // add this link to the support for the joint
    _supports[jidx][link->get_index()] = true;

    // add all supports from the outer joints of this link
    list<shared_ptr<RIGIDBODY> > child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(shared_ptr<RIGIDBODY> child, child_links)
    {
      // don't process children with lower link indices (loops)
      if (child->get_index() < link->get_index())
        continue;

      // get the inner explicit joint
      shared_ptr<JOINT> child_joint = child->get_inner_joint_explicit();
      unsigned jiidx = child_joint->get_index();

      // setup the supports
      for (unsigned i=0; i< links.size(); i++)
        if (_supports[jiidx][i])
          _supports[jidx][i] = true;
    }  

    // add the parent of this link to the queue, if it is not the base and
    // if the parent's link index is lower
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (!parent->is_base() && parent->get_index() < link->get_index())
      link_queue.push(parent);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    FILE_LOG(LOG_DYNAMICS) << "supports: " << endl;
    for (unsigned i=0; i< _supports.size(); i++)
    {
      std::ostringstream str;
      for (unsigned j=0; j< _supports[i].size(); j++)
        str << (bool) _supports[i][j] << " ";
      FILE_LOG(LOG_DYNAMICS) << i << ": " << str.str() << endl;
    }
  }

  // resize H 
  H.set_zero(body->num_joint_dof_explicit(), body->num_joint_dof_explicit());

  // ************************************************************************
  // compute spatial composite inertias 
  // ************************************************************************

  // now determine the composite inertias
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // put all leaf links into a queue
  assert(link_queue.empty());
  for (unsigned i=1; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);

  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();

    // get the index for this link
    unsigned i = link->get_index();
    
    // see whether this link has already been processed
    if (body->_processed[i])
      continue;

    // verify that all children have been processed
    if (!body->all_children_processed(link))
      continue;    

    // process the parent link, if possible
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    if (parent && i > parent->get_index())
    {
      // put the parent on the queue
      link_queue.push(parent);
    
      // get the parent index
      unsigned h = parent->get_index();
    
      // add this inertia to its parent
      Ic[h] += POSE3::transform(Ic[h].pose, Ic[i]); 

      if (LOGGING(LOG_DYNAMICS))
      {
        MATRIXN X;
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (child) link " << link << ": " << std::endl << Ic[i].to_matrix(X);
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (child) link " << link << ": " << std::endl << Ic[i];
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (parent) link " << parent << ": " << std::endl << Ic[h].to_matrix(X);
        FILE_LOG(LOG_DYNAMICS) << "  composite inertia for (parent) link " << parent << ": " << std::endl << Ic[h];
      }
    }

    // indicate that the link has been processed
    body->_processed[i] = true;
  }

  // ************************************************************************
  // compute H
  // ************************************************************************

  // compute the forces
  _momenta.resize(links.size());
  for (unsigned i=0; i < ejoints.size(); i++)
  {
    shared_ptr<RIGIDBODY> outboard = ejoints[i]->get_outboard_link(); 
    unsigned oidx = outboard->get_index();
    const std::vector<SVELOCITY>& s = ejoints[i]->get_spatial_axes();
    POSE3::transform(Ic[oidx].pose, s, _sprime);
    SPARITH::mult(Ic[oidx], _sprime, _momenta[oidx]);
    if (_sprime.size() > 0)
    {
      for (unsigned j=0; j< ejoints[i]->num_dof(); j++)
      {
        FILE_LOG(LOG_DYNAMICS) << "s[ " << j << "]: " << _sprime[j] << std::endl;
        FILE_LOG(LOG_DYNAMICS) << "Is[" << j << "]: " << _momenta[oidx][j] << std::endl;
      }
    }
  } 

  // setup H
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get the number of degrees of freedom for joint i
    const unsigned NiDOF = ejoints[i]->num_dof();

    // get the starting coordinate index for this joint
    unsigned iidx = ejoints[i]->get_coord_index();

    // get the outboard link for joint i
    shared_ptr<RIGIDBODY> obi = ejoints[i]->get_outboard_link();
    unsigned oiidx = obi->get_index();

    // get the spatial axes for jointi
    const std::vector<SVELOCITY>& s = ejoints[i]->get_spatial_axes();

    // get the appropriate submatrix of H
    SHAREDMATRIXN subi = H.block(iidx, iidx+NiDOF, iidx, iidx+NiDOF); 

    // compute the H term for i,i
    transform_and_transpose_mult(s, _momenta[oiidx], subi);

    // determine what will be the new value for m
    for (unsigned j=i+1; j< ejoints.size(); j++)
    {
      // get the number of degrees of freedom for joint j 
      const unsigned NjDOF = ejoints[j]->num_dof();

      // get the outboard link for joint j
      shared_ptr<RIGIDBODY> obj = ejoints[j]->get_outboard_link();
      unsigned ojidx = obj->get_index();

      // if link j is not supported by joint i, contribution to H is zero
      if (!_supports[ejoints[i]->get_index()][ojidx])
        continue;

      // get the starting coordinate index for joint j
      unsigned jidx = ejoints[j]->get_coord_index();

      // get the appropriate submatrices of H
      SHAREDMATRIXN subj = H.block(iidx, iidx+NiDOF, jidx, jidx+NjDOF); 
      SHAREDMATRIXN subjT = H.block(jidx, jidx+NjDOF, iidx, iidx+NiDOF); 

      // compute the appropriate submatrix of H
      transform_and_transpose_mult(s, _momenta[ojidx], subj);

      // set the transposed part
      MATRIXN::transpose(subj, subjT);
    }
  }

  FILE_LOG(LOG_DYNAMICS) << "joint space inertia: " << endl << H;
}

/// Calculates the generalized inertia matrix for the given representation
/**
 * Generic method provided for use with generalized coordinates.
 */
void CRB_ALGORITHM::calc_generalized_inertia(SHAREDMATRIXN& M)
{
  // do the precalculation
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  precalc(body);

  // get the set of links
  ReferenceFrameType rftype = body->get_computation_frame_type();
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();

// TODO: this should be able to be removed; this function should already be
//       computed using precalc(.)
  // get the joint space inertia and composite inertias
//  calc_joint_space_inertia(body, _H, _Ic);

  // get the number of base degrees-of-freedom
  const unsigned N_BASE_DOF = (body->is_floating_base()) ? 6 : 0;
  const unsigned SPATIAL_DIM = 6;

  // resize M and set H
  M.resize(N_BASE_DOF + body->num_joint_dof_explicit(), N_BASE_DOF + body->num_joint_dof_explicit());
  M.set_sub_mat(0, 0, _H);

  // look for simplest case
  if (!body->is_floating_base())
    return;

  // get the coordinates at which the base start 
  const unsigned BASE_START = M.rows() - SPATIAL_DIM;

  // ************************************************************************
  // floating base
  // ************************************************************************

  // get the number of explicit joint degrees-of-freedom
  const unsigned NjDOF = body->num_joint_dof_explicit();

  // get components of M
  SHAREDMATRIXN Ic0 = M.block(BASE_START, M.rows(), BASE_START, M.columns());
  SHAREDMATRIXN KS = M.block(0, BASE_START, BASE_START, M.columns()); 
  SHAREDMATRIXN K = M.block(BASE_START, M.rows(), 0, BASE_START); 

  // get composite inertia in desired frame
  shared_ptr<RIGIDBODY> base = body->get_base_link();
  shared_ptr<const POSE3> P = base->get_gc_pose();
  POSE3::transform(P, _Ic[base->get_index()]).to_PD_matrix(Ic0); 

  // compute K
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get the spatial axes for the joint
    shared_ptr<JOINT> joint = ejoints[i];
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes();
    if (joint->num_dof() == 0)
      continue;

    // get the index for this joint
    unsigned jidx = joint->get_coord_index();

    // get the outboard link and link index
    shared_ptr<RIGIDBODY> outboard = joint->get_outboard_link();
    unsigned oidx = outboard->get_index();

    // transform and multiply
    transform_and_mult(P, _Ic[oidx], s, _Is);

    // compute the requisite columns of K
    SHAREDMATRIXN Kb = K.block(0, SPATIAL_DIM, jidx, jidx+joint->num_dof()); 
    SHAREDMATRIXN KSb = KS.block(jidx, jidx+joint->num_dof(), 0, SPATIAL_DIM); 
    SPARITH::to_matrix(_Is, Kb); 
    MATRIXN::transpose(Kb, KSb);
  }

  FILE_LOG(LOG_DYNAMICS) << "[H K'; K Ic0] (permuted): " << std::endl << M;
}

/// Performs necessary pre-computations for computing accelerations or applying impulses
void CRB_ALGORITHM::precalc(shared_ptr<RC_ARTICULATED_BODY> body)
{
  // tolerance for not recomputing/refactorizing inertia matrix
  const double REFACTOR_TOL = 1e-4;

  // get the links and joints for the body
  const vector<shared_ptr<JOINT> >& joints = body->get_explicit_joints();

  // get the generalized coordinates
  static VECTORN gc, tmpv;
  body->get_generalized_coordinates(DYNAMIC_BODY::eEuler, gc);
  if (_gc_last.size() == 0 || ((tmpv = gc) -= _gc_last).norm_inf() > REFACTOR_TOL)
  {
    // compute spatial isolated inertias and generalized inertia matrix
    // do the calculations
    calc_generalized_inertia(body);

    // attempt to do a Cholesky factorization of M
    MATRIXN& fM = this->_fM;
    MATRIXN& M = this->_M;
    if ((_rank_deficient = !_LA->factor_chol(fM = M)))
//  if ((_rank_deficient = !factorize_cholesky(fM = M)))
    {
      std::cerr << "CRBAlgorithm::precalc() warning- Cholesky factorization of generalized inertia matrix failed" << std::endl;
      fM = M;
      _LA->svd(fM, _uM, _sM, _vM);
    }

    _gc_last = gc;
  }
}

/// Executes the composite rigid-body method
void CRB_ALGORITHM::calc_fwd_dyn()
{
  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // do necessary pre-calculations
  precalc(body);

  // execute the appropriate algorithm
  if (body->is_floating_base())
    calc_fwd_dyn_floating_base(body);
  else
    calc_fwd_dyn_fixed_base(body);
   
  // update the link accelerations
  update_link_accelerations(body);
}

/// Executes the composite rigid-body method without computing and factorizing inertia matrix
/**
 * This method is useful when the inertia matrix has already been computed-
 * considerable computation will then be avoided.
 */
void CRB_ALGORITHM::calc_fwd_dyn_special()
{
  // get the body and the reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // execute the appropriate algorithm
  if (body->is_floating_base())
    calc_fwd_dyn_floating_base(body);
  else
    calc_fwd_dyn_fixed_base(body);
   
  // update the link accelerations
  update_link_accelerations(body);
}

/// Solves for acceleration using the body inertia matrix
VECTORN& CRB_ALGORITHM::M_solve(VECTORN& xb) 
{
  // do necessary pre-calculations
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  precalc(body);

  // setup xb
  SHAREDVECTORN xb_shared = xb.segment(0, xb.rows());

  M_solve_noprecalc(xb_shared); 
  return xb;
}

/// Solves for acceleration using the body inertia matrix
SHAREDVECTORN& CRB_ALGORITHM::M_solve(SHAREDVECTORN& xb) 
{
  // do necessary pre-calculations
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  precalc(body);

  M_solve_noprecalc(xb); 
  return xb;
}

/// Solves for acceleration using the body inertia matrix
MATRIXN& CRB_ALGORITHM::M_solve(MATRIXN& XB)
{
  // do necessary pre-calculations
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  precalc(body);

  // setup XB
  SHAREDMATRIXN XB_shared = XB.block(0, XB.rows(), 0, XB.columns());

  M_solve_noprecalc(XB_shared); 
  return XB;
}

/// Solves for acceleration using the body inertia matrix
SHAREDMATRIXN& CRB_ALGORITHM::M_solve(SHAREDMATRIXN& XB)
{
  // do necessary pre-calculations
  shared_ptr<RC_ARTICULATED_BODY> body(_body);
  precalc(body);

  M_solve_noprecalc(XB); 
  return XB;
}

/// Solves for acceleration using the body inertia matrix
VECTORN& CRB_ALGORITHM::M_solve_noprecalc(VECTORN& xb)
{
  // setup xb
  SHAREDVECTORN xb_shared = xb.segment(0, xb.rows());

  // solve
  M_solve_noprecalc(xb_shared);

  return xb;
}

/// Solves for acceleration using the body inertia matrix
SHAREDVECTORN& CRB_ALGORITHM::M_solve_noprecalc(SHAREDVECTORN& xb)
{
  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
    _LA->solve_LS_fast(_uM, _sM, _vM, xb);
  else
    _LA->solve_chol_fast(_fM, xb);

  return xb;
}

/// Solves for acceleration using the body inertia matrix
MATRIXN& CRB_ALGORITHM::M_solve_noprecalc(MATRIXN& XB)
{
  // get a shared matrix
  SHAREDMATRIXN XB_shared = XB.block(0, XB.rows(), 0, XB.columns());

  // solve
  M_solve_noprecalc(XB_shared);

  return XB;
}

/// Solves for acceleration using the body inertia matrix
SHAREDMATRIXN& CRB_ALGORITHM::M_solve_noprecalc(SHAREDMATRIXN& XB)
{
  // determine whether the matrix is rank-deficient
  if (this->_rank_deficient)
    _LA->solve_LS_fast(_uM, _sM, _vM, XB);
  else
    _LA->solve_chol_fast(_fM, XB);

  return XB;
}

/// Executes the composite rigid-body method on an articulated body with a fixed base
void CRB_ALGORITHM::calc_fwd_dyn_fixed_base(shared_ptr<RC_ARTICULATED_BODY> body)
{
  // get the set of links and joints for the articulated body
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();

  // ***********************************************************************
  // first, calculate C
  // ***********************************************************************

  // call inverse dynamics to calculate C
  SFORCE f0;
  calc_generalized_forces(f0, _C);
  
  // get the number of degrees-of-freedom
  unsigned nDOF = _C.size();

  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************
  _Q.resize(nDOF);
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    unsigned j = ejoints[i]->get_coord_index();
    _Q.set_sub_vec(j, ejoints[i]->force);
  }

  FILE_LOG(LOG_DYNAMICS) << "H: " << std::endl << _M;
  FILE_LOG(LOG_DYNAMICS) << "C: " << _C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "Q: " << _Q << std::endl;

  // subtract C from Q
  _Q -= _C;

  // get the pointer to the joint-space acceleration vector
  VECTORN& qdd = this->_qdd;

  // compute joint accelerations
  M_solve_noprecalc(qdd = _Q);

  FILE_LOG(LOG_DYNAMICS) << "qdd: " << qdd << std::endl;

  // set qdd
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    unsigned j = ejoints[i]->get_coord_index();
    qdd.get_sub_vec(j, j+ejoints[i]->num_dof(), ejoints[i]->qdd);
  }

  // set spatial acceleration of the base
  this->_a0.set_zero(); 
  SACCEL abase = links.front()->get_accel();
  links.front()->set_accel(SACCEL::zero(links.front()->get_computation_frame()));
}

/// Executes the composite rigid-body method on an articulated body with a floating base
/**
 * This algorithm is taken from [Featherstone, 1987], p. 123.  This is only
 * calculated in the global frame.
 */
void CRB_ALGORITHM::calc_fwd_dyn_floating_base(shared_ptr<RC_ARTICULATED_BODY> body)
{
  const boost::shared_ptr<const POSE3> GLOBAL;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() entered" << std::endl;

  // get the set of links and explicit joints
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();

  // calculate C
  SFORCE f0;
  calc_generalized_forces(f0, _C);

  // compute generalized forces in proper frame
  SFORCE f0x = POSE3::transform(get_computation_frame(body), f0);

  // get the number of degrees-of-freedom
  unsigned nDOF = _C.size();
  
  // ***********************************************************************
  // setup vector Q of actuator forces (note: this differs from 
  // [Featherstone, 87] p. 119, and may not be correct..  be prepared to 
  // change this
  // ***********************************************************************

  _Q.resize(nDOF);
  for (unsigned i=0; i< ejoints.size(); i++)
  {
     unsigned j = ejoints[i]->get_coord_index();
    _Q.set_sub_vec(j, ejoints[i]->force);
  }

  if (LOGGING(LOG_DYNAMICS))
  {
    POSE3 base_pose = *links[0]->get_pose();
    base_pose.update_relative_pose(GLOBAL);
    FILE_LOG(LOG_DYNAMICS) << "base pose: " << base_pose << std::endl;
  }
  FILE_LOG(LOG_DYNAMICS) << "Q: " << _Q << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "C: " << _C << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "M: " << std::endl << this->_M;

  // setup the simulataneous equations to solve: [Featherstone, 1987], eq. 7.24
  SPARITH::concat(_Q -= _C, -f0x, _b);
  FILE_LOG(LOG_DYNAMICS) << "b: " << _b << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "link + external forces on base: " << f0x << std::endl;

  // solve for accelerations
  M_solve_noprecalc(_augV = _b); 
  FILE_LOG(LOG_DYNAMICS) << "b: " << _b << std::endl;

  // get pointers to a0 and qdd vectors
  SACCEL& a0 = this->_a0;
  VECTORN& qdd = this->_qdd;

  // swap components of a0
  const unsigned SPATIAL_DIM = 6;
  const unsigned BASE_START = _augV.size() - SPATIAL_DIM;
  std::swap(_augV[BASE_START+0], _augV[BASE_START+3]);
  std::swap(_augV[BASE_START+1], _augV[BASE_START+4]);
  std::swap(_augV[BASE_START+2], _augV[BASE_START+5]);

  // get out a0, qdd
  qdd = _augV.segment(0, BASE_START);
  a0 = _augV.segment(BASE_START,_augV.size());
  a0.pose = f0x.pose;

  // set the base acceleration
  links.front()->set_accel(a0);

  FILE_LOG(LOG_DYNAMICS) << "base spatial acceleration: " << a0 << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "joint accelerations:";

  // set qdd
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    unsigned j = ejoints[i]->get_coord_index();
    qdd.get_sub_vec(j, j+ejoints[i]->num_dof(), ejoints[i]->qdd);
  }

  FILE_LOG(LOG_DYNAMICS) << qdd << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_fwd_dyn_floating_base() exited" << std::endl;
}

/// Computes the vector "C", used for forward dynamics computation, using the recursive Newton-Euler algorithm
/**
 * \return the spatial vector of forces on the base, which can be ignored for forward dynamics for fixed bases
 */
void CRB_ALGORITHM::calc_generalized_forces(SFORCE& f0, VECTORN& C)
{
  const unsigned SPATIAL_DIM = 6;
  queue<shared_ptr<RIGIDBODY> > link_queue;
  SFORCE w;

  // get the body and the reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // get the set of links and joints
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();
  if (links.empty())
  {
    C.resize(0);
    f0.set_zero();
    return;
  }

  // **************************************************************************
  // first, compute forward dynamics using RNE algorithm; we copy the algorithm
  // here, because we need some data out of it
  // **************************************************************************

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() entered" << std::endl;

  // ** STEP 1: compute accelerations

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // setup the map of link accelerations
  _a.resize(links.size());

  // setup the acceleration for the base
  _a[base->get_index()].set_zero(base->get_velocity().pose);
  
  // add all child links of the base to the processing queue
  list<shared_ptr<RIGIDBODY> > children;
  base->get_child_links(std::back_inserter(children));
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, children)
    link_queue.push(rb);
    
  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
  body->_processed[base->get_index()] = true;

  // process all links
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();  
    unsigned i = link->get_index();
    body->_processed[i] = true;

    // push all children of the link onto the queue, unless they were already
    // processed  
    list<shared_ptr<RIGIDBODY> > child_links;
    link->get_child_links(std::back_inserter(child_links)); 
    BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, child_links)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the link's parent
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    unsigned h = parent->get_index();

    // get the joint for this link
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());

    // get the spatial link velocity
    const SVELOCITY& vx = link->get_velocity(); 

    // get spatial axes and derivative for this link's inner joint
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes();
    const std::vector<SVELOCITY>& sdot = joint->get_spatial_axes_dot();

    // get the current joint velocity
    const VECTORN& qd = joint->qd;

    // **** compute acceleration

    // add this link's contribution
    POSE3::transform(vx.pose, s, _sprime);
    if (_sprime.empty())
      _a[i].set_zero(vx.pose);
    else
    {
      SVELOCITY sqd = SPARITH::mult(_sprime, qd);
      _a[i] = vx.cross(sqd);
    }
    if (!sdot.empty())
      _a[i] += SACCEL(SPARITH::mult(POSE3::transform(_a[i].pose, sdot, _sprime), qd)); 

    // now add parent's contribution
    _a[i] += SPARITH::transform_accel(_a[i].pose, _a[h]);
    FILE_LOG(LOG_DYNAMICS) << " computing link velocity / acceleration; processing link " << link->body_id << std::endl;
    if (s.size() > 0)
      FILE_LOG(LOG_DYNAMICS) << "  spatial joint velocity: " << (SPARITH::mult(s,qd)) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link velocity: " << link->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  link accel: " << _a[i] << std::endl;
  }
  
  // ** STEP 2: compute link forces -- backward recursion
  // use a map to determine which links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // setup a map of link forces, all set to zero initially
  _w.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    _w[i].set_zero();
    _w[i].pose = links[i]->get_computation_frame();
  }

  // add all leaf links to the queue
  for (unsigned i=0; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // if this link has already been processed, do not process it again
    if (body->_processed[i])
      continue;

    // verify all children have been processed
    if (!body->all_children_processed(link))
      continue;
   
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->body_id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << _w[i] << std::endl;    
    if (LOGGING(LOG_DYNAMICS) && link != body->get_base_link())
      FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (link->get_inertia() * _a[i]) << std::endl;

    // add I*a to the link force and Euler torque components 
    const SVELOCITY& vx = link->get_velocity(); 
    FILE_LOG(LOG_DYNAMICS) << "  I * v = " << (link->get_inertia() * vx) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  v x I * v = " << vx.cross(link->get_inertia() * vx) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force (before inertial forces): " << _w[i] << std::endl;
    // I*(v x s*qd) = momentum x momentum
    
    _w[i] += link->get_inertia() * _a[i] + vx.cross(link->get_inertia() * vx);
    FILE_LOG(LOG_DYNAMICS) << "  inertial force: " << (link->get_inertia() * _a[i] + vx.cross(link->get_inertia() * vx)) << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force (+ I*a + v x Iv): " << _w[i] << std::endl;

    // subtract external forces
    SFORCE wext = link->sum_forces(); 
    _w[i] -= wext;
    FILE_LOG(LOG_DYNAMICS) << "  external forces: " << wext << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << _w[i] << std::endl;

    // update the parent force and add parent for processing (if parent)
    shared_ptr<RIGIDBODY> parent = link->get_parent_link();
    if (parent)
    {
      unsigned h = parent->get_index();
      _w[h] += POSE3::transform(_w[h].pose, _w[i]);
      link_queue.push(parent);
    }

    // indicate that this link has been processed
    body->_processed[i] = true;
  }
  
  // ** STEP 3: compute centrifugal/Coriolis/gravity forces (C)

  // determine the length of the C vector
  const unsigned nDOF = body->num_joint_dof_explicit();

  // compute C
  C.resize(nDOF);
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get links, joints, etc.
    shared_ptr<JOINT> joint = ejoints[i];
    shared_ptr<RIGIDBODY> ob = joint->get_outboard_link();

    // get indices
    unsigned jidx = joint->get_coord_index();
    unsigned oidx = ob->get_index();

    // compute appropriate components of C
    SHAREDVECTORN Csub = C.segment(jidx, jidx+joint->num_dof()); 
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes();
    transform_and_transpose_mult(s, _w[oidx], Csub);

    FILE_LOG(LOG_DYNAMICS) << " -- computing C for link " << ob << std::endl;
    if (LOGGING(LOG_DYNAMICS))
    {
      SFORCE w_mixed = POSE3::transform(ob->get_mixed_pose(), _w[oidx]);
      std::vector<SVELOCITY> s_mixed;
      POSE3::transform(ob->get_mixed_pose(), s, s_mixed);
      if (!s.empty())
      {
        FILE_LOG(LOG_DYNAMICS) << " -- s' (mixed pose): " << s_mixed[0] << std::endl;
        FILE_LOG(LOG_DYNAMICS) << " -- force (mixed pose): " << w_mixed << std::endl;
      }
    }
    FILE_LOG(LOG_DYNAMICS) << "   -- forces: " << _w[oidx] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "   -- component of C: " << Csub << std::endl;
  }  

  FILE_LOG(LOG_DYNAMICS) << "------------------------------------------------" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "forces on base: " << _w[0] << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() exited" << std::endl;

  // store forces on base
  f0 = _w[0]; 
}

/// Computes the vector "C", w/o inertial forces, using the recursive Newton-Euler algorithm
/**
 * \return the spatial vector of forces on the base, which can be ignored for forward dynamics for fixed bases
 */
void CRB_ALGORITHM::calc_generalized_forces_noinertial(SFORCE& f0, VECTORN& C)
{
  const unsigned SPATIAL_DIM = 6;
  queue<shared_ptr<RIGIDBODY> > link_queue;
  SFORCE w;

  // get the body and the reference frame
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // get the set of links and joints
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();
  if (links.empty())
  {
    C.resize(0);
    f0.set_zero();
    return;
  }

  // **************************************************************************
  // first, compute forward dynamics using RNE algorithm; we copy the algorithm
  // here, because we need some data out of it
  // **************************************************************************

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() entered" << std::endl;

  // ** STEP 1: compute accelerations

  // get the base link
  shared_ptr<RIGIDBODY> base = links.front();

  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;
  body->_processed[base->get_index()] = true;


  // ** STEP 1: compute link forces -- backward recursion
  // use a map to determine which links have been processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // setup a map of link forces, all set to zero initially
  _w.resize(links.size());
  for (unsigned i=0; i< links.size(); i++)
  {
    _w[i].set_zero();
    _w[i].pose = links[i]->get_computation_frame();
  }

  // add all leaf links to the queue
  for (unsigned i=0; i< links.size(); i++)
    if (body->treat_link_as_leaf(links[i]))
      link_queue.push(links[i]);
      
  // process all links up to, but not including, the base
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();

    // if this link has already been processed, do not process it again
    if (body->_processed[i])
      continue;

    // verify all children have been processed
    if (!body->all_children_processed(link))
      continue;
   
    FILE_LOG(LOG_DYNAMICS) << " computing necessary force; processing link " << link->body_id << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  currently determined link force: " << _w[i] << std::endl;    
    if (LOGGING(LOG_DYNAMICS) && link != body->get_base_link())
      FILE_LOG(LOG_DYNAMICS) << "  I * a = " << (link->get_inertia() * _a[i]) << std::endl;

    // subtract external forces
    SFORCE wext = link->sum_forces(); 
    _w[i] -= wext;
    FILE_LOG(LOG_DYNAMICS) << "  external forces: " << wext << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  force on link after subtracting external force: " << _w[i] << std::endl;

    // update the parent force and add parent for processing (if parent)
    shared_ptr<RIGIDBODY> parent = link->get_parent_link();
    if (parent)
    {
      unsigned h = parent->get_index();
      _w[h] += POSE3::transform(_w[h].pose, _w[i]);
      link_queue.push(parent);
    }

    // indicate that this link has been processed
    body->_processed[i] = true;
  }
  
  // ** STEP 2: compute C

  // determine the length of the C vector
  const unsigned nDOF = body->num_joint_dof_explicit();

  // compute C
  C.resize(nDOF);
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    // get links, joints, etc.
    shared_ptr<JOINT> joint = ejoints[i];
    shared_ptr<RIGIDBODY> ob = joint->get_outboard_link();

    // get indices
    unsigned jidx = joint->get_coord_index();
    unsigned oidx = ob->get_index();

    // compute appropriate components of C
    SHAREDVECTORN Csub = C.segment(jidx, jidx+joint->num_dof()); 
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes();
    transform_and_transpose_mult(s, _w[oidx], Csub);

    FILE_LOG(LOG_DYNAMICS) << " -- computing C for link " << ob << std::endl;
    if (LOGGING(LOG_DYNAMICS))
    {
      SFORCE w_mixed = POSE3::transform(ob->get_mixed_pose(), _w[oidx]);
      std::vector<SVELOCITY> s_mixed;
      POSE3::transform(ob->get_mixed_pose(), s, s_mixed);
      if (!s.empty())
      {
        FILE_LOG(LOG_DYNAMICS) << " -- s' (mixed pose): " << s_mixed[0] << std::endl;
        FILE_LOG(LOG_DYNAMICS) << " -- force (mixed pose): " << w_mixed << std::endl;
      }
    }
    FILE_LOG(LOG_DYNAMICS) << "   -- forces: " << _w[oidx] << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "   -- component of C: " << Csub << std::endl;
  }  

  FILE_LOG(LOG_DYNAMICS) << "------------------------------------------------" << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "forces on base: " << _w[0] << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::calc_generalized_forces() exited" << std::endl;

  // store forces on base
  f0 = _w[0]; 
}

/// Updates all link accelerations (except the base)
void CRB_ALGORITHM::update_link_accelerations(shared_ptr<RC_ARTICULATED_BODY> body)
{
  queue<shared_ptr<RIGIDBODY> > link_queue;

  // get the set of links and their velocities
  const vector<shared_ptr<RIGIDBODY> >& links = body->get_links();

  // if there are no links, there is nothing to do
  if (links.empty())
    return;

  // mark all links as not processed
  for (unsigned i=0; i< links.size(); i++)
    body->_processed[i] = false;

  // get the base link and mark it as processed
  shared_ptr<RIGIDBODY> base = links.front();
  body->_processed[base->get_index()] = true;
  
  // get the spatial acceleration of the base link (should have already been
  // computed)
  base->set_accel(this->_a0);

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::update_link_accelerations() entered" << std::endl;
  
  // add all children of the base to the link queue
  list<shared_ptr<RIGIDBODY> > children;
  base->get_child_links(std::back_inserter(children));
  BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, children)
    link_queue.push(rb);
  
  // propagate link accelerations 
  while (!link_queue.empty())
  {
    // get the link off of the front of the queue 
    shared_ptr<RIGIDBODY> link = link_queue.front();
    link_queue.pop();    
    unsigned i = link->get_index();
    body->_processed[i] = true;

    // push all unprocessed children of the link onto the queue
    list<shared_ptr<RIGIDBODY> > children;
    link->get_child_links(std::back_inserter(children));
    BOOST_FOREACH(shared_ptr<RIGIDBODY> rb, children)
      if (!body->_processed[rb->get_index()])
        link_queue.push(rb);

    // get the inner joint and the parent link
    shared_ptr<JOINT> joint(link->get_inner_joint_explicit());
    shared_ptr<RIGIDBODY> parent(link->get_parent_link());
    unsigned h = parent->get_index();
 
    // set link acceleration
    const SACCEL& ah = parent->get_accel();
    SACCEL ai = SPARITH::transform_accel(link->get_accel().pose, ah);

    // get the link spatial axis
    const std::vector<SVELOCITY>& s = joint->get_spatial_axes(); 

    // determine the link accel
    POSE3::transform(ai.pose, s, _sprime);
    if (!_sprime.empty())
    {
      SVELOCITY sqd = SPARITH::mult(_sprime, joint->qd);
      ai += SACCEL(link->get_velocity().cross(sqd));
      ai += SACCEL(SPARITH::mult(_sprime, joint->qdd)); 
    }
    link->set_accel(ai);

    FILE_LOG(LOG_DYNAMICS) << "    -- updating link " << link << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- parent acceleration: " << ah << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- velocity: " << link->get_velocity() << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qd: " << joint->qd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- qdd: " << joint->qdd << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "      -- acceleration: " << ai << std::endl;
  }

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::update_link_accelerations() exited" << std::endl;
}

/*
/// Implements RCArticulatedBodyFwdDynAlgo::apply_impulse()
void CRB_ALGORITHM::apply_impulse(const SFORCE& w, shared_ptr<RIGIDBODY> link)
{
  // An alternative method for applying impulses using generalized coordinates
  // below...
  SMatrix6 Xcp = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, link->get_position());
  SMatrix6 X0p = SMatrix6::calc_spatial_transform(IDENTITY_3x3, point, IDENTITY_3x3, ZEROS_3);
  SVector6 jx(jj[0], jj[1], jj[2], jk[0], jk[1], jk[2]);
 
  // convert the impulse to a generalized impulse
  VECTORN gf;
  body->convert_to_generalized_force(link, point, jj, jk, gf);

  // get the generalized inertia and invert it
  MATRIXN M;
  body->get_generalized_inertia(M);
  LinAlg::inverse_PD(M);

  // get the generalized velocity
  VECTORN gv;
  body->get_generalized_velocity(gv);

  // compute the change in velocity
  VECTORN dv;
  M.mult(gf, dv);

  SMatrix6 Xc0 = SMatrix6::calc_spatial_transform(IDENTITY_3x3, ZEROS_3, IDENTITY_3x3, link->get_position());
  SMatrix6 X0c = SMatrix6::calc_spatial_transform(IDENTITY_3x3, link->get_position(), IDENTITY_3x3, ZEROS_3);
  MATRIXN Mx;
  body->get_generalized_inertia(Mx, eGlobal);
  SMatrix6 M2(Mx.begin());

  FILE_LOG(LOG_DYNAMICS) << "generalized inertia (global frame): " << std::endl << Mx;
  FILE_LOG(LOG_DYNAMICS) << "spatial transform (global to centroidal): " << std::endl << Xc0;
  FILE_LOG(LOG_DYNAMICS) << "[permuted] spatial inertia (centroidal frame): " << std::endl << (Xc0 * M2 * X0c);

  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (global frame): " << (X0p * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "spatial impulse (centroidal frame): " << (Xcp * jx) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "base transform: " << std::endl << body->get_links().front()->get_pose();
  FILE_LOG(LOG_DYNAMICS) << "impulse: " << jj << " / " << (jk + Vector3::cross(jj, r)) << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "generalized impulse: " << gf << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "inverse generalized inertia: " << std::endl << M;
  FILE_LOG(LOG_DYNAMICS) << "generalized v: " << gv << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "delta v: " << dv << std::endl;
  gv += dv;
  FILE_LOG(LOG_DYNAMICS) << "new v: " << gv << std::endl;
  body->set_generalized_velocity(gv);
  FILE_LOG(LOG_DYNAMICS) << "new base linear velocity: " << body->get_links().front()->get_lvel() << std::endl;
  FILE_LOG(LOG_DYNAMICS) << "new base angular velocity: " << body->get_links().front()->get_avel() << std::endl;
}
*/

/// TODO: fix this for bodies with kinematic loops
/// Applies an impulse to an articulated body with a floating base; complexity O(n^2)
void CRB_ALGORITHM::apply_impulse(const SMOMENTUM& w, shared_ptr<RIGIDBODY> link)
{
  const boost::shared_ptr<const POSE3> GLOBAL;
  const unsigned OPSPACE_DIM = 6;
  SVELOCITY dv0;
  VECTORN& b = _b;
  VECTORN& augV = _augV;
  VECTORN& workv = _workv;

  // get the body
  shared_ptr<RC_ARTICULATED_BODY> body(_body);

  // do necessary pre-calculations
  precalc(body); 

  // does not work for bodies with kinematic loops
  assert(body->_ejoints.empty());

  // get the base link
  shared_ptr<RIGIDBODY> base = body->get_base_link();

  // compute the Jacobian for the floating base w.r.t. the contact point
  _J.clear();

  // compute the Jacobian with respect to the contact point
  shared_ptr<RIGIDBODY> l = link;
  shared_ptr<JOINT> j;
  while ((j = l->get_inner_joint_explicit()))
  {
    // compute the Jacobian column(s) for the joint
    const std::vector<SVELOCITY>& s = j->get_spatial_axes();

    // transform spatial axes to global frame
    POSE3::transform(GLOBAL, s, _sprime);

    // set the column(s) of the Jacobian
    _J.insert(_J.end(), _sprime.begin(), _sprime.end());

    // set l to its parent
    l = shared_ptr<RIGIDBODY>(l->get_parent_link());
  }

  // transform the impulse to the global frame
  SMOMENTUM w0 = POSE3::transform(GLOBAL, w); 

  // compute the impulse applied to the joints
  SPARITH::transpose_mult(_J, w0, workv);

  FILE_LOG(LOG_DYNAMICS) << "  impulse (last frame): " << w0 << std::endl;

  // special case: floating base
  if (body->_floating_base)
  {
    // determine the index where the base starts
    const unsigned BASE_START = workv.size();

    // form vector to solve for b
    SPARITH::concat(workv, w0, b);
 
    // compute changes in base and joint velocities
    M_solve_noprecalc(workv = b);

    // swap base velocity change linear and angular components 
    const unsigned BASE_A = BASE_START + 0;
    const unsigned BASE_B = BASE_START + 1;
    const unsigned BASE_G = BASE_START + 2;
    const unsigned BASE_X = BASE_START + 3;
    const unsigned BASE_Y = BASE_START + 4;
    const unsigned BASE_Z = BASE_START + 5;
    std::swap(workv[BASE_A], workv[BASE_X]);
    std::swap(workv[BASE_B], workv[BASE_Y]);
    std::swap(workv[BASE_G], workv[BASE_Z]);
 
    // get change in base and change in joint velocities
    VECTOR3 dv0_angular(workv[BASE_A], workv[BASE_B], workv[BASE_G]);
    VECTOR3 dv0_linear(workv[BASE_X], workv[BASE_Y], workv[BASE_Z]);
    SVELOCITY dv0(dv0_angular, dv0_linear, w0.pose);

    // update the base velocity
    SVELOCITY basev = base->get_velocity();
    basev += POSE3::transform(basev.pose, dv0);
    base->set_velocity(basev);

    FILE_LOG(LOG_DYNAMICS) << "  change in base velocity: " << dv0 << std::endl;
    FILE_LOG(LOG_DYNAMICS) << "  new base velocity: " << basev << std::endl;
  }

  // apply the change and update link velocities
  const vector<shared_ptr<JOINT> >& ejoints = body->get_explicit_joints();
  for (unsigned i=0; i< ejoints.size(); i++)
  {
    unsigned idx = ejoints[i]->get_coord_index();
    FILE_LOG(LOG_DYNAMICS) << " joint " << ejoints[i] << " qd: " << ejoints[i]->qd << "  dqd: " << workv.segment(idx, idx+ejoints[i]->num_dof()) << std::endl;  
    ejoints[i]->qd += workv.segment(idx, idx+ejoints[i]->num_dof());
  }  
  body->update_link_velocities();

  // reset all force and torque accumulators -- impulses drive them to zero
  body->reset_accumulators();

  FILE_LOG(LOG_DYNAMICS) << "CRBAlgorithm::apply_impulse() exited" << std::endl;
}


