/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Ravelin/UndefinedAxisException.h>

using boost::shared_ptr;
using std::vector;
using namespace Ravelin;

/// Initializes the joint
/**
 * The axes of rotation are each set to [0 0 0].  The inboard
 * and outboard links are set to NULL.
 */
UNIVERSALJOINT::UNIVERSALJOINT() : JOINT()
{
  const unsigned SPATIAL_DIM = 6;

  // init joint data
  init_data();

  // init the joint axes
  _u[eAxis1].set_zero();
  _u[eAxis2].set_zero();
  _u[eAxis1].pose = _F;
  _u[eAxis2].pose = _F;
  _h2.set_zero();

  // setup the spatial axis derivative to zero
  _s_dot.resize(num_dof());
  for (unsigned i=0; i< num_dof(); i++)
    _s_dot[i].pose = _F;
}

/// Gets the axis for this joint
VECTOR3 UNIVERSALJOINT::get_axis(Axis a) const
{
  const unsigned X = 0;

  // axis one is already set 
  if (a == eAxis1)
    return _u[DOF_1];
  else
    return _u[DOF_2];
}

/// Sets an axis of this joint
/**
 * The local axis for this joint does not take the orientation of the 
 * inboard link into account; thus, if the orientation of the inboard link 
 * changes, then the local axis remains constant.
 * \param axis a unit vector
 * \sa get_axis_global()
 * \sa set_axis_global()
 */
void UNIVERSALJOINT::set_axis(const VECTOR3& axis, Axis a) 
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // normalize the axis in case the user did not
  VECTOR3 naxis = VECTOR3::normalize(axis); 

  // set the axis
  _u[a] = POSE3::transform_vector(get_pose(), naxis); 

  // update the spatial axes
  update_spatial_axes(); 

  // set the second axis in outer link coordinates, if we just changed it
  if (a == eAxis2)
  {
    shared_ptr<const POSE3> Fi = get_inboard_pose();
    shared_ptr<const POSE3> Fo = get_outboard_pose();
    _h2 = POSE3::transform_vector(Fo, naxis);
  }
}        

/// Updates the spatial axis for this joint
void UNIVERSALJOINT::update_spatial_axes()
{
  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // call the parent method
  JOINT::update_spatial_axes();

  // update the spatial axes in joint pose 
  _s[0].set_angular(_u[DOF_1]);
  _s[0].set_linear(ZEROS_3);

  // update the spatial axes in link coordinates
  _s_dot[0].set_angular(ZEROS_3);
  _s_dot[0].set_linear(ZEROS_3);
}

/// Gets the spatial axes for this joint
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& UNIVERSALJOINT::get_spatial_axes()
{
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes() called with NULL outboard link");

  // get the transformed second axis
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 u1(R * ORIGIN3(_u[DOF_2]), get_pose());

  // update the second spatial axes
  _s[1].set_angular(u1);
  _s[1].set_linear(ZEROS_3);

  // use the JOINT function to do the rest
  return JOINT::get_spatial_axes();
}

/// Gets the derivative of the spatial-axis
/**
 * \note these spatial axes are not constant, unlike many joints.
 */
const vector<SVELOCITY>& UNIVERSALJOINT::get_spatial_axes_dot()
{
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // set zeros
  VECTOR3 ZEROS_3(0.0, 0.0, 0.0, get_pose());

  // get the inboard and outboard links
  shared_ptr<const POSE3> inboard = get_inboard_pose();
  shared_ptr<const POSE3> outboard = get_outboard_pose();
  if (!inboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes_dot() called with NULL inboard link");
  if (!outboard)
    throw std::runtime_error("UNIVERSALJOINT::get_spatial_axes_dot() called with NULL outboard link");

  // get the transformed second axis
  const double S1 = std::sin(q[DOF_1]+q_tare[DOF_1]);
  const double C1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  REAL qd1 = qd[DOF_1];
  VECTOR3 u2(0.0, -S1*qd1, C1*qd1, get_pose());

  // get the second spatial axis time derivative. This is:
  // R * \dot{u2} + \dot{R} * u
  // note that \dot{u2} is zero
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 omega = _u[DOF_1] * qd1;
  ORIGIN3 u1 = ORIGIN3::cross(ORIGIN3(omega), R * ORIGIN3(_u[DOF_2]));
  _s_dot[1].set_angular(VECTOR3(u1, get_pose()));
  _s_dot[1].set_linear(ZEROS_3);

/*
  // get q and qd
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;
  const VECTORN& qd = this->qd;

  // compute some needed quantities
  const REAL c1 = std::cos(q[DOF_1]+q_tare[DOF_1]);
  const REAL s1 = std::sin(q[DOF_1]+q_tare[DOF_1]);

  // form the time derivative of the spatial axis for the second DOF; note that spatial
  // axis for first DOF is constant, so time-derivative is zero 
  MATRIX3 Rx = MATRIX3::zero();
  REAL qd1 = qd[DOF_1];
  ORIGIN3 v(0.0, -s1*qd1, c1*qd1); 
  Rx.set_column(Y, v);
  ORIGIN3 cross = ORIGIN3::cross(ORIGIN3(1.0, 0.0, 0.0), v); 
  Rx.set_column(Z, cross);

  // compute u
  ORIGIN3 ux = (_R * Rx * MATRIX3::transpose(_R)).get_column(Y);
  VECTOR3 u(ux, get_pose());

  // get the time derivative of second spatial axis
  REAL qd1 = qd[DOF_1];
  VECTOR3 u2(0.0, -s1*qd1, c1*qd1, get_pose());

  // get the second spatial axis time derivative. This is:
  // R * \dot{u2} + \dot{R} * u
  // note that \dot{u2} is zero
  MATRIX3 R = AANGLE(_u[DOF_1], q[DOF_1]+q_tare[DOF_1]);
  VECTOR3 omega = _u[DOF_1] * qd1;
  VECTOR3 dotRu = VECTOR3::cross(omega, VECTOR3(R * ORIGIN3(_u[DOF_2]), get_pose()));

  // update the spatial axis in link coordinates; note that axis 1 is always
  // set to zero (init'd in constructor)
  _s_dot[1].set_upper(dotRu);
  _s_dot[1].set_lower(ZEROS_3);
*/

  return _s_dot;
}

/// Determines (and sets) the value of Q from the axes and the inboard link and outboard link transforms
void UNIVERSALJOINT::determine_q(VECTORN& q)
{
  std::cerr << "SPHERICALJOINT::determine_q() warning- determine_q(.) is not currently functional" << std::endl;
}

/// Gets the (local) transform for this joint
MATRIX3 UNIVERSALJOINT::get_rotation() const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get q and _q_tare
  const VECTORN& q = this->q;
  const VECTORN& q_tare = this->_q_tare;

  // compute some needed quantities
  const REAL Q1 = q[DOF_1]+q_tare[DOF_1];
  const REAL Q2 = q[DOF_2]+q_tare[DOF_2];

  // get the first rotation matrix
  ORIGIN3 u1(_u[DOF_1]);
  MATRIX3 R = AANGLE(u1, Q1);

  // get the second axis 
  ORIGIN3 u2 = R * ORIGIN3(_u[DOF_2]);

  // transform using orientation matrix
  return R * AANGLE(u2, Q2);
}

/// Gets the transform induced by this joint
shared_ptr<const POSE3> UNIVERSALJOINT::get_induced_pose()
{
  // get the rotation
  _Fprime->q = get_rotation();
  return _Fprime;
}

/// Computes the constraint jacobian
/*
void UNIVERSALJOINT::calc_constraint_jacobian(RigidBodyPtr body, unsigned index, REAL Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  shared_ptr<const POSE3> inner = get_inboard_pose();
  shared_ptr<const POSE3> outer = get_outboard_pose();

  // make sure that _u (and by extension _h2) is set
  if (_u[eAxis1].norm_sq() < std::numeric_limits<REAL>::epsilon() ||
      _u[eAxis2].norm_sq() < std::numeric_limits<REAL>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (REAL) 0.0;
    return;
  }

  // get the information necessary to compute the constraint equations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const VECTOR3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const VECTOR3& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
  const REAL q1x = q1.x;
  const REAL q1y = q1.y;
  const REAL q1z = q1.z;
  const REAL q1w = q1.w;
  const REAL p1x = p1[X];
  const REAL p1y = p1[Y];
  const REAL p1z = p1[Z];
  const REAL q2x = q2.x;
  const REAL q2y = q2.y;
  const REAL q2z = q2.z;
  const REAL q2w = q2.w;
  const REAL p2x = -p2[X];
  const REAL p2y = -p2[Y];
  const REAL p2z = -p2[Z];
  const REAL u0x = _u[DOF_1][X];
  const REAL u0y = _u[DOF_1][Y];
  const REAL u0z = _u[DOF_1][Z];
  const REAL h2x = _h2[X];
  const REAL h2y = _h2[Y];
  const REAL h2z = _h2[Z];

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
  {
    switch (index)
    {
      case 0:
        Cq[0] = 1.0;    
        Cq[1] = 0.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*p1x*q1w + 2*p1z*q1y - 2*p1y*q1z;
        Cq[4] = 4*q1x*p1x + 2*q1y*p1y + 2*q1z*p1z;
        Cq[5] = 2*p1z*q1w + 2*p1y*q1x;
        Cq[6] = 2*p1z*q1x - 2*p1y*q1w;
        break;

      case 1:
        Cq[0] = 0.0;    
        Cq[1] = 1.0;    
        Cq[2] = 0.0;    
        Cq[3] = 4*p1y*q1w - 2*p1z*q1x + 2*p1x*q1z;
        Cq[4] = 2*q1y*p1x - 2*q1w*p1z;
        Cq[5] = 2*p1x*q1x + 4*p1y*q1y + 2*p1z*q1z;
        Cq[6] = 2*p1x*q1w + 2*p1z*q1y;
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 1.0;
        Cq[3] = 4*p1z*q1w + 2*p1y*q1x - 2*p1x*q1y;
        Cq[4] = 2*q1z*p1x + 2*q1w*p1y;
        Cq[5] = 2*p1y*q1z - 2*p1x*q1w;
        Cq[6] = 4*p1z*q1z + 2*p1y*q1y + 2*p1x*q1x;
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                2*(q2x*q2y + q2w*q2z)*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                (-1 + 2*q2w*q2w + 2*q2x*q2x)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z)) + 
                h2y*(2*(q2w*q2x + q2y*q2z)*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                (-1 + 2*(q2w*q2w + q2y*q2y))*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                2*(q2x*q2y - q2w*q2z)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z)) + 
                h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
                (-2*q1y*u0x + 2*q1x*u0y + 4*q1w*u0z) + 
                2*(-(q2w*q2x) + q2y*q2z)*
                (2*q1z*u0x + 4*q1w*u0y - 2*q1x*u0z) + 
                2*(q2w*q2y + q2x*q2z)*
                (4*q1w*u0x - 2*q1z*u0y + 2*q1y*u0z));
        Cq[4] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
                (2*q1z*u0x + 2*q1w*u0y) + 
                2*(q2x*q2y + q2w*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
                (-1 + 2*q2w*q2w + 2*q2x*q2x)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
                h2y*(2*(q2w*q2x + q2y*q2z)*(2*q1z*u0x + 2*q1w*u0y) + 
                (-1 + 2*(q2w*q2w + q2y*q2y))*
                (2*q1y*u0x - 2*q1w*u0z) + 
                2*(q2x*q2y - q2w*q2z)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
                h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
                (2*q1z*u0x + 2*q1w*u0y) + 
                2*(-(q2w*q2x) + q2y*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
                2*(q2w*q2y + q2x*q2z)*
                (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z));
        Cq[5] = h2x*(2*(-(q2w*q2y) + q2x*q2z)*
               (2*q1z*u0x + 2*q1w*u0y) + 
               2*(q2x*q2y + q2w*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
               (-1 + 2*q2w*q2w + 2*q2x*q2x)*
               (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
               h2y*(2*(q2w*q2x + q2y*q2z)*(2*q1z*u0x + 2*q1w*u0y) + 
               (-1 + 2*(q2w*q2w + q2y*q2y))*
               (2*q1y*u0x - 2*q1w*u0z) + 
               2*(q2x*q2y - q2w*q2z)*
               (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z)) + 
               h2z*((-1 + 2*(q2w*q2w + q2z*q2z))*
              (2*q1z*u0x + 2*q1w*u0y) + 
              2*(-(q2w*q2x) + q2y*q2z)*(2*q1y*u0x - 2*q1w*u0z) + 
              2*(q2w*q2y + q2x*q2z)*
              (4*q1x*u0x + 2*q1y*u0y + 2*q1z*u0z));
        Cq[6] = h2x*((-1 + 2*q2w*q2w + 2*q2x*q2x)*
              (-2*q1w*u0y + 2*q1x*u0z) + 
              2*(q2x*q2y + q2w*q2z)*(2*q1w*u0x + 2*q1y*u0z) + 
              2*(-(q2w*q2y) + q2x*q2z)*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z)) + 
              h2y*(2*(q2x*q2y - q2w*q2z)*(-2*q1w*u0y + 2*q1x*u0z) + 
              (-1 + 2*(q2w*q2w + q2y*q2y))*
              (2*q1w*u0x + 2*q1y*u0z) + 
              2*(q2w*q2x + q2y*q2z)*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z)) + 
              h2z*(2*(q2w*q2y + q2x*q2z)*(-2*q1w*u0y + 2*q1x*u0z) + 
              2*(-(q2w*q2x) + q2y*q2z)*(2*q1w*u0x + 2*q1y*u0z) + 
              (-1 + 2*(q2w*q2w + q2z*q2z))*
              (2*q1x*u0x + 2*q1y*u0y + 4*q1z*u0z));

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    switch (index)
    {
      case 0:
        Cq[0] = -1.0;     
        Cq[1] = 0.0;      
        Cq[2] = 0.0;      
        Cq[3] = -(4*p2x*q2w + 2*p2z*q2y - 2*p2y*q2z);
        Cq[4] = -(4*q2x*p2x + 2*q2y*p2y + 2*q2z*p2z);
        Cq[5] = -(2*p2z*q2w + 2*p2y*q2x);
        Cq[6] = -(2*p2z*q2x - 2*p2y*q2w);
        break;

      case 1:
        Cq[0] = 0.0;      
        Cq[1] = -1.0;     
        Cq[2] = 0.0;      
        Cq[3] = -(4*p2y*q2w - 2*p2z*q2x + 2*p2x*q2z);
        Cq[4] = -(2*q2y*p2x - 2*q2w*p2z);
        Cq[5] = -(2*p2x*q2x + 4*p2y*q2y + 2*p2z*q2z);
        Cq[6] = -(2*p2x*q2w + 2*p2z*q2y);
        break;

      case 2:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = -1.0;
        Cq[3] = -(4*p2z*q2w + 2*p2y*q2x - 2*p2x*q2y);
        Cq[4] = -(2*q2z*p2x + 2*q2w*p2y);
        Cq[5] = -(2*p2y*q2z - 2*p2x*q2w);
        Cq[6] = -(4*p2z*q2z + 2*p2y*q2y + 2*p2x*q2x);
        break;

      case 3:
        Cq[0] = 0.0;
        Cq[1] = 0.0;
        Cq[2] = 0.0;
        Cq[3] = h2z*(2*q2y*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) - 
                2*q2x*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                4*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(-2*q2z*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                4*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2x*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2x*(4*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2z*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) - 
                2*q2y*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[4] = h2z*(2*q2z*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) - 
                2*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z)) + 
                h2y*(2*q2y*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2x*(4*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[5] = h2z*(2*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*
                u0x + 2*(q1x*q1y - q1w*q1z)*u0y + 
                2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2z*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z)) + 
                h2x*(2*q2x*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) - 
                2*q2w*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(2*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                4*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));
        Cq[6] = h2x*(2*q2w*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                2*q2x*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2y*(-2*q2w*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z)) + 
                h2z*(2*q2x*((-1 + 2*q1w*q1w + 2*q1x*q1x)*u0x + 
                2*(q1x*q1y - q1w*q1z)*u0y + 2*(q1w*q1y + q1x*q1z)*u0z) + 
                2*q2y*(2*(q1x*q1y + q1w*q1z)*u0x + 
                (-1 + 2*(q1w*q1w + q1y*q1y))*u0y + 
                2*(-(q1w*q1x) + q1y*q1z)*u0z) + 
                4*q2z*(2*(-(q1w*q1y) + q1x*q1z)*u0x + 
                2*(q1w*q1x + q1y*q1z)*u0y + 
                (-1 + 2*(q1w*q1w + q1z*q1z))*u0z));

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}

/// Computes the constraint jacobian
void UNIVERSALJOINT::calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, REAL Cq[7])
{
  const unsigned X = 0, Y = 1, Z = 2, SPATIAL_DIM = 7;

  // get the two links
  shared_ptr<const POSE3> inner = get_inboard_pose();
  shared_ptr<const POSE3> outer = get_outboard_pose();

  // make sure that _u (and by extension _h2) is set
  if (_u[eAxis1].norm_sq() < std::numeric_limits<REAL>::epsilon() ||
      _u[eAxis2].norm_sq() < std::numeric_limits<REAL>::epsilon())
    throw UndefinedAxisException(); 

  // mke sure that body is one of the links
  if (inner != body && outer != body)
  {
    for (unsigned i=0; i< SPATIAL_DIM; i++)
      Cq[i] = (REAL) 0.0;
    return;
  }

  // get the information necessary to compute the constraint equations
  const Quat& q1 = inner->get_orientation();
  const Quat& q2 = outer->get_orientation();
  const Quat qd1 = Quat::deriv(q1, inner->get_avel());
  const Quat qd2 = Quat::deriv(q2, outer->get_avel());
  const VECTOR3& p1 = inner->get_outer_joint_data(outer).com_to_joint_vec;
  const VECTOR3& p2 = body->get_inner_joint_data(inner).joint_to_com_vec_of;
  const REAL qx1 = q1.x;
  const REAL qy1 = q1.y;
  const REAL qz1 = q1.z;
  const REAL qw1 = q1.w;
  const REAL p1x = p1[X];
  const REAL p1y = p1[Y];
  const REAL p1z = p1[Z];
  const REAL qx2 = q2.x;
  const REAL qy2 = q2.y;
  const REAL qz2 = q2.z;
  const REAL qw2 = q2.w;
  const REAL p2x = -p2[X];
  const REAL p2y = -p2[Y];
  const REAL p2z = -p2[Z];
  const REAL ux = _u[DOF_1][X];
  const REAL uy = _u[DOF_1][Y];
  const REAL uz = _u[DOF_1][Z];
  const REAL h2x = _h2[X];
  const REAL h2y = _h2[Y];
  const REAL h2z = _h2[Z];
  const REAL dqw1 = qd1.w;
  const REAL dqx1 = qd1.x;
  const REAL dqy1 = qd1.y;
  const REAL dqz1 = qd1.z;
  const REAL dqw2 = qd2.w;
  const REAL dqx2 = qd2.x;
  const REAL dqy2 = qd2.y;
  const REAL dqz2 = qd2.z;

  // setup the constraint equations (from Shabana, p. 436), eq. 7.176
  if (body == inner)
  {
    switch (index)
    {
      case 0:
        Cq[0] = (REAL) 0.0;    
        Cq[1] = (REAL) 0.0;    
        Cq[2] = (REAL) 0.0;    
        Cq[3] = 4*p1x*dqw1 + 2*p1z*dqy1 - 2*p1y*dqz1;
        Cq[4] = 4*dqx1*p1x + 2*dqy1*p1y + 2*dqz1*p1z;
        Cq[5] = 2*p1z*dqw1 + 2*p1y*dqx1;
        Cq[6] = 2*p1z*dqx1 - 2*p1y*dqw1;
        break;

      case 1:
        Cq[0] = (REAL) 0.0;    
        Cq[1] = (REAL) 0.0;    
        Cq[2] = (REAL) 0.0;    
        Cq[3] = 4*p1y*dqw1 - 2*p1z*dqx1 + 2*p1x*dqz1;
        Cq[4] = 2*dqy1*p1x - 2*dqw1*p1z;
        Cq[5] = 2*p1x*dqx1 + 4*p1y*dqy1 + 2*p1z*dqz1;
        Cq[6] = 2*p1x*dqw1 + 2*p1z*dqy1;
        break;

      case 2:
        Cq[0] = (REAL) 0.0;
        Cq[1] = (REAL) 0.0;
        Cq[2] = (REAL) 0.0;
        Cq[3] = 4*p1z*dqw1 + 2*p1y*dqx1 - 2*p1x*dqy1;
        Cq[4] = 2*dqz1*p1x + 2*dqw1*p1y;
        Cq[5] = 2*p1y*dqz1 - 2*p1x*dqw1;
        Cq[6] = 4*p1z*dqz1 + 2*p1y*dqy1 + 2*p1x*dqx1;
        break;

      case 3:
        Cq[0] = (REAL) 0.0;
        Cq[1] = (REAL) 0.0;
        Cq[2] = (REAL) 0.0;
        Cq[3] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (-2*dqy1*ux + 2*dqx1*uy + 4*dqw1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqz1*ux + 4*dqw1*uy - 2*dqx1*uz) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (4*dqw1*ux - 2*dqz1*uy + 2*dqy1*uz) + 
   (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qy1*ux + 2*qx1*uy + 4*qw1*uz)\
    + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qz1*ux + 4*qw1*uy - 2*qx1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qw1*ux - 2*qz1*uy + 2*qy1*uz);
        Cq[4] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*(2*dqz1*ux + 2*dqw1*uy)
     + (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qz1*ux + 2*qw1*uy) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqy1*ux - 2*dqw1*uz) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (4*dqx1*ux + 2*dqy1*uy + 2*dqz1*uz) + 
   (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qy1*ux - 2*qw1*uz) + (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (4*qx1*ux + 2*qy1*uy + 2*qz1*uz);
        Cq[5] = (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (-2*dqw1*ux + 2*dqz1*uy) + 
   (2*h2x*(-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(-2*qw1*ux + 2*qz1*uy) + 
   (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (2*dqx1*uy + 2*dqw1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqx1*ux + 4*dqy1*uy + 2*dqz1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (2*qx1*uy + 2*qw1*uz) + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qx1*ux + 4*qy1*uy + 2*qz1*uz);
        Cq[6] = (h2x*(-1 + 2*qw2*qw2 + 2*qx2*qx2) + 
      2*h2y*(qx2*qy2 - qw2*qz2) + 2*h2z*(qw2*qy2 + qx2*qz2))*
    (-2*dqw1*uy + 2*dqx1*uz) + 
   (h2y*(-1 + 2*(qw2*qw2 + qy2*qy2)) + 
      2*h2x*(qx2*qy2 + qw2*qz2) + 2*h2z*(-(qw2*qx2) + qy2*qz2))*
    (2*dqw1*ux + 2*dqy1*uz) + 
   (2*h2x*(-(qw2*qy2) + qx2*qz2) + 2*h2y*(qw2*qx2 + qy2*qz2) + 
      h2z*(-1 + 2*(qw2*qw2 + qz2*qz2)))*
    (2*dqx1*ux + 2*dqy1*uy + 4*dqz1*uz) + 
   (h2x*(4*dqw2*qw2 + 4*dqx2*qx2) + 
      2*h2y*(-(dqz2*qw2) + dqy2*qx2 + dqx2*qy2 - dqw2*qz2) + 
      2*h2z*(dqy2*qw2 + dqz2*qx2 + dqw2*qy2 + dqx2*qz2))*
    (-2*qw1*uy + 2*qx1*uz) + (2*h2y*(2*dqw2*qw2 + 2*dqy2*qy2) + 
      2*h2x*(dqz2*qw2 + dqy2*qx2 + dqx2*qy2 + dqw2*qz2) + 
      2*h2z*(-(dqx2*qw2) - dqw2*qx2 + dqz2*qy2 + dqy2*qz2))*
    (2*qw1*ux + 2*qy1*uz) + (2*h2x*
       (-(dqy2*qw2) + dqz2*qx2 - dqw2*qy2 + dqx2*qz2) + 
      2*h2y*(dqx2*qw2 + dqw2*qx2 + dqz2*qy2 + dqy2*qz2) + 
      2*h2z*(2*dqw2*qw2 + 2*dqz2*qz2))*(2*qx1*ux + 2*qy1*uy + 4*qz1*uz);

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
  else
  {
    switch (index)
    {
      case 0:
        Cq[0] = (REAL) 0.0;     
        Cq[1] = (REAL) 0.0;      
        Cq[2] = (REAL) 0.0;      
        Cq[3] = -(4*p2x*dqw2 + 2*p2z*dqy2 - 2*p2y*dqz2);
        Cq[4] = -(4*dqx2*p2x + 2*dqy2*p2y + 2*dqz2*p2z);
        Cq[5] = -(2*p2z*dqw2 + 2*p2y*dqx2);
        Cq[6] = -(2*p2z*dqx2 - 2*p2y*dqw2);
        break;

      case 1:
        Cq[0] = (REAL) 0.0;      
        Cq[1] = (REAL) 0.0;     
        Cq[2] = (REAL) 0.0;      
        Cq[3] = -(4*p2y*dqw2 - 2*p2z*dqx2 + 2*p2x*dqz2);
        Cq[4] = -(2*dqy2*p2x - 2*dqw2*p2z);
        Cq[5] = -(2*p2x*dqx2 + 4*p2y*dqy2 + 2*p2z*dqz2);
        Cq[6] = -(2*p2x*dqw2 + 2*p2z*dqy2);
        break;

      case 2:
        Cq[0] = (REAL) 0.0;
        Cq[1] = (REAL) 0.0;
        Cq[2] = (REAL) 0.0;
        Cq[3] = -(4*p2z*dqw2 + 2*p2y*dqx2 - 2*p2x*dqy2);
        Cq[4] = -(2*dqz2*p2x + 2*dqw2*p2y);
        Cq[5] = -(2*p2y*dqz2 - 2*p2x*dqw2);
        Cq[6] = -(4*p2z*dqz2 + 2*p2y*dqy2 + 2*p2x*dqx2);
        break;

      case 3:
        Cq[0] = (REAL) 0.0;
        Cq[1] = (REAL) 0.0;
        Cq[2] = (REAL) 0.0;
        Cq[3] = (4*h2x*qw2 + 2*h2z*qy2 - 2*h2y*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (4*h2y*qw2 - 2*h2z*qx2 + 2*h2x*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (4*h2z*qw2 + 2*h2y*qx2 - 2*h2x*qy2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (4*dqw2*h2x - 2*dqz2*h2y + 2*dqy2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqz2*h2x + 4*dqw2*h2y - 2*dqx2*h2z)*
    (2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (-2*dqy2*h2x + 2*dqx2*h2y + 4*dqw2*h2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[4] = (4*h2x*qx2 + 2*h2y*qy2 + 2*h2z*qz2)*
    ((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (-2*h2z*qw2 + 2*h2x*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (2*h2y*qw2 + 2*h2x*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (4*dqx2*h2x + 2*dqy2*h2y + 2*dqz2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqy2*h2x - 2*dqw2*h2z)*(2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (2*dqz2*h2x + 2*dqw2*h2y)*(2*(-(qw1*qy1) + qx1*qz1)*ux + 
      2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[5] = (2*h2z*qw2 + 2*h2y*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (2*h2x*qx2 + 4*h2y*qy2 + 2*h2z*qz2)*
    (2*(dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (-2*h2x*qw2 + 2*h2y*qz2)*(2*
       (-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (2*dqx2*h2y + 2*dqw2*h2z)*((-1 + 2*qw1*qw1 + 2*qx1*qx1)*
       ux + 2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqx2*h2x + 4*dqy2*h2y + 2*dqz2*h2z)*
    (2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (-2*dqw2*h2x + 2*dqz2*h2y)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);
        Cq[6] = (-2*h2y*qw2 + 2*h2z*qx2)*((4*dqw1*qw1 + 4*dqx1*qx1)*ux + 
      2*(-(dqz1*qw1) + dqy1*qx1 + dqx1*qy1 - dqw1*qz1)*uy + 
      2*(dqy1*qw1 + dqz1*qx1 + dqw1*qy1 + dqx1*qz1)*uz) + 
   (2*h2x*qw2 + 2*h2z*qy2)*(2*
       (dqz1*qw1 + dqy1*qx1 + dqx1*qy1 + dqw1*qz1)*ux + 
      2*(2*dqw1*qw1 + 2*dqy1*qy1)*uy + 
      2*(-(dqx1*qw1) - dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uz) + 
   (2*h2x*qx2 + 2*h2y*qy2 + 4*h2z*qz2)*
    (2*(-(dqy1*qw1) + dqz1*qx1 - dqw1*qy1 + dqx1*qz1)*ux + 
      2*(dqx1*qw1 + dqw1*qx1 + dqz1*qy1 + dqy1*qz1)*uy + 
      2*(2*dqw1*qw1 + 2*dqz1*qz1)*uz) + 
   (-2*dqw2*h2y + 2*dqx2*h2z)*
    ((-1 + 2*qw1*qw1 + 2*qx1*qx1)*ux + 
      2*(qx1*qy1 - qw1*qz1)*uy + 2*(qw1*qy1 + qx1*qz1)*uz) + 
   (2*dqw2*h2x + 2*dqy2*h2z)*(2*(qx1*qy1 + qw1*qz1)*ux + 
      (-1 + 2*(qw1*qw1 + qy1*qy1))*uy + 
      2*(-(qw1*qx1) + qy1*qz1)*uz) + 
   (2*dqx2*h2x + 2*dqy2*h2y + 4*dqz2*h2z)*
    (2*(-(qw1*qy1) + qx1*qz1)*ux + 2*(qw1*qx1 + qy1*qz1)*uy + 
      (-1 + 2*(qw1*qw1 + qz1*qz1))*uz);

      default:
        throw std::runtime_error("Invalid joint constraint index!");
    }
  }
}
*/

/// Evaluates the constraint equations
void UNIVERSALJOINT::evaluate_constraints(REAL C[])
{
  shared_ptr<const POSE3> GLOBAL;
  const unsigned X = 0, Y = 1, Z = 2;

  // get the two links
  shared_ptr<const POSE3> inner = get_inboard_pose();
  shared_ptr<const POSE3> outer = get_outboard_pose();

  // This code was developed using [Shabana, 2003], p. 438; variable names
  // have been altered however

  // determine h1 and h2 in global coordinates
  VECTOR3 h1 = inner->transform_vector(GLOBAL, _u[DOF_1]);
  VECTOR3 h2 = outer->transform_vector(GLOBAL, _h2);

  // determine the global positions of the attachment points and subtract them
  VECTOR3 r1 = get_location(false);
  VECTOR3 r2 = get_location(true);
  VECTOR3 r12 = r1 - r2;

  // evaluate the constraint equations
  C[0] = r12[X];
  C[1] = r12[Y];
  C[2] = r12[Z];
  C[3] = h1.dot(h2);
}


