/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// ------------------------------------------------------------------
// Example of a double pendulum, adapted from code by James R. Taylor 
// ------------------------------------------------------------------

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "integrate.h"
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/RevoluteJointd.h>

using boost::shared_ptr;
using namespace Ravelin;

int main( void ) {
  const unsigned X = 0, Y = 1, Z = 2;
  const shared_ptr<const Pose3d> GLOBAL_3D;

  // current time
  double t = 0.0;

  // integration step size
  const double DT = 0.01;

// uncomment to log dynamics
//  Log<Ravelin::OutputToFile>::reporting_level = LOG_DYNAMICS;

  // the double pendulum has joints, so we need an articulated body 
  shared_ptr<RCArticulatedBodyd> pendulum( new RCArticulatedBodyd() );
  pendulum->body_id = "double-pendulum";

  // set the algorithm to the Composite Rigid Body Method
  pendulum->algorithm_type = RCArticulatedBodyd::eCRB;

  // set the dynamics computations to take place in link inertia frames
  // (these are frame defined at the center-of-mass of each link and changing
  //  orientation as the link changes orientation)
  pendulum->set_computation_frame_type(eLinkInertia);

  // vectors for references used in defining the body
  std::vector< shared_ptr<RigidBodyd> > links;
  std::vector< shared_ptr<Jointd > > joints;

  // create the fixed base;
  shared_ptr<RigidBodyd> base( new RigidBodyd() );
  {
    // identify the link
    base->body_id = "base";

    // this next line isn't necessary- it disables physics computation for the
    // base- because we will denote the articulated body as having a fixed
    // base, but it is good practice
    base->set_enabled( false );
 
    // set the pose of the base 
    Quatd ori(Quatd::normalize(Quatd(0,0,0,1)));
    Origin3d position(0,0,0);
    Pose3d pose( ori, position );
    base->set_pose( pose ); 

    // add the base to the set of links
    links.push_back( base );
  }

  // the pendulum body is a cylinder
  shared_ptr<RigidBodyd> link1( new RigidBodyd() );
  {
    // setup cylinder radius and height
    const double R = 0.025;
    const double H = 1.0;

    // setup inertia for cylinder- it will be defined with respect to the
    // centroid of the cylinder
    SpatialRBInertiad J;
    J.pose = link1->get_pose();
    J.m = 1.0;

    // inertia matrix for cylinder is taken from online resources; note that
    // primary axis of cylinder is aligned with y-axis 
    J.J.set_zero(3,3);
    J.J(X,X) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 
    J.J(Y,Y) = 0.5*J.m*R*R; 
    J.J(Z,Z) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 

    // assign the link1 parameters
    link1->body_id = "link1";                              // identify the link
    link1->set_inertia( J );                // set the inertia of the cylinder
    link1->set_enabled( true );                 // enable physics for the link
 
    // set the pose of the link
    // the height of the cylinder is 1.0, and the first joint will be located
    // at (0,0,0). In the zero configuration, the centroid of the cylinder will
    // therefore be located at (0,-0.5,0)
    Quatd ori(0,0,0,1);
    Origin3d position(0,-H/2.0,0);
    Pose3d pose( ori, position );
    link1->set_pose( pose ); 

    // add the link to the set of links
    links.push_back( link1 ); 
  }

  // the pendulum body is a cylinder
  shared_ptr<RigidBodyd> link2( new RigidBodyd() );
  {
    // setup cylinder radius and height
    const double R = 0.025;
    const double H = 1.0;

    // setup inertia for cylinder- it will be defined with respect to the
    // centroid of the cylinder
    SpatialRBInertiad J;
    J.pose = link2->get_pose();
    J.m = 1.0;

    // inertia matrix for cylinder is taken from online resources; note that
    // primary axis of cylinder is aligned with y-axis 
    J.J.set_zero(3,3);
    J.J(X,X) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 
    J.J(Y,Y) = 0.5*J.m*R*R; 
    J.J(Z,Z) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 

    // assign the link2 parameters
    link2->body_id = "link2";                              // identify the link
    link2->set_inertia( J );                // set the inertia of the cylinder
    link2->set_enabled( true );                 // enable physics for the link
 
    // set the pose of the link
    // the height of the cylinder is 1.0, and the second joint will be located
    // at (0,-1,0) at the zero configuration. In the zero configuration for
    // this joint, the centroid of the cylinder will therefore be located at 
    // (0,-1.5,0)
    Quatd ori(0,0,0,1);
    Origin3d position(0,-1.0-H/2.0,0);
    Pose3d pose( ori, position );
    link2->set_pose( pose ); 

    // add the link to the set of links
    links.push_back( link2 ); 
  }

  // create the first revolute joint
  boost::shared_ptr<RevoluteJointd> j1( new RevoluteJointd() );
  {
    // compute the position of the joint... center of the base w.r.t global frame
    Vector3d position(0.0, 0.0, 0.0, GLOBAL_3D);
    j1->set_location( position, base, link1 );

    // Note: set_location(...) must precede set_axis(...)
    // compute the axis of rotation
    Vector3d axis(0,0,1.0,GLOBAL_3D);   // actuates around global x-axis
    j1->set_axis( axis ); 

    // set the joint name 
    j1->joint_id = "j1";

    // add the j1 to the set of joints
    joints.push_back( j1 );
  }

  // create the second revolute joint
  boost::shared_ptr<RevoluteJointd> j2( new RevoluteJointd() );
  {
    // compute the position of the joint... center of the base w.r.t global frame
    Vector3d position(0.0, -1.0, 0.0, GLOBAL_3D);
    j2->set_location( position, link1, link2 );

    // compute the axis of rotation
    Vector3d axis(0,0,1.0,GLOBAL_3D);   // actuates around global x-axis
    j2->set_axis( axis ); 

    // set the joint name 
    j2->joint_id = "j1";

    // add j2 to the set of joints
    joints.push_back( j2 );
  }

  // construct the pendulum articulated body from the set of links and joints created above
  pendulum->set_links_and_joints( links, joints );

  // pendulum has a fixed base
  pendulum->set_floating_base(false);

  // set the initial configuration of the pendulum to give it potential energy
  VectorNd gc(2);
  gc[0] = M_PI_2;
  gc[1] = 0.0;
  pendulum->set_generalized_coordinates_euler(gc);
 
  // start the main simulation loop
  while(true) {
    // integrate the body forward
    integrate_euler(pendulum, DT); 

    // get and output the pendulum generalized coordinates
    pendulum->get_generalized_coordinates_euler(gc);
    std::cout << "t: " << t << " q: " << gc << std::endl;

    // update the time
    t += DT;
  }

  return 0;
}

