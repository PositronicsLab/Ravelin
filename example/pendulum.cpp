/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// ------------------------------------------------------------------
// Example of a pendulum, adapted from code by James R. Taylor 
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

  // the pendulum has a joint, so we need an articulated body 
  shared_ptr<RCArticulatedBodyd> pendulum( new RCArticulatedBodyd() );
  pendulum->body_id = "pendulum";

  // set the algorithm to the Composite Rigid Body Method
  pendulum->algorithm_type = RCArticulatedBodyd::eCRB;

  // set the dynamics computations to take place in link inertia frames
  // (these are frame defined at the center-of-mass of each link and changing
  //  orientation as the link changes orientation)
  pendulum->set_computation_frame_type(eLinkInertia);

  // vectors of links and joints used to define the body
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
  shared_ptr<RigidBodyd> link( new RigidBodyd() );
  {
    // setup cylinder radius and height
    const double R = 0.025;
    const double H = 1.0;

    // setup inertia for cylinder- it will be defined with respect to the
    // centroid of the cylinder
    SpatialRBInertiad J;
    J.pose = link->get_pose();
    J.m = 1.0;

    // inertia matrix for cylinder is taken from online resources; note that
    // primary axis of cylinder is aligned with y-axis 
    J.J.set_zero(3,3);
    J.J(X,X) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 
    J.J(Y,Y) = 0.5*J.m*R*R; 
    J.J(Z,Z) = 1.0/12*J.m*H*H + 0.25*J.m*R*R; 

    // assign the link parameters
    link->body_id = "link";                              // identify the link
    link->set_inertia( J );                // set the inertia of the cylinder
    link->set_enabled( true );                 // enable physics for the link
 
    // set the pose of the link 
    Quatd ori(0,0,0,1);
    Origin3d position(0,-0.5,0);
    Pose3d pose( ori, position );
    link->set_pose( pose ); 

    // add the link to the set of links
    links.push_back( link ); 
  }

  // create the revolute joint
  boost::shared_ptr<RevoluteJointd> revolute( new RevoluteJointd() );
  {
    // set the name of the joint
    revolute->joint_id = "q0";

    // set the position of the joint... center of the base w.r.t global frame
    Vector3d position(0.0, 0.0, 0.0, GLOBAL_3D);
    revolute->set_location( position, base, link );

    // set the axis of ori
    Vector3d axis(0,0,1.0,GLOBAL_3D);   // actuates around global z-axis

    // assign the revolute parameters
    // Note: set_location(...) must precede set_axis(...)
    revolute->set_axis( axis ); 

    // add the revolute to the set of joints
    joints.push_back( revolute );
  }

  // construct the pendulum articulated body from the sets of links and joints 
  pendulum->set_links_and_joints( links, joints );

  // pendulum has a fixed base
  pendulum->set_floating_base(false);

  // set the pendulum initial configuration so that it has potential energy
  VectorNd gc(1);
  gc[0] = M_PI_2;
  pendulum->set_generalized_coordinates_euler(gc);
 
  // start the main simulation loop
  while(true) {
    // integrate the body forward
    integrate_euler(pendulum, DT); 

    // get the pose of the link 
    boost::shared_ptr<const Pose3d> pose = link->get_pose();

    // copy the pose, b/c we want to alter its relative pose
    Pose3d P = *pose;

    // update the relative pose to be in the global frame
    P.update_relative_pose(GLOBAL_3D);

    // print the pose to the console 
    std::cout << "t: " << t << ", pose: " << P << std::endl;

    // update the time
    t += DT;
  }

  return 0;
}

