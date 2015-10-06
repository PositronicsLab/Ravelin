/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// ------------------------------------------------------------------
// Example of a block falling under the influence of gravity 
// integration. Applies a gravitational force along the -y axis.
// Adapted from code by James R. Taylor.
// ------------------------------------------------------------------


#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Ravelin/RigidBodyd.h>
#include <Ravelin/ArticulatedBodyd.h>
#include "integrate.h"

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;

int main( void ) {

  const shared_ptr<const Pose3d> GLOBAL_3D;
  const unsigned X = 0, Y = 1, Z = 2;

  // integration step size
  const double DT = 1e-3;

  // current time
  double t = 0.0;

  // create the rigid body
  shared_ptr<RigidBodyd> rb( new RigidBodyd() );
  rb->body_id = "block";
  rb->set_enabled( true );
  rb->set_pose( Pose3d( Quatd(0,0,0,1), Origin3d(0,0,0) ) );

  // setup the inertia for the block; inertia is defined at the geometric
  // centroid of the block
  const double XSZ = 1.0, YSZ = 1.0, ZSZ = 1.0;  // box dimensions
  SpatialRBInertiad J;
  J.m = 100.0;
  J.J.set_zero();

  // inertia matrix definition comes from standard texts
  J.J(X,X) = (YSZ*YSZ + ZSZ*ZSZ); 
  J.J(Y,Y) = (XSZ*XSZ + ZSZ*ZSZ); 
  J.J(Z,Z) = (XSZ*XSZ + YSZ*YSZ); 
  J.J *= (2.0*J.m/12.0);

  // set the pose that the inertia is defined w.r.t. 
  J.pose = rb->get_pose();

  // set the inertia
  rb->set_inertia( J );

  while(true) {
    // integrate the equations of motion of the box forward in time
    integrate_euler( rb, DT );

    // get the pose of the box
    boost::shared_ptr<const Pose3d> pose = rb->get_pose();

    // copy the pose of the box, b/c we want to alter its relative pose
    Pose3d P = *pose;

    // update the relative pose to be in the global frame
    // (this is for pedagogical purposes- the relative pose of the box will 
    //  already be defined w.r.t. the global frame)
    P.update_relative_pose(GLOBAL_3D);
    std::cout << "t: " << t << " x: " << P.x << std::endl;

    // update current time
    t += DT;
  }

  return 0;
}

