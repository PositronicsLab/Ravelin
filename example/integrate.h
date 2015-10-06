/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// -------------------------------------------------------------
// Function for integrating a body forward using explicit Euler
// integration. Applies a gravitational force along the -y axis.
// -------------------------------------------------------------

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <Ravelin/RigidBodyd.h>
#include <Ravelin/ArticulatedBodyd.h>
#include <Ravelin/VectorNd.h>

// integrates a dynamic body forward by applying gravitational forces along the 
void integrate_euler(boost::shared_ptr<Ravelin::DynamicBodyd> body, double dt)
{
  const double G = 9.8;  // acceleration due to gravity
  Ravelin::VectorNd gc, gv, gve, ga;

  // clear all forces on the body 
  body->reset_accumulators();

  // attempt to get the body as an articulated body
  boost::shared_ptr<Ravelin::ArticulatedBodyd> ab = boost::dynamic_pointer_cast<Ravelin::ArticulatedBodyd>(body);
  if (ab)
  {
    // add gravitational forces to each link
    for (unsigned i=0; i< ab->get_links().size(); i++)
    {
      boost::shared_ptr<Ravelin::RigidBodyd> link = ab->get_links()[i];
      const double mg = G*link->get_inertia().m;

      // the force is applied at the center of mass and along the y-axis,
      // so we use the mixed pose for this (the mixed pose is defined at the
      // center of mass of the link and is aligned with the global frame)
      Ravelin::SForced f(0.0, -mg, 0.0, 0.0, 0.0, 0.0, link->get_mixed_pose());

      // add the force
      link->add_force(f);
    }
  }
  else
  {
    boost::shared_ptr<Ravelin::RigidBodyd> rb = boost::dynamic_pointer_cast<Ravelin::RigidBodyd>(body);
    assert(rb);

    // setup gravitational force 
    const double mg = G*rb->get_inertia().m;

    // the force is applied at the center of mass and along the y-axis,
    // so we use the mixed pose for this (the mixed pose is defined at the
    // center of mass of the body and is aligned with the global frame)
    Ravelin::SForced f(0.0, -mg, 0.0, 0.0, 0.0, 0.0, rb->get_mixed_pose());
    rb->add_force(f);
  }

  // compute forward dynamics (calculate accelerations)
  body->calc_fwd_dyn();

  // get the current velocity using Euler parameters (if relevant)
  body->get_generalized_velocity(Ravelin::DynamicBodyd::eEuler, gve);

  // integrate the generalized velocity (using spatial parameters) forward
  body->get_generalized_acceleration(ga);
  body->get_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);
  ga *= dt;
  gv += ga;

  // integrate the generalized position forward
  body->get_generalized_coordinates_euler(gc);
  gve *= dt;
  gc += gve;
  body->set_generalized_coordinates_euler(gc);

  // *now* set the generalized velocity- velocity may be defined using 
  // the body coordinates, so we want to set it last
  body->set_generalized_velocity(Ravelin::DynamicBodyd::eSpatial, gv);
} 


#endif

