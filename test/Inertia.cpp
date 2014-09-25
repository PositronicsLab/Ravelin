#include <Ravelin/Matrix3d.h>
#include <Ravelin/MovingTransform3d.h>
#include <Ravelin/Pose3d.h>
#include "gtest/gtest.h"

using boost::shared_ptr;
using namespace Ravelin;

double rand_double()
{
  return (double) rand() / RAND_MAX * 2.0 - 1.0;
}

// verifies that J*a = f and inv(J)*f = a
TEST(InertiaTest, BackwardForward)
{
  // setup a spatial rigid body inertia matrix
  SpatialRBInertiad J;
  J.m = 1.0;
  J.h = Origin3d(rand_double(), rand_double(), rand_double());
  J.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // setup a random spatial force
  SForced f;
  f.set_upper(Vector3d(rand_double(), rand_double(), rand_double()));
  f.set_lower(Vector3d(rand_double(), rand_double(), rand_double()));

  // compute the spatial acceleration
  SAcceld a = J.inverse_mult(f);

  // make sure forward and inverse match
  SForced fprime = J * a;
  for (unsigned i=0; i< fprime.size(); i++)
    EXPECT_NEAR(f[i], fprime[i], 1e-6);
}

// verifies that
TEST(InertiaTest, Transform)
{
  // setup the pose for the inertia
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = Origin3d(rand_double(), rand_double(), rand_double());
  P->q = Quatd(rand_double(), rand_double(), rand_double(), rand_double());
  P->q.normalize();

  // setup the inertia in the pose frame
  SpatialRBInertiad J(P);
  J.m = 1.0;
  J.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // compute the inertia in the global frame
  SpatialRBInertiad J0 = Pose3d::transform(shared_ptr<Pose3d>(), J);

  // setup a random force 
  SForced f0;
  f0.set_upper(Vector3d(rand_double(), rand_double(), rand_double()));
  f0.set_lower(Vector3d(rand_double(), rand_double(), rand_double()));
  
  // compute the spatial acceleration in the global frame
  SAcceld a0 = J0.inverse_mult(f0);

  // make sure forward and inverse match
  SForced fprime0 = J0 * a0;
  for (unsigned i=0; i< fprime0.size(); i++)
    EXPECT_NEAR(f0[i], fprime0[i], 1e-6);

  // now compute everything in the other frame
  SForced f = Pose3d::transform(P, f0);

  // compute the spatial acceleration
  SAcceld a = J.inverse_mult(f);

  // convert a to global frame
  SAcceld aprime0 = Pose3d::transform(shared_ptr<Pose3d>(), a);

  // compare a and a0prime
  for (unsigned i=0; i< aprime0.size(); i++)
    EXPECT_NEAR(a0[i], aprime0[i], 1e-6);
}

