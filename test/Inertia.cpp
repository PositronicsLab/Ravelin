#include <Ravelin/Matrix3d.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/SpatialABInertiad.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
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

// verifies that inertia transformation does not alter calculations
TEST(InertiaTest, Transform2)
{
  // setup the pose for the inertia
  shared_ptr<Pose3d> P(new Pose3d);
  P->x = Origin3d(rand_double(), rand_double(), rand_double());
  P->q = Quatd(rand_double(), rand_double(), rand_double(), rand_double());
  P->q.normalize();

  // setup the inertia in the pose frame
  SpatialRBInertiad J(P);
  J.m = 2.0;
  J.h = Origin3d(rand_double(), rand_double(), rand_double());
  J.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // compute the inertia in the global frame
  SpatialRBInertiad J0 = Pose3d::transform(shared_ptr<Pose3d>(), J);

  // setup a random acceleration 
  SAcceld a0;
  a0.set_upper(Vector3d(rand_double(), rand_double(), rand_double()));
  a0.set_lower(Vector3d(rand_double(), rand_double(), rand_double()));
  
  // compute the spatial force 
  SForced f0 = J0 * a0;

  // transform a0 to P's frame
  SAcceld a = Pose3d::transform(P, a0);

  // compute the difference between f and f0
  SForced f0prime = Pose3d::transform(shared_ptr<Pose3d>(), J*a);
 
  for (unsigned i=0; i< f0prime.size(); i++)
    EXPECT_NEAR(f0[i], f0prime[i], 1e-6);
}

// verifies that added inertias are equivalent
TEST(InertiaTest, AddInertias)
{
  // setup the first inertia
  SpatialRBInertiad J1;
  J1.m = (double) rand() / RAND_MAX;
  J1.h = Origin3d(rand_double(), rand_double(), rand_double());
  J1.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // setup the second inertia
  SpatialRBInertiad J2; 
  J2.m = (double) rand() / RAND_MAX;
  J2.h = Origin3d(rand_double(), rand_double(), rand_double());
  J2.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // add the two inertias together
  SpatialRBInertiad J = J1 + J2;
  
  // make a matrices out of J1, J2, and J
  MatrixNd J1m(6,6), J2m(6,6), Jm(6,6), Jam(6,6);
  J1.to_matrix(J1m);
  J2.to_matrix(J2m);
  J.to_matrix(Jm);

  // add J1m and J2m
  J1m += J2m;

  // make a spatial articulated inertia matrix out of J
  SpatialABInertiad Ja = J;
  Ja.to_matrix(Jam);

  // make spatial articulated inertia matrices out of J1 and J2
  MatrixNd Ja1m(6,6);
  SpatialABInertiad Ja1 = J1;
  SpatialABInertiad Ja2 = J2;
  Ja1 += Ja2;
  Ja1.to_matrix(Ja1m);

  // all should be equivalent
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), J1m(i,j), 1e-6);

  // should be equivalent to spatial ab inertia too
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), Jam(i,j), 1e-6);

  // should be equivalent to spatial ab inertia too
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), Ja1m(i,j), 1e-6);
}

// verifies that subtracted inertias are equivalent
TEST(InertiaTest, SubInertias)
{
  // setup the first inertia
  SpatialRBInertiad J1;
  J1.m = (double) rand() / RAND_MAX;
  J1.h = Origin3d(rand_double(), rand_double(), rand_double());
  J1.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // setup the second inertia
  SpatialRBInertiad J2; 
  J2.m = (double) rand() / RAND_MAX;
  J2.h = Origin3d(rand_double(), rand_double(), rand_double());
  J2.J = Matrix3d(1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0);

  // subtract the two inertias 
  SpatialRBInertiad J = J1 - J2;
  
  // make matrices out of J1, J2, and J
  MatrixNd J1m(6,6), J2m(6,6), Jm(6,6), Jam(6,6);
  J1.to_matrix(J1m);
  J2.to_matrix(J2m);
  J.to_matrix(Jm);

  // subtract J1m and J2m
  J1m -= J2m;

  // make a spatial articulated inertia matrix out of J
  SpatialABInertiad Ja = J;
  Ja.to_matrix(Jam);

  // make spatial articulated inertia matrices out of J1 and J2
  MatrixNd Ja1m(6,6);
  SpatialABInertiad Ja1 = J1;
  SpatialABInertiad Ja2 = J2;
  Ja1 -= Ja2;
  Ja1.to_matrix(Ja1m);

  // all should be equivalent
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), J1m(i,j), 1e-6);

  // should be equivalent to spatial ab inertia too
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), Jam(i,j), 1e-6);

  // should be equivalent to spatial ab inertia too
  for (unsigned i=0; i< Jm.rows(); i++)
    for (unsigned j=0; j< Jm.columns(); j++)
      EXPECT_NEAR(Jm(i,j), Ja1m(i,j), 1e-6);
}

