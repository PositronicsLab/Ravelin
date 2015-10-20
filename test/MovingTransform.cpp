#include <Ravelin/Matrix3d.h>
#include <Ravelin/MovingTransform3d.h>
#include <Ravelin/Pose3d.h>
#include "gtest/gtest.h"

using boost::shared_ptr;
using namespace Ravelin;

TEST(MovingFrameTest, SourceToGlobal)
{
  SVelocityd dummy, vs, vt;
  shared_ptr<Pose3d> GLOBAL, Ts, Tt;

  // setup source pose
  Ts = shared_ptr<Pose3d>(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  Tt = shared_ptr<Pose3d>(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // setup velocities
  vs = SVelocityd(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  vt = SVelocityd(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);

  // do the conversion, testing source to global
  MovingTransform3d M1 = MovingTransform3d::calc_transform(Ts, GLOBAL, vs, dummy);

  EXPECT_NEAR(M1.Edot(0,0), +1.2566, 5e-3);
  EXPECT_NEAR(M1.Edot(0,1), +0.1688, 5e-3);
  EXPECT_NEAR(M1.Edot(0,2), -0.1205, 5e-3);
  EXPECT_NEAR(M1.Edot(1,0), -0.8512, 5e-3);
  EXPECT_NEAR(M1.Edot(1,1), -0.6260, 5e-3);
  EXPECT_NEAR(M1.Edot(1,2), -0.3564, 5e-3);
  EXPECT_NEAR(M1.Edot(2,0), -0.4612, 5e-3);
  EXPECT_NEAR(M1.Edot(2,1), +1.0470, 5e-3);
  EXPECT_NEAR(M1.Edot(2,2), +0.9936, 5e-3);
  EXPECT_NEAR(M1.rdot[0], -1.0596, 5e-3);
  EXPECT_NEAR(M1.rdot[1], +3.1335, 5e-3);
  EXPECT_NEAR(M1.rdot[2], +2.4343, 5e-3);
}

TEST(MovingFrameTest, GlobalToTarget)
{
  SVelocityd dummy, vs, vt;
  shared_ptr<Pose3d> GLOBAL, Ts, Tt;

  // setup source pose
  Ts = shared_ptr<Pose3d>(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  Tt = shared_ptr<Pose3d>(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // setup velocities
  vs = SVelocityd(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  vt = SVelocityd(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);

  // do the conversion, testing global to target
  MovingTransform3d M2 = MovingTransform3d::calc_transform(GLOBAL, Tt, dummy, vt);

  EXPECT_NEAR(M2.Edot(0,0), +0.3847, 5e-3);
  EXPECT_NEAR(M2.Edot(0,1), -0.6955, 5e-3);
  EXPECT_NEAR(M2.Edot(0,2), +0.0206, 5e-3);
  EXPECT_NEAR(M2.Edot(1,0), -0.5566, 5e-3);
  EXPECT_NEAR(M2.Edot(1,1), -0.2475, 5e-3);
  EXPECT_NEAR(M2.Edot(1,2), +0.0248, 5e-3);
  EXPECT_NEAR(M2.Edot(2,0), -0.4248, 5e-3);
  EXPECT_NEAR(M2.Edot(2,1), -0.3040, 5e-3);
  EXPECT_NEAR(M2.Edot(2,2), +0.0239, 5e-3);
  EXPECT_NEAR(M2.rdot[0], +1.0187, 5e-3);
  EXPECT_NEAR(M2.rdot[1], -0.1332, 5e-3);
  EXPECT_NEAR(M2.rdot[2], -0.7145, 5e-3);
}

TEST(MovingFrameTest, SourceToTargetRelToGlobal)
{
  SVelocityd dummy, vs, vt;
  shared_ptr<Pose3d> GLOBAL, Ts, Tt;

  // setup source pose
  Ts = shared_ptr<Pose3d>(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  Tt = shared_ptr<Pose3d>(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // setup velocities
  vs = SVelocityd(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  vt = SVelocityd(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);

  // do the conversion, testing source to target
  MovingTransform3d M3 = MovingTransform3d::calc_transform(Ts, Tt, vs, vt);

  EXPECT_NEAR(M3.Edot(0,0), 0.3664, 5e-3);
  EXPECT_NEAR(M3.Edot(0,1), 0.2136, 5e-3);
  EXPECT_NEAR(M3.Edot(0,2), 0.4513, 5e-3);
  EXPECT_NEAR(M3.Edot(1,0), 0.6395, 5e-3);
  EXPECT_NEAR(M3.Edot(1,1), 0.5950, 5e-3);
  EXPECT_NEAR(M3.Edot(1,2), 1.0534, 5e-3);
  EXPECT_NEAR(M3.Edot(2,0), 1.2118, 5e-3);
  EXPECT_NEAR(M3.Edot(2,1), -0.8153, 5e-3);
  EXPECT_NEAR(M3.Edot(2,2), -0.3272, 5e-3);
  EXPECT_NEAR(M3.rdot[0], -3.3583, 5e-3);
  EXPECT_NEAR(M3.rdot[1], +3.6356, 5e-3);
  EXPECT_NEAR(M3.rdot[2], +2.9654, 5e-3);
}

TEST(MovingFrameTest, SourceToTargetRelR)
{
  SVelocityd dummy, vs, vt;
  shared_ptr<Pose3d> GLOBAL, Ts, Tt;

  // setup pose r
  shared_ptr<Pose3d> r(new Pose3d);
  r->rpose = GLOBAL;
  r->x = Origin3d(-0.2938, -0.8479, -1.1201);
  r->q = Matrix3d(0.9062, 0.1872, -0.3791, -0.1507, -0.6946, -0.7034, -0.3950,  0.6946, -0.6013);

  // setup source pose
  Ts = shared_ptr<Pose3d>(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  Tt = shared_ptr<Pose3d>(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // update poses to point to r
  Ts->update_relative_pose(r);
  Tt->update_relative_pose(r);

  // setup velocities
  vs = SVelocityd(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  vt = SVelocityd(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);

  // transform velocities to frame r
  SVelocityd vsr = Pose3d::transform(r, vs);
  SVelocityd vtr = Pose3d::transform(r, vt);

  // do the conversion, testing source to target
  MovingTransform3d M4 = MovingTransform3d::calc_transform(Ts, Tt, vsr, vtr);

  EXPECT_NEAR(M4.Edot(0,0), 0.3664, 5e-3);
  EXPECT_NEAR(M4.Edot(0,1), 0.2136, 5e-3);
  EXPECT_NEAR(M4.Edot(0,2), 0.4513, 5e-3);
  EXPECT_NEAR(M4.Edot(1,0), 0.6395, 5e-3);
  EXPECT_NEAR(M4.Edot(1,1), 0.5950, 5e-3);
  EXPECT_NEAR(M4.Edot(1,2), 1.0534, 5e-3);
  EXPECT_NEAR(M4.Edot(2,0), 1.2118, 5e-3);
  EXPECT_NEAR(M4.Edot(2,1), -0.8153, 5e-3);
  EXPECT_NEAR(M4.Edot(2,2), -0.3272, 5e-3);
  EXPECT_NEAR(M4.rdot[0], -3.3583, 5e-3);
  EXPECT_NEAR(M4.rdot[1], +3.6356, 5e-3);
  EXPECT_NEAR(M4.rdot[2], +2.9654, 5e-3);
}

TEST(MovingFrameTest, SourceToTargetRelMAndN)
{
  SVelocityd dummy, vs, vt;
  shared_ptr<Pose3d> GLOBAL, Ts, Tt;

  // setup pose m
  shared_ptr<Pose3d> m(new Pose3d);
  m->rpose = GLOBAL;
  m->x = Origin3d(-.1952, -.2176, -.3031);
  m->q = Matrix3d(0.8849, 0.378, 0.2722, 0.4634, -0.6546, -0.5973, -0.0476, 0.6546, -0.7544);

  // setup pose n
  shared_ptr<Pose3d> n(new Pose3d);
  n->rpose = GLOBAL;
  n->x = Origin3d(-0.2938, -0.8479, -1.1201);
  n->q = Matrix3d(0.9062, 0.1872, -0.3791, -0.1507, -0.6946, -0.7034, -0.3950,  0.6946, -0.6013);

  // setup source pose
  Ts = shared_ptr<Pose3d>(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  Tt = shared_ptr<Pose3d>(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // update poses to point to m and n
  Ts->update_relative_pose(m);
  Tt->update_relative_pose(n);

  // setup velocities
  vs = SVelocityd(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  vt = SVelocityd(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);

  // transform velocities to frames m and n
  SVelocityd vsm = Pose3d::transform(m, vs);
  SVelocityd vtn = Pose3d::transform(n, vt);

  // do the conversion, testing source to target
  MovingTransform3d M5 = MovingTransform3d::calc_transform(Ts, Tt, vsm, vtn);

  EXPECT_NEAR(M5.Edot(0,0), 0.3664, 5e-3);
  EXPECT_NEAR(M5.Edot(0,1), 0.2136, 5e-3);
  EXPECT_NEAR(M5.Edot(0,2), 0.4513, 5e-3);
  EXPECT_NEAR(M5.Edot(1,0), 0.6395, 5e-3);
  EXPECT_NEAR(M5.Edot(1,1), 0.5950, 5e-3);
  EXPECT_NEAR(M5.Edot(1,2), 1.0534, 5e-3);
  EXPECT_NEAR(M5.Edot(2,0), 1.2118, 5e-3);
  EXPECT_NEAR(M5.Edot(2,1), -0.8153, 5e-3);
  EXPECT_NEAR(M5.Edot(2,2), -0.3272, 5e-3);
  EXPECT_NEAR(M5.rdot[0], -3.3583, 5e-3);
  EXPECT_NEAR(M5.rdot[1], +3.6356, 5e-3);
  EXPECT_NEAR(M5.rdot[2], +2.9654, 5e-3);
}

TEST(MovingFrameTest, AccelConversion)
{
  shared_ptr<Pose3d> GLOBAL;
  const double DT = 1e-4;

  // setup a rigid body pose 
  shared_ptr<Pose3d> body_frame(new Pose3d); 
  body_frame->rpose = GLOBAL;
  body_frame->x = Origin3d(-.1952, -.2176, -.3031);
  body_frame->q = Matrix3d(0.8849, 0.378, 0.2722, 0.4634, -0.6546, -0.5973, -0.0476, 0.6546, -0.7544);

  // setup a frame defined in the global system 
  shared_ptr<Pose3d> pose(new Pose3d); 
  pose->rpose = GLOBAL;
  pose->x = Origin3d(1, 2, 3);

  // set the velocity of the body
  SVelocityd body_v(1, 2, 3, 4, 5, 6, body_frame);
  SVelocityd old_body_v = body_v;
  SVelocityd pose_v(5.0, 7.0, 9.0, 1.0, 2.0, 3.0, GLOBAL);; 

  // set the acceleration for the body
  SAcceld body_a(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, body_frame);

  // do the conversion, testing source to target
  MovingTransform3d M5 = MovingTransform3d::calc_transform(body_frame, pose, body_v, pose_v);

  // now compute the transformation of the acceleration to pose
  SAcceld pose_a = M5.transform(body_a);
std::cerr << "body_v: " << body_v << std::endl;
std::cerr << "body_a: " << body_a << std::endl;
std::cerr << "r: " << M5.r << std::endl;
std::cerr << "E: " << std::endl << M5.E;
std::cerr << "\\dot{r}: " << M5.rdot << std::endl;
std::cerr << "\\dot{E}: " << std::endl << M5.Edot;

  // *** now numerically integrate body_frame forward by dt
  // update the body frame origin
  body_frame->x += Origin3d(body_v.get_linear())*DT;
  Quatd qd = Quatd::deriv(body_frame->q, body_v.get_angular());
  body_frame->q += qd*DT;
  body_frame->q.normalize();

  // update the velocity
  for (unsigned i=0; i< body_v.size(); i++)
    body_v[i] += body_a[i]*DT;

  // given the new poses and new velocities, compute acceleration numerically
  SVelocityd vdiff = Pose3d::transform(pose, body_v) - Pose3d::transform(pose, old_body_v);
  vdiff /= DT;

  // test numerically
  for (unsigned i=0; i< vdiff.size(); i++)
    EXPECT_NEAR(vdiff[i], pose_a[i], 1e-3);
}


