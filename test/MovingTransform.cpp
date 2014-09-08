#include <Ravelin/Matrix3d.h>
#include <Ravelin/MovingTransform3d.h>
#include <Ravelin/Pose3d.h>

using boost::shared_ptr;
using namespace Ravelin;

int main()
{
  shared_ptr<Pose3d> GLOBAL;

  // setup source pose
  shared_ptr<Pose3d> Ts(new Pose3d);
  Ts->rpose = GLOBAL;
  Ts->x = Origin3d(-0.3538, -0.8236, -1.5771);
  Ts->q = Matrix3d(-0.1271, 0.9900, 0.0617, 0.3209, 0.0999, -0.9418, -0.9385,   -0.0999, -0.330);

  shared_ptr<Pose3d> Tt(new Pose3d);
  Tt->rpose = GLOBAL;
  Tt->x = Origin3d(-1.3337, 1.1275, 0.3502);
  Tt->q = Matrix3d(0.8729, 0.3269, 0.3622, 0.4845, -0.6683, -0.5645, 0.0575, 0.6683, -0.7417);

  // setup velocities
  SVelocityd vs(-0.9792, -1.1564, -0.5336, -2.0026, 0.9642, 0.5201, GLOBAL);
  SVelocityd vt(-0.0200, -0.0348, -0.7982,  1.0187, -0.1332, -0.7145, GLOBAL);
  SVelocityd dummy;

  // do the conversion, testing source to global
  MovingTransform3d M1 = MovingTransform3d::calc_transform(Ts, GLOBAL, vs, dummy);

  // check Edot and rdot
  std::cout << "Edot: " << std::endl << M1.Edot;
  std::cout << "rdot: " << M1.rdot<< std::endl;

  // do the conversion, testing global to target
  MovingTransform3d M2 = MovingTransform3d::calc_transform(GLOBAL, Tt, dummy, vt);

  // answers should be:

  // check Edot and rdot
  std::cout << "Edot: " << std::endl << M2.Edot;
  std::cout << "rdot: " << M2.rdot<< std::endl;

  // do the conversion, testing source to target
  MovingTransform3d M3 = MovingTransform3d::calc_transform(Ts, Tt, vs, vt);

  // check Edot and rdot
  std::cout << "Edot: " << std::endl << M3.Edot;
  std::cout << "rdot: " << M3.rdot<< std::endl;
}

