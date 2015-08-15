#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

void integrate(shared_ptr<RCArticulatedBodyd> body, double t, double dt)
{
  const shared_ptr<const Pose3d> GLOBAL_3D;
  VectorNd gc, gv, ga;

  // get the last link on the body
  shared_ptr<RigidBodyd> last_link = body->get_links().back();

  // output the generalized coordinates of the body
  body->get_generalized_coordinates(DynamicBodyd::eEuler, gc);
  std::cout << t;
  if (gc.size() > 0)
  {
    for (unsigned i=0; i< gc.size(); i++)
      std::cout << " " << gc[i];
    std::cout << std::endl;
  }
  else
  {
    Transform3d wP = Pose3d::calc_relative_pose(last_link->get_pose(), GLOBAL_3D);
    std::cout << " " << wP.x[0] << " " << wP.x[1] << " "<< wP.x[2] << " " << wP.q.x << " " << " " << wP.q.y << " " << " " << wP.q.z << " " << " " << wP.q.w << " " << std::endl;
  } 

  // clear the force accumulator on the body
  body->reset_accumulators();

  // add a force to the center of the body
  SForced f(last_link->get_inertial_pose());
  f[0] = std::cos(t);
  f[1] = std::sin(t);
  f[2] = std::cos(2*t);
  f[3] = std::sin(-3*t);
  f[4] = std::cos(3*t);
  f[5] = std::sin(5*t);
  last_link->add_force(f);

  // compute forward dynamics
  body->calc_fwd_dyn();

  // integrate the generalized velocity forward
  body->get_generalized_acceleration(ga);
  body->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  ga *= dt;
  gv += ga;
  body->set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // integrate the generalized position forward
  body->get_generalized_velocity(DynamicBodyd::eEuler, gv);
  body->get_generalized_coordinates(DynamicBodyd::eEuler, gc);
  gv *= dt;
  gc += gv;
  body->set_generalized_coordinates(DynamicBodyd::eSpatial, gc);
} 

int main(int argc, char* argv[])
{
  if (argc < 1)
  {
    std::cerr << "syntax: TestJoint <urdf file>" << std::endl;
    return -1;
  }

  // read in the body file
  std::string fname = std::string(argv[1]);
  std::string name = "body";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  for (unsigned i=0; i< 1000; i++)
    integrate(rcab, i/1000.0, 1.0/1000);
}

