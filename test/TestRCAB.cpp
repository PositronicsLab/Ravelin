#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

void integrate(shared_ptr<RCArticulatedBodyd> pendulum, double dt)
{
  // compute forward dynamics
  pendulum->calc_fwd_dyn();

  // integrate the generalized velocity forward
  VectorNd gc, gv, ga;
  pendulum->get_generalized_acceleration(ga);
  pendulum->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  ga *= dt;
  gv += ga;
  pendulum->set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // integrate the generalized position forward
  pendulum->get_generalized_velocity(DynamicBodyd::eEuler, gv);
  pendulum->get_generalized_coordinates(DynamicBodyd::eEuler, gc);
  gv *= dt;
  gc += gv;
  pendulum->set_generalized_coordinates(DynamicBodyd::eSpatial, gc);
} 

int main()
{
  // read in the pendulum file
  std::string fname = "pendulum.urdf";
  std::string name = "pendulum";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  for (unsigned i=0; i< 1000; i++)
    integrate(rcab, 1.0/1000);
}

