#include <Ravelin/SForced.h>
#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/Log.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

void integrate(shared_ptr<RCArticulatedBodyd> pendulum, double dt)
{
  // clear all forces on the pendulum
  pendulum->reset_accumulators();

  // add gravity
  const double G = 9.8;
  shared_ptr<RigidBodyd> link1 = pendulum->get_links()[1];
  shared_ptr<RigidBodyd> link2 = pendulum->get_links()[2];
  const double mg1 = G*link1->get_inertia().m;
  const double mg2 = G*link2->get_inertia().m;
  SForced f1(0.0, -mg1, 0.0, 0.0, 0.0, 0.0, link1->get_mixed_pose());
  SForced f2(0.0, -mg2, 0.0, 0.0, 0.0, 0.0, link2->get_mixed_pose());
  link1->add_force(f1); 
  link2->add_force(f2); 

  // compute forward dynamics
  pendulum->calc_fwd_dyn();

  // integrate the generalized velocity forward
  VectorNd gc, gv, gve, ga;
  pendulum->get_generalized_acceleration(ga);
  pendulum->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  pendulum->get_generalized_velocity(DynamicBodyd::eSpatial, gve);
  ga *= dt;
  gv += ga;

  // integrate the generalized position forward
  pendulum->get_generalized_coordinates_euler(gc);
  static double t = 0.0;
  std::cout << t << " " << gc << std::endl;
  gve *= dt;
  gc += gve;
  pendulum->set_generalized_coordinates_euler(gc);
  pendulum->set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  // update t
  t += dt;
} 

int main()
{
  // read in the pendulum file
  std::string fname = "pendulum.urdf";
  std::string name = "pendulum";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // setup dynamics
  Log<OutputToFile>::reporting_level = LOG_DYNAMICS;

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_floating_base(false);
  rcab->set_links_and_joints(links, joints); 
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  for (unsigned i=0; i< 1000; i++)
    integrate(rcab, 1.0/1000);

}

