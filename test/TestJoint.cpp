#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/Log.h>
#include <Ravelin/Constants.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

/// Sets the velocity for a body using a sequence
void set_velocity(shared_ptr<RCArticulatedBodyd> body)
{
  VectorNd gv;
  body->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  for (unsigned i=0; i< gv.size(); i++)
    gv[i] = (double) (i+1);
  body->set_generalized_velocity(DynamicBodyd::eSpatial, gv);
}

void integrate(shared_ptr<RCArticulatedBodyd> body, double t, double dt)
{
  const shared_ptr<const Pose3d> GLOBAL_3D;
  VectorNd gc, gv, ga;

  // get all links 
  const vector<shared_ptr<RigidBodyd> >& links = body->get_links();

  // output the generalized coordinates of the body
  body->get_generalized_coordinates_euler(gc);
  std::cout << t;
  if (gc.size() > 0)
  {
    for (unsigned i=0; i< gc.size(); i++)
      std::cout << " " << gc[i];
  }
  else
  {
//    Transform3d wP = Pose3d::calc_relative_pose(last_link->get_pose(), GLOBAL_3D);
//    std::cout << " " << wP.x[0] << " " << wP.x[1] << " "<< wP.x[2] << " " << wP.q.x << " " << " " << wP.q.y << " " << " " << wP.q.z << " " << " " << wP.q.w << " " << std::endl;
  } 

  // clear the force accumulator on the body
  body->reset_accumulators();

  // add forces to the center of each link
  for (unsigned i=0; i< links.size(); i++)
  { 
    SForced f(links[i]->get_inertial_pose());
    f[0] = std::cos(t*(i*1.01));
    f[1] = std::sin(t*(i*1.01));
    f[2] = std::cos(2*t*(i*1.01));
    f[3] = std::sin(-2*t*(i*1.01));
    f[4] = std::cos(3*t*(i*1.01));
    f[5] = std::sin(-3*t*(i*1.01));
    links[i]->add_force(f);
  }

  // compute forward dynamics
  body->calc_fwd_dyn();

  // integrate the generalized velocity forward
  body->get_generalized_acceleration(ga);
  body->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  ga *= dt;
  gv += ga;
  body->set_generalized_velocity(DynamicBodyd::eSpatial, gv);

  if (gc.size() > 0)
  {
    for (unsigned i=0; i< ga.size(); i++)
      std::cout << " " << ga[i];
    std::cout << std::endl;
  }

  // integrate the generalized position forward
  body->get_generalized_velocity(DynamicBodyd::eEuler, gv);
  body->get_generalized_coordinates_euler(gc);
  gv *= dt;
  gc += gv;
  body->set_generalized_coordinates_euler(gc);
} 

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr << "syntax: TestJoint <urdf file> <steps>" << std::endl;
    return -1;
  }

  // log dynamics
  Log<OutputToFile>::reporting_level = LOG_DYNAMICS;
  OutputToFile::stream.open("dynamics.log");

  // read in the body file
  std::string fname = std::string(argv[1]);
  std::string name = "body";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;
set_velocity(rcab);
  // get the number of steps
  const unsigned STEPS = std::atoi(argv[2]);

  for (unsigned i=0; i< STEPS; i++)
    integrate(rcab, i/1000.0, 1.0/1000);

  OutputToFile::stream.close();
}

