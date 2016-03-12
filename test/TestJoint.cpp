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

void calc_dynamics(shared_ptr<RCArticulatedBodyd> body, double t)
{
  const shared_ptr<const Pose3d> GLOBAL_3D;
  VectorNd ga;

  // get all links 
  const vector<shared_ptr<RigidBodyd> >& links = body->get_links();

  // clear the force accumulator on the body
  body->reset_accumulators();

  // add forces to the center of each link
  for (unsigned i=0; i< links.size(); i++)
  { 
    SForced f(links[i]->get_pose());
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

  // output the generalized acceleration
  body->get_generalized_acceleration(ga);
  std::cout << t;
  if (ga.size() > 0)
  {
    for (unsigned i=0; i< ga.size(); i++)
      std::cout << " " << ga[i];
    std::cout << std::endl;
  }
} 

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "syntax: TestJoint <urdf file>" << std::endl;
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

  // set velocity
  set_velocity(rcab);

  // calculate dynamics
  calc_dynamics(rcab, 0.0);
  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  OutputToFile::stream.close();
}

