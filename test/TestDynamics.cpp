#include <gtest/gtest.h>
#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/Log.h>
#include <Ravelin/Constants.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

const double DT = 1e-3;
const unsigned INV_DT = (unsigned) (1.0/DT);

void integrate(shared_ptr<RCArticulatedBodyd> body, double t, double dt)
{
  const shared_ptr<const Pose3d> GLOBAL_3D;
  VectorNd gc, gv, ga;

  // get all links 
  const vector<shared_ptr<RigidBodyd> >& links = body->get_links();

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

  // integrate the generalized position forward
  body->get_generalized_velocity(DynamicBodyd::eEuler, gv);
  body->get_generalized_coordinates_euler(gc);
  gv *= dt;
  gc += gv;
  body->set_generalized_coordinates_euler(gc);
} 

class DynamicsTest : public ::testing::Test {
  public:
  static const char* filename;
};

const char* DynamicsTest::filename;

TEST_F(DynamicsTest, DynamicsLinkxLinkCOM)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
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

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsLinkxGlobal)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
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

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsLinkxJoint)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
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

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsLink)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
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

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsGlobal)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
  std::string name = "body";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsLinkCOM)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
  std::string name = "body";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsJoint)
{
  VectorNd gc1, gc2;

  // read in the body file
  std::string fname(filename);
  std::string name = "body";
  vector<shared_ptr<RigidBodyd> > links;
  vector<shared_ptr<Jointd> > joints;
  URDFReaderd::read(fname, name, links, joints);

  // create a new articulated body
  shared_ptr<RCArticulatedBodyd> rcab(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // prepare to do everything again
  links.clear();
  joints.clear();
  URDFReaderd::read(fname, name, links, joints);
  rcab = shared_ptr<RCArticulatedBodyd>(new RCArticulatedBodyd);
  rcab->set_links_and_joints(links, joints); 

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one second at a step size of 0.001
  for (unsigned i=0; i< INV_DT; i++)
    integrate(rcab, i*DT, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  gc1 -= gc2;
  ASSERT_NEAR(gc1.norm(), 0.0, EPS_DOUBLE);
}

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "syntax: TestDynamics <urdf file>" << std::endl;
    return -1;
  }

  // set the filename
  DynamicsTest::filename = argv[1];

  // run Google tests
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

