#include <gtest/gtest.h>
#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/Log.h>
#include <Ravelin/Constants.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

const double DT = 1e-3;

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
} 

/// Sets the velocity for a body using a sequence
void set_velocity(shared_ptr<RCArticulatedBodyd> body)
{
  VectorNd gv;
  body->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  for (unsigned i=0; i< gv.size(); i++)
    gv[i] = (double) (i+1);
  body->set_generalized_velocity(DynamicBodyd::eSpatial, gv);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations 
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsLinkxLinkInertia)
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkInertia);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(DynamicsTest, DynamicsInertia)
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
  rcab->set_computation_frame_type(eLinkInertia);
  rcab->algorithm_type = RCArticulatedBodyd::eCRB;

  // set generalized velocity using sequence
  set_velocity(rcab);

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc1);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkInertia);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // calculate accelerations
  calc_dynamics(rcab, 0.0);

  // get values out
  rcab->get_generalized_acceleration(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
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

