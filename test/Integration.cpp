#include <gtest/gtest.h>
#include <Ravelin/URDFReaderd.h>
#include <Ravelin/RCArticulatedBodyd.h>
#include <Ravelin/Log.h>
#include <Ravelin/Constants.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;

const double DT = 1e-3;

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

/// Sets the velocity for a body using a sequence
void set_velocity(shared_ptr<RCArticulatedBodyd> body)
{
  VectorNd gv;
  body->get_generalized_velocity(DynamicBodyd::eSpatial, gv);
  for (unsigned i=0; i< gv.size(); i++)
    gv[i] = (double) (i+1);
  body->set_generalized_velocity(DynamicBodyd::eSpatial, gv);
}

class IntegrationTest : public ::testing::Test {
  public:
  static const char* filename;
};

const char* IntegrationTest::filename;

TEST_F(IntegrationTest, DynamicsLinkxGlobal)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsLinkxJoint)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsLinkxLinkInertia)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsLink)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsGlobal)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eGlobal);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsLinkCOM)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLinkCOM);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsJoint)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eJoint);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

TEST_F(IntegrationTest, DynamicsInertia)
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

  // save the generalized coordinates
  VectorNd gc;
  rcab->get_generalized_coordinates_euler(gc);

  // set generalized velocity using sequence
  set_velocity(rcab);

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc1);

  // restore the generalized coordinates and velocity
  rcab->set_generalized_coordinates_euler(gc);
  set_velocity(rcab);

  // set dynamics algorithm and frame
  rcab->set_computation_frame_type(eLink);
  rcab->algorithm_type = RCArticulatedBodyd::eFeatherstone;

  // integrate for one step at a step size of 0.001
  integrate(rcab, 0.0, DT);

  // get values out
  rcab->get_generalized_coordinates_euler(gc2);

  // compare values
  for (unsigned i=0; i< gc1.size(); i++) 
    ASSERT_NEAR(gc1[i], gc2[i], EPS_DOUBLE);
}

int main(int argc, char* argv[])
{
  // set the filename
  IntegrationTest::filename = "../test/pendulum-fixed.urdf";

  // run Google tests
  ::testing::InitGoogleTest(&argc, argv);
  int test_result = RUN_ALL_TESTS();
  if (test_result != 0)
    return test_result;

  // set the filename
  IntegrationTest::filename = "../test/rmp_440SE.urdf";

  // run Google tests
  ::testing::InitGoogleTest(&argc, argv);
  test_result = RUN_ALL_TESTS();
  if (test_result != 0)
    return test_result;

  // set the filename
  IntegrationTest::filename = "../test/07-physics.urdf";
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

