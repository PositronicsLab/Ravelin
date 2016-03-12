/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

class RIGIDBODY;
class RC_ARTICULATED_BODY;
class JOINT;

/// Used to read the simulator state from URDF
class URDFREADER
{
  public:
    static bool read(const std::string& fname, std::string& name, std::vector<boost::shared_ptr<RIGIDBODY> >& links, std::vector<boost::shared_ptr<JOINT> >& joints);
    
  private:
    class URDFData
    {
      public:
        ~URDFData()
        {
        }

        std::map<boost::shared_ptr<JOINT>, boost::shared_ptr<RIGIDBODY> > joint_parent, joint_child;
        std::map<boost::shared_ptr<RIGIDBODY>, boost::shared_ptr<POSE3> > inertial_poses;
        std::map<std::string, std::pair<VectorNd, std::string> > materials;
    };

    static void find_outboards(const URDFData& data, boost::shared_ptr<RIGIDBODY> link, std::vector<std::pair<boost::shared_ptr<JOINT>, boost::shared_ptr<RIGIDBODY> > >& outboards, std::map<boost::shared_ptr<RIGIDBODY>, boost::shared_ptr<RIGIDBODY> >& parents);
    static void output_data(const URDFData& data, boost::shared_ptr<RIGIDBODY> link);
    static boost::shared_ptr<JOINT> find_joint(const URDFData& data, boost::shared_ptr<RIGIDBODY> outboard_link);
    static void find_children(const URDFData& data, boost::shared_ptr<RIGIDBODY> link, std::queue<boost::shared_ptr<RIGIDBODY> >& q, std::map<boost::shared_ptr<RIGIDBODY>, boost::shared_ptr<RIGIDBODY> >& parents);
    static MATRIX3 read_inertia(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static double read_mass(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static POSE3 read_origin(boost::shared_ptr<const XMLTree> node, URDFData& data);
    static void read_inertial(boost::shared_ptr<const XMLTree> node, URDFData& data, boost::shared_ptr<RIGIDBODY> link);
    static void read_axis(boost::shared_ptr<const XMLTree> node, URDFData& data, boost::shared_ptr<JOINT> joint); 
    static boost::shared_ptr<RIGIDBODY> read_parent(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<boost::shared_ptr<RIGIDBODY> >& links); 
    static boost::shared_ptr<RIGIDBODY> read_child(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<boost::shared_ptr<RIGIDBODY> >& links); 
    static void read_joint(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<boost::shared_ptr<RIGIDBODY> >& links, std::vector<boost::shared_ptr<JOINT> >& joints); 
    static void read_joints(boost::shared_ptr<const XMLTree> node, URDFData& data, const std::vector<boost::shared_ptr<RIGIDBODY> >& links, std::vector<boost::shared_ptr<JOINT> >& joints); 
    static void read_links(boost::shared_ptr<const XMLTree> node, URDFData& data, std::vector<boost::shared_ptr<RIGIDBODY> >& links); 
    static void read_link(boost::shared_ptr<const XMLTree> node, URDFData& data, std::vector<boost::shared_ptr<RIGIDBODY> >& links); 
    static bool read_robot(boost::shared_ptr<const XMLTree> node, URDFData& data, std::string& name, std::vector<boost::shared_ptr<RIGIDBODY> >& links, std::vector<boost::shared_ptr<JOINT> >& joints); 
}; // end class


