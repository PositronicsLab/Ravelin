/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using std::set;
using std::make_pair;
using std::queue;
using std::vector;
using std::pair;
using std::map;
using std::string;
using std::list;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// Reads an XML file and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
bool URDFREADER::read(const string& fname, std::string& name, vector<shared_ptr<RIGIDBODY> >& links, vector<shared_ptr<JOINT> >& joints)
{
  // *************************************************************
  // going to remove any path from the argument and change to that
  // path; this is done so that all files referenced from the
  // local path of the XML file are found
  // *************************************************************

  // set the filename to use as the argument, by default
  string filename = fname;

  // get the current pathname
  size_t BUFSIZE = 8192;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "URDFReader::read() - unable to allocate sufficient memory!" << std::endl;
      return false;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != string::npos)
  {
    // get the new working path
    string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,string::npos);
  }

  // read the XML Tree 
  shared_ptr<const XMLTree> tree = XMLTree::read_from_xml(filename);
  if (!tree)
  {
    std::cerr << "URDFReader::read() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    chdir(cwd.get());
    return false;
  }
  
  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RIGIDBODY is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // read and construct (single) robot 
  if (strcasecmp(tree->name.c_str(), "Robot") == 0)
  {
    URDFData data;
    if (!read_robot(tree, data, name, links, joints))
      return false;
  }
  else
  {
    std::cerr << "URDFReader::read() error - root element of URDF file is not a 'Robot' tag" << std::endl;
    return false;
  }

  // change back to the initial working directory
  chdir(cwd.get());

  return true;
}

/// Reads an XML string and constructs all read objects
/**
 * \return a map of IDs to read objects
 */
bool URDFREADER::read_from_string(const string& content, std::string& name, vector<shared_ptr<RIGIDBODY> >& links, vector<shared_ptr<JOINT> >& joints)
{
  // read the XML Tree 
  shared_ptr<const XMLTree> tree = XMLTree::read_from_string(content);
  if (!tree)
  {
    std::cerr << "URDFReader::read() - unable to read xml content " << std::endl;
    return false;
  }
  
  // ********************************************************************
  // NOTE: read_from_xml() (via process_tag()) treats all nodes at the
  // same level; it is irrelevant to it whether a RIGIDBODY is
  // inside or outside of its encapsulating body.  It will construct the
  // objects properly; nodes that rely on hierarchies in the XML file must
  // provide this processing themselves (see RCArticulatedBody for an example)
  // ********************************************************************

  // read and construct (single) robot 
  if (strcasecmp(tree->name.c_str(), "Robot") == 0)
  {
    URDFData data;
    if (!read_robot(tree, data, name, links, joints))
      return false;
  }
  else
  {
    std::cerr << "URDFReader::read() error - root element of URDF is not a 'Robot' tag" << std::endl;
    return false;
  }

  return true;
}

/// Reads and constructs a robot object
bool URDFREADER::read_robot(shared_ptr<const XMLTree> node, URDFData& data, string& name, vector<shared_ptr<RIGIDBODY> >& links, vector<shared_ptr<JOINT> >& joints)
{
  // read the robot name
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (name_attrib)
  {
    name = name_attrib->get_string_value();
    read_links(node, data, links);
    read_joints(node, data, links, joints);
  }
  else
  {
    std::cerr << "URDFReader::read_robot() - robot name not specified! not processing further..." << std::endl;
    return false;
  }

  return true;
}

/// Reads robot links
void URDFREADER::read_links(shared_ptr<const XMLTree> node, URDFData& data, vector<shared_ptr<RIGIDBODY> >& links)
{
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_link(*i, data, links);
}

/// Reads robot joints 
void URDFREADER::read_joints(shared_ptr<const XMLTree> node, URDFData& data, const vector<shared_ptr<RIGIDBODY> >& links, vector<shared_ptr<JOINT> >& joints)
{
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    read_joint(*i, data, links, joints);
}

/// Attempts to read a robot link from the given node
void URDFREADER::read_link(shared_ptr<const XMLTree> node, URDFData& data, vector<shared_ptr<RIGIDBODY> >& links)
{
  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Link") != 0)
    return;

  // link must have the name attribute
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (!name_attrib)
    std::cerr << "URDFReader::read_link() - link name not specified! not processing further..." << std::endl;

  // read and construct the link
  shared_ptr<RIGIDBODY> link(new RIGIDBODY);
  link->body_id = name_attrib->get_string_value();

  // read link properties
  read_inertial(node, data, link);

  // add the link to the set of links
  links.push_back(link);
}

/// Finds all children of the given link
void URDFREADER::find_children(const URDFData& data, shared_ptr<RIGIDBODY> link, queue<shared_ptr<RIGIDBODY> >& q, map<shared_ptr<RIGIDBODY>, shared_ptr<RIGIDBODY> >& parents)
{
  for (map<shared_ptr<JOINT>, shared_ptr<RIGIDBODY> >::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
    if (i->second == link)
    {
      assert(data.joint_child.find(i->first) != data.joint_child.end());
      shared_ptr<RIGIDBODY> child = data.joint_child.find(i->first)->second;
      q.push(child);
      parents[child] = link;
    } 
}

/// Finds the inner joint for a link
shared_ptr<JOINT> URDFREADER::find_joint(const URDFData& data, shared_ptr<RIGIDBODY> outboard)
{
  for (map<shared_ptr<JOINT>, shared_ptr<RIGIDBODY> >::const_iterator i = data.joint_child.begin(); i != data.joint_child.end(); i++)
    if (i->second == outboard)
      return i->first;

  return shared_ptr<JOINT>();
}

#ifdef USE_OSG
/// Copies this matrix to an OpenSceneGraph Matrixd object
static void to_osg_matrix(const POSE3& src, osg::Matrixd& tgt)
{
  // get the rotation matrix
  MATRIX3 M = src.q;

  // setup the rotation components of tgt
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = M(i,j);

  // setup the translation components of tgt
  for (unsigned i=X; i<= Z; i++)
    tgt(W,i) = src.x[i];

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (double) 0.0;
  tgt(W,W) = (double) 1.0;
}

/// Copies an OpenSceneGraph Matrixd object to this matrix 
static void from_osg_matrix(const osg::Matrixd& src, POSE3& tgt)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  MATRIX3 R;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      R(j,i) = src(i,j);

  tgt.q = R;
  tgt.x = ORIGIN3(src(W,X), src(W, Y), src(W, Z));
}
#endif

void URDFREADER::find_outboards(const URDFData& data, shared_ptr<RIGIDBODY> link, vector<pair<shared_ptr<JOINT>, shared_ptr<RIGIDBODY> > >& outboards, map<shared_ptr<RIGIDBODY>, shared_ptr<RIGIDBODY> >& parents)
{
  for (map<shared_ptr<JOINT>, shared_ptr<RIGIDBODY> >::const_iterator i = data.joint_parent.begin(); i != data.joint_parent.end(); i++)
    if (i->second == link)
    {
      assert(data.joint_child.find(i->first) != data.joint_child.end());
      shared_ptr<RIGIDBODY> outboard = data.joint_child.find(i->first)->second;
      outboards.push_back(make_pair(i->first, outboard));
      parents[outboard] = link;
    } 
}

void URDFREADER::output_data(const URDFData& data, shared_ptr<RIGIDBODY> link)
{
  #ifdef DEBUG_URDF
  std::cout << "link id: " << link->id << std::endl;
  shared_ptr<JOINT> joint = find_joint(data, link);
  if (joint)
  {
    std::cout << "  Ravelin joint position: " << joint->get_local() << std::endl;
    shared_ptr<REVOLUTEJOINT> rj = dynamic_pointer_cast<REVOLUTEJOINT>(joint);
    shared_ptr<PRISMATICJOINT> pj = dynamic_pointer_cast<PRISMATICJOINT>(joint);
    if (rj)
      std::cout << "  Ravelin joint axis: " << rj->get_axis() << std::endl;
    else if (pj)
      std::cout << "  Ravelin joint axis: " << pj->get_axis() << std::endl;
  }
  std::cout << "  Ravelin pose: " << std::endl << *link->get_pose();
  if (!link->geometries.empty())
  {
    std::cout << "  Ravelin collision pose: " << std::endl << *link->geometries.front()->get_pose();
  }
  #endif
}

/// Attempts to read a robot joint from the given node
void URDFREADER::read_joint(shared_ptr<const XMLTree> node, URDFData& data, const vector<shared_ptr<RIGIDBODY> >& links, vector<shared_ptr<JOINT> >& joints)
{
  const shared_ptr<const POSE3> GLOBAL;
  shared_ptr<JOINT> joint;
  shared_ptr<RIGIDBODY> inboard, outboard;

  // see whether the node name is correct
  if (strcasecmp(node->name.c_str(), "Joint") != 0)
    return;

  // link must have the name attribute
  XMLAttrib* name_attrib = node->get_attrib("name");
  if (!name_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint name not specified! not processing further..." << std::endl;
    return;
  }

  // link must have the name attribute
  XMLAttrib* type_attrib = node->get_attrib("type");
  if (!type_attrib)
  {
    std::cerr << "URDFReader::read_joint() - joint type not specified! not processing further..." << std::endl;
    return;
  }

  // read and construct the joint
  if (strcasecmp(type_attrib->get_string_value().c_str(), "revolute") == 0)
  {
    shared_ptr<REVOLUTEJOINT> rj(new REVOLUTEJOINT);
    joint = rj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "continuous") == 0)
  {
    shared_ptr<REVOLUTEJOINT> rj(new REVOLUTEJOINT);
    joint = rj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "prismatic") == 0)
  {
    shared_ptr<PRISMATICJOINT> pj(new PRISMATICJOINT);
    joint = pj;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "fixed") == 0)
    joint = shared_ptr<FIXEDJOINT>(new FIXEDJOINT);
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "floating") == 0)
  {
    std::cerr << "URDFReader::read_joint() - [deprecated] floating joint type specified! not processing further..." << std::endl;
    return;
  }
  else if (strcasecmp(type_attrib->get_string_value().c_str(), "planar") == 0)
  {
    std::cerr << "URDFReader::read_joint() - planar joint type currently unsupported in Ravelin! not processing further..." << std::endl;
    return;
  }
  else
  {
    std::cerr << "URDFReader::read_joint() - invalid joint type specified! not processing further..." << std::endl;
    return;
  }

  // read and verify required properties
  joint->joint_id = name_attrib->get_string_value();
  if (!(inboard = read_parent(node, data, links)))
  {
    std::cerr << "URDFReader::read_joint() - failed to properly read parent link! not processing further..." << std::endl;
    return;
  }
  if (!(outboard = read_child(node, data, links)))
  {
    std::cerr << "URDFReader::read_joint() - failed to properly read child link! not processing further..." << std::endl;
    return;
  }

  // setup the appropriate pointers
  data.joint_parent[joint] = inboard;
  data.joint_child[joint] = outboard;

  // joint frame is defined relative to the parent link frame
  shared_ptr<POSE3> origin(new POSE3(read_origin(node, data)));
  origin->rpose = inboard->get_pose(); 
  VECTOR3 location_origin(0.0, 0.0, 0.0, origin);
  VECTOR3 location = POSE3::transform_point(GLOBAL, location_origin);

  // setup a second pose, which is the inertial frame
  assert(data.inertial_poses.find(outboard) != data.inertial_poses.end());
  shared_ptr<POSE3> inertial_frame(new POSE3(*data.inertial_poses[outboard]));
  inertial_frame->rpose = origin;

  // update the outboard link pose
  inertial_frame->update_relative_pose(outboard->get_pose()->rpose);
  outboard->set_pose(*inertial_frame);

  // setup the inboard and outboard links for the joint
  joint->set_location(location, inboard, outboard);

  // read optional properties
  read_axis(node, data, joint);

  // add the joint to the set of joints 
  joints.push_back(joint);
}

/// Attempts to read the parent for the joint
shared_ptr<RIGIDBODY> URDFREADER::read_parent(shared_ptr<const XMLTree> node, URDFData& data, const vector<shared_ptr<RIGIDBODY> >& links)
{
  // look for the tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "parent") == 0)
    {
      // read the link attribute
      XMLAttrib* link_attrib = (*i)->get_attrib("link");
      if (!link_attrib)
        continue;
      string link_id = link_attrib->get_string_value();

      // find parent link
      for (unsigned j=0; j< links.size(); j++)
        if (links[j]->body_id == link_id)
          return links[j];
    }
  }

  return shared_ptr<RIGIDBODY>();
}

/// Attempts to read the child for the joint
shared_ptr<RIGIDBODY> URDFREADER::read_child(shared_ptr<const XMLTree> node, URDFData& data, const vector<shared_ptr<RIGIDBODY> >& links)
{
  // look for the tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "child") == 0)
    {
      // read the link attribute
      XMLAttrib* link_attrib = (*i)->get_attrib("link");
      if (!link_attrib)
        continue;
      string link_id = link_attrib->get_string_value();

      // find child link
      for (unsigned j=0; j< links.size(); j++)
        if (links[j]->body_id == link_id)
          return links[j];
    }
  }

  return shared_ptr<RIGIDBODY>();
}

/// Attempts to read axis for the joint
void URDFREADER::read_axis(shared_ptr<const XMLTree> node, URDFData& data, shared_ptr<JOINT> joint)
{
  VECTOR3 axis(1,0,0);
  bool axis_specified = false;

  // look for the axis tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "axis") == 0)
    {
      // read the attributes first
      XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (!xyz_attrib)
        continue;
      xyz_attrib->get_vector_value(axis);
      axis_specified = true;
    }
  }

  // get the outboard link - it's pose is identical to the joint pose
  shared_ptr<RIGIDBODY> outboard = joint->get_outboard_link();

  // setup the axis frame
  axis.pose = joint->get_pose();

  // verify that joint is of the proper type
  shared_ptr<REVOLUTEJOINT> rj = dynamic_pointer_cast<REVOLUTEJOINT>(joint);
  shared_ptr<PRISMATICJOINT> pj = dynamic_pointer_cast<PRISMATICJOINT>(joint);
  if (rj || pj)
  {
    if (rj)
      rj->set_axis(axis);
    else if (pj)
      pj->set_axis(axis);
    else
      assert(false);
  }
  else if (axis_specified)
    std::cerr << "URDFReader::read_axis() - joint axis specified for joint w/o axis!" << std::endl;
}

/// Attempts to read and set link inertial properties 
void URDFREADER::read_inertial(shared_ptr<const XMLTree> node, URDFData& data, shared_ptr<RIGIDBODY> link)
{
  // setup linear algebra object
  LINALG LA;

  // look for the inertial tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertial") == 0)
    {
      REAL mass = read_mass(*i, data);
      MATRIX3 inertia = read_inertia(*i, data);

      // verify that inertial properties are good
      MATRIX3 inertia_copy = inertia;
      if (mass <= 0.0 || !LA.is_SPD(inertia_copy, -1.0))
        link->set_enabled(false); 

      // read the inertial frame
      shared_ptr<POSE3> origin(new POSE3(read_origin(*i, data)));

      // set the inertial frame relative to the link frame
      origin->rpose = link->get_pose();

      // set inertial properties
      SPATIAL_RB_INERTIA J(link->get_pose());
      J.m = mass;
      J.J = inertia;
      link->set_inertia(J);

      // add the inertial pose data
      data.inertial_poses[link] = origin;

      // reading inertial was a success, attempt to read no further...
      // (multiple inertial tags not supported)
      return;
    }
  }
}

/// Attempts to read an "origin" tag
POSE3 URDFREADER::read_origin(shared_ptr<const XMLTree> node, URDFData& data)
{
  ORIGIN3 xyz;
  VECTOR3 rpy;

  // set both to zero
  xyz.set_zero();
  rpy.set_zero();

  // look for the tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "origin") == 0)
    {
      // look for xyz attribute 
      XMLAttrib* xyz_attrib = (*i)->get_attrib("xyz");
      if (xyz_attrib)
        xyz_attrib->get_origin_value(xyz);

      // look for rpy attribute
      XMLAttrib* rpy_attrib = (*i)->get_attrib("rpy");
      if (rpy_attrib)
        rpy_attrib->get_vector_value(rpy);

      // reading tag was a success, attempt to read no further...
      // (multiple such tags not supported)
      break;
    }
  }

  QUAT rpy_quat = QUAT::rpy(rpy[0], rpy[1], rpy[2]);
  return POSE3(rpy_quat, xyz);
}

/// Attempts to read a "mass" tag
double URDFREADER::read_mass(shared_ptr<const XMLTree> node, URDFData& data)
{
  // look for the tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "mass") == 0)
    {
      // look for the "value" attribute
      XMLAttrib* value_attrib = (*i)->get_attrib("value");
      if (value_attrib)
      {
        // reading tag was a success, attempt to read no further...
        // (multiple such tags not supported)
        double value;
        value_attrib->get_real_value(value);
        return value; 
      }
    }
  }

  // couldn't find the tag.. return 0
  return 0.0;
}

/// Attempts to read an "inertia" tag
MATRIX3 URDFREADER::read_inertia(shared_ptr<const XMLTree> node, URDFData& data)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup J to zero initially
  MATRIX3 J = MATRIX3::zero();

  // look for the tag
  const list<shared_ptr<XMLTree> >& child_nodes = node->children;
  for (list<shared_ptr<XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    if (strcasecmp((*i)->name.c_str(), "inertia") == 0)
    {
      // look for the six attributes
      XMLAttrib* ixx_attrib = (*i)->get_attrib("ixx");
      XMLAttrib* ixy_attrib = (*i)->get_attrib("ixy");
      XMLAttrib* ixz_attrib = (*i)->get_attrib("ixz");
      XMLAttrib* iyy_attrib = (*i)->get_attrib("iyy");
      XMLAttrib* iyz_attrib = (*i)->get_attrib("iyz");
      XMLAttrib* izz_attrib = (*i)->get_attrib("izz");

      // set values from present attributes
      if (ixx_attrib)
         ixx_attrib->get_real_value(J(X,X));
      if (iyy_attrib)
        iyy_attrib->get_real_value(J(Y,Y));
      if (izz_attrib)
        izz_attrib->get_real_value(J(Z,Z));
      if (ixy_attrib)
      {
        ixy_attrib->get_real_value(J(X,Y));
        J(Y,X) = J(X,Y);
      }
      if (ixz_attrib)
      {
        ixz_attrib->get_real_value(J(X,Z));
        J(Z,X) = J(X,Z);
      } 
      if (iyz_attrib)
      {
        iyz_attrib->get_real_value(J(Y,Z));
        J(Z,Y) = J(Y,Z);
      }

      // reading tag was a success, attempt to read no further...
      // (multiple such tags not supported)
      return J;
    }
  }

  // no inertia read.. return default J (0 matrix)
  return J;
}



