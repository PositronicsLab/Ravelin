/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RAVELIN_XML_TREE_H
#define _RAVELIN_XML_TREE_H

#include <boost/enable_shared_from_this.hpp>
#include <list>
#include <string>
#include <set>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <Ravelin/Origin3d.h>
#include <Ravelin/Origin3f.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/VectorNf.h>
#include <Ravelin/SVector6d.h>
#include <Ravelin/SVector6f.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/MatrixNf.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/Pose3d.h>
#include <Ravelin/Pose3f.h>

namespace Ravelin {

/// Attributes used for XML nodes
class XMLAttrib
{
  public:
    XMLAttrib(const std::string& name, const std::string& string_value);
    XMLAttrib(const std::string& name, double real_value);
    XMLAttrib(const std::string& name, float real_value);
    XMLAttrib(const std::string& name, int int_value);
    XMLAttrib(const std::string& name, unsigned unsigned_value);
    XMLAttrib(const std::string& name, const std::vector<SVelocityd>& velocity_values);
    XMLAttrib(const std::string& name, const std::vector<SVelocityf>& velocity_values);
    XMLAttrib(const std::string& name, const std::vector<SAcceld>& accel_values);
    XMLAttrib(const std::string& name, const std::vector<SAccelf>& accel_values);
    XMLAttrib(const std::string& name, double roll, double pitch, double yaw);
    XMLAttrib(const std::string& name, float roll, float pitch, float yaw);
    XMLAttrib(const std::string& name, const Vector2d& vector_value);
    XMLAttrib(const std::string& name, const Vector2f& vector_value);
    XMLAttrib(const std::string& name, const Vector3d& vector_value);
    XMLAttrib(const std::string& name, const Vector3f& vector_value);
    XMLAttrib(const std::string& name, const VectorNd& vector_value);
    XMLAttrib(const std::string& name, const VectorNf& vector_value);
    XMLAttrib(const std::string& name, const SVector6d& vector_value);
    XMLAttrib(const std::string& name, const SVector6f& vector_value);
    XMLAttrib(const std::string& name, const MatrixNd& matrix_value);
    XMLAttrib(const std::string& name, const MatrixNf& matrix_value);
    XMLAttrib(const std::string& name, const Matrix3d& matrix_value);
    XMLAttrib(const std::string& name, const Matrix3f& matrix_value);
    XMLAttrib(const std::string& name, const Quatd& quat_value);
    XMLAttrib(const std::string& name, const Quatf& quat_value);
    XMLAttrib(const std::string& name, const Origin3d& origin_value);
    XMLAttrib(const std::string& name, const Origin3f& origin_value);
    XMLAttrib(const std::string& name, bool bool_value);
    XMLAttrib(const std::string& name, long long_value);
    const std::string& get_string_value() { processed = true; return value; }
    static std::string str(double value);
    static std::string str(float value);
    void get_real_value(double& value);
    void get_real_value(float& value);
    void get_origin_value(Origin3d& origin);
    void get_origin_value(Origin3f& origin);
    void get_velocity_values(std::vector<SVelocityd>& values);
    void get_velocity_values(std::vector<SVelocityf>& values);
    void get_accel_values(std::vector<SAcceld>& values);
    void get_accel_values(std::vector<SAccelf>& values);
    int get_int_value() { processed = true; return std::atoi(value.c_str()); }
    unsigned get_unsigned_value() { processed = true; return (unsigned) std::atoi(value.c_str()); }
    bool get_bool_value();
    long get_long_value() { return std::atol(value.c_str()); }
    std::list<std::string> get_strings_value();
    void get_quat_value(Quatd& q);
    void get_axis_angle_value(Quatd& q);
    void get_rpy_value(Quatd& q);
    void get_quat_value(Quatf& q);
    void get_axis_angle_value(Quatf& q);
    void get_rpy_value(Quatf& q);
    void get_vector_value(VectorNd& v);
    void get_vector_value(Vector2d& v);
    void get_vector_value(Vector3d& v);
    void get_vector_value(SVector6d& v);
    void get_matrix_value(Matrix3d& m);
    void get_matrix_value(MatrixNd& m);
    void get_vector_value(VectorNf& v);
    void get_vector_value(Vector2f& v);
    void get_vector_value(Vector3f& v);
    void get_vector_value(SVector6f& v);
    void get_matrix_value(Matrix3f& m);
    void get_matrix_value(MatrixNf& m);
    bool operator==(const XMLAttrib& a) const { return name == a.name; }
    bool operator<(const XMLAttrib& a) const { return name < a.name; }

    /// The name of the attribute
    std::string name;

    /// The value in string form
    std::string value;

    /// Indicates whether this attribute has been processed
    bool processed;
};

/// An XML tree used for serialization 
class XMLTree : public boost::enable_shared_from_this<XMLTree>
{
  public:
    XMLTree(const std::string& name);
    XMLTree(const std::string& name, const std::list<XMLAttrib>& attributes);
    static boost::shared_ptr<const XMLTree> read_from_xml(const std::string& name);
    XMLAttrib* get_attrib(const std::string& attrib_name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_child_nodes(const std::string& name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_child_nodes(const std::list<std::string>& name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_descendant_nodes(const std::string& name) const;

    /// Adds a child tree to this tree; also sets the parent node
    void add_child(boost::shared_ptr<XMLTree> child) { children.push_back(child); child->set_parent(shared_from_this()); }

    /// Sets the parent of this tree (if any)
    void set_parent(boost::shared_ptr<XMLTree> parent) { _parent = parent; }

    /// Gets the parent of this tree (if any)
    boost::weak_ptr<XMLTree> get_parent() const { return _parent; }

    /// The set of attributes of this node
    std::set<XMLAttrib> attribs;  

    /// The list of children of this node
    std::list<boost::shared_ptr<XMLTree> > children;

    /// The name of this node
    std::string name;

    /// The ID of this node
    std::string id;

    /// Any 'content' of this node
    std::string content;

    /// The object (if any) represented by this node
    boost::shared_ptr<void> object;

    /// Indicates whether this tag has been processed
    bool processed;

  private:
    boost::weak_ptr<XMLTree> _parent;
    static boost::shared_ptr<const XMLTree> construct_xml_tree(xmlNode* root);
}; // end class

std::ostream& operator<<(std::ostream& out, const XMLTree& tree);
std::ostream& operator<<(std::ostream& out, const XMLAttrib& attr);

} // end namespace

#endif

