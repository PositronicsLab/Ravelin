/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef ARTICULATED_BODY 
#error This class is not to be included by the user directly. Use ArticulatedBodyd.h or ArticulatedBodyf.h instead.
#endif

class RIGIDBODY;
class JOINT;

/// Abstract class for articulated bodies
class ARTICULATED_BODY : public virtual DYNAMIC_BODY
{
  public:
    ARTICULATED_BODY();
    virtual ~ARTICULATED_BODY() {}
    virtual bool is_floating_base() const = 0;
    virtual boost::shared_ptr<RIGIDBODY> get_base_link() const = 0;
    unsigned num_constraint_eqns_explicit() const;
    unsigned num_constraint_eqns_implicit() const;
    virtual void rotate(const QUAT& q);
    virtual void translate(const ORIGIN3& o);
    virtual REAL calc_kinetic_energy(boost::shared_ptr<const POSE3> P = boost::shared_ptr<const POSE3>());
    boost::shared_ptr<RIGIDBODY> find_link(const std::string& id) const; 
    boost::shared_ptr<JOINT> find_joint(const std::string& id) const; 
    void get_adjacent_links(std::list<sorted_pair<boost::shared_ptr<RIGIDBODY> > >& links) const;
    virtual void set_links_and_joints(const std::vector<boost::shared_ptr<RIGIDBODY> >& links, const std::vector<boost::shared_ptr<JOINT> >& joints);
    virtual unsigned num_joint_dof() const;
    void find_loops(std::vector<unsigned>& loop_indices, std::vector<std::vector<unsigned> >& loop_links) const;
    virtual MATRIXN& calc_jacobian(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);
    virtual MATRIXN& calc_jacobian_dot(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);
    virtual MATRIXN& calc_jacobian(boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);
    virtual MATRIXN& calc_jacobian_dot(boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);

    /// Articulated bodies are always enabled
    virtual bool is_enabled() const { return true; }

    /// Gets the number of degrees-of-freedom permitted by explicit constraints
    virtual unsigned num_joint_dof_explicit() const = 0;

    /// Gets the number of degrees-of-freedom permitted by implicit constraints
    virtual unsigned num_joint_dof_implicit() const = 0;

    /// Gets the set of links
    virtual const std::vector<boost::shared_ptr<RIGIDBODY> >& get_links() const { return _links; }

    /// Gets the set of joints
    virtual const std::vector<boost::shared_ptr<JOINT> >& get_joints() const { return _joints; }

    /// Gets the set of explicit joints
    virtual const std::vector<boost::shared_ptr<JOINT> >& get_explicit_joints() const = 0; 

    /// Gets the set of implicit joints
    virtual const std::vector<boost::shared_ptr<JOINT> >& get_implicit_joints() const = 0; 

    /// Gets shared pointer to this object as type ARTICULATED_BODY
    boost::shared_ptr<ARTICULATED_BODY> get_this() { return boost::dynamic_pointer_cast<ARTICULATED_BODY>(shared_from_this()); }

    /// Gets shared pointer to this object as type const ArticulateBody
    boost::shared_ptr<const ARTICULATED_BODY> get_this() const { return boost::dynamic_pointer_cast<const ARTICULATED_BODY>(shared_from_this()); }

    /// Abstract method for applying an impulse to this articulated body
    /**
     * \param w the impulsive force 
     * \param link link in the articulated body where the impulse is applied
     */
    virtual void apply_impulse(const SMOMENTUM& w, boost::shared_ptr<RIGIDBODY> link) = 0;
      
    /// Method for resetting the force and torque accumulators on all links
    virtual void reset_accumulators() = 0;

  protected:
    /// Vector for processing links
    std::vector<unsigned> _processed;

    /// Method for "compiling" the body
    /**
     * Compilation is necessary anytime the structure or topolgy of the body
     * changes.  Generally, this will only be necessary once - after a call
     * to set_links() or set_joints().
     */
    virtual void compile();

    /// The set of links for this articulated body
    std::vector<boost::shared_ptr<RIGIDBODY> > _links;

    /// The set of joints for this articulated body
    std::vector<boost::shared_ptr<JOINT> > _joints;

  private:
    ARTICULATED_BODY(const ARTICULATED_BODY& ab) {}

    // joint constraint violation
    std::vector<REAL> _cvio;

    // joint velocity tolerances (for joints at constraints)
    std::vector<REAL> _cvel_vio;

    // temporary variables
    VECTORN _dq;

}; // end class


