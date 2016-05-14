/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#ifndef JOINT
#error This class is not to be included by the user directly. Use Jointd.h or Jointf.h instead.
#endif

class ARTICULATED_BODY;
class RIGIDBODY;

/// Defines a bilateral constraint (a joint)
class JOINT : public virtual_enable_shared_from_this<JOINT>
{
  public:
    enum ConstraintType { eUnknown, eExplicit, eImplicit };
    enum DOFs { DOF_1=0, DOF_2=1, DOF_3=2, DOF_4=3, DOF_5=4, DOF_6=5 };

    JOINT();
    virtual const std::vector<SVELOCITY>& get_spatial_axes();
    void add_force(const VECTORN& force);
    void reset_force(); 
    ConstraintType get_constraint_type() const { return _constraint_type; }
    boost::shared_ptr<ARTICULATED_BODY> get_articulated_body();

    /// The ID of this joint
    std::string joint_id;

    /// Sets whether this constraint is implicit or explicit (or unknown)
    void set_constraint_type(ConstraintType type) { _constraint_type = type; }

    virtual void set_inboard_link(boost::shared_ptr<RIGIDBODY> link, bool update_pose);
    virtual void set_outboard_link(boost::shared_ptr<RIGIDBODY> link, bool update_pose);
    void set_location(const VECTOR3& p, boost::shared_ptr<RIGIDBODY> inboard, boost::shared_ptr<RIGIDBODY> outboard);
    VECTOR3 get_location(bool use_outboard = false) const;
    virtual void set_inboard_pose(boost::shared_ptr<const POSE3> inboard_pose, bool update_joint_pose);
    virtual void set_outboard_pose(boost::shared_ptr<POSE3> outboard_pose, bool update_joint_pose);
    virtual void update_spatial_axes();
    virtual void evaluate_constraints_dot(REAL C[]);
    virtual void set_q_tare(const VECTORN& tare) { _q_tare = tare; }
    virtual const VECTORN& get_q_tare() const { return _q_tare; }  

    /// Gets the inboard link for this joint
    boost::shared_ptr<RIGIDBODY> get_inboard_link() const { return (_inboard_link.expired()) ? boost::shared_ptr<RIGIDBODY>() : boost::shared_ptr<RIGIDBODY>(_inboard_link); }

    /// Gets the outboard link for this joint
    boost::shared_ptr<RIGIDBODY> get_outboard_link() const { return (_outboard_link.expired()) ? boost::shared_ptr<RIGIDBODY>() : boost::shared_ptr<RIGIDBODY>(_outboard_link); }

    /// The acceleration of this joint
    VECTORN qdd;

    /// The actuator force (user/controller sets this)
    VECTORN force;

    /// Constraint forces calculated by forward dynamics
    VECTORN lambda;

    /// Gets the joint index (returns UINT_MAX if not set)
    unsigned get_index() const { return _joint_idx; }

    /// Sets the joint index
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_index(unsigned index) { _joint_idx = index; }

    /// Gets the constraint index for this joint
    unsigned get_constraint_index() const { return _constraint_idx; }

    /// Sets the constraint index for this joint
    void set_constraint_index(unsigned idx) { _constraint_idx = idx; } 

    /// Sets the coordinate index for this joint
    /**
     * This is set automatically by the articulated body. Users should not
     * change this index or unknown behavior will result.
     */
    void set_coord_index(unsigned index) { _coord_idx = index; }

    /// Gets the starting coordinate index for this joint
    unsigned get_coord_index() const { return _coord_idx; }

    /// Gets the pose of this joint (relative to the inboard pose instead of the outboard pose as returned by get_pose_from_outboard())
    boost::shared_ptr<const POSE3> get_pose() const { return _F; };

    /// Gets the pose of this joint relative to the outboard pose (rather than the inboard pose as returned by get_pose())
    boost::shared_ptr<const POSE3> get_pose_from_outboard() const { return _Fb; };

    /// Gets whether this joint is in a singular configuration
    /**
     * \note only used by reduced-coordinate articulated bodies
     */
    virtual bool is_singular_config() const = 0;

    /// Gets the number of constraint equations for this joint
    virtual unsigned num_constraint_eqns() const { return 6 - num_dof(); }

    /// Evaluates the constraint equations for this joint
    /**
     * When the joint constraints are exactly satisfied, the result will be
     * the zero vector.
     * \param C a vector of size num_constraint_eqns(); contains the evaluation
     *        (on return)
     * \note only used by maximal-coordinate articulated bodies
     */
    virtual void evaluate_constraints(REAL C[]) = 0;

    /// Abstract method to get the spatial axes derivatives for this joint
    /**
     * Only applicable for reduced-coordinate articulated bodies
     */
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot() = 0;

    /// Abstract method to get the local transform for this joint
    /**
     * The local transform for the joint transforms the coordinate frame
     * attached to the joint center and aligned with the inner link frame.
     */
    virtual boost::shared_ptr<const POSE3> get_induced_pose() = 0;

    /// Abstract method to determine the value of Q (joint position) from current transforms
    virtual void determine_q(VECTORN& q) = 0;

    /// Gets the number of degrees-of-freedom for this joint
    virtual unsigned num_dof() const = 0;
  
    /// The position of this joint
    VECTORN q;

    /// The velocity of this joint
    VECTORN qd;

    /// Computes the constraint Jacobian for this joint with respect to the given body
    /**
     * \param inboard 'true' if the Jacobian is to be computed w.r.t. the
     *        inboard link; 'false' for the outboard
     * \param Cq a 6 x ndof matrix for the given body (on return) 
     */
    virtual void calc_constraint_jacobian(bool inboard, MATRIXN& Cq) = 0;
 
     /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body
    /**
     * \param inboard 'true' if the Jacobian is to be computed w.r.t. the
     *        inboard link; 'false' for the outboard
     * \param Cq a 6 x ndof matrix for the given body (on return)
     */
    virtual void calc_constraint_jacobian_dot(bool inboard, MATRIXN& Cq) = 0;

  protected:
    void calc_constraint_jacobian_numeric(bool inboard, MATRIXN& Cq);
    bool transform_jacobian(MATRIXN& J, bool use_inboard, MATRIXN& output);
    void invalidate_pose_vectors() { get_outboard_link()->invalidate_pose_vectors(); }
    boost::shared_ptr<const POSE3> get_inboard_pose() { if (_inboard_link.expired()) return boost::shared_ptr<const POSE3>(); return get_inboard_link()->get_pose(); }
    boost::shared_ptr<const POSE3> get_outboard_pose() { if (_outboard_link.expired()) return boost::shared_ptr<const POSE3>(); return get_outboard_link()->get_pose(); }

    /// The frame induced by the joint 
    boost::shared_ptr<POSE3> _Fprime;

    /// The frame of this joint
    boost::shared_ptr<POSE3> _F;

    /// The frame of this joint _backward_ from the outboard link
    boost::shared_ptr<POSE3> _Fb;

    /// Method for initializing all variables in the joint
    /**
     * This method should be called at the beginning of all constructors of
     * all derived classes.
     */
    virtual void init_data();

    /// The spatial axes (in joint frame) for the joint
    /**
     * Spatial axes are used in the dynamics equations for reduced-coordinate
     * articulated bodies only.
     */
    std::vector<SVELOCITY> _s;

    /// The stored "tare" value for the initial joint configuration
    /**
     * The tare value is the value that the joint assumes in the an
     * articulated body's initial configuration. This value is necessary
     * so that- when the body's joints are set to the zero vector- the body
     * re-enters the initial configuration.
     */
    VECTORN _q_tare;

    boost::weak_ptr<RIGIDBODY> _inboard_link;
    boost::weak_ptr<RIGIDBODY> _outboard_link;
    ConstraintType _constraint_type;
    unsigned _joint_idx;
    unsigned _coord_idx;
    unsigned _constraint_idx;
}; // end class


