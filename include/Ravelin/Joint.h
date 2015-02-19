/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 license 
 ****************************************************************************/

#ifndef JOINT
#error This class is not to be included by the user directly. Use Jointd.h or Jointf.h instead.
#endif

/// Defines a bilateral constraint (a joint)
class JOINT
{
  public:
    enum DOFs { DOF_1=0, DOF_2=1, DOF_3=2, DOF_4=3, DOF_5=4, DOF_6=5 };

    JOINT();
    virtual const std::vector<SVELOCITY>& get_spatial_axes();
    void set_location(const VECTOR3& p);
    VECTOR3 get_location(bool use_outboard = false) const;
    virtual void set_inboard_pose(boost::shared_ptr<const POSE3> inboard_pose, bool update_joint_pose);
    virtual void set_outboard_pose(boost::shared_ptr<POSE3> outboard_pose, bool update_joint_pose);
    virtual void update_spatial_axes();
    void evaluate_constraints_dot(REAL C[6]);
    virtual void determine_q_tare();

    /// Gets the pose of this joint (relative to the inboard pose instead of the outboard pose as returned by get_pose_from_outboard())
    boost::shared_ptr<const POSE3> get_pose() const { return _F; };

    // Gets the inboard pose 
    boost::shared_ptr<const POSE3> get_inboard_pose() const { return _F->rpose; } 

    // Gets the outboard pose 
    boost::shared_ptr<const POSE3> get_outboard_pose() const { return _Fb->rpose; } 

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

  protected:

    /// The frame induced by the joint 
    boost::shared_ptr<POSE3> _Fprime;

    /// The frame of this joint
    boost::shared_ptr<POSE3> _F;

    /// The frame of this joint _backward_ from the outboard link
    boost::shared_ptr<POSE3> _Fb;

    /// Computes the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
//    virtual void calc_constraint_jacobian(RigidBodyPtr body, unsigned index, REAL Cq[]) = 0;
 
     /// Computes the time derivative of the constraint Jacobian for this joint with respect to the given body in Rodrigues parameters
    /**
     * \param body the body with which the constraint Jacobian will be 
     *        calculated; if the body is not either the inner or outer link,
     *        the constraint Jacobian will be a zero matrix
     * \param index the index of the constraint equation for which to calculate
     * \param Cq a vector that contains the corresponding column of the
     *        constraint Jacobian on return
     */
//    virtual void calc_constraint_jacobian_dot(RigidBodyPtr body, unsigned index, REAL Cq[]) = 0;

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
}; // end class


