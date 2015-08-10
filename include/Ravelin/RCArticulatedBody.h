/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef RCARTICULATEDBODY 
#error This class is not to be included by the user directly. Use RCArticulatedBodyd.h or RCArticulatedBodyf.h instead.
#endif

class JOINT;

/// Defines an articulated body for use with reduced-coordinate dynamics algorithms
/**
 * Reduced-coordinate articulated bodies cannot rely upon the integrator to automatically update
 * the states (i.e., positions, velocities) of the links, as is done with maximal-coordinate 
 * articulated bodies.  Rather, the integrator updates the joint positions and velocities; the
 * states are obtained from this reduced-coordinate representation.
 * Notes about concurrency: <br /><br />
 *
 * It is generally desirable to be able to run forward dynamics and inverse 
 * dynamics algorithms concurrently to simulate actual robotic systems.   In 
 * general, derived classes should not operate on state variables
 * (joint positions, velocities, accelerations and floating base positions, 
 * velocites, and accelerations) directly during execution of the algorithm.  
 * Rather, derived classes should operate on copies of the state
 * variables, updating the state variables on conclusion of the algorithms.  
 */
class RCARTICULATEDBODY : public ARTICULATEDBODY
{
  friend class CRBALGORITHM;
  friend class FSABALGORITHM;

  public:
    enum ForwardDynamicsAlgorithmType { eFeatherstone, eCRB }; 
    RCARTICULATEDBODY();
    virtual ~RCARTICULATEDBODY() {}
    virtual void reset_accumulators();
    virtual void update_link_poses();    
    virtual void update_link_velocities();
    virtual void apply_impulse(const SMOMENTUM& w, boost::shared_ptr<RIGIDBODY> link);
    virtual void calc_fwd_dyn();
    boost::shared_ptr<RCARTICULATEDBODY> get_this() { return boost::dynamic_pointer_cast<RCARTICULATEDBODY>(shared_from_this()); }
    boost::shared_ptr<const RCARTICULATEDBODY> get_this() const { return boost::dynamic_pointer_cast<const RCARTICULATEDBODY>(shared_from_this()); }
    virtual void set_generalized_forces(const SHAREDVECTORN& gf);
    virtual void set_generalized_forces(const VECTORN& gf) { DYNAMICBODY::set_generalized_forces(gf); }
    virtual void add_generalized_force(const SHAREDVECTORN& gf);
    virtual void add_generalized_force(const VECTORN& gf) { DYNAMICBODY::add_generalized_force(gf); }
    virtual void apply_generalized_impulse(const SHAREDVECTORN& gj);
    virtual void apply_generalized_impulse(const VECTORN& gj) { DYNAMICBODY::apply_generalized_impulse(gj); }
    virtual void set_generalized_coordinates(DYNAMICBODY::GeneralizedCoordinateType gctype, const SHAREDVECTORN& gc);
    virtual void set_generalized_coordinates(DYNAMICBODY::GeneralizedCoordinateType gctype, const VECTORN& gc) { DYNAMICBODY::set_generalized_coordinates(gctype, gc); }
    virtual void set_generalized_velocity(DYNAMICBODY::GeneralizedCoordinateType gctype, const SHAREDVECTORN& gv);
    virtual void set_generalized_velocity(DYNAMICBODY::GeneralizedCoordinateType gctype, const VECTORN& gv) { DYNAMICBODY::set_generalized_velocity(gctype, gv); }
    virtual SHAREDMATRIXN& get_generalized_inertia(SHAREDMATRIXN& M);
    virtual MATRIXN& get_generalized_inertia(MATRIXN& M) { return DYNAMICBODY::get_generalized_inertia(M); }
    virtual SHAREDVECTORN& get_generalized_forces(SHAREDVECTORN& f);
    virtual VECTORN& get_generalized_forces(VECTORN& f) { return DYNAMICBODY::get_generalized_forces(f); }
    virtual SHAREDVECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLEBODY> body, const SFORCE& w, SHAREDVECTORN& gf);
    virtual VECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLEBODY> body, const SFORCE& w, VECTORN& gf) { return DYNAMICBODY::convert_to_generalized_force(body, w, gf); }
    virtual unsigned num_generalized_coordinates(DYNAMICBODY::GeneralizedCoordinateType gctype) const;
    virtual void set_links_and_joints(const std::vector<boost::shared_ptr<RIGIDBODY> >& links, const std::vector<boost::shared_ptr<JOINT> >& joints);
    virtual unsigned num_joint_dof_implicit() const;
    virtual unsigned num_joint_dof_explicit() const { return _n_joint_DOF_explicit; }
    void set_floating_base(bool flag);
    virtual void set_computation_frame_type(ReferenceFrameType rftype);
    virtual VECTORN& solve_generalized_inertia(const VECTORN& b, VECTORN& x) { return DYNAMICBODY::solve_generalized_inertia(b, x); }
    virtual SHAREDMATRIXN& transpose_solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X);
    virtual SHAREDVECTORN& solve_generalized_inertia(const SHAREDVECTORN& v, SHAREDVECTORN& result);
    virtual SHAREDMATRIXN& solve_generalized_inertia(const SHAREDMATRIXN& m, SHAREDMATRIXN& result);
    virtual boost::shared_ptr<const POSE3> get_gc_pose() const; 
    virtual void validate_position_variables();
    virtual SHAREDVECTORN& get_generalized_coordinates(DYNAMICBODY::GeneralizedCoordinateType gctype, SHAREDVECTORN& gc);
    virtual VECTORN& get_generalized_coordinates(DYNAMICBODY::GeneralizedCoordinateType gctype, VECTORN& gc) { return DYNAMICBODY::get_generalized_coordinates(gctype, gc); }
    virtual SHAREDVECTORN& get_generalized_velocity(DYNAMICBODY::GeneralizedCoordinateType gctype, SHAREDVECTORN& gv);
    virtual VECTORN& get_generalized_velocity(DYNAMICBODY::GeneralizedCoordinateType gctype, VECTORN& gv) { return DYNAMICBODY::get_generalized_velocity(gctype, gv); }
    virtual SHAREDVECTORN& get_generalized_acceleration(SHAREDVECTORN& ga);
    void set_generalized_acceleration(const SHAREDVECTORN& a);
    virtual VECTORN& get_generalized_acceleration(VECTORN& ga) { return DYNAMICBODY::get_generalized_acceleration(ga); }

    template <class V>
    void get_generalized_acceleration_generic(V& ga);

    template <class V>
    void get_generalized_coordinates_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_coordinates_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_velocity_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gv);

    template <class V>
    void get_generalized_velocity_generic(DYNAMICBODY::GeneralizedCoordinateType gctype, V& gv);

    /// Gets whether the base of this body is fixed or "floating"
    virtual bool is_floating_base() const { return _floating_base; }

    /// Gets the number of DOF of the explicit joints in the body, not including floating base DOF
    virtual unsigned num_joint_dof() const { return _n_joint_DOF_explicit + num_joint_dof_implicit(); }

    /// Gets the base link
    virtual boost::shared_ptr<RIGIDBODY> get_base_link() const { return (!_links.empty()) ? _links.front() : boost::shared_ptr<RIGIDBODY>(); }

    /// The forward dynamics algorithm
    ForwardDynamicsAlgorithmType algorithm_type;

  protected:
    /// Whether this body uses a floating base
    bool _floating_base;
  
     virtual void compile();

    /// The number of DOF of the explicit joint constraints in the body (does not include floating base DOF!)
    unsigned _n_joint_DOF_explicit;

    /// Gets the vector of explicit joint constraints
    const std::vector<boost::shared_ptr<JOINT> >& get_explicit_joints() const { return _ejoints; }

  private:
    RCARTICULATEDBODY(const RCARTICULATEDBODY& rcab) {}
    virtual MATRIXN& calc_jacobian_column(boost::shared_ptr<JOINT> joint, const VECTOR3& point, MATRIXN& Jc);
/*
    virtual MATRIXN& calc_jacobian_floating_base(const VECTOR3& point, MATRIXN& J);
*/
    bool all_children_processed(boost::shared_ptr<RIGIDBODY> link) const;

    /// The vector of explicit joint constraints
    std::vector<boost::shared_ptr<JOINT> > _ejoints;

    /// The vector of implicit joint constraints
    std::vector<boost::shared_ptr<JOINT> > _ijoints;

    /// Indicates when position data has been invalidated
    bool _position_invalidated;

    /// The CRB algorithm
    CRBALGORITHM _crb;

    /// The FSAB algorithm
    FSABALGORITHM _fsab;

    /// Linear algebra object
    boost::shared_ptr<LINALG> _LA;

    static REAL sgn(REAL x);
    bool treat_link_as_leaf(boost::shared_ptr<RIGIDBODY> link) const;
    void update_factorized_generalized_inertia();
    static bool supports(boost::shared_ptr<JOINT> joint, boost::shared_ptr<RIGIDBODY> link);
    void determine_generalized_forces(VECTORN& gf) const;
    void determine_generalized_accelerations(VECTORN& xdd) const;
    void determine_constraint_force_transform(MATRIXN& K) const;
    void determine_implicit_constraint_movement_jacobian(MATRIXN& D);
    void determine_implicit_constraint_jacobian(MATRIXN& J);
    void determine_implicit_constraint_jacobian_dot(MATRIXN& J);
    void set_implicit_constraint_forces(const VECTORN& lambda);
}; // end class

#include "RCArticulatedBody.inl"


