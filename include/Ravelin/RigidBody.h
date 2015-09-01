/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef RIGIDBODY 
#error This class is not to be included by the user directly. Use RigidBodyd.h or RigidBodyf.h instead.
#endif

class JOINT;
class ARTICULATED_BODY;

/// Represents a single rigid body
/**
 *  Contains information needed to represent a rigid body, including position
 *  and velocity (both linear and angular), mass, center of
 *  mass, inertia matrix, collision data, and visualization data.  This class
 *  is used for both non-articulated and articulated rigid bodies, though not
 *  all member data may be used in the latter.
 * \todo implement rest matrix
 */
class RIGIDBODY : public virtual SINGLE_BODY 
{
  friend class ARTICULATED_BODY;
  friend class RC_ARTICULATED_BODY;
  friend class MCARTICULATED_BODY;
  friend class JOINT;

  public:
    enum Compliance { eRigid, eCompliant};
    RIGIDBODY();
    virtual ~RIGIDBODY() {}
    void add_force(const SFORCE& w);
    void set_pose(const POSE3& pose);
    virtual void set_inertial_pose(const POSE3& pose);
    void apply_impulse(const SMOMENTUM& w);
    virtual void rotate(const QUAT& q);
    virtual void translate(const ORIGIN3& o);
    virtual void calc_fwd_dyn();
    const SPATIAL_RB_INERTIA& get_inertia();
    void set_inertia(const SPATIAL_RB_INERTIA& J);
    boost::shared_ptr<const POSE3> get_inertial_pose() const { return _jF; }

    bool is_child_link(boost::shared_ptr<const RIGIDBODY> query) const;
    bool is_descendant_link(boost::shared_ptr<const RIGIDBODY> query) const;
    const SVELOCITY& get_velocity();
    void set_velocity(const SVELOCITY& xd);
    void set_accel(const SACCEL& xdd);
    virtual const SACCEL& get_accel();
    void set_velocity(const SACCEL& xdd);
    virtual void set_generalized_forces(const SHAREDVECTORN& gf);
    virtual void set_generalized_forces(const VECTORN& gf) { DYNAMIC_BODY::set_generalized_forces(gf); }
    virtual SHAREDMATRIXN& get_generalized_inertia(SHAREDMATRIXN& M);
    virtual MATRIXN& get_generalized_inertia(MATRIXN& M) { return DYNAMIC_BODY::get_generalized_inertia(M); }
    virtual SHAREDVECTORN& get_generalized_coordinates(DYNAMIC_BODY::GeneralizedCoordinateType gctype, SHAREDVECTORN& gc);
    virtual SHAREDVECTORN& get_generalized_velocity(DYNAMIC_BODY::GeneralizedCoordinateType gctype, SHAREDVECTORN& gv);
    virtual SHAREDVECTORN& get_generalized_acceleration(SHAREDVECTORN& ga);
    virtual void add_generalized_force(const SHAREDVECTORN& gf);
    virtual void add_generalized_force(const VECTORN& gf) { DYNAMIC_BODY::add_generalized_force(gf); }
    virtual void apply_generalized_impulse(const SHAREDVECTORN& gf);
    virtual void apply_generalized_impulse(const VECTORN& gj) { DYNAMIC_BODY::apply_generalized_impulse(gj); }
    virtual void set_generalized_coordinates(DYNAMIC_BODY::GeneralizedCoordinateType gctype, const SHAREDVECTORN& gc);
    virtual void set_generalized_coordinates(DYNAMIC_BODY::GeneralizedCoordinateType gctype, const VECTORN& gc) { DYNAMIC_BODY::set_generalized_coordinates(gctype, gc); }
    virtual void set_generalized_acceleration(const SHAREDVECTORN& ga);
    virtual void set_generalized_velocity(DYNAMIC_BODY::GeneralizedCoordinateType gctype, const SHAREDVECTORN& gv);
    virtual void set_generalized_velocity(DYNAMIC_BODY::GeneralizedCoordinateType gctype, const VECTORN& gv) { DYNAMIC_BODY::set_generalized_velocity(gctype, gv); }
    virtual SHAREDVECTORN& get_generalized_forces(SHAREDVECTORN& f);
    virtual VECTORN& get_generalized_forces(VECTORN& f) { return DYNAMIC_BODY::get_generalized_forces(f); }
    virtual SHAREDVECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf);
    virtual VECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLE_BODY> body, const SFORCE& w, VECTORN& gf) { return DYNAMIC_BODY::convert_to_generalized_force(body, w, gf); }
    virtual unsigned num_generalized_coordinates(DYNAMIC_BODY::GeneralizedCoordinateType gctype) const;
    virtual SHAREDMATRIXN& transpose_solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X);
    SHAREDMATRIXN& transpose_solve_generalized_inertia_single(const SHAREDMATRIXN& B, SHAREDMATRIXN& X);
    virtual SHAREDMATRIXN& solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X);
    virtual SHAREDVECTORN& solve_generalized_inertia(const SHAREDVECTORN& b, SHAREDVECTORN& x);
    boost::shared_ptr<RIGIDBODY> get_parent_link() const;
    boost::shared_ptr<JOINT> get_inner_joint_explicit() const;
    void add_inner_joint(boost::shared_ptr<JOINT> j);
    void add_outer_joint(boost::shared_ptr<JOINT> j);
    void remove_inner_joint(boost::shared_ptr<JOINT> joint);
    void remove_outer_joint(boost::shared_ptr<JOINT> joint);
    virtual REAL calc_kinetic_energy();
    virtual VECTOR3 calc_point_vel(const VECTOR3& p) const;
    bool is_base() const;
    bool is_ground() const;
    virtual boost::shared_ptr<const POSE3> get_computation_frame() const;
    virtual void set_computation_frame_type(ReferenceFrameType rftype);
    virtual MATRIXN& calc_jacobian(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);
    virtual MATRIXN& calc_jacobian_dot(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J);
    const SFORCE& sum_forces();
    void reset_accumulators();
    SFORCE calc_euler_torques();
    void set_enabled(bool flag);
    void invalidate_pose_vectors();

    template <class OutputIterator>
    OutputIterator get_parent_links(OutputIterator begin) const;

    template <class OutputIterator>
    OutputIterator get_child_links(OutputIterator begin) const;

    /// Gets the shared pointer for <b>this</b>
    boost::shared_ptr<RIGIDBODY> get_this() { return boost::dynamic_pointer_cast<RIGIDBODY>(shared_from_this()); }

    /// Gets the shared const pointer for <b>this</b>
    boost::shared_ptr<const RIGIDBODY> get_this() const { return boost::dynamic_pointer_cast<const RIGIDBODY>(shared_from_this()); }

    /// Gets the current pose of this body
    boost::shared_ptr<const POSE3> get_pose() const { return _F; }

    // Gets the pose used in generalized coordinates calculations
    boost::shared_ptr<const POSE3> get_gc_pose() const { return _F2; }


    // Gets the "mixed" pose of this body (pose origin at the body's reference point but pose aligned with global frame)
    boost::shared_ptr<const POSE3> get_mixed_pose() const { return _F2; }

    /// Synonym for get_mass() (implements SingleBody::calc_mass())
    REAL calc_mass() const { return _Jm.m; }

    /// Gets the mass of this body
    virtual REAL get_mass() const { return _Jm.m; }

    /// Gets whether this body is enabled
    bool is_enabled() const { return _enabled; }

    /// Gets the articulated body corresponding to this body
    /**
     * \return a pointer to the articulated body, or NULL if this body is not
     *         a link an articulated body
     */
    boost::shared_ptr<ARTICULATED_BODY> get_articulated_body() const { return (_abody.expired()) ? boost::shared_ptr<ARTICULATED_BODY>() : boost::shared_ptr<ARTICULATED_BODY>(_abody); }

    /// Sets the articulated body corresponding to this body
    /**
     * \param body a pointer to the articulated body or NULL if this body is
     *        not a link in an articulated body
     */
    virtual void set_articulated_body(boost::shared_ptr<ARTICULATED_BODY> body) { _abody = body; }

    /// Gets the number of child links of this link
    unsigned num_child_links() const { return _outer_joints.size(); }

    /// Gets whether this body is an end-effector (i.e., the number of child links is zero) in an articulated body
    bool is_end_effector() const { assert (!_abody.expired()); return !is_base() && _outer_joints.empty(); }

    /// Removes all inner joints from this link
    void clear_inner_joints() { _inner_joints.clear(); }

    /// Removes all outer joints from this link
    void clear_outer_joints() { _outer_joints.clear(); }

    /// Gets the link index (returns std::numeric_limits<unsigned>::max() if not set)
    unsigned get_index() const { return _link_idx; }

    /// Sets the link index
    /**
     * This is set automatically by the articulated body.  Users should not
     * change this index or unknown behavior will result.
     */
    void set_index(unsigned index) { _link_idx = index; }

    /// Gets the set of inner joints for this link
    const std::set<boost::shared_ptr<JOINT> >& get_inner_joints() const { return _inner_joints; }

    /// Gets the list of outer joints for this link
    const std::set<boost::shared_ptr<JOINT> >& get_outer_joints() const { return _outer_joints; }

    /// Compliance value, determines event type
    Compliance compliance;

  private:
    template <class V>
    void get_generalized_coordinates_generic(DYNAMIC_BODY::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void get_generalized_velocity_generic(DYNAMIC_BODY::GeneralizedCoordinateType gctype, V& gv);

    template <class V>
    void get_generalized_acceleration_generic(V& ga);

    template <class V>
    void set_generalized_coordinates_generic(DYNAMIC_BODY::GeneralizedCoordinateType gctype, V& gc);

    template <class V>
    void set_generalized_velocity_generic(DYNAMIC_BODY::GeneralizedCoordinateType gctype, V& gv);

    template <class V>
    void set_generalized_acceleration_generic(V& ga);

    void set_force(const SFORCE& w);
    void apply_generalized_impulse_single(const SHAREDVECTORN& gf);
    SHAREDMATRIXN& get_generalized_inertia_single(SHAREDMATRIXN& M);
    virtual SHAREDMATRIXN& get_generalized_inertia_inverse(SHAREDMATRIXN& M) const;
    SHAREDVECTORN& get_generalized_forces_single(SHAREDVECTORN& f);
    SHAREDVECTORN& convert_to_generalized_force_single(boost::shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf);
    unsigned num_generalized_coordinates_single(DYNAMIC_BODY::GeneralizedCoordinateType gctype) const;
    SHAREDMATRIXN& solve_generalized_inertia_single(const SHAREDMATRIXN& B, SHAREDMATRIXN& X);
    SHAREDVECTORN& solve_generalized_inertia_single(const SHAREDVECTORN& b, SHAREDVECTORN& x);
    boost::shared_ptr<RIGIDBODY> get_parent_link(boost::shared_ptr<JOINT> j) const;
    boost::shared_ptr<RIGIDBODY> get_child_link(boost::shared_ptr<JOINT> j) const;

    /// Indicates whether link frame velocity is valid (up-to-date)
    bool _xdi_valid;

    /// Indicates whether inner joint frame velocity is valid (up-to-date)
    bool _xdj_valid;

    /// Indicates whether inertial frame velocity is valid (up-to-date)
    bool _xdm_valid;

    /// Indicates whether global frame velocity is valid (up-to-date)
    bool _xd0_valid;

    /// Indicates whether link frame acceleration is valid (up-to-date)
    bool _xddi_valid;

    /// Indicates whether global frame acceleration is valid (up-to-date)
    bool _xdd0_valid;

    /// Indicates whether inertial frame acceleration is valid (up-to-date)
    bool _xddm_valid;

    /// Indicates whether link frame force is valid (up-to-date)
    bool _forcei_valid;

    /// Indicates whether inner joint frame force is valid (up-to-date)
    bool _forcej_valid;

    /// Indicates whether inertial frame force is valid (up-to-date)
    bool _forcem_valid;

    /// Indicates whether global frame force is valid (up-to-date)
    bool _force0_valid;

    /// Indicates whether the global frame inertia matrix is valid
    bool _J0_valid;

    /// Indicates whether the link frame inertia matrix is valid
    bool _Ji_valid;

    /// Indicates whether the link com frame inertia matrix is valid
    bool _Jcom_valid;

    /// Indicates whether the inner joint frame inertia matix is valid
    bool _Jj_valid;

    /// Spatial rigid body inertia matrix (global frame)
    SPATIAL_RB_INERTIA _J0;

    /// Velocity (global frame)
    SVELOCITY _xd0;

    /// Acceleration (global frame)
    SACCEL _xdd0;

    /// Cumulative force on the body (global frame)
    SFORCE _force0;

    /// Spatial rigid body inertia matrix (inertial frame)
    SPATIAL_RB_INERTIA _Jm;

    /// Velocity (inertial frame)
    SVELOCITY _xdm;

    /// Acceleration (inertial frame)
    SACCEL _xddm;

    /// Cumulative force on the body (inertial frame)
    SFORCE _forcem;

    /// Spatial rigid body inertia matrix (link frame)
    SPATIAL_RB_INERTIA _Ji;

    /// Velocity (link frame)
    SVELOCITY _xdi;

    /// Acceleration (link frame)
    SACCEL _xddi;

    /// Cumulative force on the body (link frame)
    SFORCE _forcei;

    /// Spatial rigid body inertia matrix (link COM frame)
    SPATIAL_RB_INERTIA _Jcom;

    /// Velocity (link com frame)
    SVELOCITY _xdcom;

    /// Acceleration (link com frame)
    SACCEL _xddcom;

    /// Cumulative force on the body (com frame)
    SFORCE _forcecom;

    /// The link index (if a link in an articulated body)
    unsigned _link_idx;

    /// Flag for determining whether or not the body is physically enabled
    bool _enabled;

  protected:
    void update_mixed_pose();

    /// Indicates whether inner joint frame acceleration is valid (up-to-date)
    bool _xddj_valid;

    /// Spatial rigid body inertia matrix (inner joint frame)
    SPATIAL_RB_INERTIA _Jj;

    /// Velocity (inner joint frame)
    SVELOCITY _xdj;

    /// Cumulative force on the body (inner joint frame)
    SFORCE _forcej;

    /// Acceleration (inner joint frame)
    SACCEL _xddj;

    /// reference pose for this body
    boost::shared_ptr<POSE3> _F;

    /// secondary pose for this body
    boost::shared_ptr<POSE3> _F2;

    /// inertial pose for this body
    boost::shared_ptr<POSE3> _jF;

    /// Pointer to articulated body (if this body is a link)
    boost::weak_ptr<ARTICULATED_BODY> _abody;

    /// Inner joints and associated data
    std::set<boost::shared_ptr<JOINT> > _inner_joints;

    /// Outer joints and associated data
    std::set<boost::shared_ptr<JOINT> > _outer_joints;

}; // end class

// incline inline functions
#include "RigidBody.inl"


