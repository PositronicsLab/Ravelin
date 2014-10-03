/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef UNIVERSALJOINT
#error This class is not to be included by the user directly. Use UniversalJointd.h or UniversalJointf.h instead.
#endif

/// Defines a joint for purely rotational motion 
class UNIVERSALJOINT : public virtual JOINT
{
  public:
    enum Axis { eAxis1, eAxis2 };
    UNIVERSALJOINT();
    VECTOR3 get_axis(Axis a) const;
    void set_axis(const VECTOR3& axis, Axis a);
    virtual void update_spatial_axes();    
    virtual void determine_q(VECTORN& q);
    virtual boost::shared_ptr<const POSE3> get_induced_pose();
    virtual const std::vector<SVELOCITY>& get_spatial_axes();
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot();
    virtual unsigned num_dof() const { return 2; }
    virtual void evaluate_constraints(REAL C[]);

    /// Universal joint is never singular 
    virtual bool is_singular_config() const { return false; }

    /// The tolerance to which a joint configuration is considered singular
    /**
     * \note if this tolerance is too low, then dynamics computation may
     * become unstable; if the tolerance is too high, then dynamics computation
     * will be slower.  A safe value is 1e-2.
     */
    REAL SINGULAR_TOL;

  protected:
    bool assign_axes();
    static bool rel_equal(REAL x, REAL y);
    MATRIX3 get_rotation() const;

    /// The local joint axes
    VECTOR3 _u[2];

    /// The second joint axis in outboard pose frame
    VECTOR3 _h2;

    /// The derivative of the spatial axis
    std::vector<SVELOCITY> _s_dot;

//    virtual void calc_constraint_jacobian(RigidBodyPtr, unsigned index, REAL Cq[7]);
//    virtual void calc_constraint_jacobian_dot(RigidBodyPtr, unsigned index, REAL Cq[7]);
    void setup_joint();

}; // end class
