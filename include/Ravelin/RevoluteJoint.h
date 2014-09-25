/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef REVOLUTEJOINT
#error This class is not to be included by the user directly. Use RevoluteJointd.h or RevoluteJointf.h instead.
#endif

/// Defines a rotary joint 
class REVOLUTEJOINT : public virtual JOINT
{
  public:
    REVOLUTEJOINT();
    virtual void update_spatial_axes();    
    virtual void determine_q(VECTORN& q);
    virtual boost::shared_ptr<const POSE3> get_induced_pose();
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot();
    virtual unsigned num_dof() const { return 1; }
    virtual void evaluate_constraints(REAL C[]);
    VECTOR3 get_axis() const;
    void set_axis(const VECTOR3& axis);

    /// Revolute joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

  protected:
//    virtual void calc_constraint_jacobian(RigidBodyPtr, unsigned index, REAL Cq[7]);
//    virtual void calc_constraint_jacobian_dot(RigidBodyPtr, unsigned index, REAL Cq[7]);

    /// The joint axis (defined in inner relative pose coordinates)
    VECTOR3 _u;

    /// Two unit vectors that make a orthonormal basis with _u
    VECTOR3 _ui, _uj;

    /// The joint axis (defined in outer relative pose coordinates)
    VECTOR3 _v2;

    /// The time derivative of the spatial axis -- should be zero vector 6x1
    std::vector<SVELOCITY> _s_dot;
}; // end class

