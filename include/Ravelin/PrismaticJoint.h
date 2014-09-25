/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef PRISMATICJOINT
#error This class is not to be included by the user directly. Use PrismaticJointd.h or PrismaticJointf.h instead.
#endif

/// Defines a sliding joint 
class PRISMATICJOINT : public virtual JOINT
{
  public:
    PRISMATICJOINT();
    virtual void update_spatial_axes();    
    virtual void determine_q(VECTORN& q);
    virtual boost::shared_ptr<const POSE3> get_induced_pose();
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot();
    virtual unsigned num_dof() const { return 1; }
    virtual void evaluate_constraints(REAL C[]);
    VECTOR3 get_axis() const { return _u; }
    void set_axis(const VECTOR3& axis);

    /// Prismatic joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

  protected:
//    virtual void calc_constraint_jacobian(RigidBodyPtr, unsigned index, REAL Cq[7]);
//    virtual void calc_constraint_jacobian_dot(RigidBodyPtr, unsigned index, REAL Cq[7]);

    /// The axis of the joint (inboard pose frame)
    VECTOR3 _u;

    /// The vector from the inboard pose to the outboard pose in inboard pose frame
    VECTOR3 _ui;

    /// Vector attached to outboard pose and initially orthogonal to joint axis; vector specified in outer pose frame
    VECTOR3 _uj;

    /// The joint axis on the outboard pose; vector specified in outboard pose frame 
    VECTOR3 _v2;

    /// The time derivative of the spatial axis -- should be zero vector 6x1
    std::vector<SVELOCITY> _s_dot;
}; // end class

