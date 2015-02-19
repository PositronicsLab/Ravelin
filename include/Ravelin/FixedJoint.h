/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef FIXEDJOINT
#error This class is not to be included by the user directly. Use FixedJointd.h or FixedJointf.h instead.
#endif

/// Defines a joint for fixing two bodies together or fixing one body to the ground
class FIXEDJOINT : public virtual JOINT
{
  public:
    FIXEDJOINT();
    virtual void update_spatial_axes();    
    virtual void determine_q(VECTORN& q) { }
    virtual boost::shared_ptr<const POSE3> get_induced_pose();
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot();
    virtual unsigned num_dof() const { return 0; }
    virtual void evaluate_constraints(REAL C[]);
    virtual void set_inboard_pose(boost::shared_ptr<const POSE3> inboard_pose, bool update_joint_pose);
    virtual void set_outboard_pose(boost::shared_ptr<POSE3> outboard_pose, bool update_joint_pose);

    /// Fixed joint can never be in a singular configuration
    virtual bool is_singular_config() const { return false; }

  private:
/*
    virtual void calc_constraint_jacobian(RigidBodyPtr, unsigned index, REAL Cq[7]);
    virtual void calc_constraint_jacobian_dot(RigidBodyPtr, unsigned index, REAL Cq[7]);
*/
    void setup_joint();

    /// The relative transform from the inboard link to the outboard link
    boost::shared_ptr<POSE3> _T;

    /// The orientation constant that we attempt to maintain
    VECTOR3 _rconst;

    /// The vector from the inner link to the outer link in inner link frame
    VECTOR3 _ui;

    /// The time derivative of the spatial axis -- should be zero vector 6x1
    std::vector<SVELOCITY> _s_dot;

    /// Temporaries for absolute coordinate calculations
    boost::shared_ptr<POSE3> _F1, _F2;
}; // end class

