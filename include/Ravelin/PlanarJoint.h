/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef PLANARJOINT
#error This class is not to be included by the user directly. Use PlanarJointd.h or PlanarJointf.h instead.
#endif

/// Defines a joint that constrains motion to a plane 
class PLANARJOINT : public virtual JOINT
{
  public:
    PLANARJOINT();
    virtual void update_spatial_axes();    
    virtual void determine_q(VECTORN& q);
    virtual boost::shared_ptr<const POSE3> get_induced_pose();
    virtual unsigned num_dof() const { return 3; }
    virtual void evaluate_constraints(REAL C[]);
    virtual void calc_constraint_jacobian(bool inboard, MATRIXN& Cq);
    virtual void calc_constraint_jacobian_dot(bool inboard, MATRIXN& Cq);
    virtual bool is_singular_config() const { return false; }
    void set_normal(const VECTOR3& normal);
    virtual const std::vector<SVELOCITY>& get_spatial_axes_dot() { return _s_dot; }

  protected:
    void update_offset();

    /// Vectors orthogonal to the normal vector in the outboard link frame
    VECTOR3 _vi, _vj;

    /// The plane normal and tangents (global frame)
    VECTOR3 _normal, _tan1, _tan2;

    /// The plane offset such that n'*x = offset
    REAL _offset;

    /// The derivative of the spatial axis
    std::vector<SVELOCITY> _s_dot;

}; // end class

