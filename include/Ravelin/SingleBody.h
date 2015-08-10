/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SINGLEBODY 
#error This class is not to be included by the user directly. Use SingleBodyd.h or SingleBodyf.h instead.
#endif

class ARTICULATEDBODY;

/// Superclass for both rigid and deformable bodies 
class SINGLEBODY : public DYNAMICBODY 
{
  public:
    virtual ~SINGLEBODY() {}
    virtual boost::shared_ptr<DYNAMICBODY> get_super_body() const;

    /// Gets the computation frame for the body
    virtual boost::shared_ptr<const POSE3> get_computation_frame() const = 0;

    /// Gets the mass of the body (for gravity calculation)
    virtual REAL get_mass() const = 0;

    /// Gets the pose of the body
    virtual boost::shared_ptr<const POSE3> get_pose() const = 0;

    /// Gets the acceleration of the body
    virtual const SACCEL& get_accel() = 0;

    /// Gets the velocity of the body
    virtual const SVELOCITY& get_velocity() = 0;

    /// Applies an impulse at a point on the body
    virtual void apply_impulse(const SMOMENTUM& w) = 0;

    /// Calculates the mass of the body
    virtual REAL calc_mass() const = 0;

    /// Gets the articulated body that this body is a part of (if any)
    virtual boost::shared_ptr<ARTICULATEDBODY> get_articulated_body() const = 0;

    /// Determines whether the body is enabled
    virtual bool is_enabled() const = 0;

    /// Calculates the velocity at a point on the body *in the body frame*
    virtual VECTOR3 calc_point_vel(const VECTOR3& point) const = 0;

    /// Calculates the velocity at a point on the body in a given direction
    REAL calc_point_vel(const VECTOR3& point, const VECTOR3& dir);

}; // end class

