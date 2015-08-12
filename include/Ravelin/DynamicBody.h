/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef DYNAMIC_BODY 
#error This class is not to be included by the user directly. Use DynamicBodyd.h or DynamicBodyf.h instead.
#endif

class SINGLE_BODY;

/// Superclass for deformable bodies and single and multi-rigid bodies  
class DYNAMIC_BODY : public virtual boost::enable_shared_from_this<DYNAMIC_BODY>
{
  public:
    enum GeneralizedCoordinateType { eEuler, eSpatial };

    virtual ~DYNAMIC_BODY() {}

    /// The identifier for this body
    std::string body_id;

    /// The Jacobian transforms from the generalized coordinate from to the given frame
    virtual MATRIXN& calc_jacobian(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J) = 0;
    virtual MATRIXN& calc_jacobian_dot(boost::shared_ptr<const POSE3> source_pose, boost::shared_ptr<const POSE3> target_pose, boost::shared_ptr<DYNAMIC_BODY> body, MATRIXN& J) = 0;

    /// Validates position-based variables (potentially dangerous for a user to call)
    virtual void validate_position_variables() { };

    /// Validates velocity-based variables (potentially dangerous for a user to call)
    virtual void validate_velocity_variables() { };

    /// Sets the computation frame type for this body
    virtual void set_computation_frame_type(ReferenceFrameType rftype) = 0;

    /// Gets the computation frame type for this body
    ReferenceFrameType get_computation_frame_type() const { return _rftype; }

    /// Forces a recalculation of forward dynamics
    virtual void calc_fwd_dyn() = 0;

    /// Resets the force and torque accumulators on the dynamic body
    virtual void reset_accumulators() = 0;

    /// Rotates the dynamic body by the given orientation 
    virtual void rotate(const QUAT& q) = 0;

    /// Translates the dynamic body by the given translation 
    virtual void translate(const ORIGIN3& o) = 0;

    /// Calculates the kinetic energy of the body in an arbitrary frame
    virtual REAL calc_kinetic_energy() = 0;

    /// Gets the frame for generalized coordinates
    virtual boost::shared_ptr<const POSE3> get_gc_pose() const = 0;

    /// Gets the number of generalized coordinates
    virtual unsigned num_generalized_coordinates(GeneralizedCoordinateType gctype) const = 0;

    /// Sets the generalized forces on the body
    virtual void set_generalized_forces(const SHAREDVECTORN& gf) = 0;

    /// Sets the generalized forces on the body
    virtual void set_generalized_forces(const VECTORN& gf)
    {
      const SHAREDVECTORN gf_shared = gf.segment(0, gf.size()).get();
      set_generalized_forces(gf_shared);
    }

    /// Adds a generalized force to the body
    virtual void add_generalized_force(const SHAREDVECTORN& gf) = 0;

    /// Adds a generalized force to the body
    virtual void add_generalized_force(const VECTORN& gf)
    {
      const SHAREDVECTORN gf_shared = gf.segment(0, gf.size()).get();
      add_generalized_force(gf_shared);
    }

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(const SHAREDVECTORN& gj) = 0;

    /// Applies a generalized impulse to the body
    virtual void apply_generalized_impulse(const VECTORN& gj)
    {
      const SHAREDVECTORN gj_shared = gj.segment(0, gj.size()).get();
      apply_generalized_impulse(gj_shared);
    }

    /// Gets the generalized coordinates of this body
    virtual SHAREDVECTORN& get_generalized_coordinates(GeneralizedCoordinateType gctype, SHAREDVECTORN& gc) = 0;

    /// Gets the generalized coordinates of this body
    virtual VECTORN& get_generalized_coordinates(GeneralizedCoordinateType gctype, VECTORN& gc)
    {
      const unsigned NGC = num_generalized_coordinates(gctype);
      gc.resize(NGC);
      SHAREDVECTORN gc_shared = gc.segment(0, gc.size());
      get_generalized_coordinates(gctype, gc_shared);
      return gc;
    }

    /// Gets the generalized velocity of this body
    virtual SHAREDVECTORN& get_generalized_velocity(GeneralizedCoordinateType gctype, SHAREDVECTORN& gv) = 0;

    /// Gets the generalized velocity of this body
    virtual VECTORN& get_generalized_velocity(GeneralizedCoordinateType gctype, VECTORN& gv)
    {
      const unsigned NGC = num_generalized_coordinates(gctype);
      gv.resize(NGC);
      SHAREDVECTORN gv_shared = gv.segment(0, gv.size());
      get_generalized_velocity(gctype, gv_shared);
      return gv;
    }

    /// Sets the generalized velocity of this body
    virtual void set_generalized_acceleration(const SHAREDVECTORN& ga) = 0;

    /// Sets the generalized acceleration of this body
    virtual void set_generalized_acceleration(const VECTORN& ga)
    {
      const SHAREDVECTORN ga_shared = ga.segment(0, ga.size()).get();
      set_generalized_acceleration(ga_shared);
    }

    /// Gets the generalized velocity of this body
    virtual SHAREDVECTORN& get_generalized_acceleration(SHAREDVECTORN& ga) = 0;

    /// Gets the generalized acceleration of this body
    virtual VECTORN& get_generalized_acceleration(VECTORN& ga) 
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      ga.resize(NGC);
      SHAREDVECTORN ga_shared = ga.segment(0, ga.size());
      get_generalized_acceleration(ga_shared);
      return ga; 
    }

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const SHAREDVECTORN& gc) = 0;

    /// Sets the generalized coordinates of this body
    virtual void set_generalized_coordinates(GeneralizedCoordinateType gctype, const VECTORN& gc) 
    {
      const SHAREDVECTORN gc_shared = gc.segment(0, gc.size()).get();
      set_generalized_coordinates(gctype, gc_shared);
    }

    /// Sets the generalized velocity of this body
    /**
      * \param gv the generalized velocity
      * \note uses the current generalized coordinates
      */
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const SHAREDVECTORN& gv) = 0;

    /// Sets the generalized velocity of this body
    virtual void set_generalized_velocity(GeneralizedCoordinateType gctype, const VECTORN& gv) 
    {
      const SHAREDVECTORN gv_shared = gv.segment(0, gv.size()).get();
      set_generalized_velocity(gctype, gv_shared);
    }

    /// Gets the generalized inertia of this body
    MATRIXN& get_generalized_inertia(MATRIXN& M)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      M.resize(NGC, NGC);
      SHAREDMATRIXN X = M.block(0, NGC, 0, NGC);
      get_generalized_inertia(X);
      return M;
    }

    /// Gets the generalized inertia of this body
    virtual SHAREDMATRIXN& get_generalized_inertia(SHAREDMATRIXN& M) = 0;

    /// Solves using the inverse generalized inertia
    virtual SHAREDMATRIXN& solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X) = 0;

    /// Solves using the inverse generalized inertia
    virtual MATRIXN& solve_generalized_inertia(const MATRIXN& B, MATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDMATRIXN Bx = B.block(0, B.rows(), 0, B.columns()).get();
      X.resize(NGC, B.columns());
      SHAREDMATRIXN Xx = X.block(0, X.rows(), 0, X.columns());
      solve_generalized_inertia(Bx, Xx);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual MATRIXN& solve_generalized_inertia(const SHAREDMATRIXN& B, MATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      X.resize(NGC, B.columns());
      SHAREDMATRIXN Xx = X.block(0, X.rows(), 0, X.columns());
      solve_generalized_inertia(B, Xx);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual SHAREDMATRIXN& solve_generalized_inertia(const MATRIXN& B, SHAREDMATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDMATRIXN Bx = B.block(0, B.rows(), 0, B.columns()).get();
      solve_generalized_inertia(Bx, X);
      return X;
    }

    /// Solves using the inverse generalized inertia
    virtual VECTORN& solve_generalized_inertia(const VECTORN& b, VECTORN& x)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDVECTORN bx = b.segment(0, b.rows()).get();
      x.resize(NGC);
      SHAREDVECTORN xx = x.segment(0, x.rows());
      solve_generalized_inertia(bx, xx);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual VECTORN& solve_generalized_inertia(const SHAREDVECTORN& b, VECTORN& x)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      x.resize(NGC);
      SHAREDVECTORN xx = x.segment(0, x.rows());
      solve_generalized_inertia(b, xx);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual SHAREDVECTORN& solve_generalized_inertia(const VECTORN& b, SHAREDVECTORN& x)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDVECTORN bx = b.segment(0, b.rows()).get();
      solve_generalized_inertia(bx, x);
      return x;
    }

    /// Solves using the inverse generalized inertia
    virtual SHAREDVECTORN& solve_generalized_inertia(const SHAREDVECTORN& b, SHAREDVECTORN& x) = 0;

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual MATRIXN& transpose_solve_generalized_inertia(const MATRIXN& B, MATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDMATRIXN Bx = B.block(0, B.rows(), 0, B.columns()).get();
      X.resize(NGC, B.rows());
      SHAREDMATRIXN Xx = X.block(0, X.rows(), 0, X.columns());
      transpose_solve_generalized_inertia(Bx, Xx);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual MATRIXN& transpose_solve_generalized_inertia(const SHAREDMATRIXN& B, MATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      X.resize(NGC, B.rows());
      SHAREDMATRIXN Xx = X.block(0, X.rows(), 0, X.columns());
      transpose_solve_generalized_inertia(B, Xx);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual SHAREDMATRIXN& transpose_solve_generalized_inertia(const MATRIXN& B, SHAREDMATRIXN& X)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      const SHAREDMATRIXN Bx = B.block(0, B.rows(), 0, B.columns()).get();
      transpose_solve_generalized_inertia(Bx, X);
      return X;
    }

    /// Solves the transpose matrix using the inverse generalized inertia
    virtual SHAREDMATRIXN& transpose_solve_generalized_inertia(const SHAREDMATRIXN& B, SHAREDMATRIXN& X) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual SHAREDVECTORN& get_generalized_forces(SHAREDVECTORN& f) = 0;

    /// Gets the external forces on this body
    /**
     * \note uses the current generalized coordinates
     */
    virtual VECTORN& get_generalized_forces(VECTORN& f)
    {
      f.resize(num_generalized_coordinates(DYNAMIC_BODY::eSpatial));
      SHAREDVECTORN f_shared = f.segment(0, f.size());
      get_generalized_forces(f_shared);
      return f;
    }

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the force 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual SHAREDVECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLE_BODY> body, const SFORCE& w, SHAREDVECTORN& gf) = 0;

    /// Converts a force to a generalized force
    /**
     * \param body the actual rigid body to which the force/torque is applied 
     *               (at the center-of-mass)
     * \param w the force 
     * \param gf the generalized force, on return
     * \note uses the current generalized coordinates
     */
    virtual VECTORN& convert_to_generalized_force(boost::shared_ptr<SINGLE_BODY> body, const SFORCE& w, VECTORN& gf)
    {
      const unsigned NGC = num_generalized_coordinates(DYNAMIC_BODY::eSpatial);
      gf.resize(NGC);
      SHAREDVECTORN gf_shared = gf.segment(0, gf.size());
      convert_to_generalized_force(body, w, gf_shared);
      return gf;
    }

  protected:

    /// The computation frame type
    ReferenceFrameType _rftype;

    /// Temporaries for use with integration
    VECTORN gc, gv, gcgv, xp, xv, xa;

}; // end class


