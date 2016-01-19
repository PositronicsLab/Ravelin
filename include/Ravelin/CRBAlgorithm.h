/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef CRB_ALGORITHM 
#error This class is not to be included by the user directly. Use CRBAlgorithmd.h or CRBAlgorithmf.h instead.
#endif

/// Computes forward dynamics using composite-rigid body method
class CRB_ALGORITHM
{
  friend class RC_ARTICULATED_BODY;

  public:
    CRB_ALGORITHM();
    ~CRB_ALGORITHM() {}
    boost::shared_ptr<RC_ARTICULATED_BODY> get_body() const { return boost::shared_ptr<RC_ARTICULATED_BODY>(_body); }
    void set_body(boost::shared_ptr<RC_ARTICULATED_BODY> body) { _body = body; setup_parent_array(); }
    void calc_fwd_dyn();
    void apply_impulse(const SMOMENTUM& w, boost::shared_ptr<RIGIDBODY> link);
    void calc_generalized_inertia(SHAREDMATRIXN& M);
    void calc_generalized_inertia(SHAREDMATRIXN& M, boost::shared_ptr<const POSE3> P);
    void calc_generalized_forces(SFORCE& f0, VECTORN& C);
    void calc_generalized_forces_noinertial(SFORCE& f0, VECTORN& C);
    bool factorize_cholesky(MATRIXN& M);
    VECTORN& M_solve(VECTORN& xb);
    SHAREDVECTORN& M_solve(SHAREDVECTORN& xb);
    MATRIXN& M_solve(MATRIXN& XB);
    SHAREDMATRIXN& M_solve(SHAREDMATRIXN& XB);

  private:
    void calc_fwd_dyn_special();
    static boost::shared_ptr<const POSE3> get_computation_frame(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    std::vector<unsigned> _lambda;
    void setup_parent_array();

    /// The body that this algorithm operates on
    boost::weak_ptr<RC_ARTICULATED_BODY> _body;

    /// The spatial acceleration of the base computed on the last call to calc_fwd_dyn()
    SACCEL _a0;

    /// The vector of joint accelerations computed on the last call to calc_fwd_dyn()
    VECTORN _qdd;

    /// The joint space inertia matrix H (fixed base) or augmented matrix [I_0^c K; K^s H] (floating base, see [Featherstone 1987], p. 123) used to compute forward dynamics for floating bases
    MATRIXN _M;

    /// A factorization (or possibly inverse) of the matrix M; note that we compute this b/c we generally may need to solve multiple systems of linear equations using this matrix as a LHS at different times -- always in global frame
    MATRIXN _fM;

    /// Determines whether the system of equations for forward dynamics is rank-deficient
     bool _rank_deficient;

    void calc_joint_space_inertia(boost::shared_ptr<RC_ARTICULATED_BODY> body, MATRIXN& H, std::vector<SPATIAL_RB_INERTIA>& Ic);
    void apply_coulomb_joint_friction(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void precalc(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_generalized_inertia(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_fwd_dyn_fixed_base(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_fwd_dyn_floating_base(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void update_link_accelerations(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    static void to_spatial7_inertia(const SPATIAL_RB_INERTIA& I, const QUAT& q, MATRIXN& I7);
    VECTORN& M_solve_noprecalc(VECTORN& xb);
    MATRIXN& M_solve_noprecalc(MATRIXN& XB);
    SHAREDVECTORN& M_solve_noprecalc(SHAREDVECTORN& xb);
    SHAREDMATRIXN& M_solve_noprecalc(SHAREDMATRIXN& XB);
    void transform_and_mult(boost::shared_ptr<const POSE3> P, const SPATIAL_RB_INERTIA& I, const std::vector<SVELOCITY>& s, std::vector<SMOMENTUM>& Is);

  private:
    // temporaries for transform_and_transpose_mult() functions
    std::vector<SFORCE> _tandt_fx;
    std::vector<SMOMENTUM> _tandt_wx, _Isprime;
    std::vector<SVELOCITY> _tandt_tx;

    // temporary for calc_fwd_dyn() 
    std::vector<SACCEL> _a;

    // temporary for calc_generalized_forces() 
    std::vector<SFORCE> _w;

    // temporary spatial axes
    std::vector<SVELOCITY> _sprime;

    // temporaries for solving and linear algebra
    boost::shared_ptr<LINALG> _LA;
    MATRIXN _uM, _vM;
    VECTORN _sM;

    // temporaries for calc_generalized_inertia()
    MATRIXN _H;
    VECTORN _rowi;
    std::vector<SPATIAL_RB_INERTIA> _Ic;
    std::vector<SMOMENTUM> _Is;

    // temporaries for calc_joint_space_inertia()
    MATRIXN _workM, _sub;
    std::vector<std::vector<bool> > _supports;
    std::vector<std::vector<SMOMENTUM> > _momenta;

    // temporaries for calc_fwd_dyn_fixed_base(), calc_fwd_dyn_floating_base()
    VECTORN _C, _Q, _b, _augV;

    // temporaries for applying impulse
    VECTORN _workv;
    std::vector<SVELOCITY> _J;

    // precalc
    VECTORN _gc_last;

    #include "CRBAlgorithm.inl"
}; // end class

