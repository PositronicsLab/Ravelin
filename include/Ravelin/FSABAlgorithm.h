/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef FSAB_ALGORITHM 
#error This class is not to be included by the user directly. Use FSABAlgorithmd.h or FSABAlgorithmf.h instead.
#endif

class RC_ARTICULATED_BODY;

/// Implements Featherstone's algorithm for forward dynamics
/**
 * Implements Featherstone's algorithm for articulated bodies.  Featherstone's 
 * algorithm runs in O(n) time [n = # of joints].  This implementation is based
 * on Brian Mirtich's Ph. D. thesis, and remains pretty consistent with it. 
 * There are a couple of changes, to produce a nice implementation.  The user 
 * need not be concerned with these issues, but they are useful to know for 
 * debugging.
 * <ol>
 * <li>Mirtich labels his links from 1..n, and considers the base to be link 0; 
 * the total number of links is considered to be n, rather than n+1.  I make 
 * the total number of links n+1 and treat the links the same as the base.  I 
 * do this so that the user thinks of the base as a link for purposes of link 
 * connectivity.</li>
 * <li>Mirtich labels his joints from 0..n-1</li>.  When labeling the link in 
 * Mirtich's style, link i and joint i match up (joint i is link i's inner 
 * joint).  When labeling the link in my style, joint i-1 is the corresponding 
 * joint for link i.</li>
 * </ol>
 * Note that one critical note for manipulator setup is that the base is the 
 * first link in the list of links.
 */
class FSAB_ALGORITHM 
{
  friend class RC_ARTICULATED_BODY;

  public:
    FSAB_ALGORITHM();
    ~FSAB_ALGORITHM() {}
    boost::shared_ptr<RC_ARTICULATED_BODY> get_body() const { return boost::shared_ptr<RC_ARTICULATED_BODY>(_body); }
    void set_body(boost::shared_ptr<RC_ARTICULATED_BODY> body) { _body = body; }
    void calc_fwd_dyn();
    void calc_inverse_generalized_inertia_noprecalc(MATRIXN& iM);
    void solve_generalized_inertia_noprecalc(SHAREDVECTORN& v);
    void solve_generalized_inertia_noprecalc(SHAREDMATRIXN& Y);
    void apply_generalized_impulse(const VECTORN& gj);
    void apply_impulse(const SMOMENTUM& j, boost::shared_ptr<RIGIDBODY> link);
    void calc_spatial_inertias(boost::shared_ptr<RC_ARTICULATED_BODY> body);

    /// The body that this algorithm operates on
    boost::weak_ptr<RC_ARTICULATED_BODY> _body;

    /// The spatial accelerations
    std::vector<SACCEL> _a;

    /// The articulated body inertias
    std::vector<SPATIAL_AB_INERTIA> _I;

    /// The articulated body spatial zero accelerations
    std::vector<SFORCE> _Z;

    /// Vector of link velocity updates
    std::vector<SVELOCITY> _dv;

    /// The spatial coriolis vectors
    std::vector<SACCEL> _c;

    /// The expressions I*s
    std::vector<std::vector<SMOMENTUM> > _Is;

    /// Cholesky factorizations sIs
    std::vector<MATRIXN> _sIs;

    /// SVDs of sIs
    std::vector<MATRIXN> _usIs, _vsIs;
    std::vector<VECTORN> _ssIs;

    /// Determines whether the equations for a joint are rank deficient 
    std::vector<bool> _rank_deficient;

    /// The temporary expression Q - I*s'*c - s'*Z
    std::vector<VECTORN> _mu;

  private:
    // pointer to the linear algebra routines
    boost::shared_ptr<LINALG> _LA;

    /// work variables 
    VECTORN _workv, _workv2, _sTY, _qd_delta, _sIsmu, _Qi, _Q;
    MATRIXN _sIss, _workM;
    std::vector<SMOMENTUM> _Y;

    /// processed vector
    std::vector<bool> _processed;

    void calc_fwd_dyn_special();
    static REAL sgn(REAL x);
    static void push_children(boost::shared_ptr<RIGIDBODY> link, std::queue<boost::shared_ptr<RIGIDBODY> >& q);
    void apply_coulomb_joint_friction(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void apply_generalized_impulse(unsigned index, VECTORN& vgj);
    void set_spatial_velocities(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_spatial_accelerations(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_spatial_zero_accelerations(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    void calc_spatial_coriolis_vectors(boost::shared_ptr<RC_ARTICULATED_BODY> body);
    VECTORN& solve_sIs(unsigned idx, const VECTORN& v, VECTORN& result) const;
    MATRIXN& solve_sIs(unsigned idx, const MATRIXN& v, MATRIXN& result) const;
    MATRIXN& transpose_solve_sIs(unsigned idx, const std::vector<SVELOCITY>& m, MATRIXN& result) const;
}; // end class


