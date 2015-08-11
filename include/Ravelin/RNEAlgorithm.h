/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef RNE_ALGORITHM 
#error This class is not to be included by the user directly. Use RNEAlgorithmd.h or RNEAlgorithmf.h instead.
#endif

/// Wrapper class for passing data to and from inverse dynamics algorithms
struct RCArticulatedBodyInvDynData
{
  public:

    /// External force applied to this link 
    SFORCE wext;

    /// Inner joint acceleration (desired)
    VECTORN qdd;
};
  
/// Implementation of the Recursive Newton-Euler algorithm for inverse dynamics
/**
 * Algorithm taken from Featherstone, 1987.
 */ 
class RNE_ALGORITHM
{
  public:
    std::map<boost::shared_ptr<JOINT>, VECTORN> calc_inv_dyn(boost::shared_ptr<RC_ARTICULATED_BODY> body, const std::map<boost::shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data);
    void calc_constraint_forces(boost::shared_ptr<RC_ARTICULATED_BODY> body);

  private:
    std::map<boost::shared_ptr<JOINT>, VECTORN> calc_inv_dyn_fixed_base(boost::shared_ptr<RC_ARTICULATED_BODY> body, const std::map<boost::shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data) const;
    std::map<boost::shared_ptr<JOINT>, VECTORN> calc_inv_dyn_floating_base(boost::shared_ptr<RC_ARTICULATED_BODY> body, const std::map<boost::shared_ptr<RIGIDBODY>, RCArticulatedBodyInvDynData>& inv_dyn_data) const;
};


