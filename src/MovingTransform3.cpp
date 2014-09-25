/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;

/// Sets up a default moving transform
MOVINGTRANSFORM3::MOVINGTRANSFORM3()
{
  source = boost::shared_ptr<const POSE3>();
  target = boost::shared_ptr<const POSE3>();
  r.set_zero();
  E.set_identity();
  rdot.set_zero();
  Edot.set_zero();
  v.set_zero();
}

/// Copies a moving transform to another
MOVINGTRANSFORM3& MOVINGTRANSFORM3::operator=(const MOVINGTRANSFORM3& source)
{
  this->source = source.source;
  this->target = source.target;
  this->r = source.r;
  this->E = source.E;
  this->rdot = source.rdot;
  this->Edot = source.Edot;
  this->v = source.v;
  return *this;
}

/// Transforms a spatial acceleration (defined relative to the source frame's relative pose) to a spatial acceleration (relative to the target frame's relative pose)
/**
 * \param a an acceleration defined in the source frame
 */
SACCEL MOVINGTRANSFORM3::transform(const SACCEL& a) const
{
  #ifndef NEXCEPT
  if (a.pose != source)
    throw FrameException();
  #endif

  // setup r and rdot as vectors
  VECTOR3 rv(r, source);
  VECTOR3 rdotv(rdot, source);

  // get the components of a 
  VECTOR3 atop = a.get_upper();
  VECTOR3 abot = a.get_lower();

  // do the calculations for X*a
  VECTOR3 Etop(E * ORIGIN3(atop), target);
  VECTOR3 cross = VECTOR3::cross(rv, atop);
  VECTOR3 Xa_bot(E * ORIGIN3(abot - cross), target);

  // get the components of v
  VECTOR3 vtop = v.get_upper();
  VECTOR3 vbot = v.get_lower();

  // now do the calculations for dX/dt*a
  // d/dt | E    0 |  =  | d/dt E                0      | 
  //      | -Erx E |     | -d/dt E*rx - E*rdotx  d/dt E |
  VECTOR3 Edottop(Edot * ORIGIN3(vtop), target);
  VECTOR3 cross1 = VECTOR3::cross(rv, vtop);
  VECTOR3 cross2 = VECTOR3::cross(rdotv, vtop);
  VECTOR3 Xdotv_bot(-Edot*ORIGIN3(cross1 - vbot) - E*ORIGIN3(cross2), target); 

  // finally, set the results
  SACCEL result;
  result.set_upper(Etop + Edottop);
  result.set_lower(Xa_bot + Xdotv_bot);

  // set the pose
  result.pose = target;

  return result;
}

/// Computes the relative transformation between two moving poses 
/**
 * \param source the source frame
 * \param target the target frame
 * \param vs the velocity of the source frame
 * \param vt the velocity of the target frame
 * \return a moving transformation
 * \note If source/target are not defined with respect to the GLOBAL frame,
 *       the caller must take care that the intermediate frames are either
 *       stationary, or that their velocity has been accounted for in defining
 *       vs/vt.
 */
MOVINGTRANSFORM3 MOVINGTRANSFORM3::calc_transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const SVELOCITY& vs, const SVELOCITY& vt)
{
  const boost::shared_ptr<const POSE3> GLOBAL; 
  boost::shared_ptr<const POSE3> r, s; 
  MOVINGTRANSFORM3 result;

  // setup the source and targets
  result.source = source;
  result.target = target;

  // check for special case: no transformation 
  if (source == target)
  {
    result.r.set_zero();
    result.rdot.set_zero();
    result.E.set_identity();
    result.Edot.set_zero();
    result.v = POSE3::transform(source, vs);
    return result;
  }

  // check for special case: transformation to global frame
  if (!target)
  {
    // make sure that the target is zero
    #ifndef NEXCEPT
    if (vt.get_linear().norm() > EPS || vt.get_angular().norm() > EPS)
      throw std::runtime_error("Target frame is global frame but velocity is non-zero");
    #endif

    // combine transforms from this to i: this will give aTl
    QUAT q = source->q;
    ORIGIN3 x = source->x;
    s = source;
    while (s)
    {
      s = s->rpose;
      if (!s)
        break;
      x = s->x + s->q * x;
      q = s->q * q;
    }

    // setup result r and E
    result.E = q;
    result.r = result.E.transpose_mult(-x);

    // store the velocity
    result.v = POSE3::transform(source, vs);

    // get the velocity of s in the target frame
    SVELOCITY vs_target = POSE3::transform(GLOBAL, vs);

    // setup result Edot
    MATRIX3 vs_target_hat = MATRIX3::skew_symmetric(vs_target.get_angular());
    result.Edot = vs_target_hat * result.E;

    // setup result rdot
    MATRIX3 vs_hat = MATRIX3::skew_symmetric(vs.get_angular());
    ORIGIN3 vs_target_o(vs_target.get_linear());
    result.rdot = result.E.transpose_mult(vs_hat * x - vs_target_o);

    // setup the transform 
    return result; 
  }

  // check for special case: transformation from global frame
  if (!source)
  {
    // make sure that the source is zero
    #ifndef NEXCEPT
    if (vs.get_linear().norm() > EPS || vs.get_angular().norm() > EPS)
      throw std::runtime_error("Source frame is global frame but velocity is non-zero");
    #endif

    // combine transforms from target to q
    QUAT q = target->q;
    ORIGIN3 x = target->x;
    r = target;
    while (r)
    {
      r = r->rpose;
      if (!r)
        break;
      x = r->x + r->q * x; 
      q = r->q * q;
    }

    // compute the inverse pose of the right 
    q = QUAT::invert(q);
    x = q * (-x);

    // setup result r and E
    result.E = q;
    result.r = result.E.transpose_mult(-x);

    // set the velocity of the source
    result.v.set_zero(source);

    // get the velocity of t in the global frame
    SVELOCITY vt_global = POSE3::transform(GLOBAL, vt);

    // setup result Edot
    MATRIX3 vt_hat = MATRIX3::skew_symmetric(-vt_global.get_angular());
    result.Edot = result.E.mult(vt_hat);

    // setup result rdot
    result.rdot = ORIGIN3(vt_global.get_linear());

    return result;
  }

  // if both transforms are defined relative to the same frame, this is easy
  if (source->rpose == target->rpose)
  {
    // store the velocity of the source
    result.v = POSE3::transform(source, vs);

    if (!source->rpose)
    {
      // compute the inverse pose of p 
      QUAT inv_target_q = QUAT::invert(target->q);
      QUAT q = inv_target_q * source->q;
      ORIGIN3 x = inv_target_q * (source->x - target->x);

      // setup the result
      result.E = q;
      result.r = result.E.transpose_mult(-x);

      // setup rotation matrices for calculating Edot
      const QUAT& wQs = source->q;
      QUAT sQw = wQs.inverse();
      const QUAT& tQw = inv_target_q;
      QUAT wQt = tQw.inverse();
      MATRIX3 tRw = tQw;
      MATRIX3 wRs = wQs;

      // E = tRw * wRs
      MATRIX3 vs_hat = MATRIX3::skew_symmetric(vs.get_angular());
      MATRIX3 vt_hat = MATRIX3::skew_symmetric(vt.get_angular());
      
      // E = tRw * wRs
      // d/dt E = tRw * \dot{wRs} + \dot{tRw} * wRs
      result.Edot = tRw * vs_hat * wRs - tRw * vt_hat * wRs;   

      // r = sRw * (xw^target - xw^source)
      // d/dt r =  d/dt sRw * (xw^target - xw^source) +
      //           wRs' * d/dt (xw^target - dw^source)
      // d/dt sRw = (d/dt wRs)' = (skew(omegas) * wRs)' = wRs' * -skew(omegas)
      ORIGIN3 vt_o(vt.get_linear());
      ORIGIN3 vs_o(vs.get_linear());
      result.rdot = -(sQw * (vs_hat * (target->x - source->x)))
                    +sQw * (vt_o - vs_o);

      return result;
    }
    else
    {
      // compute the inverse pose of p 
      QUAT inv_target_q = QUAT::invert(target->q);
      QUAT q = inv_target_q * source->q;
      ORIGIN3 x = inv_target_q * (source->x - target->x);

      // setup the result
      result.E = q;
      result.r = result.E.transpose_mult(-x);

      // setup rotation matrices for calculating Edot
      const QUAT& rQs = source->q;
      const QUAT& tQr = inv_target_q;

      // get transform from r to global
      r = source->rpose;
      TRANSFORM3 wTr = POSE3::calc_relative_pose(r, GLOBAL);

      // setup tRw and wRs
      QUAT tQw = tQr * wTr.q.inverse();
      QUAT wQt = tQw.inverse();
      QUAT wQs = wTr.q * rQs; 
      QUAT sQw = wQs.inverse();

      // get the velocities in the global frame
      SVELOCITY vsw = wTr.transform(vs);
      SVELOCITY vtw = wTr.transform(vt);

      // E = tRw * wRs
      MATRIX3 vs_hat = MATRIX3::skew_symmetric(vsw.get_angular());
      MATRIX3 vt_hat = MATRIX3::skew_symmetric(vtw.get_angular());
      
      // E = tRw * wRs
      // d/dt E = tRw * \dot{wRs} + \dot{tRw} * wRs
      // \dot{tRw} = \dot{wRt}' = (vt_hat * wRt)' = -tRw * vt_hat
      MATRIX3 tRw = tQw;
      MATRIX3 wRs = wQs;
      result.Edot = tRw * vs_hat * wRs - tRw * vt_hat * wRs;   

      // get the difference between the two frames in the world frame
      VECTOR3 diff = wTr.transform_vector(VECTOR3(target->x - source->x, r));

      // r = wRs' * (xw^target - xw^source)
      // d/dt r =  d/dt wRs' * (xw^target - xw^source) +
      //           wRs' * d/dt (xw^target - dw^source)
      ORIGIN3 vt_o(vtw.get_linear());
      ORIGIN3 vs_o(vsw.get_linear());
      result.rdot = -(sQw * (vs_hat * ORIGIN3(diff))) + 
                    sQw * (vt_o - vs_o);

      return result;
    }
  }
  else
  {
    // store the velocity of the source
    result.v = POSE3::transform(source, vs);

    // search for the common link - we arbitrary move up the target while
    // one step at a time while searching through all levels of the source 
    unsigned i = std::numeric_limits<unsigned>::max();
    r = target;
    while (true)
    {
      if (POSE3::is_common(source, r, i))
        break;
      else
      {
        assert(r);
        r = r->rpose;
      } 
    } 
    
    // combine transforms from this to i: this will give rTs, where r is
    // the common frame and s is the source
    QUAT left_q = source->q;
    ORIGIN3 left_x = source->x;
    s = source;
    for (unsigned j=0; j < i; j++)
    {
      s = s->rpose;
      if (!s)
        break;
      left_x = s->x + s->q * left_x;
      left_q = s->q * left_q;
    }

    // combine transforms from target to r, this will give tTr
    QUAT right_q = target->q;
    ORIGIN3 right_x = target->x;
    while (target != r)
    {
      target = target->rpose;
      if (!target)
        break;
      right_x = target->x + target->q * right_x; 
      right_q = target->q * right_q;
    }

    // compute the inverse pose of the right 
    QUAT inv_right_q = QUAT::invert(right_q);

    // multiply the inverse pose of p by this 
    QUAT q = inv_right_q * left_q;      
    ORIGIN3 x = inv_right_q * (left_x - right_x);

    // get the transform from r to the global frame
    TRANSFORM3 wTr = POSE3::calc_relative_pose(r, GLOBAL);

    // get the velocities of s and t in the global frame
    SVELOCITY vsw = POSE3::transform(GLOBAL, vs);
    SVELOCITY vtw = POSE3::transform(GLOBAL, vt);

    // setup result r and E
    result.E = q;
    result.r = result.E.transpose_mult(-x);

    // setup rotation matrices for calculating Edot
    MATRIX3 wRs = wTr.q * left_q;
    MATRIX3 tRw = inv_right_q * wTr.q.inverse();

    // E = tRw * wRs
    // d/dt E = tRw * \dot{wRs} + \dot{tRw} * wRs
    MATRIX3 vs_hat = MATRIX3::skew_symmetric(vsw.get_angular());
    MATRIX3 vt_hat = MATRIX3::skew_symmetric(vtw.get_angular());
    result.Edot = tRw * vs_hat * wRs - tRw * vt_hat * wRs;   

    // get the positions of the origins in the global frame
    ORIGIN3 target_x = right_x + wTr.x;
    ORIGIN3 source_x = left_x + wTr.x;

    // r = wRs' * (x^target - x^source)
    // d/dt r =  d/dt wRs' * (x^target - x^source) +
    //           wRs' * d/dt (x^target - x^source)
    ORIGIN3 vt_o(vtw.get_linear());
    ORIGIN3 vs_o(vsw.get_linear());
    result.rdot = -wRs.transpose_mult(vs_hat) * (target_x - source_x) + 
                  wRs.transpose_mult(vt_o - vs_o);

    return result;
  }
}


