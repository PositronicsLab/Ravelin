/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef MOVINGTRANSFORM3
#error This class is not to be included by the user directly. Use MovingTransform3d.h or MovingTransform3f.h instead. 
#endif

class POSE3;

/// A transformation between two rigid body poses 
class MOVINGTRANSFORM3
{
  public:
    MOVINGTRANSFORM3();
    SACCEL transform(const SACCEL& a) const;
/*
    POSE3 transform(const POSE3& p) const;
    POSE3 inverse_transform(const POSE3& p) const;
    SACCEL inverse_transform(const SACCEL& t) const;
    TRANSFORM3& invert();
    TRANSFORM3 inverse() const { return invert(*this); }
    static TRANSFORM3 invert(const TRANSFORM3& m);
*/
    MOVINGTRANSFORM3& operator=(const MOVINGTRANSFORM3& source);
    static MOVINGTRANSFORM3 calc_transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target, const SVELOCITY& vs, const SVELOCITY& vt);

    /// the velocity of the original frame
    SVELOCITY v;

    /// the r vector
    ORIGIN3 r;

    /// the time derivative of the r vector
    ORIGIN3 rdot;

    /// the E matrix 
    MATRIX3 E;

    /// the time derivative of the E matrix 
    MATRIX3 Edot;

    /// the "source" pose
    boost::shared_ptr<const POSE3> source; 

    /// the "target" pose
    boost::shared_ptr<const POSE3> target; 

  private:
    void transform_spatial(const SVECTOR6& w, SVECTOR6& result) const;
    void inverse_transform_spatial(const SVECTOR6& w, SVECTOR6& result) const;
}; // end class

std::ostream& operator<<(std::ostream& out, const MOVINGTRANSFORM3& m);

