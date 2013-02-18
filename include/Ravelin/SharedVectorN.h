/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef SHAREDVECTORN
#error This class is not to be included by the user directly. Use SharedVectorNff.h or SharedVectorNd.h instead.
#endif

/// A generic N-dimensional floating point vector
class SHAREDVECTORN 
{
  friend class VECTORN;
  friend class MATRIXN;
  friend class SHAREDMATRIXN;

  public:
    SHAREDVECTORN();
    SHAREDVECTORN(const SHAREDVECTORN& source);
    virtual ~SHAREDVECTORN() {}
    SHAREDVECTORN& normalize() { assert(norm() > EPS); operator*=((REAL) 1.0/norm()); return *this; }
    unsigned size() const { return _len; }
    static REAL norm(const SHAREDVECTORN& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const SHAREDVECTORN& v) { return SHAREDVECTORN::dot(v, v); }
    REAL norm_inf() const { return norm_inf(*this); }
    REAL norm1() const { return norm1(*this); }
    REAL norm() const { return norm(*this); }
    REAL norm_sq() const { return norm_sq(*this); }
    SHAREDVECTORN get_sub_vec(unsigned start, unsigned end);
    SHAREDVECTORN& resize(unsigned N, bool preserve = false);
    SHAREDVECTORN& operator=(const VECTOR3& source);
    SHAREDVECTORN& operator=(const SHAREDVECTORN& source);
    SHAREDVECTORN& operator=(const VECTORN& source);
    SHAREDVECTORN& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    REAL* data() { assert(_data); return _data.get() + _start; }
    const REAL* data() const { assert(_data); return _data.get() + _start; }
    SHAREDVECTORN& resize(unsigned m, unsigned n, bool preserve = false) { if (n != 1) throw MissizeException(); return resize(m, preserve); }

    unsigned rows() const { return _len; }
    unsigned columns() const { return 1; }
    unsigned leading_dim() const { return _len; }
    unsigned inc() const { return _inc; }

    // inline code specific to SharedVectorN
    #include "SharedVectorN.inl"

    // inline code for VectorN/SharedVectorN
    #define XVECTORN SHAREDVECTORN
    #include "XVectorN.inl"
    #undef XVECTORN

  protected:
    boost::shared_array<REAL> _data;
    unsigned _start;
    unsigned _inc;
    unsigned _len;
}; // end class

std::ostream& operator<<(std::ostream& out, const SHAREDVECTORN& v);
std::istream& operator>>(std::istream& in, SHAREDVECTORN& v);

