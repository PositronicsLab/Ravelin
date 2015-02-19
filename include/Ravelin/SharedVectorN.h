/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef SHAREDVECTORN
#error This class is not to be included by the user directly. Use SharedVectorNff.h or SharedVectorNd.h instead.
#endif

class CONST_SHAREDVECTORN;
class VECTORN;

/// A generic N-dimensional floating point vector
class SHAREDVECTORN 
{
  friend class VECTORN;
  friend class MATRIXN;
  friend class MATRIX2;
  friend class MATRIX3;
  friend class CONST_SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDMATRIXN;

  public:
    SHAREDVECTORN();
    SHAREDVECTORN(const SHAREDVECTORN& source) { reset_from(source); }
    SHAREDVECTORN(unsigned len, unsigned inc, unsigned start, SharedResizable<REAL> data);
    void reset_from(const SHAREDVECTORN& source);
    virtual ~SHAREDVECTORN() {}
    SHAREDVECTORN& normalize() { assert(norm() > EPS); operator*=((REAL) 1.0/norm()); return *this; }
    unsigned size() const { return _len; }
    static REAL norm(const SHAREDVECTORN& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const SHAREDVECTORN& v) { return SHAREDVECTORN::dot(v, v); }
    REAL norm_inf() const { return norm_inf(*this); }
    REAL norm1() const { return norm1(*this); }
    REAL norm() const { return norm(*this); }
    REAL norm_sq() const { return norm_sq(*this); }
    SHAREDVECTORN segment(unsigned start, unsigned end);
    CONST_SHAREDVECTORN segment(unsigned start, unsigned end) const;
    SHAREDVECTORN& resize(unsigned N, bool preserve = false);
    SHAREDVECTORN& operator=(const VECTOR3& source);
    SHAREDVECTORN& operator=(const SHAREDVECTORN& source);
    SHAREDVECTORN& operator=(const CONST_SHAREDVECTORN& source);
    SHAREDVECTORN& operator=(const VECTORN& source);
    SHAREDVECTORN& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    REAL* data() { return _data.get() + _start; }
    const REAL* data() const { return _data.get() + _start; }
    SHAREDVECTORN& resize(unsigned m, unsigned n, bool preserve = false) { if (n != 1) throw MissizeException(); return resize(m, preserve); }

    /// Resets the shared structure
    void reset() { _data.reset(); _start = _inc = _len = 0; }

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
    SharedResizable<REAL> _data;
    unsigned _start;
    unsigned _inc;
    unsigned _len;
}; // end class

/// A generic N-dimensional floating point vector
class CONST_SHAREDVECTORN 
{
  friend class VECTORN;
  friend class MATRIXN;
  friend class MATRIX2;
  friend class MATRIX3;
  friend class SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDMATRIXN;

  public:
    CONST_SHAREDVECTORN();
    CONST_SHAREDVECTORN(unsigned len, unsigned inc, unsigned start, SharedResizable<REAL> data);
    CONST_SHAREDVECTORN(const SHAREDVECTORN& source) { reset_from(source); }
    CONST_SHAREDVECTORN(const CONST_SHAREDVECTORN& source) { reset_from(source); }
    const SHAREDVECTORN get() const; 
    void reset_from(const SHAREDVECTORN& source);
    void reset_from(const CONST_SHAREDVECTORN& source);
    virtual ~CONST_SHAREDVECTORN() {}
    unsigned size() const { return _len; }
    static REAL norm(const CONST_SHAREDVECTORN& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const CONST_SHAREDVECTORN& v) { return CONST_SHAREDVECTORN::dot(v, v); }
    REAL norm_inf() const { return norm_inf(*this); }
    REAL norm1() const { return norm1(*this); }
    REAL norm() const { return norm(*this); }
    REAL norm_sq() const { return norm_sq(*this); }
    CONST_SHAREDVECTORN segment(unsigned start, unsigned end) const;
    CONST_SHAREDVECTORN& resize(unsigned N, bool preserve = false);
    const REAL* data() const { return _data.get() + _start; }
    CONST_SHAREDVECTORN& resize(unsigned m, unsigned n, bool preserve = false) 
    { 
      #ifndef NEXCEPT
      if (n != 1) 
        throw MissizeException(); 
      #endif
      return resize(m, preserve); 
    }

    /// Resets the shared structure
    void reset() { _data.reset(); _start = _inc = _len = 0; }
    unsigned rows() const { return _len; }
    unsigned columns() const { return 1; }
    unsigned leading_dim() const { return _len; }
    unsigned inc() const { return _inc; }

    // inline code specific to SharedVectorN
    #include "ConstSharedVectorN.inl"

  protected:
    SharedResizable<REAL> _data;
    unsigned _start;
    unsigned _inc;
    unsigned _len;
}; // end class

std::ostream& operator<<(std::ostream& out, const SHAREDVECTORN& v);
std::istream& operator>>(std::istream& in, SHAREDVECTORN& v);
std::ostream& operator<<(std::ostream& out, const CONST_SHAREDVECTORN& v);

