/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef VECTORN
#error This class is not to be included by the user directly. Use VectorNf.h or VectorNd.h instead. 
#endif

class MATRIXN;
class SHAREDMATRIXN;
class CONST_SHAREDMATRIXN;

/// A generic N-dimensional floating point vector
class VECTORN 
{
  public:
    template <class ForwardIterator>
    VECTORN(ForwardIterator begin, ForwardIterator end);

    VECTORN();
    VECTORN(unsigned N);
    VECTORN(const VECTORN& source);
    VECTORN(const SHAREDVECTORN& source);
    VECTORN(const CONST_SHAREDVECTORN& source);
    VECTORN(const VECTOR2& v);
    VECTORN(const VECTOR3& v);
    VECTORN(const MATRIXN& v);
    VECTORN(const SHAREDMATRIXN& v);
    VECTORN(const CONST_SHAREDMATRIXN& v);
    VECTORN(unsigned N, const REAL* array);
    static VECTORN construct_variable(unsigned N, ...);
    virtual ~VECTORN() {}
    VECTORN& normalize() { assert(norm() > EPS); operator*=((REAL) 1.0/norm()); return *this; }
    unsigned size() const { return _data.size(); }
    static REAL norm_sq(const VECTORN& v) { return VECTORN::dot(v, v); }
    static REAL norm(const VECTORN& v) { return std::sqrt(norm_sq(v)); } 
    REAL norm_inf() const { return norm_inf(*this); }
    REAL norm1() const { return norm1(*this); }
    REAL norm() const { return norm(*this); }
    REAL norm_sq() const { return norm_sq(*this); }
    static VECTORN one(unsigned N);
    static VECTORN& concat(const VECTORN& v1, const VECTORN& v2, VECTORN& result);
    VECTORN& augment(const VECTORN& v);
    VECTORN& resize(unsigned N, bool preserve = false) { _data.resize(N, preserve); return *this; }
    SHAREDVECTORN segment(unsigned start_idx, unsigned end_idx);
    CONST_SHAREDVECTORN segment(unsigned start_idx, unsigned end_idx) const;
    static VECTORN zero(unsigned n);
    VECTORN& operator=(REAL r);
    VECTORN& operator=(const VECTOR2& source);
    VECTORN& operator=(const VECTOR3& source);
    VECTORN& operator=(const VECTORN& source);
    VECTORN& operator=(const SHAREDVECTORN& source);
    VECTORN& operator=(const CONST_SHAREDVECTORN& source);
    VECTORN& operator=(const MATRIXN& source);
    VECTORN& operator=(const SHAREDMATRIXN& source);
    VECTORN& operator=(const CONST_SHAREDMATRIXN& source);
    VECTORN& operator/=(REAL scalar) { return operator*=((REAL) 1.0/scalar); }
    REAL* data() { return _data.get(); }
    const REAL* data() const { return _data.get(); }
    static VECTORN& parse(const std::string& s, VECTORN& v);
    static VECTORN parse(const std::string& s) { VECTORN v; parse(s, v); return v; }
    VECTORN& resize(unsigned m, unsigned n, bool preserve = false);

    void free_memory() { resize(0); compress(); }
    void compress() { _data.compress(); }
    unsigned rows() const { return _data.size(); }
    unsigned columns() const { return 1; }
    unsigned leading_dim() const { return _data.size(); }
    unsigned inc() const { return 1; }

    // inline code specific to VectorN
    #include "VectorN.inl"

    // inline code for VectorN/SharedVectorN
    #define XVECTORN VECTORN
    #include "XVectorN.inl"
    #undef XVECTORN

  protected:
    SharedResizable<REAL> _data;
}; // end class

std::ostream& operator<<(std::ostream& out, const VECTORN& v);
std::istream& operator>>(std::istream& in, VECTORN& v);

