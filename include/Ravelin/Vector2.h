/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef VECTOR2
#error This class is not to be included by the user directly. Use VECTOR2.h or Vector2f.h instead.
#endif

class ORIGIN2;
class POSE2;

/// A two-dimensional floating point vector used for computational geometry calculations
class VECTOR2
{
  public:
    VECTOR2(boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>()) { this->pose = pose; }
    VECTOR2(boost::shared_ptr<POSE2> pose) { this->pose = boost::const_pointer_cast<const POSE2>(pose); }
    VECTOR2(REAL x, REAL y, boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>());
    VECTOR2(REAL x, REAL y, boost::shared_ptr<POSE2> pose);
    VECTOR2(const REAL* array, boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>());
    VECTOR2(const REAL* array, boost::shared_ptr<POSE2> pose);
    VECTOR2(const ORIGIN2& o, boost::shared_ptr<POSE2> pose) { this->pose = boost::const_pointer_cast<const POSE2>(pose); operator=(o); } 
    VECTOR2(const ORIGIN2& o, boost::shared_ptr<const POSE2> pose) { this->pose = pose; operator=(o); } 
    REAL dot(const VECTOR2& v) const { return v[0]*_data[0] + v[1]*_data[1]; }
    static REAL dot(const VECTOR2& v1, const VECTOR2& v2) { return v1[0]*v2[0] + v1[1]*v2[1]; }
    REAL norm_inf() const { return std::max(std::fabs(_data[0]), std::fabs(_data[1])); }
    REAL norm() const { return std::sqrt(norm_sq()); }
    REAL norm_sq() const { return dot(*this, *this); }
    void normalize() { assert(norm() > std::numeric_limits<REAL>::epsilon()); operator/=(norm()); }
    static VECTOR2 normalize(const VECTOR2& v) { VECTOR2 w = v; w.normalize(); return w; }
    unsigned size() const { return 2; }
    static REAL norm(const VECTOR2& v) { return std::sqrt(norm_sq(v)); }
    static REAL norm_sq(const VECTOR2& v) { return v.dot(v); }
    VECTOR2& set_zero() { _data[0] = _data[1] = 0.0; return *this; }
    VECTOR2& set_one() { _data[0] = _data[1] = 1.0; return *this; }
    VECTOR2& set_zero(boost::shared_ptr<const POSE2> pose) { _data[0] = _data[1] = 0.0; this->pose = pose; return *this; }
    VECTOR2& set_one(boost::shared_ptr<const POSE2> pose) { _data[0] = _data[1] = 1.0; this->pose = pose; return *this; }
    static VECTOR2 zero(boost::shared_ptr<const POSE2> pose = boost::shared_ptr<const POSE2>()) { return VECTOR2(0.0, 0.0, pose); }
    VECTOR2& operator=(const VECTOR2& v) { _data[0] = v[0]; _data[1] = v[1]; pose = v.pose; return *this; }
    VECTOR2& operator=(const ORIGIN2& o);
    VECTOR2 operator+(const VECTOR2& v) const { VECTOR2 w = *this; return w += v; }
    VECTOR2 operator-(const VECTOR2& v) const { VECTOR2 w = *this; return w -= v; }
    VECTOR2 operator+(const ORIGIN2& o) const { VECTOR2 v = *this; return v += o; }
    VECTOR2 operator-(const ORIGIN2& o) const { VECTOR2 v = *this; return v -= o; }
    VECTOR2& operator+=(const VECTOR2& v);
    VECTOR2& operator-=(const VECTOR2& v);
    VECTOR2& operator+=(const ORIGIN2& o);
    VECTOR2& operator-=(const ORIGIN2& o);
    VECTOR2& operator*=(REAL scalar) { _data[0] *= scalar; _data[1] *= scalar; return *this; }
    VECTOR2& operator/=(REAL scalar) { _data[0] /= scalar; _data[1] /= scalar; return *this; }
    VECTOR2 operator*(REAL scalar) const { VECTOR2 v = *this; v *= scalar; return v; }
    VECTOR2 operator/(REAL scalar) const { VECTOR2 v = *this; v /= scalar; return v; }
    VECTOR2 operator-() const { return VECTOR2(-_data[0], -_data[1]); }
    REAL* data() { return _data; }
    const REAL* data() const { return _data; }
    unsigned rows() const { return 2; }
    unsigned columns() const { return 1; }
    unsigned inc() const { return 1; }
    unsigned leading_dim() const { return 2; }
    REAL& operator[](const unsigned i);
    const REAL& operator[](const unsigned i) const;
    VECTOR2& resize(unsigned m, unsigned n, bool preserve = false);
    COLUMN_ITERATOR begin() { return column_iterator_begin(); }
    CONST_COLUMN_ITERATOR begin() const { return column_iterator_begin(); }
    COLUMN_ITERATOR end() { return column_iterator_end(); }
    CONST_COLUMN_ITERATOR end() const { return column_iterator_end(); }
    COLUMN_ITERATOR column_iterator_begin();
    CONST_COLUMN_ITERATOR column_iterator_begin() const;
    COLUMN_ITERATOR column_iterator_end();
    CONST_COLUMN_ITERATOR column_iterator_end() const;
    ROW_ITERATOR row_iterator_begin();
    CONST_ROW_ITERATOR row_iterator_begin() const;
    ROW_ITERATOR row_iterator_end();
    CONST_ROW_ITERATOR row_iterator_end() const;
    REAL* data(unsigned i);
    const REAL* data(unsigned i) const;
    const REAL& x() const { return _data[0]; }
    const REAL& y() const { return _data[1]; }
    REAL& x() { return _data[0]; }
    REAL& y() { return _data[1]; }

    /// Computes the "perp" operator
    VECTOR2 perp() const { return VECTOR2(_data[1], -_data[0]); }

    /// Computes the dot "perp" product
    REAL dot_perp(const VECTOR2& v) const { return dot(v.perp()); }

    /// The pose that this vector is defined in
    boost::shared_ptr<const POSE2> pose;

  private:
    REAL _data[2];
}; // end class

inline VECTOR2 operator*(REAL scalar, const VECTOR2& v) { return v * scalar; }

/// Writes a VECTOR2 to the specified stream
inline std::ostream& operator<<(std::ostream& out, const VECTOR2& v)
{
  out << '[' << v[0] << ',' << ' ' << v[1] << ']' << ' ';
  return out;
};

