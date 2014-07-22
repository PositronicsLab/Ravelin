/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// ***********************************************************************
// template specializations BEGIN 
// ***********************************************************************
namespace Ravelin {

/// Gets the specified column of the matrix
/*
 * \param i the 0-index of the column
 */
template <>
ORIGIN2& MATRIX2::get_column(unsigned i, ORIGIN2& result) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i*SZ;
  result.x() = _data[st++];
  result.y() = _data[st++];
  return result;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
template <>
ORIGIN2& MATRIX2::get_row(unsigned i, ORIGIN2& result) const
{
  const unsigned SZ = 2;
  assert(i < SZ);

  unsigned st = i;
  result.x() = _data[st];
  st += SZ;
  result.y() = _data[st];
  return result;
}

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param o a 2D vector 
 */
template <>
MATRIX2& MATRIX2::set_column(unsigned i, const ORIGIN2& o)
{
  const unsigned SZ = 2;
  assert(i < SZ);

  const REAL* odata = o.data();
  unsigned st = i*SZ;
  _data[st++] = odata[0];
  _data[st] = odata[1];
  return *this;
}

/// Sets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param o a 2D vector
 * \note the number of columns of this must be three!
 */
template <>
MATRIX2& MATRIX2::set_row(unsigned i, const ORIGIN2& o)
{
  const unsigned SZ = 2;
  assert(i < SZ);

  const REAL* odata = o.data();
  unsigned st = i;
  _data[st] = odata[0];
  st += SZ;
  _data[st] = odata[1]; 
  return *this;
}


} // end namespace

// ***********************************************************************
// template specializations END 
// ***********************************************************************

/// Constructs a matrix from an array
/**
 * \param array an array of 4 REAL values in column-major format
 */
MATRIX2::MATRIX2(const REAL* array)
{
  for (unsigned i=0; i< 4; i++)
    _data[i] = array[i];
}

/// Constructs a matrix from 4 values
/**
 * The resulting matrix will appear as follows: <br />
 * m00 m01 <br />
 * m10 m11 
 */
MATRIX2::MATRIX2(REAL m00, REAL m01, REAL m10, REAL m11)
{
  const unsigned X = 0, Y = 1;
  _data[0] = m00;
  _data[1] = m10;
  _data[2] = m01;
  _data[3] = m11;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
ORIGIN2 MATRIX2::mult(const ORIGIN2& o) const
{
  ORIGIN2 result;
  REAL* rdata = result.data();
  const REAL* odata = o.data();
  rdata[0] = _data[0]*odata[0] + _data[2]*odata[1];
  rdata[1] = _data[1]*odata[0] + _data[3]*odata[1];
  return result;
}

/// Multiplies the transpose of this matrix by a vector and returns the result in a new vector
ORIGIN2 MATRIX2::transpose_mult(const ORIGIN2& o) const
{
  ORIGIN2 result;
  REAL* rdata = result.data();
  const REAL* odata = o.data();
  rdata[0] = _data[0]*odata[0] + _data[1]*odata[1];
  rdata[1] = _data[2]*odata[0] + _data[3]*odata[1];
  return result;
}

/// Multiplies the transpose of this matrix by a matrix and returns the result in a new matrix 
MATRIX2 MATRIX2::transpose_mult(const MATRIX2& m) const
{
  MATRIX2 result;
  result._data[0] = _data[0]*m._data[0] + _data[1]*m._data[1];
  result._data[1] = _data[2]*m._data[0] + _data[3]*m._data[1];
  result._data[2] = _data[0]*m._data[2] + _data[1]*m._data[3];
  result._data[3] = _data[2]*m._data[2] + _data[3]*m._data[3];
  return result;
}

/// Multiplies the transpose of this matrix by the transpose of a matrix and returns the result in a new matrix 
MATRIX2 MATRIX2::transpose_mult_transpose(const MATRIX2& m) const
{
  MATRIX2 result;
  result._data[0] = _data[0]*m._data[0] + _data[1]*m._data[2];
  result._data[1] = _data[2]*m._data[0] + _data[3]*m._data[2];
  result._data[2] = _data[0]*m._data[1] + _data[1]*m._data[3];
  result._data[3] = _data[2]*m._data[1] + _data[3]*m._data[3];
  return result;
}

/// Multiplies this matrix by the transpose of a matrix and returns the result in a new matrix 
MATRIX2 MATRIX2::mult_transpose(const MATRIX2& m) const
{
  MATRIX2 result;
  result._data[0] = _data[0]*m._data[0] + _data[2]*m._data[2];
  result._data[1] = _data[1]*m._data[0] + _data[3]*m._data[2];
  result._data[2] = _data[0]*m._data[1] + _data[2]*m._data[3];
  result._data[3] = _data[1]*m._data[1] + _data[3]*m._data[3];
  return result;
}

/// Multiplies this matrix by another 2x2 matrix
MATRIX2 MATRIX2::mult(const MATRIX2& m) const
{
  MATRIX2 result;
  result._data[0] = _data[0]*m._data[0] + _data[2]*m._data[1];
  result._data[1] = _data[1]*m._data[0] + _data[3]*m._data[1];
  result._data[2] = _data[0]*m._data[2] + _data[2]*m._data[3];
  result._data[3] = _data[1]*m._data[2] + _data[3]*m._data[3];
  return result;
}

/// Inverts this matrix 
MATRIX2& MATRIX2::invert()
{
  *this = invert(*this);
  return *this;
}

/// Calculates the determinant for a 2x2 matrix
REAL MATRIX2::det() const
{
  return _data[0]*_data[3] - _data[1]*_data[2];
}

/// Determines whether this is an orthonormal matrix
bool MATRIX2::is_orthonormal() const
{
  REAL determinant = det();
  return (std::fabs(determinant - 1.0) < EPS || std::fabs(determinant + (REAL) 1.0) < EPS);
}

/// Determines the inverse of the given matrix
MATRIX2 MATRIX2::invert(const MATRIX2& m)
{
  // compute the determinant
  REAL dt = (REAL) 1.0/m.det();
  return MATRIX2(m._data[3]*dt, -m._data[2]*dt, -m._data[1]*dt, m._data[0]*dt);
}

/// Sets this matrix to the rotation matrix of the specified angle around the Z axis
void MATRIX2::set_rot_Z(REAL angle)
{
  *this = rot_Z(angle);
}

/// Returns the rotation matrix of the specified angle around the Z axis
MATRIX2 MATRIX2::rot_Z(REAL angle)
{
  REAL sina = std::sin(angle);
  REAL cosa = std::cos(angle);
  MATRIX2 r;
  r._data[0] = cosa;
  r._data[1] = sina;
  r._data[2] = -sina;
  r._data[3] = cosa;
  return r;
}

/// Graham-Schmidt orthonormalization
bool MATRIX2::orthonormalize(VECTOR2& a, VECTOR2& b)
{
  // Gram-Schmidt orthonormalization produces vectors u0, u1, and u2 as
  // follows:
  //
  //   u0 = v0/|v0|
  //   u1 = (v1 - Dot(u0,v1)*u0)/|v1 - Dot(u0,v1)*u0|

  // Compute u0.
  REAL flen = a.norm();
  if (flen == (REAL) 0.0)
    return false;

  // Compute u1.
  REAL fd0 = VECTOR2::dot(a, b);
  b -= fd0*a;
  flen = b.norm();
  return (flen != (REAL) 0.0);
}

/// Makes the matrix orthonormal using Gram-Schmidt orthogonalization
bool MATRIX2::orthonormalize()
{
  // The algorithm uses Gram-Schmidt orthogonalization applied to the
  // columns of M.
  VECTOR2 a(xx(),xy());
  VECTOR2 b(yx(),yy());
  if (orthonormalize(a,b))
  {
    xx() = a.x();
    xy() = a.y();
    yx() = b.x();
    yy() = b.y();
    return true;
  }
  
  return false;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
ORIGIN2 MATRIX2::get_row(unsigned i) const
{
  ORIGIN2 o;
  get_row(i, o);
  return o;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 */
ORIGIN2 MATRIX2::get_column(unsigned i) const
{
  ORIGIN2 o;
  get_column(i, o);
  return o;
}

/// Sets this matrix to its transpose
void MATRIX2::transpose()
{
  const unsigned YX = 1, XY = 3;
  std::swap(_data[XY], _data[YX]);
}

/// Determines the transpose of this matrix
MATRIX2 MATRIX2::transpose(const MATRIX2& m)
{
  MATRIX2 n = m;
  n.transpose();
  return n;
}

/// Copies a matrix to this one
MATRIX2& MATRIX2::operator=(const MATRIX2& m)
{
  const unsigned SZ = 4;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m._data[i];
  return *this;
}

/// Copies a matrix to this one
MATRIX2& MATRIX2::operator=(const MATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 2 || m.columns() != 2)
    throw MissizeException();
  #endif
  const REAL* data = m.data();
  const unsigned SZ = 4;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = data[i];
  return *this;
}

/// Copies a matrix to this one
MATRIX2& MATRIX2::operator=(const SHAREDMATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 2 || m.columns() != 2)
    throw MissizeException();
  #endif
  const unsigned SZ = 4;
  CONST_COLUMN_ITERATOR miter = m.column_iterator_begin();
  for (unsigned i=0; i< SZ; i++)
    _data[i] = *miter++;
  return *this;
}

/// Copies a matrix to this one
MATRIX2& MATRIX2::operator=(const CONST_SHAREDMATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 2 || m.columns() != 2)
    throw MissizeException();
  #endif
  const unsigned SZ = 4;
  CONST_COLUMN_ITERATOR miter = m.column_iterator_begin();
  for (unsigned i=0; i< SZ; i++)
    _data[i] = *miter++;
  return *this;
}

/// Multiplies this matrix by a scalar in place
MATRIX2& MATRIX2::operator*=(REAL scalar)
{
  _data[0] *= scalar;
  _data[1] *= scalar;
  _data[2] *= scalar;
  _data[3] *= scalar;
  return *this;
}

/// Returns the negation of this matrix
MATRIX2 MATRIX2::operator-() const
{
  MATRIX2 m;
  std::transform(_data, _data+4, m.data(), std::negate<REAL>());
  return m;
}

/// Adds m to this in place
MATRIX2& MATRIX2::operator+=(const MATRIX2& m)
{
  std::transform(_data, _data+4, m.data(), _data, std::plus<REAL>());
  return *this;
}

/// Subtracts m from this in place
MATRIX2& MATRIX2::operator-=(const MATRIX2& m)
{
  std::transform(_data, _data+4, m.data(), _data, std::minus<REAL>());
  return *this;
}

/// Sets this matrix to identity
void MATRIX2::set_identity()
{
  xx() = yy() = (REAL) 1.0;
  xy() = yx() = (REAL) 0.0;
}

/// Sets this matrix to zero
void MATRIX2::set_zero()
{
  xx() = xy() = yx() = yy() = (REAL) 0.0;
}

/// Checks whether this matrix is symmetric
bool MATRIX2::is_symmetric(REAL tolerance) const
{
  // make sure tolerance is non-negative
  tolerance = std::fabs(tolerance);

  // check symmetric
  if (!OPS::rel_equal(xy(), yx(), tolerance))
    return false;

  return true;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const MATRIX2& m)
{
  const unsigned SZ = 2; 
  for (unsigned i=0; i< SZ; i++)
  {
    for (unsigned j=0; j< SZ-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,SZ-1) << std::endl;
  }
   
  return out;
}

REAL& MATRIX2::operator()(const unsigned i, unsigned j)
{
  #ifndef NEXCEPT
  if (i >= 2 || j >= 2)
    throw InvalidIndexException();
  #endif
  return _data[j*2+i];
}

const REAL& MATRIX2::operator()(const unsigned i, unsigned j) const
{
  #ifndef NEXCEPT
  if (i >= 2 || j >= 2)
    throw InvalidIndexException();
  #endif
  return _data[j*2+i];
}

REAL* MATRIX2::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 4)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

const REAL* MATRIX2::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 4)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

MATRIX2& MATRIX2::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 2 || n != 2)
    throw std::runtime_error("Attempt to resize fixed-length matrix!");
  #endif

  return *this;
}
 

