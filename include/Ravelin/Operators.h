#ifndef _RAVELIN_OPERATORS_H_
#define _RAVELIN_OPERATORS_H_

#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/MissizeException.h>


namespace Ravelin {

/// Does a transposition operation
template <class M1, class M2>
M2& transpose(const M1& m1, M2& m2)
{
  if (sizeof(*m1.data()) != sizeof(*m2.data()))
    throw DataMismatchException();
 
  // resize m2
  m2.resize(m1.columns(), m2.rows());

  // setup necessary constants
  const unsigned LD1 = m1.leading_dim();
  const unsigned LD2 = m2.leading_dim();

  // case 1: double type
  if (sizeof(*m1.data()) == sizeof(double))
  {
    const double* d1 = (const double*) m1.data(); 
    double* d2 = (double*) m2.data(); 

    // do the transposition
    for (unsigned i=0, i1=0; i< m1.rows(); i++)
    {
      for (unsigned j=0, i2=i; j< m1.columns(); j++, i2 += LD2)
        d2[i2] = d1[i1+j]; 
      i1 += LD1;
    }
  }
  else
  {
    assert(sizeof(*m1.data()) == sizeof(float));
    const float* d1 = (const float*) m1.data(); 
    float* d2 = (float*) m2.data(); 

    // do the transposition
    for (unsigned i=0, i1=0; i< m1.rows(); i++)
    {
      for (unsigned j=0, i2=i; j< m1.columns(); j++, i2 += LD2)
        d2[i2] = d1[i1+j];
      i1 += LD1;
    }
  }

  return m2;
}

/// Determines whether two floats are equal
inline bool rel_equal(float x, float y)
{
  return (std::fabs(x-y) <= EPS_FLOAT * std::max(std::fabs(x), std::max(std::fabs(y), 1.0f)));
}

/// Determines whether two floats are equal
inline bool rel_equal(float x, float y, float tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), 1.0f)));
}

/// Determines whether two doubles are equal
inline bool rel_equal(double x, double y)
{
  return (std::fabs(x-y) <= EPS_DOUBLE * std::max(std::fabs(x), std::max(std::fabs(y), 1.0)));
}

/// Determines whether two doubles are equal
inline bool rel_equal(double x, double y, double tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), 1.0)));
}

/// Multiples a matrix by the transpose of a matrix or vector
template <class T, class U, class V>
V& mult_transpose(const T& x, const U& y, V& z)
{
  // verify that we can multiply these
  if (x.columns() != y.columns())
    throw MissizeException();
  if (sizeof(x.data()) != sizeof(y.data()))
    throw DataMismatchException();
  if (sizeof(x.data()) != sizeof(z.data()))
    throw DataMismatchException();

  // resize the result
  z.resize(x.rows(), y.rows());

  // do the multiplication
  CBLAS::gemm(CblasNoTrans, CblasTrans, x.rows(), y.rows(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

  return z;
}

/// Multiples a matrix by a matrix or vector
template <class T, class U, class V>
V& mult(const T& x, const U& y, V& z)
{
  // verify that we can multiply these
  if (x.columns() != y.rows())
    throw MissizeException();
  if (sizeof(x.data()) != sizeof(y.data()))
    throw DataMismatchException();
  if (sizeof(x.data()) != sizeof(z.data()))
    throw DataMismatchException();

  // resize the result
  z.resize(x.rows(), y.columns());

  // do the multiplication
  if (y.columns() > 1)
    // matrix
    CBLAS::gemm(CblasNoTrans, CblasNoTrans, x.rows(), y.columns(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);
  else
    // vector
    CBLAS::gemv(CblasColMajor, CblasNoTrans, x.rows(), x.columns(), 1.0, x.data(), x.leading_dim(), y.data(), y.inc(), 0.0, z.data(), z.inc());

  return z;
}

/// Multiples the transpose of a matrix or vector by a matrix or vector
template <class T, class U, class V>
V& transpose_mult(const T& x, const U& y, V& z)
{
  // verify that we can multiply these
  if (x.rows() != y.rows())
    throw MissizeException();
  if (sizeof(x.data()) != sizeof(y.data()))
    throw DataMismatchException();
  if (sizeof(x.data()) != sizeof(z.data()))
    throw DataMismatchException();

  // resize the result
  z.resize(x.columns(), y.columns());

  // multiply
  CBLAS::gemm(CblasTrans, CblasNoTrans, x.columns(), y.columns(), x.rows(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

  return z;
}

/// Multiples a matrix by a matrix or vector
template <class T, class U, class V>
V& transpose_mult_transpose(const T& x, const U& y, V& z)
{
  // verify that we can multiply these
  if (x.rows() != y.columns())
    throw MissizeException();
  if (sizeof(x.data()) != sizeof(y.data()))
    throw DataMismatchException();
  if (sizeof(x.data()) != sizeof(z.data()))
    throw DataMismatchException();

  // resize the result
  z.resize(x.columns(), y.rows());

  // do the multiplication
  CBLAS::gemm(CblasTrans, CblasTrans, x.columns(), y.rows(), x.rows(), 1.0, x.data(), x.leading_dim(), y.data(), y.leading_dim(), 0.0, z.data(), z.leading_dim);

  return z;
}

template <class U, class V, class M>
M& outer_prod(const U& x, const V& y, M& z)
{
  // make sure both are vectors
  if (x.columns() != 1 || y.columns() != 1)
    throw MissizeException();
  if (sizeof(x.data()) != sizeof(y.data()))
    throw DataMismatchException();
  if (sizeof(x.data()) != sizeof(z.data()))
    throw DataMismatchException();

  // get n and m
  unsigned m = x.rows();
  unsigned n = y.rows();

  // resize the matrix and set to zero
  z.set_zero(m,n);

  // call BLAS outer product routine
  if (m > 0 && n > 0)
    CBLAS::ger(CblasColMajor, m, n, x.data(), x.inc(), y.data(), y.inc(), z.data(), z.leading_dim());

  return z;
}

} // end namespace

#endif

