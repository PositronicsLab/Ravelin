/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <boost/algorithm/minmax.hpp>
#include <cstring>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <boost/algorithm/minmax.hpp>
#include <Ravelin/cblas.h>
#include "clapack.h"
#include <Ravelin/MissizeException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/NumericalException.h>
#include <Ravelin/SingularException.h>
#include <Ravelin/Opsf.h>
#include <Ravelin/LinAlgf.h>

using namespace Ravelin;
using std::vector;
using boost::shared_array;
using std::endl;

// ***********************************************************************
// single precision routines start here
// ***********************************************************************

float LinAlgf::log2(float x)
{
  return std::log(x)/std::log(2.0f);
}

/// Calls LAPACK function for Givens rotation
void LinAlgf::lartg_(float* F, float* G, float* CS, float* SN, float* R)
{
  slartg_(F, G, CS, SN, R);
}

/// Calls LAPACK function for solving tridiagonal systems of linear equations
void LinAlgf::gtsv_(INTEGER* N, INTEGER* NRHS, SINGLE* DL, SINGLE* D, SINGLE* DU, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sgtsv_(N, NRHS, DL, D, DU, B, LDB, INFO);
}

/// Calls LAPACK function for solving triangular systems of linear equations
void LinAlgf::trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, SINGLE* AP, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  char DIAG = 'N';
  strtrs_(UPLO, TRANS, &DIAG, N, NRHS, AP, LDA, B, LDB, INFO);
}

/// Calls LAPACK function *ormqr_
void LinAlgf::ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, SINGLE* A, INTEGER* LDA, SINGLE* TAU, SINGLE* C, INTEGER* LDC, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  SINGLE WORK_QUERY;
  sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK_QUERY, &LWORK, INFO);
  assert(INFO == 0);

  // declare memory
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // do the real call now
  sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, workv().data(), &LWORK, INFO);
  assert(INFO == 0);
}

/// Calls LAPACK function for doing LDL' factorization of a matrix
void LinAlgf::sptrf_(char* UPLO, INTEGER* N, SINGLE* AP, INTEGER* IPIV, INTEGER* INFO)
{
  ssptrf_(UPLO, N, AP, IPIV, INFO);
}

/// Calls LAPACK function for solution to Ax=b using LDL' factorization
void LinAlgf::sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* AP, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  ssptrs_(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for least-squares solution to Ax=b using SVD
void LinAlgf::gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, SINGLE* RCOND, INTEGER* INFO)
{
  // create array to hold singular values
  INTEGER min_mn = std::min(*M, *N);
  workv().resize(min_mn);

   // necessary for computing WORK sizes
  INTEGER ISPEC = (INTEGER) 9;
  INTEGER TMP = (INTEGER) 0;
  const char* NAME = "SGELSD";
  const char* OPTS = " ";
  INTEGER smlsiz = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, NRHS, &TMP, strlen(NAME), strlen(OPTS));
  assert(smlsiz > (INTEGER) 0);
  INTEGER NLVL = std::max((INTEGER) 0, (INTEGER) (log2((float) min_mn/(float) (smlsiz+1))+1));
  INTEGER LIWORK = std::max((INTEGER) 1, 3*min_mn*NLVL + 11*min_mn);
  iworkv().resize(LIWORK);
  INTEGER RANK;

  // first do WORK query
  float WORK_QUERY;
  INTEGER LWORK = -1;
  sgelsd_(M, N, NRHS, A, M, B, LDB, workv().data(), RCOND, &RANK, &WORK_QUERY, &LWORK, &iworkv().front(), INFO);
  assert(*INFO == 0);

  // setup LWORK 
  LWORK = (INTEGER) WORK_QUERY;

  // setup second WORK array
  workv2().resize(std::max((INTEGER) 1, LWORK));

  // compute
  sgelsd_(M, N, NRHS, A, M, B, LDB, workv().data(), RCOND, &RANK, workv2().data(), &LWORK, &iworkv().front(), INFO);
}

/// Calls appropriate LAPACK function for solving systems of linear equations Ax=b, where A is symmetric indefinite
void LinAlgf::sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{ 
  INTEGER LWORK = -1;
  float WORK_QUERY;

  // first, determine workspace size
  ssysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return;

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // call LAPACK
  ssysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for svd (divide and conquer)
void LinAlgf::gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* S, SINGLE* U, INTEGER* LDU, SINGLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  float WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;
  iworkv().resize(8*minmn);

  // call LAPACK to determine the optimal workspace size
  sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, &iworkv().front(), INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK); 

  // call LAPACK once again
  sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, workv().data(), &LWORK, &iworkv().front(), INFO);
}

/// Calls LAPACK function for svd 
void LinAlgf::gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* S, SINGLE* U, INTEGER* LDU, SINGLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  float WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;

  // call LAPACK to determine the optimal workspace size
  sgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK); 

  // call LAPACK once again
  sgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for computing eigenvalues and eigenvectors
void LinAlgf::syevd_(char* JOBZ, char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* EVALS, INTEGER* INFO)
{
  // do work query first
  float WORK_QUERY;
  INTEGER LWORK = -1;
  INTEGER IWORK_QUERY;
  INTEGER LIWORK = -1;
  ssyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK_QUERY, &LWORK, &IWORK_QUERY, &LIWORK, INFO);

  // set array sizes
  LWORK = (INTEGER) WORK_QUERY;
  LIWORK = IWORK_QUERY;
  workv().resize(LWORK);
  iworkv().resize(LIWORK);
  ssyevd_(JOBZ, UPLO, N, A, LDA, EVALS, workv().data(), &LWORK, &iworkv().front(), &LIWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations using LU factorization
void LinAlgf::gesv_(INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* X, INTEGER* LDX, INTEGER* INFO)
{
  sgesv_(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

/// Calls LAPACK function for solving a system of linear equations from a Cholesky factorization
void LinAlgf::potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  spotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
}

/// Calls LAPACK function for computing matrix inverse using Cholesky factorization
void LinAlgf::potri_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* INFO)
{
  spotri_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for Cholesky factorization
void LinAlgf::potrf_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* INFO)
{
  spotrf_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for solving system of equations Ax=b, where A is PSD
void LinAlgf::posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sposv_(UPLO, N, NRHS, A, N, B, LDB, INFO);
}

/// Calls LAPACK function for inverting symmetric indefinite matrix using factorization
void LinAlgf::sytri_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  workv().resize(*N);
  ssytri_(UPLO, N, A, LDA, IPIV, workv().data(), INFO);
}

/// Calls LAPACK function for factorizing symmetric, indefinite matrix 
void LinAlgf::sytrf_(char* UPLO, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO) 
{
  // perform workspace query for factorization
  float WORK_QUERY;
  INTEGER LWORK = -1;
  ssytrf_(UPLO, N, A, LDA, IPIV, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup WORK array
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // perform the necessary factorization
  ssytrf_(UPLO, N, A, LDA, IPIV, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
void LinAlgf::geqp3_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* JPVT, SINGLE* TAU, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  SINGLE work_query;
  sgeqp3_(M, N, A, LDA, JPVT, TAU, &work_query, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) work_query;
  workv().resize(LWORK);

  // do QR factorization
  sgeqp3_(M, N, A, LDA, JPVT, TAU, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
void LinAlgf::geqrf_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, SINGLE* TAU, INTEGER* INFO)
{
  // determine LWORK
  const char* NAME = "SGEQRF";
  const char* OPTS = " ";
  INTEGER ISPEC = 1L;
  INTEGER TMP = -1;
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, &TMP, &TMP, strlen(NAME), strlen(OPTS));
  assert(NB >= 0);

  // setup WORK vectors
  INTEGER LWORK = NB*(*N);
  workv().resize(LWORK);
  sgeqrf_(M, N, A, LDA, TAU, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for LU factorization 
void LinAlgf::getrf_(INTEGER* M, INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  sgetrf_(M, N, A, LDA, IPIV, INFO);
}

/// Calls LAPACK function for matrix inversion using LU factorization
void LinAlgf::getri_(INTEGER* N, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{ 
  // compute block size
  INTEGER TMP = -1;
  INTEGER ISPEC = 1;
  const char* NAME = "SGETRI";
  const char* OPTS = " ";
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, N, &TMP, &TMP, &TMP, strlen(NAME), strlen(OPTS));

  // setup LAPACK parameters
  INTEGER LWORK = (*N)*NB;

  // setup work vector
  workv().resize(LWORK);
  
  sgetri_(N, A, LDA, IPIV, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations with a given LU factorization
void LinAlgf::getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, SINGLE* A, INTEGER* LDA, INTEGER* IPIV, SINGLE* B, INTEGER* LDB, INTEGER* INFO)
{
  sgetrs_(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for forming Q from a QR factorization
void LinAlgf::orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, SINGLE* A, INTEGER* LDA, SINGLE* TAU, INTEGER* INFO)
{
  // do a workspace query
  INTEGER LWORK = -1;
  SINGLE WORK_QUERY;
  sorgqr_(M, N, K, A, LDA, TAU, &WORK_QUERY, &LWORK, INFO);

  // initialize the work array
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // call the function for real
  sorgqr_(M, N, K, A, LDA, TAU, workv().data(), &LWORK, INFO);
}

#include <Ravelin/fdefs.h>
#include "LinAlg.cpp"
#include <Ravelin/undefs.h>

// external methods
extern int solve_superlu(bool notrans, bool CSR, int m, int n, int nrhs, int nnz, int* col_indices, int* row_ptr, float* A_nz, float* x, float* b);

/// Does a LU factorization of a sparse matrix
VectorNf& LinAlgf::solve_sparse_direct(const SparseMatrixNf& A, const VectorNf& b, Transposition trans, VectorNf& x)
{
  #ifdef USE_SUPERLU
  // verify sizes
  if (A.rows() != A.columns() || A.rows() != b.rows())
    throw MissizeException();  

  // resize x
  x.resize(b.size());

  // setup constants
  const int nrhs = b.columns();
  int m = (int) A.rows();
  int n = (int) A.columns();

  // setup A 
  int nnz = (int) A.get_nnz();
  float* nz_val = (float*) A.get_data();
  int* col_indices = (int*) A.get_indices();
  int* row_ptr = (int*) A.get_ptr();

  // check info
  int info = solve_superlu(trans == eNoTranspose, A.get_storage_type() == SparseMatrixNf::eCSR, m, n, nrhs, nnz, col_indices, row_ptr, nz_val, x.data(), (float*) b.data());
  if (info > 0)
    throw SingularException();
  #else
  throw std::runtime_error("Ravelin not built with SuperLU support!");
  #endif

  return x;
}

/// Does a LU factorization of a sparse matrix
MatrixNf& LinAlgf::solve_sparse_direct(const SparseMatrixNf& A, const MatrixNf& B, Transposition trans, MatrixNf& X)
{
  #ifdef USE_SUPERLU
  // verify sizes
  if (A.rows() != A.columns() || A.rows() != B.rows())
    throw MissizeException();  

  // resize x
  X.resize(A.columns(), B.columns());

  // setup constants
  const int nrhs = B.columns();
  int m = (int) A.rows();
  int n = (int) A.columns();

  // setup A 
  int nnz = (int) A.get_nnz();
  float* nz_val = (float*) A.get_data();
  int* col_indices = (int*) A.get_indices();
  int* row_ptr = (int*) A.get_ptr();

  // check info
  int info = solve_superlu(trans == eNoTranspose, A.get_storage_type() == SparseMatrixNf::eCSR, m, n, nrhs, nnz, col_indices, row_ptr, nz_val, X.data(), (float*) X.data());
  if (info > 0)
    throw SingularException();
  #else
  throw std::runtime_error("Ravelin not built with SuperLU support!");
  #endif

  return X;
}

