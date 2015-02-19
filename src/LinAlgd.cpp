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
#include <stdexcept>
#include <boost/algorithm/minmax.hpp>
#include <Ravelin/cblas.h>
#include "clapack.h"
#include <Ravelin/MissizeException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/NumericalException.h>
#include <Ravelin/SingularException.h>
#include <Ravelin/SparseMatrixNd.h>
#include <Ravelin/LinAlgd.h>

using namespace Ravelin;
using std::vector;
using boost::shared_array;
using std::endl;

double LinAlgd::log2(double x)
{
  return std::log(x)/std::log(2.0);
}

/// Calls LAPACK function for Givens rotation
void LinAlgd::lartg_(double* F, double* G, double* CS, double* SN, double* R)
{
  dlartg_(F, G, CS, SN, R);
}

/// Calls LAPACK function for solving tridiagonal systems of linear equations
void LinAlgd::gtsv_(INTEGER* N, INTEGER* NRHS, DOUBLE* DL, DOUBLE* D, DOUBLE* DU, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dgtsv_(N, NRHS, DL, D, DU, B, LDB, INFO);
}

/// Calls LAPACK function for solving triangular systems of linear equations
void LinAlgd::trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, DOUBLE* AP, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  char DIAG = 'N';
  dtrtrs_(UPLO, TRANS, &DIAG, N, NRHS, AP, LDA, B, LDB, INFO);
}

/// Calls LAPACK function *ormqr_
void LinAlgd::ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, DOUBLE* C, INTEGER* LDC, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  DOUBLE WORK_QUERY;
  dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO == 0);

  // declare memory
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // do the real call now
  dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, workv().data(), &LWORK, INFO);
  assert(*INFO == 0);
}

/// Calls LAPACK function for doing LDL' factorization of a matrix
void LinAlgd::sptrf_(char* UPLO, INTEGER* N, DOUBLE* AP, INTEGER* IPIV, INTEGER* INFO)
{
  dsptrf_(UPLO, N, AP, IPIV, INFO);
}

/// Calls LAPACK function for solution to Ax=b using LDL' factorization
void LinAlgd::sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* AP, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dsptrs_(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for least-squares solution to Ax=b using SVD
void LinAlgd::gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, DOUBLE* RCOND, INTEGER* INFO)
{
  // create array to hold singular values
  INTEGER min_mn = std::min(*M, *N);
  workv2().resize(min_mn);

   // necessary for computing WORK sizes
  INTEGER ISPEC = (INTEGER) 9;
  INTEGER TMP = (INTEGER) 0;
  const char* NAME = "DGELSD";
  const char* OPTS = " ";
  INTEGER smlsiz = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, NRHS, &TMP, strlen(NAME), strlen(OPTS));
  assert(smlsiz > (INTEGER) 0);
  INTEGER NLVL = std::max((INTEGER) 0, (INTEGER) (log2((double) min_mn/(double) (smlsiz+1))+1));
  INTEGER LIWORK = std::max((INTEGER) 1, 3*min_mn*NLVL + 11*min_mn);
  iworkv().resize(LIWORK);
  INTEGER RANK;

  // first do WORK query
  double WORK_QUERY;
  INTEGER LWORK = -1;
  dgelsd_(M, N, NRHS, A, M, B, LDB, workv2().data(), RCOND, &RANK, &WORK_QUERY, &LWORK, &LIWORK, INFO);
  assert(*INFO == 0);

  // setup LWORK 
  LWORK = (INTEGER) WORK_QUERY;

  // setup WORK array
  workv().resize(std::max((INTEGER) 1, LWORK));

  // compute
  dgelsd_(M, N, NRHS, A, M, B, LDB, workv2().data(), RCOND, &RANK, workv().data(), &LWORK, &iworkv().front(), INFO);
}

/// Calls appropriate LAPACK function for solving systems of linear equations Ax=b, where A is symmetric indefinite
void LinAlgd::sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{ 
  INTEGER LWORK = -1;
  double WORK_QUERY;

  // first, determine workspace size
  dsysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDA, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // call LAPACK
  dsysv_(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for svd (divide and conquer)
void LinAlgd::gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* S, DOUBLE* U, INTEGER* LDU, DOUBLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  double WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;
  iworkv().resize(8*minmn);

  // call LAPACK to determine the optimal workspace size
  dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, &iworkv().front(), INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK); 

  // call LAPACK once again
  dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, V, LDVT, workv().data(), &LWORK, &iworkv().front(), INFO);
}

/// Calls LAPACK function for svd
void LinAlgd::gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* S, DOUBLE* U, INTEGER* LDU, DOUBLE* V, INTEGER* LDVT, INTEGER* INFO)
{
  double WORK_QUERY;
  INTEGER minmn = std::min(*M, *N);
  INTEGER LWORK = -1;

  // call LAPACK to determine the optimal workspace size
  dgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, &WORK_QUERY, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK); 

  // call LAPACK once again
  dgesvd_(JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDVT, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for computing eigenvalues and eigenvectors
void LinAlgd::syevd_(char* JOBZ, char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* EVALS, INTEGER* INFO)
{
  // do work query first
  double WORK_QUERY;
  INTEGER LWORK = -1;
  INTEGER IWORK_QUERY;
  INTEGER LIWORK = -1;
  dsyevd_(JOBZ, UPLO, N, A, LDA, EVALS, &WORK_QUERY, &LWORK, &IWORK_QUERY, &LIWORK, INFO);

  // set array sizes
  LWORK = (INTEGER) WORK_QUERY;
  LIWORK = IWORK_QUERY;
  workv().resize(LWORK);
  iworkv().resize(LIWORK);

  dsyevd_(JOBZ, UPLO, N, A, LDA, EVALS, workv().data(), &LWORK, &iworkv().front(), &LIWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations using LU factorization
void LinAlgd::gesv_(INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* X, INTEGER* LDX, INTEGER* INFO)
{
  dgesv_(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

/// Calls LAPACK function for solving a system of linear equation from a Cholesky factorization
void LinAlgd::potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dpotrs_(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
}

/// Calls LAPACK function for computing matrix inverse using Cholesky factorization
void LinAlgd::potri_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* INFO)
{
  dpotri_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for Cholesky factorization
void LinAlgd::potrf_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* INFO)
{
  dpotrf_(UPLO, N, A, LDA, INFO);
}

/// Calls LAPACK function for solving system of equations Ax=b, where A is PSD
void LinAlgd::posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dposv_(UPLO, N, NRHS, A, N, B, LDB, INFO);
}

/// Calls LAPACK function for inverting symmetric indefinite matrix using factorization
void LinAlgd::sytri_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  workv().resize(*N);
  dsytri_(UPLO, N, A, LDA, IPIV, workv().data(), INFO);
}

/// Calls LAPACK function for factorizing symmetric, indefinite matrix 
void LinAlgd::sytrf_(char* UPLO, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER *INFO) 
{
  // perform workspace query for factorization
  double WORK_QUERY;
  INTEGER LWORK = -1;
  dsytrf_(UPLO, N, A, LDA, IPIV, &WORK_QUERY, &LWORK, INFO);
  assert(*INFO >= 0);
  if (*INFO > 0)
    return; 

  // setup WORK array
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // perform the necessary factorization
  dsytrf_(UPLO, N, A, LDA, IPIV, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
void LinAlgd::geqp3_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* JPVT, DOUBLE* TAU, INTEGER* INFO)
{
  // determine workspace size
  INTEGER LWORK = -1;
  DOUBLE work_query;
  dgeqp3_(M, N, A, LDA, JPVT, TAU, &work_query, &LWORK, INFO);

  // setup workspace
  LWORK = (INTEGER) work_query;
  workv().resize(LWORK);

  // do QR factorization
  dgeqp3_(M, N, A, LDA, JPVT, TAU, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for QR factorization
void LinAlgd::geqrf_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, INTEGER* INFO)
{
  // determine LWORK
  const char* NAME = "DGEQRF";
  const char* OPTS = " ";
  INTEGER ISPEC = 1L;
  INTEGER TMP = -1;
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, M, N, &TMP, &TMP, strlen(NAME), strlen(OPTS));
  assert(NB > 0);

  // setup WORK vectors
  INTEGER LWORK = NB*(*N);
  workv().resize(LWORK);
  dgeqrf_(M, N, A, LDA, TAU, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for LU factorization 
void LinAlgd::getrf_(INTEGER* M, INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{
  dgetrf_(M, N, A, LDA, IPIV, INFO);
}

/// Calls LAPACK function for matrix inversion using LU factorization
void LinAlgd::getri_(INTEGER* N, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO)
{ 
  // compute block size
  INTEGER TMP = -1;
  INTEGER ISPEC = 1;
  const char* NAME = "DGETRI";
  const char* OPTS = " ";
  INTEGER NB = ilaenv_(&ISPEC, (char*) NAME, (char*) OPTS, N, &TMP, &TMP, &TMP, strlen(NAME), strlen(OPTS));

  // setup LAPACK parameters
  INTEGER LWORK = (*N)*NB;

  // setup work vector 
  workv().resize(LWORK);
  
  dgetri_(N, A, LDA, IPIV, workv().data(), &LWORK, INFO);
}

/// Calls LAPACK function for solving system of linear equations with a given LU factorization
void LinAlgd::getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, DOUBLE* A, INTEGER* LDA, INTEGER* IPIV, DOUBLE* B, INTEGER* LDB, INTEGER* INFO)
{
  dgetrs_(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
}

/// Calls LAPACK function for forming Q from a QR factorization
void LinAlgd::orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, DOUBLE* A, INTEGER* LDA, DOUBLE* TAU, INTEGER* INFO)
{
  // do a workspace query
  INTEGER LWORK = -1;
  DOUBLE WORK_QUERY;
  dorgqr_(M, N, K, A, LDA, TAU, &WORK_QUERY, &LWORK, INFO);

  // initialize the work array
  LWORK = (INTEGER) WORK_QUERY;
  workv().resize(LWORK);

  // call the function for real
  dorgqr_(M, N, K, A, LDA, TAU, workv().data(), &LWORK, INFO);
}

#include <Ravelin/ddefs.h>
#include "LinAlg.cpp"
#include <Ravelin/undefs.h>

// external methods
extern int solve_superlu(bool notrans, bool CRS, int m, int n, int nrhs, int nnz, int* col_indices, int* row_ptr, double* A_nz, double* x, double* b);

/// Does a LU factorization of a sparse matrix
VectorNd& LinAlgd::solve_sparse_direct(const SparseMatrixNd& A, const VectorNd& b, Transposition trans, VectorNd& x)
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
  double* nz_val = (double*) A.get_data();
  int* col_indices = (int*) A.get_indices();
  int* row_ptr = (int*) A.get_ptr();

  // check info
  int info = solve_superlu(trans == eNoTranspose, A.get_storage_type() == SparseMatrixNd::eCSR, m, n, nrhs, nnz, col_indices, row_ptr, nz_val, x.data(), (double*) b.data());
  if (info > 0)
    throw SingularException();
  #else
  throw std::runtime_error("Ravelin not built with SuperLU support!");
  #endif

  return x;
}

/// Does a LU factorization of a sparse matrix
MatrixNd& LinAlgd::solve_sparse_direct(const SparseMatrixNd& A, const MatrixNd& B, Transposition trans, MatrixNd& X)
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
  double* nz_val = (double*) A.get_data();
  int* col_indices = (int*) A.get_indices();
  int* row_ptr = (int*) A.get_ptr();

  // check info
  int info = solve_superlu(trans == eNoTranspose, A.get_storage_type() == SparseMatrixNd::eCSR, m, n, nrhs, nnz, col_indices, row_ptr, nz_val, X.data(), (double*) X.data());
  if (info > 0)
    throw SingularException();
  #else
  throw std::runtime_error("Ravelin not built with SuperLU support!");
  #endif

  return X;
}

