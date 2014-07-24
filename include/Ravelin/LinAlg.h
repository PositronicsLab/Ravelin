/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef LINALG
#error This class is not to be included by the user directly. Use LinAlgf.h or LinAlgd.h instead.
#endif

/// Linear algebra routines
/**
LinAlg is a set of static routines that interface to LAPACK.  I have included only very few routines here, however they should be some of the most utilized: SVD, (SVD-based) pseudo-inverse, linear equation solving, and matrix inverse.
*/
class LINALG
{
  public:
    enum SVD { eSVD1, eSVD2 };

  private:
    static REAL log2(REAL x);
    static void lartg_(REAL* F, REAL* G, REAL* CS, REAL* SN, REAL* R);
    static void gtsv_(INTEGER* N, INTEGER* NRHS, REAL* DL, REAL* D, REAL* DU, REAL* B, INTEGER* LDB, INTEGER* INFO);
    static void trtrs_(char* UPLO, char* TRANS, INTEGER* N, INTEGER* NRHS, REAL* AP, INTEGER* LDA, REAL* B, INTEGER* LDB, INTEGER* INFO);
    void ormqr_(char* SIDE, char* TRANS, INTEGER* M, INTEGER* N, INTEGER* K, REAL* A, INTEGER* LDA, REAL* TAU, REAL* C, INTEGER* LDC, INTEGER* INFO);
    static void sptrf_(char* UPLO, INTEGER* N, REAL* AP, INTEGER* IPIV, INTEGER* INFO);
    static void sptrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, REAL* AP, INTEGER* IPIV, REAL* B, INTEGER* LDB, INTEGER* INFO);
    void gelsd_(INTEGER* M, INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER*
LDA, REAL* B, INTEGER* LDB, REAL* RCOND, INTEGER* INFO);
    void sysv_(char* UPLO, INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER* LDA, INTEGER* IPIV, REAL* B, INTEGER* LDB, INTEGER* INFO);
    void gesdd_(char* JOBZ, INTEGER* M, INTEGER* N, REAL* A, INTEGER* LDA, REAL* S, REAL* U, INTEGER* LDU, REAL* V, INTEGER* LDVT, INTEGER* INFO);
    void gesvd_(char* JOBU, char* JOBV, INTEGER* M, INTEGER* N, REAL* A, INTEGER* LDA, REAL* S, REAL* U, INTEGER* LDU, REAL* V, INTEGER* LDVT, INTEGER* INFO);
    void syevd_(char* JOBZ, char* UPLO, INTEGER* N, REAL* A, INTEGER* LDA, REAL* EVALS, INTEGER* INFO);
    static void gesv_(INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER* LDA, INTEGER* IPIV, REAL* X, INTEGER* LDX, INTEGER* INFO);
    static void potrs_(char* UPLO, INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER* LDA, REAL* B, INTEGER* LDB, INTEGER* INFO);
    static void potri_(char* UPLO, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* INFO);
    static void potrf_(char* UPLO, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* INFO);
    static void posv_(char* UPLO, INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER* LDA, REAL* B, INTEGER* LDB, INTEGER* INFO);
    void sytri_(char* UPLO, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO);
    void sytrf_(char* UPLO, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO);
    void geqp3_(INTEGER* M, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* JPVT, REAL* TAU, INTEGER* INFO);
    void geqrf_(INTEGER* M, INTEGER* N, REAL* A, INTEGER* LDA, REAL* TAU, INTEGER* INFO);
    static void getrf_(INTEGER* M, INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO);
    void getri_(INTEGER* N, REAL* A, INTEGER* LDA, INTEGER* IPIV, INTEGER* INFO);
    static void getrs_(char* TRANS, INTEGER* N, INTEGER* NRHS, REAL* A, INTEGER* LDA, INTEGER* IPIV, REAL* B, INTEGER* LDB, INTEGER* INFO);
    void orgqr_(INTEGER* M, INTEGER* N, INTEGER* K, REAL* A, INTEGER* LDA, REAL* TAU, INTEGER* INFO);

  public:
    void compress();
    void free_memory();
    static void factor_LDL(MATRIXN& M, std::vector<int>& IPIV);
    MATRIXN& pseudo_invert(MATRIXN& A, REAL tol=(REAL) -1.0);
    static void givens(REAL a, REAL b, REAL& c, REAL& s);
    static MATRIX2 givens(REAL c, REAL s);
    static void householder(REAL alpha, const VECTORN& x, REAL& tau, VECTORN& v);
    static VECTORN& solve_sparse_direct(const SPARSEMATRIXN& A, const VECTORN& b, Transposition trans, VECTORN& x);
    static MATRIXN& solve_sparse_direct(const SPARSEMATRIXN& A, const MATRIXN& B, Transposition trans, MATRIXN& X);
    void update_QR_rank1(MATRIXN& Q, MATRIXN& R, const VECTORN& u, const VECTORN& v);
    void update_QR_delete_cols(MATRIXN& Q, MATRIXN& R, unsigned k, unsigned p);
    void update_QR_insert_cols(MATRIXN& Q, MATRIXN& R, MATRIXN& U, unsigned k);
    void update_QR_insert_rows(MATRIXN& Q, MATRIXN& R, MATRIXN& U, unsigned k);
    void update_QR_delete_rows(MATRIXN& Q, MATRIXN& R, unsigned k, unsigned p);

    /// work matrices
    FastThreadable<MATRIXN> workM, workM2;

    /// work matrix (for SVD)
    FastThreadable<MATRIXN> U;

    /// work STL integer vector (pivoting)
    FastThreadable<std::vector<INTEGER> > pivwork;

    /// work matrix (for SVD)
    FastThreadable<MATRIXN> V;

    /// work vector (for SVD/eigenvalues)
    FastThreadable<VECTORN> S;

    /// work vectors (LAPACK routines)
    FastThreadable<VECTORN> workv, workv2;

    /// work STL integer vector (LAPACK routines)
    FastThreadable<std::vector<INTEGER> > iworkv;

    // include templated routines here...
    #include "LinAlg.inl"

}; // end class

