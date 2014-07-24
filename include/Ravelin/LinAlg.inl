/// Solves a tridiagonal system
/**
 * \param dl the (n-1) elements on the subdiagonal (destroyed on return)
 * \param d  the n elements on the diagonal (destroyed on return)
 * \param du the n elements on the superdiagonal (destroyed on return)
 * \param XB the right hand side on input, the solution on return
 * \return the solution
 * throws SingularException if the matrix is singular
 */
template <class X>
static X& solve_tridiagonal_fast(VECTORN& dl, VECTORN& d, VECTORN& du, X& XB)
{
  // make sure everything is the proper size
  #ifndef NEXCEPT
  if (sizeof(dl.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  if (sizeof(d.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  if (sizeof(du.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // determine N
  INTEGER N = d.size();
  #ifndef NEXCEPT
  if (XB.rows() != N)
    throw MissizeException();
  #endif

  // call the tridiagonal solver
  INTEGER LDB = XB.leading_dim();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;
  gtsv_(&N, &NRHS, dl.data(), d.data(), du.data(), XB.data(), &LDB, &INFO);

  // see whether solution successful
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Performs the Cholesky factorization of a matrix
/**
 * \param A the matrix A on input; the factorized (upper triangular) matrix on output
 * \return <b>true</b> if matrix factored successfully, <b>false</b> otherwise
 */
template <class X>
static bool factor_chol(X& A)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
    return true;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER INFO;

  // perform the Cholesky factorization
  potrf_(&UPLO, &N, A.data(), &LDA, &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    return false;

  // make the matrix upper triangular
  A.zero_lower_triangle();
  return true;
}

/// Performs the LU factorization of a matrix
/**
 * \param A the m x n matrix to be factored; on exit, the factors L and U
 *        from the factorization M = P*L*U (unit diagonal elements of L are
 *        not stored)
 * \param pivwork on output, contains the pivot indices (for 1 <= i <= min(m,n),
 *        row i of A was interchanged with row pivwork[i])
 * \return <b>false</b> if A is singular, <b>true</b> otherwise
 */
template <class X>
static bool factor_LU(X& A, std::vector<int>& pivwork)
{
  if (A.rows() == 0 || A.columns() == 0)
    return true;

  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER LDA = A.leading_dim();  

  // resize pivwork if necessary
  if (pivwork.size() < std::min(A.rows(), A.columns()))
    pivwork.resize(std::min(A.rows(),A.columns()));

  // call LAPACK
  INTEGER INFO;
  getrf_(&M, &N, A.data(), &LDA, &pivwork.front(), &INFO);
  return INFO == 0;
}

/// Computes a matrix inverse using the factorization determined via factor_LU()
/**
 * \param M the LU-factored matrix
 * \param pivwork the pivots determined by factor_LU()
 * \note throws SingularException if matrix is singular
 */
template <class X>
X& inverse_LU(X& M, const std::vector<int>& pivwork)
{
  if (M.rows() == 0 || M.columns() == 0)
    return M;

  #ifndef NEXCEPT
  if (M.rows() != M.columns())
    throw NonsquareMatrixException();
  #endif

  // call lapack
  INTEGER INFO;
  INTEGER N = M.rows();
  INTEGER LDM = M.leading_dim();
  getri_(&N, M.data(), &LDM, (INTEGER*) &pivwork.front(), &INFO);
  if (INFO > 0)
    throw SingularException();

  return M;
}

/// Solves a triangular system of linear equations
/**
 * \param A the matrix
 * \param utri if <b>true</b> A is upper triangular (lower triangular otherwise)
 * \param transpose_A if <b>true</b>, solves A'*x = b
 * \param XB contains b on entry, x on return
 * \return reference to XB
 */
template <class X, class Y>
static X& solve_tri_fast(Y& A, bool utri, bool transpose_A, X& XB)
{
  #ifndef NEXCEPT
  if (A.rows() != XB.rows())
    throw MissizeException();

  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  if (A.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // setup parameters for LAPACK
  char TRANS = (transpose_A) ? 'T' : 'N';
  char UPLO = (utri) ? 'U' : 'L';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER LDB = XB.leading_dim();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;
  trtrs_(&UPLO, &TRANS, &N, &NRHS, A.data(), &LDA, XB.data(), &LDB, &INFO);

  return XB;
}

/// Solves systems of linear equations using the factorization determined via factor_LDL()
/**
 * \param M the factorization performed by factor_LDL()
 * \param XB the right hand sides on input, the vectors X on return
 */
template <class X>
static X& solve_LDL_fast(const MATRIXN& M, const std::vector<int>& pivwork, X& XB)
{
  #ifndef NEXCEPT
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();
  #endif

  // check for empty matrix
  if (M.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(M.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.rows();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = XB.leading_dim();
  INTEGER INFO;

  // call the solver routine
  sptrs_(&UPLO, &N, &NRHS, M.data(), (int*) &pivwork.front(), XB.data(), &LDB, &INFO);
  assert(INFO == 0);

  return XB;
}

/// Solves a system of linear equations using the factorization determined via factor_chol()
/**
 * \param M the Cholesky decomposition performed by factor_chol()
 * \param XB the right hand sides on input, the vectors X on return
 */
template <class X, class Y>
static X& solve_chol_fast(const Y& M, X& XB)
{
  #ifndef NEXCEPT
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();
  #endif

  // check for empty matrices
  if (M.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(M.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // setup parameters for LAPACK
  char UPLO = 'U';
  INTEGER N = M.rows();
  INTEGER NRHS = XB.columns();
  INTEGER LDA = M.leading_dim();
  INTEGER LDB = XB.leading_dim();
  INTEGER INFO;

  // call the solver routine
  potrs_(&UPLO, &N, &NRHS, (REAL*) M.data(), &LDA, XB.data(), &LDB, &INFO);
  assert(INFO == 0);

  return XB;
}

/// Solves a system of linear equations using the factorization determined via factor_LU()
/**
 * \param M the LU factorization performed by factor_LU()
 * \param XB the right hand side on input, the matrix X on return
 * \param transpose if <b>true</b>, solves M'X = B
 * \param pivwork pivots computed by factor_LU()
 */
template <class Y, class X>
static X& solve_LU_fast(const Y& M, bool transpose, const std::vector<int>& pivwork, X& XB)
{
  #ifndef NEXCEPT
  if (M.rows() != XB.rows())
    throw MissizeException();

  if (M.rows() != M.columns())
    throw NonsquareMatrixException();
  #endif

  // check for empty matrix
  if (M.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(M.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // setup parameters to LAPACK
  INTEGER N = M.rows();
  INTEGER LDM = M.leading_dim();
  INTEGER LDB = XB.leading_dim();
  INTEGER NRHS = XB.columns();
  char TRANS = (transpose) ? 'T' : 'N';
  
  // call LAPACK
  INTEGER INFO;
  getrs_(&TRANS, &N, &NRHS, (REAL*) M.data(), &LDM, (INTEGER*) &pivwork.front(), XB.data(), &LDB, &INFO);

  return XB;
}

/// Calculates the rank of a matrix
template <class X>
unsigned calc_rank(X& A, REAL tol = (REAL) -1.0)
{
  // look for easy out
  if (A.rows() == 0 || A.columns() == 0)
    return 0;

  // compute the SVD of A
  MATRIXN& Ux = U();
  MATRIXN& Vx = V();
  VECTORN& Sx = S();
  svd(A, Ux, Sx, Vx);

  // get the dimensions of A
  unsigned m = A.rows();
  unsigned n = A.columns();
  unsigned maxmn = std::max(m,n);
  unsigned minmn = std::min(m,n);

  // get the # of singular values
  REAL* Sx_data = Sx.data();
  REAL tolerance = (tol < 0.0) ? Sx_data[0] * maxmn * std::numeric_limits<REAL>::epsilon() : tol;
  unsigned ns = 0;
  for (unsigned i=Sx.size()-1; i > 0; i--, ns++)
    if (Sx_data[i] > tolerance)
      break;

  assert(minmn > ns);
  return minmn - ns;
}

/// Computes the nullspace of a matrix
/**
 * \note A is destroyed on return
 */
template <class Y>
MATRIXN& nullspace(Y& A, MATRIXN& nullspace, REAL tol = -1.0)
{
  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(nullspace.data()))
    throw DataMismatchException();
  #endif

  // look for fast way out
  if (A.rows() == 0)
  {
    nullspace.set_zero(A.columns(), A.columns());
    for (unsigned i=0; i< A.columns(); i++)
      nullspace(i,i) = (REAL) 1.0;
    return nullspace;
  }

  // compute the SVD of A
  MATRIXN& Ux = U();
  VECTORN& Sx = S();
  svd(A, Ux, Sx, nullspace);

  // get the dimensions of A
  unsigned m = A.rows();
  unsigned n = A.columns();
  boost::tuple<unsigned, unsigned> min_max = boost::minmax(m, n);
  unsigned minmn = min_max.get<0>();
  unsigned maxmn = min_max.get<1>();

  // get the # of singular values
  unsigned ns = 0;
  REAL* Sx_data = Sx.data();
  REAL tolerance = (tol < 0.0) ? Sx_data[0] * maxmn * std::numeric_limits<REAL>::epsilon() : tol;
  for (unsigned i=Sx.size()-1; i > 0; i--, ns++)
    if (Sx_data[i] > tolerance)
      break;

  // add in # of singular values from non-square matrices
  for (unsigned i=minmn; i < n; i++)
    ns++;

  // shift nullspace (if necesary)
  if (ns == nullspace.columns())
    return nullspace;
  else
  {
    COLUMN_ITERATOR bi = nullspace.block_column_iterator_begin(0, nullspace.rows(), nullspace.columns()-ns, nullspace.columns());
    std::copy(bi, bi+ns*nullspace.rows(), nullspace.column_iterator_begin());
    nullspace.resize(nullspace.rows(), ns, true);
    return nullspace;
  }
}

/// Computes the condition number of a matrix
template <class X>
REAL cond(X& A)
{
  svd(A, U(), S(), V());
  REAL* S_data = S().data();
  return S_data[0] / S_data[S().size()-1];
}

/// Determines whether a symmetric matrix is positive semi-definite
template <class X>
bool is_SPSD(X& m, REAL tol)
{
  VECTORN& evals = S();
  eig_symm(m, evals);

  REAL* ev = evals.data();
  // make tolerance positive, if it is not already
  if (tol < (REAL) 0.0)
    tol = std::fabs(ev[evals.size()-1]) * m.rows() * std::numeric_limits<REAL>::max();

  // check whether all eigenvalues are non-negative to numerical tolerance
  for (unsigned i=0; i< evals.size(); i++)
    if (ev[i] < -tol)
      return false;

  return true;
}

/// Determines whether a matrix is positive-definite
template <class X>
bool is_SPD(X& m, REAL tol)
{
  // get the eigenvalues of the matrix
  VECTORN& evals = S();
  eig_symm(m, evals);

  // make tolerance positive, if it is not already
  REAL* ev = evals.data();
  if (tol < (REAL) 0.0)
    tol = std::fabs(ev[evals.size()-1]) * m.rows() * std::numeric_limits<REAL>::max();

  // check whether all eigenvalues are positive to numerical tolerance
  for (unsigned i=0; i< evals.size(); i++)
    if (ev[i] < tol)
      return false;

  return true;
}

/// Computes the eigenvalues of the matrix A
/**
 * \param A a matrix
 * \param evals on return, the eigenvalues will be stored here in ascending order
 */
template <class X, class Y>
void eig_symm(X& A, Y& evals)
{
  // make sure A is not zero sized
  if (A.rows() == 0)
  {
    evals.resize(0);
    return;
  }

  // verify that A is square
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that the matrix is symmetric
  #ifndef NEXCEPT
  REAL tol = std::numeric_limits<REAL>::epsilon() * A.rows() * A.norm_inf();
  if (!A.is_symmetric(tol))
    std::cerr << "LinAlg::eig_symm() - matrix does not appear to be symmetric!" << std::endl;
  #endif

  // make sure that the eigenvalues array is the proper size
  evals.resize(A.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'N';
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, A.data(), &LDA, evals.data(), &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

/// Computes the eigenvalues and eigenvectors of the matrix A
/**
 * \param A a square symmetric matrix on input, eigenvectors corresponding to eigenvalues on return
 * \param evals on return, the eigenvalues will be stored here in ascending order
 */
template <class X, class Y>
void eig_symm_plus(X& A_evecs, Y& evals)
{
  // make sure that A is not zero sized
  if (A_evecs.rows() == 0 || A_evecs.columns() == 0)
  {
    evals.resize(0);
    A_evecs.resize(0,0);
    return;
  }

  #ifndef NEXCEPT
  if (A_evecs.rows() != A_evecs.columns())
    throw NonsquareMatrixException();

  if (sizeof(A_evecs.data()) != sizeof(evals.data()))
    throw DataMismatchException();

  // verify that the matrix is symmetric
  const REAL tol = std::numeric_limits<REAL>::epsilon() * A_evecs.norm_inf() * A_evecs.rows();
  if (!A_evecs.is_symmetric(tol))
    std::cerr << "LinAlg::eig_symm() - matrix does not appear to be symmetric!" << std::endl;
  #endif

  // make sure that the eigenvalues array is the proper size
  evals.resize(A_evecs.rows()); 

  // form inputs to LAPACK
  char JOBZ = 'V';
  char UPLO = 'U';
  INTEGER N = A_evecs.rows();
  INTEGER LDA = A_evecs.leading_dim();
  INTEGER INFO;

  // call LAPACK
  syevd_(&JOBZ, &UPLO, &N, A_evecs.data(), &LDA, evals.data(), &INFO);
  assert(INFO == 0);
  
  if (INFO > 0)
    throw NumericalException("Eigenvalue/eigenvector determination did not converge");
}

template <class X, class MatU, class VecS, class MatV>
void svd(X& A, MatU& U, VecS& S, MatV& V)
{
  MATRIXN& A_backup = workM();

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(U.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(S.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(V.data()))
    throw DataMismatchException();
  #endif

  // copy A
  A_backup = A;

  try
  {
    svd1(A, U, S, V);
  }
  catch (NumericalException e)
  {
    A = A_backup;
    svd2(A, U, S, V); 
  }
}

/// Does an 'in place' SVD (destroying A)
/**
 * The singular value decomposition of A is U*S*V' (' is the transpose 
 * operator); to recompose A, it will be necessary to transpose V before
 * multiplication (i.e., V is returned by the algorithm, not V').
 * Note: passed matrices and vectors U, S, and V are resized as necessary. 
 * \param A the matrix on which the SVD will be performed (destroyed on return)
 * \param U on output, a A.rows() x A.rows() orthogonal matrix
 * \param S on output, a min(A.rows(), A.columns()) length vector of singular values
 * \param V on output, a A.columns() x A.columns() orthogonal matrix
 */
template <class X, class MatU, class VecS, class MatV>
void svd1(X& A, MatU& U, VecS& S, MatV& V)
{
  // make sure that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
  {
    U.set_zero(A.rows(), A.rows());
    S.resize(0);
    V.set_zero(A.columns(), A.columns());
    return;
  } 

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(U.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(S.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(V.data()))
    throw DataMismatchException();
  #endif

  // setup U
  if (U.rows() != A.rows() || U.columns() != A.rows())
    U.resize(A.rows(), A.rows());

  // setup S
  unsigned minmn = std::min(A.rows(), A.columns());
  if (S.size() != minmn)
    S.resize(minmn);

  // setup V
  if (V.rows() != A.columns() || V.columns() != A.columns())
    V.resize(A.columns(), A.columns());

  // setup call to LAPACK
  char JOB = 'A';
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER LDA = A.leading_dim();
  INTEGER LDU = U.leading_dim();
  INTEGER LDVT = V.leading_dim();
  INTEGER INFO;

  // call LAPACK 
  gesvd_(&JOB, &JOB, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Singular value decomposition failed to converge");

  // transpose V
  V.transpose();
}

/// Does an 'in place' SVD (destroying A), using divide and conquer algorithm
/**
 * The singular value decomposition of A is U*S*V' (' is the transpose 
 * operator); to recompose A, it will be necessary to transpose V before
 * multiplication (i.e., V is returned by the algorithm, not V').
 * Note: passed matrices and vectors U, S, and V are resized as necessary. 
 * \param A the matrix on which the SVD will be performed (destroyed on return)
 * \param U on output, a A.rows() x A.rows() orthogonal matrix
 * \param S on output, a min(A.rows(), A.columns()) length vector of singular values
 * \param V on output, a A.columns() x A.columns() orthogonal matrix
 */
template <class X, class MatU, class VecS, class MatV>
void svd2(X& A, MatU& U, VecS& S, MatV& V)
{
  // make sure that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
  {
    U.set_zero(A.rows(), A.rows());
    S.resize(0);
    V.set_zero(A.columns(), A.columns());
    return;
  } 

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(U.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(S.data()))
    throw DataMismatchException();
  if (sizeof(A.data()) != sizeof(V.data()))
    throw DataMismatchException();
  #endif

  // setup U
  if (U.rows() != A.rows() || U.columns() != A.rows())
    U.resize(A.rows(), A.rows());

  // setup S
  unsigned minmn = std::min(A.rows(), A.columns());
  if (S.size() != minmn)
    S.resize(minmn);

  // setup V
  if (V.rows() != A.columns() || V.columns() != A.columns())
    V.resize(A.columns(), A.columns());

  // setup call to LAPACK
  char JOBZ = 'A';
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER LDA = A.leading_dim();
  INTEGER LDU = U.leading_dim();
  INTEGER LDVT = V.leading_dim();
  INTEGER INFO;

  // call LAPACK 
  gesdd_(&JOBZ, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, V.data(), &LDVT, &INFO);
  assert(INFO >= 0);

  if (INFO > 0)
    throw NumericalException("Singular value decomposition failed to converge");

  // transpose V
  V.transpose();
}

/// Solves a symmetric, indefinite square matrix
/**
 * \param A the matrix to be solved; the matrix is destroyed on return
 * \param XB the RHS B (A*X = B) on input; the solution, X, on return
 */
template <class X>
X& solve_symmetric_fast(MATRIXN& A, X& XB)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify A and b are compatible 
  if (A.columns() != XB.rows())
    throw MissizeException();
  #endif

  // make sure that A is not zero sized
  if (A.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = XB.leading_dim();
  pivwork().resize(N);
  INTEGER INFO;

  // call LAPACK
  sysv_(&UPLO, &N, &NRHS, A.data(), &LDA, &pivwork().front(), XB.data(), &LDB, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Inverts the symmetric, indefinite matrix A
/**
 * \param A a square, symmetric indefinite matrix; inverse will be contained
 *        here on return
 */
template <class X>
X& inverse_symmetric(X& A)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero size
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER INFO;

  // setup pivot array
  pivwork().resize(N);

  // perform the necessary factorization
  sytrf_(&UPLO, &N, A.data(), &LDA, &pivwork().front(), &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // perform the inversion
  sytri_(&UPLO, &N, A.data(), &LDA, &pivwork().front(), &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // make the matrix symmetric
  REAL* data = A.data();
  for (unsigned i=1, ii=A.leading_dim(); i< N; i++, ii += A.leading_dim())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.leading_dim())
      data[jj+i] = data[ii+j];

  return A;
}

/// Solves a system of equations A*X = B using a symmetric, positive-definite square matrix
/**
 * \param A the matrix coefficients; this matrix will be destroyed on return
 * \param XB on input B, on output X
 */
template <class Y, class X>
static X& solve_SPD_fast(Y& A, X& XB)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();

  // verify that A and b are proper size
  if (A.columns() != XB.rows())
    throw MissizeException();
  #endif

  // verify that A is not zero size
  if (A.rows() == 0 || XB.columns() == 0)
    return XB.set_zero();

  #ifndef NEXCEPT
  if (sizeof(A.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // form inputs to LAPACK
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER LDB = XB.leading_dim();
  INTEGER NRHS = XB.columns();
  INTEGER INFO;

  // call LAPACK
  posv_(&UPLO, &N, &NRHS, A.data(), &LDA, XB.data(), &LDB, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

/// Inverts the symmetric, positive-definite matrix A using Cholesky factorization
/**
 * \param A the Cholesky factorization of a matrix; contains the inverse
 *        on return
 */
template <class X>
static X& inverse_chol(X& A)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero sized
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER INFO;

  // perform the inverse
  potri_(&UPLO, &N, A.data(), &LDA, &INFO);
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  // make the matrix symmetric
  REAL* data = A.data();
  for (unsigned i=1, ii=A.leading_dim(); i< N; i++, ii += A.leading_dim())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.leading_dim())
      data[jj+i] = data[ii+j];

  return A;
}

/// Inverts the symmetric, positive-definite matrix A using Cholesky factorization
/**
 * \param A a square, symmetric positive-definite matrix; contains the inverse
 *        on return
 */
template <class X>
static X& inverse_SPD(X& A)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero sized
  if (A.rows() == 0)
    return A;

  // setup LAPACK args for Cholesky factorization
  char UPLO = 'U';
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER INFO;

  // perform the Cholesky factorization
  potrf_(&UPLO, &N, A.data(), &LDA, &INFO);
  assert(INFO >= 0);

  // perform the inverse
  potri_(&UPLO, &N, A.data(), &LDA, &INFO);
  if (INFO > 0)
    throw SingularException();

  // now, make the matrix symmetric
  REAL* data = A.data();
  for (unsigned i=1, ii=A.leading_dim(); i< N; i++, ii += A.leading_dim())
    for (unsigned j=0, jj=0; j< i; j++, jj += A.leading_dim())
      data[jj+i] = data[ii+j];

  return A;
}

/// Inverts the matrix A using LU factorization
/**
 * \param A a square matrix; contains the inverse on return
 */
template <class X>
X& invert(X& A)
{
  factor_LU(A, pivwork());
  inverse_LU(A, pivwork());

  return A;
}

/// Most robust system of linear equations solver (solves AX = B)
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param A the coefficient matrix (destroyed on return)
 * \param XB the matrix B on input, the matrix X on return
 * \param tol the tolerance for determining the rank of A; if tol < 0.0,
 *        tol is computed using machine epsilon
 */
template <class X, class Y, class Vec, class Z>
X& solve_LS_fast(const Y& U, const Vec& S, const Z& V, X& XB, REAL tol = (REAL) -1.0)
{
  // verify that U, S, V and B are appropriate sizes
  #ifndef NEXCEPT
  if (U.rows() != XB.rows())
    throw MissizeException();
  if (S.columns() > 1 || (S.size() != U.columns() && S.size() != V.rows()))
    throw MissizeException();
  if (sizeof(U.data()) != sizeof(XB.data()) || 
      sizeof(U.data()) != sizeof(S.data()) ||
      sizeof(U.data()) != sizeof(V.data()))
    throw DataMismatchException();
  #endif

  // get the dimensionality of A
  const unsigned m = U.rows();
  const unsigned n = V.columns();
  const unsigned k = XB.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy outs
  if (m == 0 || n == 0)
  {
    XB.set_zero(m,k);
    return XB;
  }

  // compute the svd
  VECTORN& Sx = workv();
  MATRIXN& workMx = workM();
  MATRIXN& workM2x = workM2();

  // determine new tolerance based on first std::singular value if necessary
  Sx.resize(S.rows());
  std::copy(S.begin(), S.end(), Sx.begin());
  REAL* Sx_data = Sx.data();
  if (tol < (REAL) 0.0)
    tol = Sx_data[0] * std::max(m,n) * std::numeric_limits<REAL>::epsilon();

  // compute 1/S
  unsigned S_len = Sx.size();

  // A is m x n, B is m x k
  // (L -> R, scaling V)    n^2 + n*min(n,m)*m + nmk [n < m < k, n < k < m]
  // (L -> R, scaling U')   m^2 + n*min(n,m)*m + nmk [m < n < k]
  // (R -> L, scaling U')   m^2 + m^2k + nmk + n^2k  [k < n < m]
  // (R -> L, scaling U'*B) m^2k + min(n,m)*k + n*min(m,n)*k [k < m < n, m < k < n]

  // compute inv(s) 
  for (unsigned i=0; i< S_len; i++)
    Sx_data[i] = (std::fabs(Sx_data[i]) > tol) ? (REAL) 1.0/Sx_data[i] : (REAL) 0.0;

  // check cases
  // case 1: n is smallest
  if (n <= m && n <= k)
  {
    // scale n columns of V
    workM2x = V;
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx_data[i], workM2x.data()+workM2x.leading_dim()*i, 1);

    // multiply scaled V by U' = workM
    workMx.resize(n, m);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, n, (REAL) 1.0, workM2x.data(), workM2x.leading_dim(), U.data(), U.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply workM * XB
    workM2x.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workMx.data(), workMx.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, workM2x.data(), workM2x.leading_dim());
    XB.resize(workM2x.rows(), workM2x.columns());
    std::copy(workM2x.row_iterator_begin(), workM2x.row_iterator_end(), XB.row_iterator_begin());
  }
  // case 2: m < n < k
  else if (m <= n && n <= k)
  {
    // scale columns of U
    workM2x = U;
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(m, Sx_data[i], workM2x.data()+workM2x.leading_dim()*i, 1);

    // multiply V by scaled U' = workM
    workMx.resize(n,m);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, m, (REAL) 1.0, V.data(), V.leading_dim(), workM2x.data(), workM2x.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply workM * XB
    workM2x.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workMx.data(), workMx.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, workM2x.data(), workM2x.leading_dim());
    XB.resize(workM2x.rows(), workM2x.columns());
    std::copy(workM2x.row_iterator_begin(), workM2x.row_iterator_end(), XB.row_iterator_begin());
  }
  // case 3: k < n < m
  else if (k <= n && n <= m)
  {
    // scale columns of U
    workM2x = U;
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx_data[i], workM2x.data()+workM2x.leading_dim()*i, 1);

    // multiply U' * XB (resulting in n x k matrix)
    workMx.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workM2x.data(), workM2x.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply V * workM
    V.mult(workMx, XB);
  }
  // case 4: n is largest
  else
  {
    assert(n >= m && n >= k);

    // scale m columns of V (n x m)
    workM2x = V;
    for (unsigned i=0; i< m; i++)
      workM2x.column(i) *= Sx_data[i];
    SHAREDMATRIXN sharedV = workM2x.block(0, n, 0, m);

    // multiply U' * XB (resulting in m x k matrix)
    U.transpose_mult(XB, workMx);

    // do the multiplication
    sharedV.mult(workMx, XB);

/*
    // scale m columns of V
    workM2x = V;
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(n, Sx_data[i], workM2x.data()+workM2x.leading_dim()*i, 1);

    // multiply U' * XB (resulting in m x k matrix)
    U.transpose_mult(XB, workMx);

    // multiply V * workM
    XB.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workM2x.data(), workM2x.leading_dim(), workMx.data(), workMx.leading_dim(), (REAL) 0.0, XB.data(), XB.leading_dim());
*/
  }

  return XB;

// NOTE: this is disabled b/c it does not work as well...

/*
  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = std::max(M,N);
  INTEGER INFO;

  // allocate storage for solving the RHS
  const INTEGER XB_rows = std::max(M,NRHS);
  const INTEGER XB_cols = std::max(N,NRHS);
  REAL* rhs = new REAL[XB_rows*XB_cols];
  std::copy(XB.data(), XB.end(), rhs);

  // compute
  gelsd_(&M, &N, &NRHS, A.data(), &M, rhs, &LDB, &stol, &INFO);

  // don't check success - SVD algorithm is very stable.. 

  // copy solution into XB
  XB.resize(N, NRHS);
  std::copy(rhs, rhs+(M*NRHS), XB.data());

  // mark A as destroyed and free memory
  A.resize(0,0);
  delete [] rhs;

  return XB;
*/
}


/// Most robust system of linear equations solver (solves AX = B)
/**
 * Solves rank-deficient and underdetermined (minimum norm solution) systems.
 * Computes least-squares solution to overdetermined systems.
 * \param A the coefficient matrix (destroyed on return)
 * \param XB the matrix B on input, the matrix X on return
 * \param tol the tolerance for determining the rank of A; if tol < 0.0,
 *        tol is computed using machine epsilon
 */
template <class X, class Y>
X& solve_LS_fast(Y& A, X& XB, SVD svd_algo, REAL tol)
{
  // verify that A and B are appropriate sizes
  #ifndef NEXCEPT
  if (A.rows() != XB.rows())
    throw MissizeException();

  if (sizeof(A.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // get the dimensionality of A
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  const unsigned k = XB.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy out
  if (m == 0 || n == 0)
  {
    XB.set_zero(n, XB.columns());
    return XB;
  }

  // compute the svd
  MATRIXN& Ux = U();
  MATRIXN& Vx = V();
  VECTORN& Sx = S();
  MATRIXN& workMx = workM();
  if (svd_algo == eSVD1)
    svd1(A, Ux, Sx, Vx);
  else
    svd2(A, Ux, Sx, Vx);  

  // determine new tolerance based on first std::singular value if necessary
  REAL* Sx_data = Sx.data();
  if (tol < (REAL) 0.0)
    tol = Sx_data[0] * std::max(m,n) * std::numeric_limits<REAL>::epsilon();

  // compute 1/S
  unsigned S_len = Sx.size();

  // A is m x n, B is m x k
  // (L -> R, scaling V)    n^2 + n*min(n,m)*m + nmk [n < m < k, n < k < m]
  // (L -> R, scaling U')   m^2 + n*min(n,m)*m + nmk [m < n < k]
  // (R -> L, scaling U')   m^2 + m^2k + nmk + n^2k  [k < n < m]
  // (R -> L, scaling U'*B) m^2k + min(n,m)*k + n*min(m,n)*k [k < m < n, m < k < n]

  // compute inv(s) 
  for (unsigned i=0; i< S_len; i++)
    Sx_data[i] = (std::fabs(Sx_data[i]) > tol) ? (REAL) 1.0/Sx_data[i] : (REAL) 0.0;

  // check cases
  // case 1: n is smallest
  if (n < m && n < k)
  {
    // scale n columns of V
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx_data[i], Vx.data()+Vx.leading_dim()*i, 1);

    // multiply scaled V by U' = workM
    workMx.resize(n, m);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, n, (REAL) 1.0, Vx.data(), Vx.leading_dim(), Ux.data(), Ux.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply workM * XB
    Vx.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workMx.data(), workMx.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, Vx.data(), Vx.leading_dim());
    XB.resize(Vx.rows(), Vx.columns());
    std::copy(Vx.row_iterator_begin(), Vx.row_iterator_end(), XB.row_iterator_begin());
  }
  // case 2: m < n < k
  else if (m < n && n < k)
  {
    // scale columns of U
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(m, Sx_data[i], Ux.data()+Ux.leading_dim()*i, 1);

    // multiply V by scaled U' = workM
    workMx.resize(n,m);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, m, (REAL) 1.0, Vx.data(), Vx.leading_dim(), Ux.data(), Ux.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply workM * XB
    Vx.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, workMx.data(), workMx.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, Vx.data(), Vx.leading_dim());
    XB.resize(Vx.rows(), Vx.columns());
    std::copy(Vx.row_iterator_begin(), Vx.row_iterator_end(), XB.row_iterator_begin());
  }
  // case 3: k < n < m
  else if (k < n && n < m)
  {
    // scale columns of U
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx_data[i], Ux.data()+Ux.leading_dim()*i, 1);

    // multiply U' * XB (resulting in n x k matrix)
    workMx.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasTrans, CblasNoTrans, n, k, m, (REAL) 1.0, Ux.data(), Ux.leading_dim(), XB.data(), XB.leading_dim(), (REAL) 0.0, workMx.data(), workMx.leading_dim());

    // multiply V * workM
    Vx.mult(workMx, XB);
  }
  // case 4: n is largest
  else
  {
    assert(n >= m && n >= k);

    // scale m columns of V
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(n, Sx_data[i], Vx.data()+Vx.leading_dim()*i, 1);

    // multiply U' * XB (resulting in m x k matrix)
    Ux.transpose_mult(XB, workMx);

    // multiply V * workM
    XB.resize(n,k);
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, m, (REAL) 1.0, Vx.data(), Vx.leading_dim(), workMx.data(), workMx.leading_dim(), (REAL) 0.0, XB.data(), XB.leading_dim());
  }

  return XB;

// NOTE: this is disabled b/c it does not work as well...

/*
  // setup LAPACK parameters
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER NRHS = XB.columns();
  INTEGER LDB = std::max(M,N);
  INTEGER INFO;

  // allocate storage for solving the RHS
  const INTEGER XB_rows = std::max(M,NRHS);
  const INTEGER XB_cols = std::max(N,NRHS);
  REAL* rhs = new REAL[XB_rows*XB_cols];
  std::copy(XB.data(), XB.end(), rhs);

  // compute
  gelsd_(&M, &N, &NRHS, A.data(), &M, rhs, &LDB, &stol, &INFO);

  // don't check success - SVD algorithm is very stable.. 

  // copy solution into XB
  XB.resize(N, NRHS);
  std::copy(rhs, rhs+(M*NRHS), XB.data());

  // mark A as destroyed and free memory
  A.resize(0,0);
  delete [] rhs;

  return XB;
*/
}

/// Solves the general system AX = B
/**
 * \param A a square matrix (destroyed on return)
 * \param XB the matrix B on input, the matrix X on return
 */
template <class X, class Y>
Y& solve_fast(X& A, Y& XB)
{  
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero size
  if (A.rows() == 0)
    return XB;

  // verify that A and b are compatible
  #ifndef NEXCEPT
  if (A.columns() != XB.rows())
    throw MissizeException();

  if (sizeof(A.data()) != sizeof(XB.data()))
    throw DataMismatchException();
  #endif

  // setup LAPACK parameters
  INTEGER N = A.rows();
  INTEGER LDA = A.leading_dim();
  INTEGER LDB = XB.leading_dim();
  INTEGER NRHS = XB.columns();
  pivwork().resize(N);
  INTEGER INFO;

  // call LAPACK (use solving routine that uses LU factorization)
  gesv_(&N, &NRHS, A.data(), &LDA, &pivwork().front(), XB.data(), &LDB, &INFO);

  // check for singularity
  assert(INFO >= 0);
  if (INFO > 0)
    throw SingularException();

  return XB;
}

template <class X, class Y>
X& solve_LS_fast1(Y& A, X& XB, REAL tol = (REAL) -1.0) 
{ 
  return solve_LS_fast(A, XB, eSVD1, tol); 
}     

template <class X, class Y>
X& solve_LS_fast2(Y& A, X& XB, REAL tol = (REAL) -1.0) 
{ 
  return solve_LS_fast(A, XB, eSVD2, tol); 
}     

/// Performs the QR factorization of a matrix with column pivoting
/**
 * Factorizes A*P = Q*R
 * \param AQ the matrix A on input; the matrix R on output
 * \param Q the matrix Q on output
 * \param PI the column pivots on output
 */
template <class ARMat, class QMat>
void factor_QR(ARMat& AR, QMat& Q, std::vector<int>& PI)
{
  // get matrix/vector
  VECTORN& tau = workv2();

  // setup constants
  const unsigned m = AR.rows();
  const unsigned n = AR.columns();

  // determine LAPACK parameters
  INTEGER M = AR.rows();
  INTEGER N = AR.columns();
  INTEGER MINMN = std::min(M, N);
  unsigned min_mn = (unsigned) std::min(M,N);

  // setup tau vector
  workv2().resize(min_mn);

  // setup PI for entry
  PI.resize(N);
  for (int i=0; i< N; i++)
    PI[i] = 0;

  // call LAPACK
  INTEGER LDA = AR.leading_dim();
  INTEGER INFO;
  geqp3_(&M, &N, AR.data(), &LDA, PI.data(), workv2().data(), &INFO);
  assert(INFO == 0);

  // setup Q
  Q.resize(m,m);
  std::copy(AR.data(), AR.data()+m*min_mn, Q.data());
  orgqr_(&M, &MINMN, &MINMN, Q.data(), &M, tau.data(), &INFO);

  // resize AR
  AR.resize(std::min(AR.rows(), AR.columns()), AR.columns(), true);

  // make R upper triangular
  for (unsigned i=0; i< AR.columns(); i++)
  {
    RowIteratord coli = AR.block_row_iterator_begin(i+1,AR.rows(),i,i+1);
    std::fill(coli, coli.end(), 0.0);
  }
}

/// Performs the QR factorization of a matrix
/**
 * \param AQ the m x n matrix A on input; the matrix min(m,n) x n R on output
 * \param Q the m x min(m,n) matrix Q on output
 */
template <class ARMat, class QMat>
void factor_QR(ARMat& AR, QMat& Q)
{
  // get matrix/vector
  VECTORN& tau = workv2();

  // setup constants
  const unsigned m = AR.rows();
  const unsigned n = AR.columns();

  // check for zero sized matrix
  if (m == 0 || n == 0)
    return;

  // determine LAPACK parameters
  INTEGER M = AR.rows();
  INTEGER N = AR.columns();
  INTEGER MINMN = std::min(M, N);
  const unsigned min_mn = (unsigned) MINMN;

  // setup tau vector
  tau.resize(min_mn);

  // call LAPACK
  INTEGER LDA = AR.leading_dim();
  INTEGER INFO;
  geqrf_(&M, &N, AR.data(), &M, tau.data(), &INFO);
  assert(INFO == 0);

  // setup Q
  Q.resize(m,m);
  std::copy(AR.data(), AR.data()+m*min_mn, Q.data());
  orgqr_(&M, &MINMN, &MINMN, Q.data(), &M, tau.data(), &INFO);

  // make R upper triangular 

  // note: R is m x n, so we don't have to resize
  for (unsigned i=0; i< AR.columns(); i++)
  {
    RowIteratord coli = AR.block_row_iterator_begin(i+1,AR.rows(),i,i+1);
    std::fill(coli, coli.end(), 0.0);
  }
}


