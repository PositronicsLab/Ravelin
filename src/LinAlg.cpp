/// Frees all allocated memory
void LINALG::free_memory()
{
  workM().resize(0,0);
  workM2().resize(0,0);
  U().resize(0,0);
  pivwork().resize(0);
  V().resize(0,0);
  S().resize(0);
  workv().resize(0);
  workv2().resize(0);
  iworkv().resize(0);
  compress();
}

/// Compresses all memory
void LINALG::compress()
{
  workM().compress();
  workM2().compress();
  U().compress();
//  pivwork().shrink_to_fit();
  V().compress();
  S().compress();
  workv().compress();
  workv2().compress();
//  iworkv().shrink_to_fit();
}

/// Performs a LDL' factorization of a symmetric, indefinite matrix
/**
 * \param A the matrix A on input; the factorized matrix on output
 */
void LINALG::factor_LDL(MATRIXN& A, vector<int>& IPIV)
{
  #ifndef NEXCEPT
  if (A.rows() != A.columns())
    throw NonsquareMatrixException();
  #endif

  // verify that A is not zero sized
  if (A.rows() == 0 || A.columns() == 0)
    return;

  // verify that A is symmetric
  #ifndef NEXCEPT
  if (!A.is_symmetric())
    std::cerr << "LINALG::factor_LDL() warning: matrix is assymetrical!" << endl;
  #endif

  // resize IPIV
  IPIV.resize(A.rows());

  // setup LAPACK args for factorization
  char UPLO = 'L';
  INTEGER N = A.rows();
  INTEGER INFO;

  // get A's data -- we're going to modify it directly
  REAL* data = A.data();

  // alter matrix to put into packed format
  for (unsigned j=0, k=0, r=0; j< A.rows(); j++, r++)
    for (unsigned i=0, s=0; i< A.columns(); i++, s+= A.rows())
      if (i >= j)
         data[k++] = data[r+s];

  // perform the factorization
  sptrf_(&UPLO, &N, A.data(), &IPIV.front(), &INFO);
  assert(INFO >= 0);
}

//// Computes the psuedo-inverse of a matrix
MATRIXN& LINALG::pseudo_invert(MATRIXN& A, REAL tol)
{
  // get the dimensionality of A
  const unsigned m = A.rows();
  const unsigned n = A.columns();
  const unsigned minmn = std::min(m, n);

  // check for easy out
  if (m == 0)
  {
    A.resize(n,0);
    return A;
  }
  else if (n == 0)
  {
    A.resize(0, m);
    return A;
  }

  // compute the svd
  MATRIXN& Ux = U();
  MATRIXN& Vx = V();
  VECTORN& Sx = S();
  svd(A, Ux, Sx, Vx);
  REAL* Sx_data = Sx.data();

  // determine new tolerance based on first std::singular value if necessary
  if (tol < 0.0)
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
    Sx_data[i] = (std::fabs(Sx_data[i]) > tol) ? (REAL) 1.0/Sx[i] : (REAL) 0.0;

  // scale U' or V, depending on whichever is smaller, m or n
  if (m < n)
    for (unsigned i=0; i< m; i++)
      CBLAS::scal(m, Sx_data[i], Ux.data()+Ux.leading_dim()*i, 1);
  else
    for (unsigned i=0; i< n; i++)
      CBLAS::scal(n, Sx_data[i], Vx.data()+Vx.leading_dim()*i, 1);

  // size the result properly
  A.resize(n, m);

  // do the multiplication (V * U')
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, minmn, (REAL) 1.0, Vx.data(), Vx.leading_dim(), Ux.data(), Ux.leading_dim(), (REAL) 0.0, A.data(), A.leading_dim());

  return A;
}

/// Less robust least squares solver (solves Ax = b)
/**
 * \note this method does not work!
 */
/*
VECTORN& LINALG::solve_LS_fast2(MATRIXN& A, VECTORN& XB)
{
  // verify that A is not zero size
  if (A.rows() == 0 || A.columns() == 0)
    return XB;

  // verify that A and b are appropriate sizes
  if (A.rows() != XB.size())
    throw MissizeException();

  // do QR factorization
  INTEGER M = A.rows();
  INTEGER N = A.columns();
  INTEGER min_mn = std::min(A.rows(), A.columns());

  // setup tau vector
  SAFESTATIC FastThreadable<VECTORN> TAU;
  TAU().resize(min_mn);

  // call LAPACK
  INTEGER INFO;
  geqrf_(&M, &N, A.data(), &M, TAU().data(), &INFO);
  assert(INFO == 0);

  // call ORMQR for work query
  char SIDE = 'L';
  char TRANS = 'T';
  INTEGER CM = XB.size();
  INTEGER CN = 1;
  INTEGER K = TAU().size();
  INTEGER LDA = A.rows();
  INTEGER LDC = CM;
  ormqr_(&SIDE, &TRANS, &CM, &CN, &K, A.data(), &LDA, TAU().data(), XB.data(), &LDC, &INFO);

  // solve triangular system
  char UPLO = 'U';
  INTEGER LDB = XB.size();
  INTEGER NRHS = 1;
  trtrs_(&UPLO, &TRANS, &N, &NRHS, A.data(), &LDA, XB.data(), &LDB, &INFO);

  return XB;
}
*/

/// Computes a givens rotation
void LINALG::givens(REAL a, REAL b, REAL& c, REAL& s)
{
  // setup LAPACK parameters
  REAL UNUSED;

  // call LAPACK (use solving routine that uses LU factorization)
  lartg_(&a, &b, &c, &s, &UNUSED);

  // reverse s
  s = -s;
}

/// Computes the givens matrix given a c and s
MATRIX2 LINALG::givens(REAL c, REAL s)
{
  return MATRIX2(c, s, -s, c);
}

/// Computes a householder vector
void LINALG::householder(REAL alpha, const VECTORN& x, REAL& tau, VECTORN& v)
{
  REAL s = x.norm_sq();
  v = x;
  if (s < (REAL) 0.0)
    tau = (REAL) 0.0;
  else
  {
    REAL t = std::sqrt(alpha*alpha + s);
    REAL v_one = (alpha <= (REAL) 0.0) ? (alpha-t) : -s/(alpha+t);
    tau = 2*v_one*v_one/(s + v_one*v_one);
    v /= v_one;
  }
}

/// Updates a QR factorization by a rank-1 update
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 */
void LINALG::update_QR_rank1(MATRIXN& Q, MATRIXN& R, const VECTORN& u, const VECTORN& v)
{
  // apply Givens rotations
  throw std::runtime_error("This method not currently implemented/working");
}

/// Updates a QR factorization by deleting p columns starting at column idx k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param k the column index to start deleting at
 * \parma p the number of columns to delete
 */
void LINALG::update_QR_delete_cols(MATRIXN& Q, MATRIXN& R, unsigned k, unsigned p)
{
  throw std::runtime_error("This method not currently implemented/working");
/*
  VECTORN workv, workv2;
  MATRIXN workM, workM2, workM3;
  vector<REAL> c, s, tau;
  vector<unsigned> select;
  vector<VECTORN> V;

  const int m = Q.rows();
  const int n = R.columns();
  const int lim = std::min(m-1,n-(int) p);
  assert(k + p <= n);

  // make k one indexed to work with our 1-indexed algorithm
  k++;

  // simplest case
  if (k == n-p+1)
  {
    R.get_sub_mat(0,m,0,k-1,workM);
    R = workM;
    Q.get_sub_mat(0,m,0,m,workM);
    Q = workM;
    return;
  }

  // next simplest case
  if (k > std::min(m-1,(int) (n-p)))
  {
    // get relevant indices of R
    select.clear();
    for (int i=1; i<= n; i++)
      if (i < k || i >= k+p)
        select.push_back(i-1);
    R.select_columns(select.begin(), select.end(), workM);
    R = workM;
    Q.get_sub_mat(0,m,0,R.rows(),workM);
    Q = workM;
    return;
  }

  // setup c, s
  c.resize(std::max(m,n)+1);
  s.resize(std::max(m,n)+1);

  // shift R
  R.get_sub_mat(0,m,k+p-1,n,workM);
  R.set_sub_mat(0,k-1,workM);

  // third simplest case: p = 1 and m >= n
  if (p == 1 && m >= n)
  {
    for (int j=k; j<= n-1; j++)
    {
      // compute Givens rotation
      givens(R(j-1,j-1), R(j,j-1), c[j], s[j]);

      // update R
      R(j-1,j-1) = c[j]*R(j-1,j-1) - s[j]*R(j,j-1);
      R.get_sub_mat(j-1,j+1,j-1,n-1,workM);
      givens(c[j], s[j]).transpose_mult(workM, workM2);
      R.set_sub_mat(j-1,j-1,workM2);
    }

    // compute upper triangular part of R
    R.get_sub_mat(0,m,0,n-1,workM);
    for (int i=1; i<= workM.rows(); i++)
      for (int j=1; j<= std::min(i-1,(int) workM.columns()); j++)
        workM(i-1,j-1) = (REAL) 0.0;
    R = workM;

    // update Q
    for (int j=k; j<= n-1; j++)
    {
      Q.get_sub_mat(0,m,j-1,j+1,workM);
      workM.mult(givens(c[j], s[j]), workM2);
      Q.set_sub_mat(0,j-1, workM2);
    }

    return;
  }

  // householder case (p > 1)
  tau.resize(lim+1);
  V.resize(lim+1);
  for (int j=k; j<= lim; j++)
  {
    // compute householder vector/scalar
    int last = std::min(j+(int) p, m);
    R.get_sub_mat(j, last, j-1, j, workv2);
    householder(R(j-1,j-1), workv2, tau[j], V[j]);

    // update R
    R(j-1,j-1) -= tau[j]*(R(j-1,j-1) + V[j].dot(workv2));
    if (j < n-p)
    {
      R.get_sub_mat(j-1,last,j,n-p, workM);
      workv.resize(V[j].size()+1);
      workv.set_sub_vec(1, V[j]);
      workv[0] = (REAL) 1.0;
      workM.transpose_mult(workv, workv2);
      workv *= tau[j];
      outer_prod(workv, workv2, workM2);
      workM -= workM2;
      R.set_sub_mat(j-1,j,workM);
    }
  }

  // setup upper triangular R
  R.get_sub_mat(0,m,0,n-p,workM);
  R = workM;
  for (int i=1; i<= workM.rows(); i++)
    for (int j=1; j<= std::min(i-1,(int) workM.columns()); j++)
      R(i-1,j-1) = (REAL) 0.0;

  // setup Q
  for (int j=k; j<= lim; j++)
  {
    int last = std::min(j+(int) p,m);
    Q.get_sub_mat(0,m,j-1,last,workM);
    workv.resize(V[j].size()+1);
    workv.set_sub_vec(1, V[j]);
    workv[0] = (REAL) 1.0;
    workv2 = workv;
    workv2 *= tau[j];
    outer_prod(workv, workv2, workM2);
    workM.mult(workM2, workM3);
    workM -= workM3;
    Q.set_sub_mat(0,j-1,workM);
  }
//  Q.get_sub_mat(0,m,0,std::min(lim+(int) p, m), workM);
  Q.get_sub_mat(0,m,0,m, workM);
  Q = workM;
*/
}

/// Updates a QR factorization by inserting one or more columns at column idx k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param U a m x p matrix, destroyed on return
 */
void LINALG::update_QR_insert_cols(MATRIXN& Q, MATRIXN& R, MATRIXN& U, unsigned k)
{
  throw std::runtime_error("This method not currently implemented/working");
/*
  SAFESTATIC VECTORN workv;
  SAFESTATIC MATRIXN workM, workM2, Qu;
  SAFESTATIC MATRIXN c, s;

  const int m = Q.rows();
  const int n = R.columns();
  const int p = U.columns();
  assert(U.rows() == m);

  // make k one indexed to work with our 1-indexed algorithm
  k++;

  // setup c, s
  c.resize(std::max(m,n)+1,std::max(m,n)+1);
  s.resize(std::max(m,n)+1,std::max(m,n)+1);

  // setup U
  Q.transpose_mult(U, workM);
  U = workM;
  if (m > n+1)
  {
    // do a QR factorization
    U.get_sub_mat(n,m,0,p,workM);
    LINALG::factor_QR(workM, Qu);
  }

  if (k <= n)
  {
    // zero out the rest with givens, stop at last column of U or last row that
    // is reached first
    int jstop = std::min((int) p,m-(int) k-2);
    for (int j=1; j<= jstop; j++)
    {
      int istart = std::min(n+j,m);
      int upfirst = std::max(istart-j,1);
      for (int i=istart; i> j; i--)
      {
        givens(U(i-2,j-1),U(i-1,j-1), c(i,j), s(i,j));

        // update U
        U(i-2,j-1) = c(i,j)*U(i-2,j-1) - s(i,j)*U(i-1,j-1);
        if (j < p)
        {
          U.get_sub_mat(i-2,i,j,p,workM);
          givens(c(i,j), s(i,j)).transpose_mult(workM, workM2);
          U.set_sub_mat(i-2,j,workM2);
        }

        // update R
        R.get_sub_mat(i-2,i,upfirst-1,n,workM);
        givens(c(i,j), s(i,j)).transpose_mult(workM, workM2);
        R.set_sub_mat(i-2,upfirst-1,workM2);

        // update one more column next i step
        upfirst--;
      }
    }
  }

  // finish R
  workM.resize(U.rows(), U.columns()+R.columns());
  if (k == 1)
  {
    workM.set_sub_mat(0, 0, U);
    workM.set_sub_mat(0, U.columns(), R);
  }
  else if (k == n+1)
  {
    workM.set_sub_mat(0, 0, R);
    workM.set_sub_mat(0, R.columns(), U);
  }
  else
  {
    R.get_sub_mat(0,m,0,k-1,workM2);
    workM.set_sub_mat(0,0,workM2);
    workM.set_sub_mat(0,workM2.columns(),U);
    R.get_sub_mat(0,m,k-1,n,workM2);
    workM.set_sub_mat(0,U.columns()+k-1,workM2);
  }
  R = workM;

  // finally, make R upper triangular
  for (int i=1; i<= (int) R.rows(); i++)
    for (int j=1; j<= std::min((int) R.columns(),i-1); j++)
      R(i-1,j-1) = (REAL) 0.0;

  // compute Q
  if (m > n+1)
  {
    Q.get_sub_mat(0,m,n,m,workM);
    workM.mult(Qu, workM2);
    Q.set_sub_mat(0,n,workM2);
  }
  if (k <= n)
  {
    int jstop = std::min((int) p,m-(int) k-2);
    for (int j=1; j<= jstop; j++)
    {
      int istart = std::min(n+j,m);
      for (int i=istart; i>= j+1; i--)
      {
        Q.get_sub_mat(0,m,i-2,i,workM);
        workM.mult(givens(c(i,j), s(i,j)), workM2);
        Q.set_sub_mat(0,i-2,workM2);
      }
    }
  }
*/
}

/// Updates a QR factorization by inserting a block of rows, starting at index k
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param U a p x n matrix (destroyed on return)
 * \param k the index to insert at
 */
void LINALG::update_QR_insert_rows(MATRIXN& Q, MATRIXN& R, MATRIXN& U, unsigned k)
{
  throw std::runtime_error("This method not currently implemented/working");
/*
  SAFESTATIC VECTORN workv, workv2, workv3, workv4;
  SAFESTATIC MATRIXN workM, workM2, workM3, Qu, Ru;
  SAFESTATIC vector<REAL> c, s, tau;
  SAFESTATIC vector<VECTORN> V;

  const int m = Q.rows();
  const int n = R.columns();
  const int lim = std::min(m,n);
  const int p = U.rows();

  // make k one indexed
  k++;

  // verify that U is the correct size
  assert(U.columns() == n);

  // simplest case: inserting a single row, m >= n
  if (m >= n && p == 1)
  {
    // resize c and s
   c.resize(n+1);
   s.resize(n+1);

    for (int j=1; j<= n; j++)
    {
      // compute Givens rotation
      givens(R(j-1,j-1), U(0, j-1), c[j], s[j]);

      // update R
      R(j-1,j-1) = c[j]*R(j-1,j-1) - s[j]*U(0,j-1);

      // update jth row of R and u
      R.get_sub_mat(j-1, j, j, n, workv);
      U.get_sub_mat(0, 1, j, n, workv2);
      workv3 = workv;
      workv3 *= c[j];
      workv4 = workv2;
      workv4 *= s[j];
      workv3 -= workv4;
      R.set_sub_mat(j-1,j, workv3, eTranspose);
      workv3 = workv;
      workv3 *= s[j];
      workv4 = workv2;
      workv4 *= c[j];
      workv3 += workv4;
      U.set_sub_mat(0,j, workv3, eTranspose);
    }

    // setup new R
    workM.resize(R.rows() + 1, R.columns());
    workM.set_sub_mat(0, 0, R);
    ITERATOR b = workM.block_iterator(R.rows(), R.rows()+1, 0, R.columns());
    std::fill_n(b, R.columns(), (REAL) 0.0);
    R = workM;

    // compute new Q
    workM.resize(Q.rows()+1, Q.columns()+1);
    workM.set_sub_mat(0,0,Q);
    b = workM.block_iterator(Q.rows(), Q.rows()+1, 0, Q.columns());
    std::fill_n(b, Q.columns(), (REAL) 0.0);
    b = workM.block_iterator(0, Q.rows(), Q.columns(), Q.columns()+1);
    std::fill_n(b, Q.rows(), (REAL) 0.0);
    workM(Q.rows(), Q.columns()) = (REAL) 1.0;
    Q = workM;
    if (k != m+1)
    {
      // permute Q
      workM2.resize(Q.rows(), Q.columns());
      Q.get_sub_mat(0,k-1,0,m+1,workM);
      workM2.set_sub_mat(0,0,workM);
      Q.get_row(m, workv);
      workM2.set_row(workM.rows(), workv);
      Q.get_sub_mat(k-1,m,0,m+1, workM);
      workM2.set_sub_mat(k, 0, workM);
      Q = workM2;
    }
    for (int j=1; j<= n; j++)
    {
      Q.get_sub_mat(0,m+1,j-1,j,workv);
      Q.get_sub_mat(0,m+1,m,m+1,workv2);
      workv3 = workv;
      workv3 *= c[j];
      workv4 = workv2;
      workv4 *= s[j];
      workv3 -= workv4;
      Q.set_sub_mat(0,j-1,workv3);
      workv3 = workv;
      workv3 *= s[j];
      workv4 = workv2;
      workv4 *= c[j];
      workv3 += workv4;
      Q.set_sub_mat(0,m,workv3);
    }

    return;
  }

  // householder case (p > 1 or m < n)
  V.resize(lim+1);
  tau.resize(lim+1);

  for (int j=1; j<= lim; j++)
  {
    // compute householder
    U.get_sub_mat(0,p,j-1,j,workv);
    householder(R(j-1,j-1), workv, tau[j], V[j]);

    // remember old jth row of R
    R.get_sub_mat(j-1,j,j,n,workv);

    // update jth row of R
    R.get_sub_mat(j-1,j,j-1,n,workv2); // workv2 = R(j,j:n)
    workv2 *= ((REAL) 1.0 - tau[j]);
    U.get_sub_mat(0,p,j-1,n,workM);   // workM = U(1:p,j:n)
    workM.transpose_mult(V[j], workv3) *= tau[j];
    workv2 -= workv3;
    R.set_sub_mat(j-1,j-1,workv2,eTranspose);

    // update trailing part if U
    if (j < n)
    {
      U.get_sub_mat(0,p,j,n,workM);                 // get X = U(1:p,j+1:n)
      workv2 = V[j];
      workv2 *= tau[j];             // tau*V
      outer_prod(workv2, V[j], workM2);   // tau*V*V'
      workM2.mult(workM, workM3);                   // tau*V*V'*X
      workM -= workM3;                              // X -= tau*V*V'*X
      outer_prod(workv2, workv, workM2);  // Y = tau*V*Rj'
      workM -= workM2;                              // X -= Y
      U.set_sub_mat(0,j,workM);
    }
  }

  // update R
  workM.resize(R.rows()+p,n);
  workM.set_sub_mat(0, 0, R);
  ITERATOR b = workM.block_iterator(R.rows(), workM.rows(), 0, n);
  std::fill_n(b, p*n, (REAL) 0.0);
  R = workM;
  if (m < n)
  {
    U.get_sub_mat(0, U.rows(), m,n, Ru);
    factor_QR(Ru, Qu);
    R.set_sub_mat(m,m,Ru);
  }

  // update Q
  workM.resize(Q.rows()+p, Q.columns()+p);
  workM.set_sub_mat(0,0,Q);
  b = workM.block_iterator(0, m, Q.columns(), Q.columns()+p);
  std::fill_n(b, m*p, (REAL) 0.0);
  b = workM.block_iterator(Q.rows(), Q.rows()+p, 0, Q.columns());
  std::fill_n(b, p*Q.columns(), (REAL) 0.0);
  b = workM.block_iterator(Q.rows(), Q.rows()+p, Q.columns(), Q.columns()+p);
  std::fill_n(b, p*p, (REAL) 0.0);
  for (unsigned i=0, jj=m, kk=p; i< p; i++, jj++, kk++)
    workM(jj-1,kk-1) = (REAL) 1.0;
  Q = workM;
  if (k != m+1)
  {
    // permute Q
    Q.get_sub_mat(0,k-1,0,m+p,workM);
    Q.get_sub_mat(m,m+p,0,m+p,workM2);
    Q.get_sub_mat(k-1,m,0,m+p,workM3);
    Q.set_sub_mat(0,0,workM);
    Q.set_sub_mat(workM.rows(),0,workM2);
    Q.set_sub_mat(workM.rows()+workM2.rows(),0,workM3);
  }
  for (int j=1; j<= lim; j++)
  {
    // remember jth column of Q
    Q.get_sub_mat(0,m+p,j-1,j,workv);

    // update jth column
    workv2 = workv;
    workv2 *= ((REAL) 1.0 - tau[j]);
    Q.get_sub_mat(0,m+p,m,m+p,workM);
    workM.mult(V[j], workv3) *= tau[j];
    workv2 -= workv3;
    Q.set_sub_mat(0,j-1,workv2);

    // update m+1:p columns of Qhat
    workv2 = V[j];
    workv2 *= tau[j];                 // s = V*tau
    outer_prod(workv, workv2, workM);       // Y = Qhatk * V' * tau
    Q.get_sub_mat(0,m+p,m,m+p, workM2);
    workM -= workM2;
    workM2.mult(V[j], workv);                         // r = X*V
    outer_prod(workv, workv2, workM2);      // Y = X*V*V'*tau
    workM += workM2;
    workM.negate();
    Q.set_sub_mat(0,m,workM);
    if (m < n)
    {
      Q.get_sub_mat(0,m+p,m,m+p,workM);
      workM.mult(Qu, workM2);
      Q.set_sub_mat(0,m,workM2);
    }
  }
*/
}

/// Updates a QR factorization by deleting a block of rows
/**
 * \param Q a m x min(m,n) matrix
 * \param R a min(m,n) x n matrix
 * \param k the index to start deleting at
 * \param p the number of rows to delete
 */
void LINALG::update_QR_delete_rows(MATRIXN& Q, MATRIXN& R, unsigned k, unsigned p)
{
/*
  throw std::runtime_error("This method not currently implemented/working");
  SAFESTATIC vector<REAL> c, s;
  SAFESTATIC VECTORN workv, workv2;
  SAFESTATIC MATRIXN W, workM, workM2, cc, ss;

  const int m = Q.rows();
  const int n = R.columns();

  assert(k+p <= m);

  // make k 1-indexed
  k++;

  // simplest case: p = 1
  if (p == 1)
  {
    // resize c and s
    c.resize(m);
    s.resize(m);

    Q.get_sub_mat(k-1, k, 0, m, workv);
    for (int j=m-1; j>= 1; j--)
    {
      // compute givens rotation and update q
      givens(workv[j-1], workv[j], c[j], s[j]);
      workv[j-1] = c[j]*workv[j-1] - s[j]*workv[j];

      // update R if there is a nonzero row
      if (j <= n)
      {
        R.get_sub_mat(j-1,j+1,j-1,n,workM);
        givens(c[j], s[j]).transpose_mult(workM, workM2);
        R.set_sub_mat(j-1,j-1,workM2);
      }
    }
    R.get_sub_mat(1,m,0,n,workM);
    R = workM;

    // compute Q
    if (k != 1)
    {
      Q.get_sub_mat(0,k-1,0,m,workM);
      Q.set_sub_mat(1,0,workM);
    }

    for (int j=m-1; j>= 2; j--)
    {
      Q.get_sub_mat(1,m,j-1,j+1,workM);
      workM.mult(givens(c[j], s[j]), workM2);
      Q.set_sub_mat(1,j-1,workM2);
    }

    // do not need to update 1st column of Q
    Q.get_sub_mat(1,m,0,1,workv);
    Q.get_sub_mat(1,m,1,2,workv2);
    workv *= s[1];
    workv2 *= c[1];
    workv += workv2;
    Q.set_sub_mat(1,1,workv);
    Q.get_sub_mat(1,m,1,m,workM);
    Q = workM;
    return;
  }

  // "standard" case
  cc.resize(p+1,m);
  ss.resize(p+1,m);
  Q.get_sub_mat(k-1,k+p-1,0,m,W);
  for (int i=1; i<= p; i++)
  {
    for (int j=m-1; j>= i; j--)
    {
      givens(W(i-1,j-1), W(i-1,j), cc(i,j), ss(i,j));
      W(i-1,j-1) = W(i-1,j-1) * cc(i,j) - W(i-1,j) * ss(i,j);
      W.get_sub_mat(i,p,j-1,j+1,workM);
      workM.mult(givens(cc(i,j), ss(i,j)), workM2);
      W.set_sub_mat(i,j-1, workM2);

      // update R if there is a nonzero row
      if (j <= n+i-1)
      {
        R.get_sub_mat(j-1,j+1,j-i,n,workM);
        givens(cc(i,j), ss(i,j)).transpose_mult(workM, workM2);
        R.set_sub_mat(j-1,j-i,workM2);
      }
    }
  }

  // compute the new R
  R.get_sub_mat(p,m,0,n,workM);
  R = workM;

  // compute the new Q
  if (k != 1)
  {
    Q.get_sub_mat(0,k-1,0,m,workM);
    Q.set_sub_mat(p,0,workM);
  }
  for (int i=1; i<= p; i++)
  {
    for (int j=m-1; j>= i+1; j--)
    {
      Q.get_sub_mat(p,m,j-1,j+1,workM);
      workM.mult(givens(cc(i,j), ss(i,j)), workM2);
      Q.set_sub_mat(p,j-1,workM2);
    }

    Q.get_sub_mat(p,m,i-1,i, workv) *= ss(i,i);
    Q.get_sub_mat(p,m,i,i+1, workv2) *= cc(i,i);
    workv += workv2;
    Q.set_sub_mat(p,i,workv);
  }

  // update Q
  Q.get_sub_mat(p,m,p,m,workM);
  Q = workM;
*/
}


