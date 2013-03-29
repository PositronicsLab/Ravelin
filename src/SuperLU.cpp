/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <vector>
#include <stdexcept>
#include <superlu/slu_ddefs.h>

struct SuperMatrices
{
  SuperMatrix A, B, X;

  SuperMatrices()
  {
    // allocate memory for A
    A.Store = (void*) SUPERLU_MALLOC(sizeof(NRformat));
    B.Store = (void*) SUPERLU_MALLOC(sizeof(DNformat));
    X.Store = (void*) SUPERLU_MALLOC(sizeof(DNformat));
    if (!A.Store || !B.Store || !X.Store)
      throw std::runtime_error("Unable to allocate memory");
  }

  ~SuperMatrices()
  {
    // free memory
    SUPERLU_FREE(A.Store);
    SUPERLU_FREE(B.Store);
    SUPERLU_FREE(X.Store);
  }
};

// creates a CRS matrix
static void create_CompRow_Matrix(SuperMatrix *A, int m, int n, int nnz,
                       double *nzval, int *colind, int *rowptr,
                       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NRformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    Astore = (NRformat*) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}

// creates a CRS matrix
static void create_CompRow_Matrix(SuperMatrix *A, int m, int n, int nnz,
                       float* nzval, int *colind, int *rowptr,
                       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NRformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    Astore = (NRformat*) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}

// creates a dense matrix
static void create_Dense_Matrix(SuperMatrix *X, int m, int n, double *x, int ldx,
                    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;

    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (double *) x;
}

// creates a dense matrix
static void create_Dense_Matrix(SuperMatrix *X, int m, int n, float *x, int ldx,
                    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;

    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (float *) x;
}

/// Does a LU factorization of a sparse matrix
int solve_superlu(bool notrans, int m, int n, int nrhs, int nnz, int* col_indices, int* row_ptr, double* A_nz, double* x, double* b)
{
  #ifdef USE_SUPERLU
  static SuperMatrices SM;
  SuperMatrix L, U;

  // temporaries
  static std::vector<int> perm_r, perm_c, etree;
  static std::vector<double> R, C, ferr, berr;

  // setup super LU options
  superlu_options_t o;
  set_default_options(&o);
  o.Fact = DOFACT;
  o.Equil = YES;
  o.IterRefine = EXTRA;
  o.SymmetricMode = NO;
  o.ConditionNumber = NO;
  o.ColPerm = NATURAL;
  o.Trans = (notrans) ? NOTRANS : TRANS;

  // setup B and X
  create_Dense_Matrix(&SM.B, m, nrhs, b, m, SLU_DN, SLU_D, SLU_GE);
  create_Dense_Matrix(&SM.X, m, nrhs, x, m, SLU_DN, SLU_D, SLU_GE);

  // setup A 
  create_CompRow_Matrix(&SM.A, m, n, nnz, A_nz, col_indices, row_ptr, 
                         SLU_NR, SLU_D, SLU_GE);  

  // resize temporaries
  etree.resize(n);
  perm_r.resize(m);
  perm_c.resize(n);
  R.resize(m);
  C.resize(n);
  ferr.resize(nrhs);
  berr.resize(nrhs);

  // do the solution
  char equed[1];
  double rpg, rcond;
  mem_usage_t mem_usage;
  SuperLUStat_t stat;
  StatInit(&stat);
  int lwork = 0;
  int info;
/*
  dgssv(&o, &SM.A, &perm_c[0], &perm_r[0], &L, &U, &SM.B, &stat, &info);
*/
  dgssvx(&o, &SM.A,     // options structure and A matrix
         &perm_c[0], &perm_r[0], &etree[0], // some temporaries 
         &equed[0], &R[0], &C[0],           // more temporaries
         &L, &U,   // L and U outputs
         NULL, lwork,      // indicate to use malloc() for data storage
         &SM.B, &SM.X,   // linear system system inputs 
         &rpg, &rcond, &ferr[0], &berr[0], &mem_usage, &stat, &info); // outputs
  StatFree(&stat);
  // check info
  assert(info >= 0);

  // free memory
  if (lwork >= 0)
  {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
  }

  return info;
  #else
  return -1;
  #endif 
}

/// Does a LU factorization of a sparse matrix
int solve_superlu(bool notrans, int m, int n, int nrhs, int nnz, int* col_indices, int* row_ptr, float* A_nz, float* x, float* b)
{
  #ifdef USE_SUPERLU
  static SuperMatrices SM;
  SuperMatrix L, U;

  // temporaries
  static std::vector<int> perm_r, perm_c, etree;
  static std::vector<double> R, C, ferr, berr;

  // setup super LU options
  superlu_options_t o;
  set_default_options(&o);
  o.Fact = DOFACT;
  o.Equil = YES;
  o.IterRefine = EXTRA;
  o.SymmetricMode = NO;
  o.ConditionNumber = NO;
  o.ColPerm = NATURAL;
  o.Trans = (notrans) ? NOTRANS : TRANS;

  // setup B and X
  create_Dense_Matrix(&SM.B, m, nrhs, b, m, SLU_DN, SLU_S, SLU_GE);
  create_Dense_Matrix(&SM.X, m, nrhs, x, m, SLU_DN, SLU_S, SLU_GE);

  // setup A 
  create_CompRow_Matrix(&SM.A, m, n, nnz, A_nz, col_indices, row_ptr, 
                         SLU_NR, SLU_S, SLU_GE);  

  // resize temporaries
  etree.resize(n);
  perm_r.resize(m);
  perm_c.resize(n);
  R.resize(m);
  C.resize(n);
  ferr.resize(nrhs);
  berr.resize(nrhs);

  // do the solution
  char equed[1];
  double rpg, rcond;
  mem_usage_t mem_usage;
  SuperLUStat_t stat;
  StatInit(&stat);
  int lwork = 0;
  int info;
/*
  sgssv(&o, &SM.A, &perm_c[0], &perm_r[0], &L, &U, &SM.B, &stat, &info);
*/
  dgssvx(&o, &SM.A,     // options structure and A matrix
         &perm_c[0], &perm_r[0], &etree[0], // some temporaries 
         &equed[0], &R[0], &C[0],           // more temporaries
         &L, &U,   // L and U outputs
         NULL, lwork,      // indicate to use malloc() for data storage
         &SM.B, &SM.X,   // linear system system inputs 
         &rpg, &rcond, &ferr[0], &berr[0], &mem_usage, &stat, &info); // outputs
  StatFree(&stat);
  // check info
  assert(info >= 0);

  // free memory
  if (lwork >= 0)
  {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
  }

  return info;
  #else
  return -1;
  #endif 
}

