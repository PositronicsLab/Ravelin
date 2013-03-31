#include <iostream>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SparseMatrixNd.h>

using namespace Ravelin;
using std::endl;
using std::cout;

MatrixNd random_sparse(unsigned SZ1, unsigned SZ2)
{
  MatrixNd m(SZ1,SZ2);
  for (unsigned i=0; i< SZ1; i++)
    for (unsigned j=0; j< SZ2; j++)
    {
      if (rand() % 2 == 0)
        m(i,j) = (double) rand() / RAND_MAX;
      else
        m(i,j) = 0.0;
    }

  return m;
}

void test_mult(const SparseMatrixNd& s1, const SparseMatrixNd& s2, const MatrixNd& d)
{
  MatrixNd rm1, rm2;
  VectorNd rv1, rv2;

  // get the first column from d
  VectorNd col1 = d.column(0);

  // test matrix/vector multiplication
  d.mult(col1, rv2);
  s1.mult(col1, rv1);
  cout << "testing matrix/vector (CSR format)  error: " << (rv1 -= rv2).norm() << endl; 
  s2.mult(col1, rv1);
  cout << "testing matrix/vector (CSC format)  error: " << (rv1 -= rv2).norm() << endl; 

  // test matrix/matrix multiplication
  d.mult(d, rm2);
  s1.mult(d, rm1);
  cout << "testing matrix/matrix (CSR format)  error: " << (rm1 -= rm2).norm_inf() << endl; 
  s2.mult(d, rm1);
  cout << "testing matrix/matrix (CSC format)  error: " << (rm1 -= rm2).norm_inf() << endl; 

  // test transpose matrix/vector multiplication
  d.transpose_mult(col1, rv2);
  s1.transpose_mult(col1, rv1);
  cout << "testing transpose matrix/vector (CSR format)  error: " << (rv1 -= rv2).norm() << endl; 
  s2.transpose_mult(col1, rv1);
  cout << "testing transpose matrix/vector (CSC format)  error: " << (rv1 -= rv2).norm() << endl; 

  // test transpose matrix/matrix multiplication
  d.transpose_mult(d, rm2);
  s1.transpose_mult(d, rm1);
  cout << "testing transpose matrix/matrix (CSR format)  error: " << (rm1 -= rm2).norm_inf() << endl; 
  s2.transpose_mult(d, rm1);
  cout << "testing transpose matrix/matrix (CSC format)  error: " << (rm1 -= rm2).norm_inf() << endl; 

  // test matrix/transpose matrix multiplication
  d.mult_transpose(d, rm2);
  s1.mult_transpose(d, rm1);
  cout << "testing matrix/transpose matrix (CSR format)  error: " << (rm1 -= rm2).norm_inf() << endl; 
  s2.mult_transpose(d, rm1);
  cout << "testing matrix/transpose matrix (CSC format)  error: " << (rm1 -= rm2).norm_inf() << endl; 

  // test transpose matrix/transpose matrix multiplication
  d.transpose_mult_transpose(d, rm2);
  s1.transpose_mult_transpose(d, rm1);
  cout << "testing transpose matrix/transpose matrix (CSR format)  error: " << (rm1 -= rm2).norm_inf() << endl; 
  s2.transpose_mult_transpose(d, rm1);
  cout << "testing transpose matrix/transpose matrix (CSC format)  error: " << (rm1 -= rm2).norm_inf() << endl; 
}

void test_plus(SparseMatrixNd& s1, const SparseMatrixNd& s2, const MatrixNd& d)
{
  // test subtraction
  s1 -= s2;
  cout << "subtraction error: " << s1.norm_inf() << std::endl;

  // test negation and addition
  s1 = SparseMatrixNd(SparseMatrixNd::eCSR, d);
  s1.negate();
  s1 += s2;
  cout << "negation and addition error: " << s1.norm_inf() << std::endl;
}

void test_to_dense(const SparseMatrixNd& s1, const SparseMatrixNd& s2, const MatrixNd& d)
{
  MatrixNd d1, d2;
  s1.to_dense(d1);
  s2.to_dense(d2);

  d1 -= d;
  d2 -= d;
  cout << "testing sparse to dense (CSR): " << d1.norm_inf() << endl;
  cout << "testing sparse to dense (CSC): " << d2.norm_inf() << endl;
}

int main()
{
  // setup a random sparse matrix in dense form
  const unsigned SZ = 4;
  MatrixNd dense = random_sparse(SZ, SZ);

  // setup the sparse matrix
  SparseMatrixNd s1(SparseMatrixNd::eCSR, dense);
  SparseMatrixNd s2(SparseMatrixNd::eCSC, dense);

  // test sparse to dense
  test_to_dense(s1, s2, dense);

  // test multiplication arithmetic 
  test_mult(s1, s2, dense);

  // test addition/subtraction arithmetic 
  test_plus(s1, s2, dense);
}

