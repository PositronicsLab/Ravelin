#ifndef REENTRANT 
FastThreadable<MATRIXN> MATRIXN::_n;
FastThreadable<VECTORN> MATRIXN::_workv;
#endif

/// Default constructor - constructs an empty matrix
MATRIXN::MATRIXN()
{
  _rows = _columns = 0;
}

/// Constructs a rows x columns dimensional (unitialized) matrix 
MATRIXN::MATRIXN(unsigned rows, unsigned columns)
{
  _rows = rows;
  _columns = columns;
  if (rows > 0 && columns > 0)
    _data.resize(_rows*_columns);
}

/// Constructs a matrix from a vector
/**
 * \param v the vector
 * \param determines whether the vector will be transposed
 */
MATRIXN::MATRIXN(const VECTORN& v, Transposition trans)
{
  _rows = _columns = 0;
  set(v, trans);
}

/// Constructs a matrix from a vector
/**
 * \param v the vector
 * \param determines whether the vector will be transposed
 */
MATRIXN::MATRIXN(const CONST_SHAREDVECTORN& v, Transposition trans)
{
  _rows = _columns = 0;
  set(v, trans);
}

/// Constructs a matrix from a vector
/**
 * \param v the vector
 * \param determines whether the vector will be transposed
 */
MATRIXN::MATRIXN(const SHAREDVECTORN& v, Transposition trans)
{
  _rows = _columns = 0;
  set(v, trans);
}

/// Constructs a matrix from an array
/**
 * \param rows the number of rows of the matrix
 * \param columns the number of columns of the matrix
 * \param array an array of rows*columns REAL values in column-major format
 */
MATRIXN::MATRIXN(unsigned rows, unsigned columns, const REAL* array)
{
  _rows = _columns = 0;
  resize(rows, columns, false);
  CBLAS::copy(rows*columns, array, 1, data(), 1);
}

/// Constructs a matrix from a MATRIX3
MATRIXN::MATRIXN(const MATRIX3& m)
{
  _rows = _columns = 0;
  operator=(m);
}

/// Constructs a matrix from a pose
/* 
MATRIXN::MATRIXN(const POSE& m)
{
  _rows = _columns = _capacity = 0;
  set(m);
}
*/

/// Copy constructor
MATRIXN::MATRIXN(const MATRIXN& source)
{
  _rows = _columns = 0;
  MATRIXN::operator=(source);
}

/// Copy constructor
MATRIXN::MATRIXN(const SHAREDMATRIXN& source)
{
  _rows = _columns = 0;
  MATRIXN::operator=(source);
}

/// Copy constructor
MATRIXN::MATRIXN(const CONST_SHAREDMATRIXN& source)
{
  _rows = _columns = 0;
  MATRIXN::operator=(source);
}

/// Sets a matrix from a MATRIX2
MATRIXN& MATRIXN::operator=(const MATRIX2& m)
{
  const unsigned SZ = 2;
  resize(SZ,SZ);
  CBLAS::copy(SZ*SZ, m.data(), 1, _data.get(), 1);
  return *this;
}

/// Sets a matrix from a MATRIX3
MATRIXN& MATRIXN::operator=(const MATRIX3& m)
{
  const unsigned SZ = 3;
  resize(SZ,SZ);
  CBLAS::copy(SZ*SZ, m.data(), 1, _data.get(), 1);
  return *this;
}

/// Sets this matrix from a vector
MATRIXN& MATRIXN::set(const VECTORN& v, Transposition trans)
{
  if (trans == eNoTranspose)
    resize(v.rows(), 1, false);
  else
    resize(1, v.columns(), false);
  std::copy(v.data(), v.data()+v.size(), _data.get());
  return *this;
}

/// Sets this matrix from a vector
MATRIXN& MATRIXN::set(const SHAREDVECTORN& v, Transposition trans)
{
  if (trans == eNoTranspose)
    resize(v.rows(), 1, false);
  else
    resize(1, v.columns(), false);
  std::copy(v.column_iterator_begin(), v.column_iterator_end(), column_iterator_begin());
  return *this;
}

/// Sets this matrix from a vector
MATRIXN& MATRIXN::set(const CONST_SHAREDVECTORN& v, Transposition trans)
{
  if (trans == eNoTranspose)
    resize(v.rows(), 1, false);
  else
    resize(1, v.columns(), false);
  std::copy(v.column_iterator_begin(), v.column_iterator_end(), column_iterator_begin());
  return *this;
}

/*
/// Sets a matrix from a pose 
MATRIXN& MATRIXN::set(const POSE& p)
{
  // resize the matrix
  const unsigned SZ = 4;
  resize(SZ,SZ);

  // copy
  const unsigned MCOL1 = 0, MCOL2 = 3, MCOL3 = 6, MCOL4 = 9;
  const unsigned COL1 = 0, COL2 = 4, COL3 = 8, COL4 = 12;
  const unsigned WX = 3, WY = 7, WZ = 11, WW = 15;
  const REAL* mdata = m.begin();
  for (unsigned i=0; i< 3; i++)
  {
    _data[COL1 + i] = mdata[MCOL1 + i];
    _data[COL2 + i] = mdata[MCOL2 + i];
    _data[COL3 + i] = mdata[MCOL3 + i];
    _data[COL4 + i] = mdata[MCOL4 + i];
  }

  // set last row
  _data[WX] = _data[WY] = _data[WZ] = (REAL) 0.0;
  _data[WW] = (REAL) 1.0;
  return *this;
}
*/

/// Constructs a MATRIXN using a variable number of double values
/**
 * \note the values are given row by row
 * \note There is no means in C++ to check the types of a list of variable
 * arguments.  If the variable arguments are not of type double, then
 * unexpected values will result. Constructing a matrix using the
 * statement MATRIXN::construct_variable(2, 2, 1.0, 1.0, 1.0, 0) is incorrect 
 * because the programmer has assumed that the integer 0 will be converted to 
 * a double type; this is not the case.
 */
MATRIXN MATRIXN::construct_variable(unsigned rows, unsigned cols, ...)
{
  MATRIXN m(rows, cols);
  std::va_list lst;
  va_start(lst, cols);
  for (unsigned i=0; i< rows; i++)
    for (unsigned j=0; j< cols; j++)
      m(i,j) = (REAL) va_arg(lst, double);
  va_end(lst);
  return m;
}

/// Removes a row from the matrix
/**
 * \note reallocates memory 
 */
MATRIXN& MATRIXN::remove_row(unsigned i)
{
  MATRIXN tmp(_rows-1, _columns);
  CONST_SHAREDMATRIXN top_block_this = block(0, i, 0, columns());
  CONST_SHAREDMATRIXN bottom_block_this = block(i+1, rows(), 0, columns());
  SHAREDMATRIXN top_block_tmp = tmp.block(0, i, 0, columns());
  SHAREDMATRIXN bottom_block_tmp = tmp.block(i, tmp.rows(), 0, columns());
  top_block_tmp = top_block_this;
  bottom_block_tmp = bottom_block_this;
 
  // now set this from the new matrix
  _data = tmp._data;
  _rows = tmp._rows;
  _columns = tmp._columns;
  return *this;
}

/// Removes a column from the matrix
/**
 * \note downsizes the matrix -- does not reallocate memory
 */
MATRIXN& MATRIXN::remove_column(unsigned i)
{
  #ifdef REENTRANT
  VECTORN workv;
  #else
  VECTORN& workv = _workv();
  #endif

  for (unsigned j=i+1; j< _columns; j++)
  {
    get_column(j, workv);
    set_column(j-1, workv);
  }

  // downsize the matrix
  _columns--;
  return *this;
}

/// Returns the zero matrix
/**
 * \param rows the number of rows
 * \param columns the number of columns of the matrix
 * \return a rows x columns zero matrix
 */
MATRIXN MATRIXN::zero(unsigned int rows, unsigned int columns)
{
  // create the new matrix
  MATRIXN m(rows, columns);

  // get the array for the matrix, and set all values to zero
  REAL* x = m.data();
  for (unsigned i=0; i< rows*columns; i++)
    x[i] = (REAL) 0.0;

  return m; 
}

/// Resizes this matrix, optionally preserving its existing elements
/**
 * \note this method keeps from reallocating memory unless absolutely
 *       necessary (i.e., if the matrix grows or preserve=true, then memory
 *       will need to be reallocated.
 */
MATRIXN& MATRIXN::resize(unsigned rows, unsigned columns, bool preserve)
{
  // if the matrix is already the proper size, exit
  if (_rows == rows && _columns == columns)
    return *this;

  // if we can downsize, do that..
  if (rows*columns <= _data.capacity() && (_rows == rows || !preserve))
  {
    _rows = rows;
    _columns = columns;
    _data.resize(_rows*_columns, false);
    return *this;
  }

  // create a new array
  _rows = rows;
  _columns = columns;
  _data.resize(_rows*_columns, true);
  return *this;
}

/// Sets the matrix to the zero matrix
MATRIXN& MATRIXN::set_zero()
{
  std::fill(_data.get(), _data.get()+_rows*_columns, (REAL) 0.0);
  return *this;
}

/// Sets this matrix to its transpose
MATRIXN& MATRIXN::transpose()
{
  #ifdef REENTRANT
  MATRIXN n;
  #else
  MATRIXN& n = _n();
  #endif

  // do fastest transpose first (if possible)
  if (_rows == 1 || _columns == 1)
  {
    std::swap(_rows, _columns);
    return *this;
  }  

  // do second fastest transpose, if possible
  if (_rows == _columns)
  {
    for (unsigned i=0; i< _rows; i++)
      for (unsigned j=i+1; j< _rows; j++)
        std::swap(_data[i*_columns+j], _data[j*_rows+i]);

    return *this;
  }

  // do slowest transpose operation
  n.resize(_columns, _rows);
  REAL* ndata = n.data();
  for (unsigned i=0; i< _rows; i++)
    for (unsigned j=0; j< _columns; j++)
      ndata[i*_columns + j] = _data[j*_rows + i];
  operator=(n);

  return *this;
}

/// Multiplies this matrix by another in place
MATRIXN& MATRIXN::operator*=(REAL scalar)
{
  // call BLAS scaling function
  if (_rows > 0 && _columns > 0)
    CBLAS::scal(_rows*_columns, scalar, _data.get(), 1);
  return *this;
}

/// Divides this matrix by a scalar in place
MATRIXN& MATRIXN::operator/=(REAL scalar)
{
  // call BLAS scaling function
  if (_rows > 0 && _columns > 0)
    CBLAS::scal(_rows*_columns, (REAL) 1.0/scalar, _data.get(), 1);
  return *this;
}

/// Negates this matrix in place
MATRIXN& MATRIXN::negate()
{
  unsigned n = _rows * _columns;
  for (unsigned i=0; i< n; i++)
    _data[i] = -_data[i];
  return *this;
}

/// Checks whether the given matrix is symmetric to the specified tolerance
bool MATRIXN::is_symmetric(REAL tolerance) const
{
  if (_rows != _columns)
    return false;

  // make sure that tolerance is positive 
  if (tolerance <= (REAL) 0.0)
    tolerance = std::numeric_limits<REAL>::epsilon() * norm_inf() * _rows;

  // check symmetry
  const unsigned LD = leading_dim();

  // loop over columns
  for (unsigned i=0, ii=0; i< _rows; i++, ii+= LD)
    for (unsigned j=0, jj=0; j< i; j++, jj+= LD)
      if (std::fabs(_data[jj+i] - _data[ii+j]) > tolerance)
        return false;

  return true;
}

/// Zeros the upper triangle of the matrix
MATRIXN& MATRIXN::zero_upper_triangle()
{
  // check for easy exit
  if (_rows == 0 || _columns == 0)
    return *this;

  // zero the upper triangle
  for (unsigned i=0, s=leading_dim(); i< _rows; i++, s+= leading_dim()+1)
    for (unsigned j=i+1, r=0; j< _columns; j++, r+= leading_dim())
      _data[r+s] = (REAL) 0.0;

  return *this;
}

/// Zeros the lower triangle of the matrix
MATRIXN& MATRIXN::zero_lower_triangle()
{
  // check for easy exit
  if (_rows == 0 || _columns == 0)
    return *this;

  // zero the lower triangle
  for (unsigned i=1, s=1; i< _rows; i++, s++)
    for (unsigned j=0, r=0; j< std::min(i, _columns); j++, r+= leading_dim())
      _data[r+s] = (REAL) 0.0;

  return *this;
}

/// Returns an identity matrix
MATRIXN MATRIXN::identity(unsigned n)
{
  MATRIXN m;
  m.set_identity(n);
  return m;
}

/// Sets this matrix to the identity matrix
MATRIXN& MATRIXN::set_identity(unsigned i)
{
  resize(i,i);
  return set_identity();
}

/// Gets the desired entry
const REAL& MATRIXN::operator()(unsigned i, unsigned j) const
{
  #ifndef NEXCEPT
  if (i >= _rows || j >= _columns)
    throw InvalidIndexException();
  #endif
  return _data[j*_rows+i];
}

/// Gets the desired entry
REAL& MATRIXN::operator()(unsigned i, unsigned j) 
{
  #ifndef NEXCEPT
  if (i >= _rows || j >= _columns)
    throw InvalidIndexException();
  #endif
  return _data[j*_rows+i];
}

/// Sets this matrix to the identity matrix
MATRIXN& MATRIXN::set_identity()
{
  #ifndef NEXCEPT
  if (_rows != _columns)
    throw MissizeException();
  #endif

  // set matrix to identity
  set_zero();

  // set diagonal entries
  for (unsigned i=0, j=0; i< _rows; i++, j+= leading_dim()+1)
    _data[j] = (REAL) 1.0;

  return *this;
}

/// Subtracts m from this 
MATRIXN& MATRIXN::operator-=(const MATRIXN& m)
{
  #ifndef NEXCEPT
  if (_rows != m._rows || _columns != m._columns)
    throw MissizeException(); 
  #endif

  const unsigned N = _rows*_columns;
  if (N > 0)
    std::transform(_data.get(), _data.get()+N, m._data.get(), _data.get(), std::minus<REAL>());
  return *this;
}

/// Adds m to this 
MATRIXN& MATRIXN::operator+=(const MATRIXN& m)
{
  #ifndef NEXCEPT
  if (_rows != m._rows || _columns != m._columns)
    throw MissizeException(); 
  #endif

  const unsigned N = _rows*_columns;
  if (N > 0)
    std::transform(_data.get(), _data.get()+N, m._data.get(), _data.get(), std::plus<REAL>());
  return *this;
}

/// Sets this to m 
MATRIXN& MATRIXN::operator=(const MATRIXN& m)
{
  // resize this (don't preserve)
  resize(m.rows(), m.columns(), false);

  const unsigned N = _rows*_columns;
  if (N > 0)
    std::copy(m._data.get(), m._data.get()+N, _data.get());
  return *this;
}

/// Sets this to m 
MATRIXN& MATRIXN::operator=(const SHAREDMATRIXN& m)
{
  // resize this (don't preserve)
  resize(m.rows(), m.columns(), false);

  if (_rows > 0 && _columns > 0)
    std::copy(m.column_iterator_begin(), m.column_iterator_end(), column_iterator_begin());
  return *this;
}

/// Sets this to m 
MATRIXN& MATRIXN::operator=(const CONST_SHAREDMATRIXN& m)
{
  // resize this (don't preserve)
  resize(m.rows(), m.columns(), false);

  if (_rows > 0 && _columns > 0)
    std::copy(m.column_iterator_begin(), m.column_iterator_end(), column_iterator_begin());
  return *this;
}

// use common routines
#define XMATRIXN MATRIXN
#include "XMatrixN.cpp"
#undef XMATRIXN

