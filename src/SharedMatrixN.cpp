SHAREDMATRIXN::SHAREDMATRIXN()
{
  _rows = _columns = _ld = _start = 0;
}

/// Constructs the matrix using the given arguments
SHAREDMATRIXN::SHAREDMATRIXN(unsigned rows, unsigned cols, unsigned leading_dim, unsigned start, SharedResizable<REAL> data)
{
  _rows = rows;
  _columns = cols;
  _ld = leading_dim;
  _start = start;
  _data = data;
}

/// Copy constructor
SHAREDMATRIXN::SHAREDMATRIXN(const SHAREDMATRIXN& source)
{
  reset_from(source);
}

/// Resets this from another shared matrix
void SHAREDMATRIXN::reset_from(const SHAREDMATRIXN& source)
{
  _rows = source.rows();
  _columns = source.columns();
  _ld = source.rows();
  _start = source._start;
  _data = source._data;
}

/// Accesses the given element
REAL& SHAREDMATRIXN::operator()(unsigned i, unsigned j)
{
  #ifndef NEXCEPT
  if (i >= _rows || j >= _columns)
    throw InvalidIndexException();
  #endif
  return _data[_start + j*_ld + i];
}

/// Accesses the given element
const REAL& SHAREDMATRIXN::operator()(unsigned i, unsigned j) const
{
  #ifndef NEXCEPT
  if (i >= _rows || j >= _columns)
    throw InvalidIndexException();
  #endif
  return _data[_start + j*_ld + i];
}

/// Sets a matrix from a MATRIX3
SHAREDMATRIXN& SHAREDMATRIXN::operator=(const MATRIX3& m)
{
  const unsigned SZ = 3;
  resize(SZ,SZ);
  set_sub_mat(0,0, m);
  return *this;
}

/// Sets a matrix from a MATRIXN
SHAREDMATRIXN& SHAREDMATRIXN::operator=(const MATRIXN& m)
{
  resize(m.rows(),m.columns());
  set_sub_mat(0,0, m);
  return *this;
}

/// Sets a matrix from a MATRIXN
SHAREDMATRIXN& SHAREDMATRIXN::operator=(const SHAREDMATRIXN& m)
{
  resize(m.rows(),m.columns());
  set_sub_mat(0,0, m);
  return *this;
}

/// Sets a matrix from a MATRIXN
SHAREDMATRIXN& SHAREDMATRIXN::operator=(const CONST_SHAREDMATRIXN& m)
{
  resize(m.rows(),m.columns());
  set_sub_mat(0,0, m);
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

/// Resizes this matrix, optionally preserving its existing elements
/**
 * \note this method keeps from reallocating memory unless absolutely
 *       necessary (i.e., if the matrix grows or preserve=true, then memory
 *       will need to be reallocated.
 */
SHAREDMATRIXN& SHAREDMATRIXN::resize(unsigned rows, unsigned columns, bool preserve)
{
  // if the matrix is already the proper size, exit
  if (_rows == rows && _columns == columns)
    return *this;
  #ifndef NEXCEPT
  else
    throw std::runtime_error("Attempt to resize shared matrix!");
  #endif
}

/// Sets the matrix to the zero matrix
SHAREDMATRIXN& SHAREDMATRIXN::set_zero()
{
  COLUMN_ITERATOR i = column_iterator_begin();
  while (i != i.end())
    *i++ = (REAL) 0.0;
  return *this;
}

/// Multiplies this matrix by another in place
SHAREDMATRIXN& SHAREDMATRIXN::operator*=(REAL scalar)
{
  // call BLAS scaling function
  if (_rows > 0)
  {
    REAL* colstart = _data.get()+_start; 
    for (unsigned i=0; i< _columns; i++)
    {
      CBLAS::scal(_rows, scalar, colstart, 1);
      colstart += leading_dim();
    }
  }

  return *this;
}

/// Divides this matrix by a scalar in place
SHAREDMATRIXN& SHAREDMATRIXN::operator/=(REAL scalar)
{
  // call BLAS scaling function
  if (_rows > 0)
  {
    REAL* colstart = _data.get()+_start; 
    for (unsigned i=0; i< _columns; i++)
    {
      CBLAS::scal(_rows, (REAL) 1.0/scalar, colstart, 1);
      colstart += leading_dim();
    }
  }

  return *this;
}

/// Negates this matrix in place
SHAREDMATRIXN& SHAREDMATRIXN::negate()
{
  return operator*=((REAL) -1.0);
}

/// Checks whether the given matrix is symmetric to the specified tolerance
bool SHAREDMATRIXN::is_symmetric(REAL tolerance) const
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
SHAREDMATRIXN& SHAREDMATRIXN::zero_upper_triangle()
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
SHAREDMATRIXN& SHAREDMATRIXN::zero_lower_triangle()
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

/// Sets this matrix to the identity matrix
SHAREDMATRIXN& SHAREDMATRIXN::set_identity()
{
  set_zero();
  for (unsigned i=0, j=0; i< _rows; i++, j+= leading_dim()+1)
    _data[j] = (REAL) 1.0;

  return *this;
}

CONST_SHAREDMATRIXN::CONST_SHAREDMATRIXN()
{
  _rows = _columns = _ld = _start = 0;
}

/// Constructs the matrix using the given arguments
CONST_SHAREDMATRIXN::CONST_SHAREDMATRIXN(unsigned rows, unsigned cols, unsigned leading_dim, unsigned start, SharedResizable<REAL> data)
{
  _rows = rows;
  _columns = cols;
  _ld = leading_dim;
  _start = start;
  _data = data;
}

/// Copy constructor
CONST_SHAREDMATRIXN::CONST_SHAREDMATRIXN(const SHAREDMATRIXN& source)
{
  reset_from(source);
}

/// Gets this object as a standard shared matrix 
/**
 * \note const-ness is not enforced by my compiler! 
 */
const SHAREDMATRIXN CONST_SHAREDMATRIXN::get() const
{
  SHAREDMATRIXN m;
  m._rows = _rows;
  m._columns = _columns;
  m._ld = _ld;
  m._start = _start;
  m._data = _data;
  return m;
}

/// Resets this from another shared matrix
void CONST_SHAREDMATRIXN::reset_from(const SHAREDMATRIXN& source)
{
  _rows = source.rows();
  _columns = source.columns();
  _ld = source.leading_dim();
  _start = source._start;
  _data = source._data;
}

/// Resets this from another shared matrix
void CONST_SHAREDMATRIXN::reset_from(const CONST_SHAREDMATRIXN& source)
{
  _rows = source.rows();
  _columns = source.columns();
  _ld = source.leading_dim();
  _start = source._start;
  _data = source._data;
}

/// Copy constructor
CONST_SHAREDMATRIXN::CONST_SHAREDMATRIXN(const CONST_SHAREDMATRIXN& source)
{
  reset_from(source);
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const CONST_SHAREDMATRIXN& m)
{
  const unsigned OUTPUT_PRECISION = 8;

  if (m.rows() == 0 || m.columns() == 0)
  {
    out << "(empty)" << std::endl;
    return out;
  }

  for (unsigned i=0; i< m.rows(); i++)
  {
    for (unsigned j=0; j< m.columns()-1; j++)
      out << std::setprecision(OUTPUT_PRECISION) << m(i,j) << " ";
    out << std::setprecision(OUTPUT_PRECISION) << m(i,m.columns()-1) << std::endl;
  }

  return out;
}

/// Accesses the given element
const REAL& CONST_SHAREDMATRIXN::operator()(unsigned i, unsigned j) const
{
  #ifndef NEXCEPT
  if (i >= _rows || j >= _columns)
    throw InvalidIndexException();
  #endif
  return _data[_start + j*_ld + i];
}

// use common routines
#define XMATRIXN SHAREDMATRIXN
#include "XMatrixN.cpp"
#undef XMATRIXN

