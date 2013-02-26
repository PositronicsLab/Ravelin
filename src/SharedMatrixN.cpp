SHAREDMATRIXN::SHAREDMATRIXN()
{
  _rows = _columns = _ld = _start = 0;
}

/// Copy constructor
SHAREDMATRIXN::SHAREDMATRIXN(const SHAREDMATRIXN& source)
{
  _rows = source.rows();
  _columns = source.columns();
  _ld = source.rows();
  _start = 0;
  _data = source._data;
}

/// Sets a matrix from a MATRIX3
SHAREDMATRIXN& SHAREDMATRIXN::operator=(const MATRIX3& m)
{
  const unsigned SZ = 3;
  resize(SZ,SZ);
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
  else
    throw std::runtime_error("Attempt to resize shared matrix!");
}

/// Sets the matrix to the zero matrix
SHAREDMATRIXN& SHAREDMATRIXN::set_zero()
{
  ITERATOR i = begin();
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

