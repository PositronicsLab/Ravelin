#include "ConstMatrixN.inl"

//////////////////////////////////////////////////////////////////////////////
// This file consists of general routines for dynamically allocated matrices
//////////////////////////////////////////////////////////////////////////////

/// Gets a shared vector for a row
SHAREDVECTORN row(unsigned i)
{
  // get the starting offset
  const unsigned OFFSET = data() - _data.get();

  SHAREDVECTORN v;
  v._data = _data;
  v._start = OFFSET + i;
  v._inc = leading_dim();
  v._len = _columns;
  return v;
}

/// Gets a shared vector for a column 
SHAREDVECTORN column(unsigned i)
{
  // get the starting offset
  const unsigned OFFSET = data() - _data.get();

  SHAREDVECTORN v;
  v._data = _data;
  v._start = OFFSET + leading_dim() * i;
  v._inc = 1;
  v._len = _rows;
  return v;
}

/// Gets a block as a shared matrix
SHAREDMATRIXN block(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > _rows || col_end < col_start || col_end > _columns)
    throw InvalidIndexException();
  #endif

  // determine the offset
  const unsigned OFFSET = data() - _data.get();

  SHAREDMATRIXN m;
  m._data = _data;
  m._rows = row_end - row_start;
  m._columns = col_end - col_start;
  m._ld = leading_dim();
  m._start = OFFSET + m._ld * col_start + row_start;
  return m;  
}

/// Sets a matrix from a vector
template <class V>
XMATRIXN& set(const V& v, Transposition trans)
{
  #ifndef NEXCEPT
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  // resize the matrix
  if (trans == eNoTranspose)
    resize(v.size(), 1);
  else
    resize(1, v.size());

  if (_rows > 0 && _columns > 0)
    CBLAS::copy(v.size(), v.data(), 1, _data.get(), 1);

  return *this;
}

/// Determines the transpose of a matrix and stores the result in a given matrix
template <class M>
static M& transpose(const XMATRIXN& m, M& result)
{
  #ifndef NEXCEPT
  if (sizeof(m.data()) != sizeof(result.data()))
    throw DataMismatchException();
  #endif

  // resize the result
  result.resize(m.columns(), m.rows());

  const REAL* mdata = m.data();
  REAL* rdata = result.data();
  for  (unsigned i=0; i< m.rows(); i++)
    for (unsigned j=0; j< m.columns(); j++)
      rdata[i*result.leading_dim()+j] = mdata[j*m.leading_dim()+i];
  
  return result;
}

/// Adds m to *this in place
template <class M>
XMATRIXN& operator+=(const M& m)
{
  #ifndef NEXCEPT
  if (_rows != m.rows() || _columns != m.columns())
    throw MissizeException(); 
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();
  #endif

  if (_rows > 0 && _columns > 0) 
    std::transform(column_iterator_begin(), column_iterator_end(), m.column_iterator_begin(), column_iterator_begin(), std::plus<REAL>());
  return *this;
}

/// Subtracts m from *this in place
template <class M>
XMATRIXN& operator-=(const M& m)
{
  #ifndef NEXCEPT
  if (_rows != m.rows() || _columns != m.columns())
    throw MissizeException(); 
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();
  #endif
  
  if (_rows > 0 && _columns > 0) 
    std::transform(column_iterator_begin(), column_iterator_end(), m.column_iterator_begin(), column_iterator_begin(), std::minus<REAL>());
  return *this;
}

/// Multiplies the diagonal matrix formed from d by the matrix m
template <class V, class W>
static W& diag_mult(const V& d, const XMATRIXN& m, W& result)
{
  #ifndef NEXCEPT
  if (d.size() != m.rows())
    throw MissizeException();
  if (sizeof(d.data()) != sizeof(m.data()))
    throw DataMismatchException();
  #endif

  result.resize(d.size(), m.columns());
  for (unsigned i=0; i< m.columns(); i++)
    std::transform(d.column_iterator_begin(), d.column_iterator_end(), m.column_iterator_begin()+m.rows()*i, result.column_iterator_begin()+result.rows()*i, std::multiplies<REAL>());

  return result;
}

/// Multiplies the diagonal matrix formed from d by the matrix transpose(m)
template <class V, class W>
static W& diag_mult_transpose(const V& d, const XMATRIXN& m, W& result)
{
  #ifndef NEXCEPT
  if (d.size() != m.columns())
    throw MissizeException();
  if (sizeof(d.data()) != sizeof(m.data()))
    throw DataMismatchException();
  #endif

  // copy the transpose of m to the result
  XMATRIXN::transpose(m, result);

  // do the specified number of times
  for (unsigned i=0; i< m.rows(); i++)
    CBLAS::scal(d.size(), d[i], result.data()+i, result.leading_dim());

  return result;
}

/// Sets this to a m x n sized zero matrix
XMATRIXN& set_zero(unsigned m, unsigned n) 
{ 
  return resize(m,n).set_zero(); 
}

