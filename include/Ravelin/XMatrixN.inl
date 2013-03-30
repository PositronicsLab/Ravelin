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

/// Gets a constant shared vector for a row
CONST_SHAREDVECTORN row(unsigned i) const
{
  // get the starting offset
  const unsigned OFFSET = data() - _data.get();

  CONST_SHAREDVECTORN v;
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

/// Gets a shared vector for a column 
CONST_SHAREDVECTORN column(unsigned i) const
{
  // get the starting offset
  const unsigned OFFSET = data() - _data.get();

  CONST_SHAREDVECTORN v;
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

/// Gets a block as a constant shared matrix
CONST_SHAREDMATRIXN block(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end) const
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > _rows || col_end < col_start || col_end > _columns)
    throw InvalidIndexException();
  #endif

  // determine the offset
  const unsigned OFFSET = data() - _data.get();

  CONST_SHAREDMATRIXN m;
  m._data = _data;
  m._rows = row_end - row_start;
  m._columns = col_end - col_start;
  m._ld = leading_dim();
  m._start = OFFSET + m._ld * col_start + row_start;
  return m;  
}


/// Get an iterator to the beginning 
ITERATOR begin()
{
  ITERATOR i;
  i._count = 0;
  i._sz = _rows*_columns;
  i._ld = leading_dim();
  i._rows = _rows;
  i._columns = _columns;
  i._data_start = i._current_data = data();
  return i;
}

/// Get an iterator to the end
ITERATOR end()
{
  ITERATOR i;
  i._sz = _rows*_columns;
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = _rows;
  i._columns = _columns;
  i._data_start = data();
  i._current_data = data() + i._ld*i._columns;
  return i;
}

/// Get an iterator to the beginning 
CONST_ITERATOR begin() const
{
  CONST_ITERATOR i;
  i._count = 0;
  i._sz = _rows*_columns;
  i._ld = leading_dim();
  i._rows = _rows;
  i._columns = _columns;
  i._data_start = i._current_data = data();
  return i;
}

/// Get an iterator to the end
CONST_ITERATOR end() const
{
  CONST_ITERATOR i;
  i._sz = _rows*_columns;
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = _rows;
  i._columns = _columns;
  i._data_start = data();
  i._current_data = data() + i._ld*i._columns;
  return i;
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
    std::transform(begin(), end(), m.begin(), begin(), std::plus<REAL>());
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
    std::transform(begin(), end(), m.begin(), begin(), std::minus<REAL>());
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
    std::transform(d.begin(), d.end(), m.begin()+m.rows()*i, result.begin()+result.rows()*i, std::multiplies<REAL>());

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

/// Multiplies the diagonal matrix formed from d by the vector v
template <class V, class U, class W>
static W& diag_mult(const V& d, const U& v, W& result)
{
  #ifndef NEXCEPT
  if (d.size() != v.size())
    throw MissizeException();
  if (sizeof(d.data()) != sizeof(v.data()))
    throw DataMismatchException();
  if (sizeof(d.data()) != sizeof(result.data()))
    throw DataMismatchException();
  #endif

  result.resize(d.size());
  std::transform(d.begin(), d.end(), v.begin(), result.begin(), std::multiplies<REAL>());
  return result;
}

