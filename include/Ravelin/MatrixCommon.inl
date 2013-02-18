/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Gets an iterator to a block
ITERATOR block_iterator(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();

  ITERATOR i;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = i._current_data = data() + i._ld*col_start + row_start;
  return i;
}

/// Gets an iterator to a block
CONST_ITERATOR block_iterator(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end) const
{
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();

  CONST_ITERATOR i;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = i._current_data = data() + i._ld*col_start + row_start;
  return i;
}

template <class V>
V& get_row(unsigned row, V& v) const
{
  if (row > rows())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // setup the iterators 
  CONST_ITERATOR source = block_iterator(row, row+1, 0, columns());
  CONST_ITERATOR source_end = source.end();  
  v.resize(source._sz);
  ITERATOR target = v.begin();

  // iterate
  while (source != source_end)
    *target++ = *source++;

  return v;
}

template <class V>
V& get_column(unsigned column, V& v) const
{
  if (column > columns())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // setup the iterators 
  CONST_ITERATOR source = block_iterator(0, rows(), column, column+1);
  CONST_ITERATOR source_end = source.end();  
  v.resize(source._sz);
  ITERATOR target = v.begin();

  // iterate
  while (source != source_end)
    *target++ = *source++;

  return v;
}

template <class V>
MATRIXX& set_row(unsigned row, const V& v)
{
  if (row > rows())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // setup the iterators 
  ITERATOR target = block_iterator(row, row+1, 0, columns());
  ITERATOR target_end = target.end();  
  CONST_ITERATOR source = v.begin();

  // iterate
  while (target != target_end)
    *target++ = *source++;

  return *this;
}

template <class V>
MATRIXX& set_column(unsigned column, const V& v)
{
  if (column > columns())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // setup the iterators 
  ITERATOR target = block_iterator(0, rows(), column, column+1);
  ITERATOR target_end = target.end();  
  CONST_ITERATOR source = v.begin();

  // iterate
  while (target != target_end)
    *target++ = *source++;

  return *this;
}

template <class T, class U>
U& transpose_mult_transpose(const T& x, U& y) const
{
  if (sizeof(data()) != sizeof(x.data()))
    throw DataMismatchException();
  if (sizeof(data()) != sizeof(y.data()))
    throw DataMismatchException();

  unsigned xrows = x.rows();
  unsigned xcols = x.columns();
  if (xcols != this->rows())
    throw MissizeException();
  y.resize(this->columns(), xrows);
  if (rows() == 0 || columns() == 0 || xrows == 0)
  {
    y.set_zero();
    return y;
  }
  CBLAS::gemm(CblasColMajor, CblasTrans, CblasTrans, columns(), xrows, rows(), (REAL) 1.0, data(), leading_dim(), x.data(), x.leading_dim(), (REAL) 0.0, y.data(), y.leading_dim()); 
  return y;
}

template <class T, class U>
U& mult_transpose(const T& x, U& y) const
{
  if (sizeof(data()) != sizeof(x.data()))
    throw DataMismatchException();
  if (sizeof(data()) != sizeof(y.data()))
    throw DataMismatchException();

  unsigned xrows = x.rows();
  unsigned xcols = x.columns();
  if (xcols != this->columns())
    throw MissizeException();
  y.resize(this->rows(), rows);
  if (rows() == 0 || columns() == 0 || xrows == 0)
  {
    y.set_zero();
    return y;
  }
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasTrans, rows(), xrows, columns(), (REAL) 1.0, data(), leading_dim(), x.data(), x.leading_dim(), (REAL) 0.0, y.data(), y.leading_dim()); 
      return y;
}

template <class T, class U>
U& transpose_mult(const T& x, U& y) const
{
  if (sizeof(data()) != sizeof(x.data()))
    throw DataMismatchException();
  if (sizeof(data()) != sizeof(y.data()))
    throw DataMismatchException();

  unsigned xrows = x.rows();
  unsigned xcols = x.columns();
  if (xrows != this->rows())
    throw MissizeException();
  y.resize(this->columns(), xcols);
  if (rows() == 0 || columns() == 0 || xcols == 0)
  {
    y.set_zero();
    return y;
  }
//  if (xcols > 1)
    CBLAS::gemm(CblasColMajor, CblasTrans, CblasNoTrans, columns(), xcols, rows(), (REAL) 1.0, data(), leading_dim(), x.data(), x.leading_dim(), (REAL) 0.0, y.data(), y.leading_dim()); 
/*
  else
    CBLAS::gemv(CblasTrans, rows(), columns(), (REAL) 1.0, data(), leading_dim(), x.data(), x.inc(), (REAL) 0.0, y.data(), y.inc());
*/
  return y;
}

template <class T, class U>
U& mult(const T& x, U& y) const
{
  if (sizeof(data()) != sizeof(x.data()))
    throw DataMismatchException();
  if (sizeof(data()) != sizeof(y.data()))
    throw DataMismatchException();

  unsigned xrows = x.rows();
  unsigned xcols = x.columns();
  if (xrows != this->columns())
    throw MissizeException();
  y.resize(rows(), xcols);
  if (rows() == 0 || columns() == 0 || xcols == 0)
  {
    y.set_zero();
    return y;
  }
//  if (xcols > 1)
    CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rows(), xcols, columns(), (REAL) 1.0, data(), leading_dim(), x.data(), x.leading_dim(), (REAL) 0.0, y.data(), y.leading_dim()); 
/*
  else
    CBLAS::gemv(CblasNoTrans, rows(), columns(), (REAL) 1.0, data(), leading_dim(), x.data(), x.inc(), (REAL) 0.0, y.data(), y.inc());
*/
  return y;
}

/// Gets the specified sub matrix
/**
 * \param row_start the row to start (inclusive)
 * \param row_end the row to end (exclusive)
 * \param col_start the column to start (inclusive)
 * \param col_end the column to end (exclusive)
 * \param transpose determines whether to store the transpose of the submatrix into m
 * \return a (row_end - row_start) x (col_end - col_start) sized matrix
 */
template <class M>
M& get_sub_mat(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end, M& m, Transposition trans = eNoTranspose) const
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  if (row_start > row_end || row_end > rows() || col_start > col_end || col_end > columns())
    throw InvalidIndexException();

  // resize the matrix
  if (trans == eNoTranspose)
    m.resize(row_end - row_start, col_end - col_start);
  else
    m.resize(col_end - col_start, row_end - row_start);

  // see whether to exit now
  if (rows() == 0 || columns() == 0)
    return m;

  // copy each column using BLAS
  if (trans == eNoTranspose)
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), data()+row_start+(col_start+i) * rows(), 1, m.data()+i*m.rows(), 1);
  else
    for (unsigned i=0; i< m.rows(); i++)
      CBLAS::copy(m.columns(), data()+row_start+(col_start+i) * rows(), 1, m.data()+i, m.rows());
  
  return m;
}

/// Sets the specified sub matrix
/**
 * \param row_start the row to start (inclusive)
 * \param col_start the column to start (inclusive)
 * \param m the source matrix
 * \param transpose determines whether the m is to be transposed 
 * \note fails assertion if m is too large to insert into this
 */
template <class M>
MATRIXX& set_sub_mat(unsigned row_start, unsigned col_start, const M& m, Transposition trans = eNoTranspose)
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  if (trans == eNoTranspose)
  {
    if (row_start + m.rows() > rows() || col_start + m.columns() > columns())
      throw MissizeException();
  }
  else if (row_start + m.columns() > rows() || col_start + m.rows() > columns())
    throw MissizeException();

  // make sure there is data to copy
  if (m.rows() == 0 || m.columns() == 0)
    return *this;

  // copy each column of m using BLAS
  if (trans == eNoTranspose)
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), m.data()+i*m.rows(), 1, data()+row_start+(col_start+i) * rows(), 1);
  else
    for (unsigned i=0; i< m.columns(); i++)
      CBLAS::copy(m.rows(), m.data()+i*m.rows(), 1, data()+row_start+col_start*rows() + i, rows());

  return *this;
}

/// Gets a submatrix of columns (not necessarily a block)
template <class ForwardIterator, class M>
M& select_columns(ForwardIterator col_start, ForwardIterator col_end, M& m) const
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  // setup vectors of selections
  const unsigned ncols = std::distance(col_start, col_end);

  // resize matrix 
  m.resize(rows(), ncols);

  // make sure there is data to copy
  if (rows() == 0 || ncols == 0)
    return m;

  // populate m
  unsigned mi;
  ForwardIterator i;
  for (i=col_start, mi=0; i != col_end; i++, mi++)
    CBLAS::copy(rows(), data()+rows()*(*i), 1, m.data()+rows()*mi, 1);

  return m;
}

/// Gets a submatrix of rows (not necessarily a block)
template <class ForwardIterator, class M>
M& select_rows(ForwardIterator row_start, ForwardIterator row_end, M& m) const
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  // setup vectors of selections
  const unsigned nrows = std::distance(row_start, row_end);

  // resize matrix 
  m.resize(nrows, columns());

  // make sure there is data to copy
  if (nrows == 0 || columns() == 0)
    return m;

  // populate m
  unsigned mi;
  ForwardIterator i;
  for (i=row_start, mi=0; i != row_end; i++, mi++)
    CBLAS::copy(columns(), data()+*i, rows(), m.data()+mi, nrows);

  return m;
}

/// Gets a submatrix (not necessarily a block)
template <class ForwardIterator1, class ForwardIterator2, class X>
X& select(ForwardIterator1 row_start, ForwardIterator1 row_end, ForwardIterator2 col_start, ForwardIterator2 col_end, X& m) const
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  // setup vectors of selections
  const unsigned nrows = std::distance(row_start, row_end);
  const unsigned ncols = std::distance(col_start, col_end);

  // resize matrix 
  m.resize(nrows, ncols);

  // make sure there is data to copy
  if (nrows == 0 || ncols == 0)
    return m;

  // get pointers to data
  REAL* mdata = m.data();
  const REAL* data = &operator()(*row_start, *col_start);

  // determine difference between first and last row
  unsigned rowbegin = *row_start;
  unsigned rowend = *row_start;
  for (ForwardIterator1 i = row_start; i != row_end; i++)
    rowend = *i;
  unsigned row_sub = rowend - rowbegin;

  // outer loop is over columns 
  for (ForwardIterator2 j=col_start; j != col_end; )
  {
    for (ForwardIterator1 i = row_start; i != row_end; )
    {
      // copy the data
      *mdata = *data;
      mdata++;

      // determine how we need to advance the rows
      unsigned row_diff = *i;
      i++;
      row_diff -= *i;
      row_diff = -row_diff;

      // if we're able, advance data_start 
      if (i != row_end)
        data += row_diff;
    }

    // determine how we need to advance the columns
    unsigned col_diff = *j;
    j++;
    col_diff -= *j;
    col_diff = -col_diff;

    // advance data
    data += (col_diff * leading_dim());
    data -= row_sub; 
  }

  return m;
}

/// Gets a submatrix (not necessarily a block)
template <class X>
X& select(const std::vector<bool>& rows, const std::vector<bool>& cols, X& m) const
{
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();
  if (rows.size() != this->rows() || cols.size() != this->columns())
    throw MissizeException();

  // setup vectors of selections
  const unsigned nrows = std::count(rows.begin(), rows.end(), true);
  const unsigned ncols = std::count(cols.begin(), cols.end(), true);

  // resize matrix 
  m.resize(nrows, ncols);

  // make sure there is data to copy
  if (nrows == 0 || ncols == 0)
    return m;

  // get pointers to data
  ITERATOR mdata = m.begin();
  CONST_ITERATOR data = begin();

  // outer loop is over columns 
  for (unsigned j=0; j < cols.size(); j++)
  {
    for (unsigned i=0; i< rows.size(); i++, data++)
    {
      if (!cols[j] || !rows[i])
        continue;

      // copy the data
      *mdata++ = *data;
    }
  }

  return m;
}

template <class M>
M& select_square(const std::vector<bool>& indices, M& result) const
{
  return select(indices, indices, result);
}

template <class ForwardIterator, class M>
M& select_square(ForwardIterator start, ForwardIterator end, M& m) const
{
  return select(start, end, start, end, m);
}

