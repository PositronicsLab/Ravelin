/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// This file includes common routines for all matrices, dynamic and fixed size, and mutable and constant (through inclusion of ConstMatrixCommon.inl) 
//////////////////////////////////////////////////////////////////////////////

// include the constant versions
#include "ConstMatrixCommon.inl"

/// Gets a column iterator to a block
COLUMN_ITERATOR block_column_iterator_begin(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();
  #endif

  COLUMN_ITERATOR i;
  i._count = 0;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = i._current_data = data() + i._ld*col_start + row_start;
  return i;
}

/// Gets a column iterator to a block
COLUMN_ITERATOR block_column_iterator_end(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();
  #endif

  COLUMN_ITERATOR i;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = data() + i._ld*col_start + row_start;
  i._current_data = i._data_start + i._sz;
  return i;
}

/// Gets a row iterator to a block
ROW_ITERATOR block_row_iterator_begin(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();
  #endif

  ROW_ITERATOR i;
  i._count = 0;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = i._current_data = data() + i._ld*col_start + row_start;
  return i;
}

/// Gets a row iterator to a block
ROW_ITERATOR block_row_iterator_end(unsigned row_start, unsigned row_end, unsigned col_start, unsigned col_end)
{
  #ifndef NEXCEPT
  if (row_end < row_start || row_end > rows() || col_end < col_start || col_end > columns())
    throw InvalidIndexException();
  #endif

  ROW_ITERATOR i;
  i._sz = (row_end - row_start)*(col_end - col_start);
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = row_end - row_start;
  i._columns = col_end - col_start;
  i._data_start = data() + i._ld*col_start + row_start;
  i._current_data = i._data_start + i._sz;
  return i;
}

template <class V>
MATRIXX& set_row(unsigned row, const V& v)
{
  #ifndef NEXCEPT
  if (row > rows())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  // setup the iterators 
  ROW_ITERATOR target = block_row_iterator_begin(row, row+1, 0, columns());
  ROW_ITERATOR target_end = target.end();  
  CONST_ROW_ITERATOR source = v.row_iterator_begin();

  // iterate
  while (target != target_end)
    *target++ = *source++;

  return *this;
}

template <class V>
MATRIXX& set_column(unsigned column, const V& v)
{
  #ifndef NEXCEPT
  if (column > columns())
    throw InvalidIndexException();

  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  // setup the iterators 
  COLUMN_ITERATOR target = block_column_iterator_begin(0, rows(), column, column+1);
  COLUMN_ITERATOR target_end = target.end();  
  CONST_COLUMN_ITERATOR source = v.column_iterator_begin();

  // iterate
  while (target != target_end)
    *target++ = *source++;

  return *this;
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
  #ifndef NEXCEPT
  if (sizeof(data()) != sizeof(m.data()))
    throw DataMismatchException();

  if (trans == eNoTranspose)
  {
    if (row_start + m.rows() > rows() || col_start + m.columns() > columns())
      throw MissizeException();
  }
  else if (row_start + m.columns() > rows() || col_start + m.rows() > columns())
    throw MissizeException();
  #endif

  // make sure there is data to copy
  if (m.rows() == 0 || m.columns() == 0)
    return *this;

  // copy each column of m using BLAS
  if (trans == eNoTranspose)
  {
    CONST_COLUMN_ITERATOR mIter = m.column_iterator_begin();
    COLUMN_ITERATOR tIter = block_column_iterator_begin(row_start, row_start+m.rows(), col_start, col_start+m.columns());
    while (mIter != mIter.end())
    {
      *tIter = *mIter;
      tIter++;
      mIter++;
    } 
  }
  else
  {
    CONST_ROW_ITERATOR mIter = m.row_iterator_begin();
    COLUMN_ITERATOR tIter = block_column_iterator_begin(row_start, row_start+m.columns(), col_start, col_start+m.rows());
    while (mIter != mIter.end())
    {
      *tIter = *mIter;
      tIter++;
      mIter++;
    } 
  }

  return *this;
}

/// Get an iterator to the beginning (iterates column by column) 
ROW_ITERATOR row_iterator_begin()
{
  ROW_ITERATOR i;
  i._count = 0;
  i._sz = rows()*columns();
  i._ld = leading_dim();
  i._rows = rows();
  i._columns = columns();
  i._data_start = i._current_data = data();
  return i;
}

/// Get an iterator to the end
ROW_ITERATOR row_iterator_end()
{
  ROW_ITERATOR i;
  i._sz = rows()*columns();
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = rows();
  i._columns = columns();
  i._data_start = data();
  i._current_data = data() + i._ld*i._columns;
  return i;
}

/// Get an iterator to the beginning (iterates column by column) 
COLUMN_ITERATOR column_iterator_begin()
{
  COLUMN_ITERATOR i;
  i._count = 0;
  i._sz = rows()*columns();
  i._ld = leading_dim();
  i._rows = rows();
  i._columns = columns();
  i._data_start = i._current_data = data();
  return i;
}

/// Get an iterator to the end
COLUMN_ITERATOR column_iterator_end()
{
  COLUMN_ITERATOR i;
  i._sz = rows()*columns();
  i._count = i._sz;
  i._ld = leading_dim();
  i._rows = rows();
  i._columns = columns();
  i._data_start = data();
  i._current_data = data() + i._ld*i._columns;
  return i;
}


