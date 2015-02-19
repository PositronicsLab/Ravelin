/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 *
 * This file contains inline code specific to SharedVectorNf/SharedVectorNd.
 ****************************************************************************/

/// Returns the desired component of this vector
const REAL& operator[](unsigned i) const
{
  #ifndef NEXCEPT
  if (i > _len)
    throw InvalidIndexException();
  #endif
  return _data[i+_start];
}

/// Gets the appropriate data element
const REAL* data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= _len)
    throw InvalidIndexException();
  #endif
  return &_data[i+_start];
}

CONST_COLUMN_ITERATOR segment_iterator(unsigned start, unsigned end) const
{
  #ifndef NEXCEPT
  if (end < start || end > _len)
    throw InvalidIndexException();
  #endif

  CONST_COLUMN_ITERATOR i;
  i._sz = end - start;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = data() + start;
  return i;
}

CONST_COLUMN_ITERATOR column_iterator_begin() const
{
  CONST_COLUMN_ITERATOR i;
  i._count = 0;
  i._sz = size();
  i._data_start = i._current_data  = data();

  // if the increment is > 1, we have a shared row vector
  if (inc() > 1)
  {
    i._rows = 1;
    i._columns = i._sz;
    i._ld = inc();
  }
  else
  {
    i._rows = i._sz;
    i._columns = 1;
    i._ld = size();
  }

  return i;
}

CONST_COLUMN_ITERATOR column_iterator_end() const
{
  CONST_COLUMN_ITERATOR i;
  i._count = size();
  i._sz = size();
  i._data_start = data();
  i._current_data  = data() + size();

  // if the increment is > 1, we have a shared row vector
  if (inc() > 1)
  {
    i._rows = 1;
    i._columns = i._sz;
    i._ld = inc();
  }
  else
  {
    i._rows = i._sz;
    i._columns = 1;
    i._ld = size();
  }

  return i;
}

CONST_COLUMN_ITERATOR begin() const { return column_iterator_begin(); }
CONST_COLUMN_ITERATOR end() const { return column_iterator_end(); }

CONST_ROW_ITERATOR row_iterator_begin() const
{
  CONST_ROW_ITERATOR i;
  i._count = 0;
  i._sz = size();
  i._data_start = i._current_data  = data();

  // if the increment is > 1, we have a shared row vector
  if (inc() > 1)
  {
    i._rows = 1;
    i._columns = i._sz;
    i._ld = inc();
  }
  else
  {
    i._rows = i._sz;
    i._columns = 1;
    i._ld = size();
  }

  return i;
}

CONST_ROW_ITERATOR row_iterator_end() const
{
  CONST_ROW_ITERATOR i;
  i._count = size();
  i._sz = size();
  i._data_start = data();
  i._current_data  = data() + size();

  // if the increment is > 1, we have a shared row vector
  if (inc() > 1)
  {
    i._rows = 1;
    i._columns = i._sz;
    i._ld = inc();
  }
  else
  {
    i._rows = i._sz;
    i._columns = 1;
    i._ld = size();
  }

  return i;
}

/// Returns true if all components of this vector are not infinite and not NaN, and false otherwise
bool is_finite() const
{
  CONST_COLUMN_ITERATOR i = column_iterator_begin();
  while (i != i.end())
  {
    if (std::isnan(*i) || std::isinf(*i))
      return false;
    i++;
  }

  return true;
}

/// Gets a sub-vector from this vector
/**
 * \param start_idx the starting index (inclusive)
 * \param end_idx the ending index (exclusive)
 */
template <class V>
V& get_sub_vec(unsigned start, unsigned end, V& v) const
{
  #ifndef NEXCEPT
  if (start > end || end > _len)
    throw InvalidIndexException();
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  // resize the vector
  const unsigned SZ = end - start;
  v.resize(SZ);

  // see whether to exit now
  if (_len == 0 || SZ == 0)
    return v;

  // copy using BLAS
  CBLAS::copy(SZ, data(start), inc(), v.data(), v.inc());

  return v;
}

/// Compares two vectors lexographically
template <class V>
bool operator<(const V& v) const
{
  #ifndef NEXCEPT 
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  CONST_COLUMN_ITERATOR data = column_iterator_begin();
  CONST_COLUMN_ITERATOR vdata = v.column_iterator_begin();

  // compare up to the shorter of the two
  unsigned shorter = std::min(size(), v.size());
  for (unsigned i=0; i< shorter; i++, data++, vdata++)
  {
    if (OPS::rel_equal(*data, *vdata))
      continue;
    return (*data < *vdata);
  }
  
  // still here?  comparison was identical to this point; if this is shorter
  // than v, return true
  if (_len < v.size())
    return true;
  
  // vectors *may* be identical
  return false;
}

/// Compares two vectors
/**
 * \note this method exists solely for convenience of use with iterators 
*/
template <class V>
bool operator==(const V& v) const
{
  #ifndef NEXCEPT
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();
  #endif

  if (_len != v.size())
    return false;

  // use iterator  
  CONST_COLUMN_ITERATOR data = column_iterator_begin();
  CONST_COLUMN_ITERATOR vdata = v.column_iterator_begin();
  for (unsigned i=0; i< _len; i++, data++, vdata++)
    if (!OPS::rel_equal(*data, *vdata))
      return false;
    
  return true;
}

/// Computes the dot-product between two vectors
template <class V1, class V2>
static REAL dot(const V1& v1, const V2& v2)
{    
  #ifndef NEXCEPT
  if (v1.size() != v2.size())
    throw MissizeException();
  if (sizeof(v1.data()) != sizeof(v2.data()))
    throw DataMismatchException();
  #endif

  if (v1.size() == 0)
    return (REAL) 0.0;

  // use the BLAS routine for dot-product
  return CBLAS::dot(v1.size(), v1.data(), v1.inc(), v2.data(), v2.inc());
}

/// Computes the dot-product between two vectors
template <class V>
REAL dot(const V& v) const
{    
  return dot(*this, v);
}

/// Gets a subvector (not necessarily contiguous)
template <class V>
V& select(const std::vector<bool>& indices, V& v) const
{
  // setup the vector
  unsigned len = std::count(indices.begin(), indices.end(), true);
  #ifndef NEXCEPT
  if (len > _len)
    throw MissizeException();
  #endif
  v.resize(len);

  // use iterator
  CONST_COLUMN_ITERATOR iter = column_iterator_begin();
  COLUMN_ITERATOR viter = v.column_iterator_begin();
  for (unsigned i = 0; i < indices.size(); i++, iter++)
    if (indices[i])
      *viter++ = *iter;

  return v;
}

/// Gets a subvector (not necessarily contiguous)
template <class ForwardIterator, class V>
V& select(ForwardIterator idx_begin, ForwardIterator idx_end, V& v) const
{
  // setup the vector
  unsigned len = std::distance(idx_begin, idx_end);
  #ifndef NEXCEPT
  if (len > _len)
    throw MissizeException();
  #endif
  v.resize(len);

  // use iterator
  COLUMN_ITERATOR viter = v.column_iterator_begin();
  CONST_COLUMN_ITERATOR iter = column_iterator_begin();
  unsigned vidx = 0;
  for (ForwardIterator i = idx_begin; i != idx_end; i++)
    *viter++ = iter[*i];

  return v;
}

/// Computes the infinity norm of this vector
static REAL norm_inf(const CONST_SHAREDVECTORN& v) 
{
  REAL nrm = (REAL) 0.0;
  CONST_COLUMN_ITERATOR i = v.column_iterator_begin();
  while (i != i.end())
    nrm = std::max(nrm , std::fabs(*i++));  
  return nrm;
}

/// Computes the l1-norm of this vector
static REAL norm1(const CONST_SHAREDVECTORN& v) 
{
  REAL nrm = (REAL) 0.0;
  CONST_COLUMN_ITERATOR i = v.column_iterator_begin();
  while (i != i.end())
    nrm += std::fabs(*i++);  
  return nrm;
}

