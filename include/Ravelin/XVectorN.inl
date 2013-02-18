/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 *
 * This file contains inline code for VectorNf/VectorNd and 
 * SharedVectorNf/SharedVectorNd.
 ****************************************************************************/

ITERATOR segment_iterator(unsigned start, unsigned end)
{
  if (end < start || end > _len)
    throw InvalidIndexException();

  ITERATOR i;
  i._sz = end - start;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = data() + start;
  return i;
}

CONST_ITERATOR segment_iterator(unsigned start, unsigned end) const
{
  if (end < start || end > _len)
    throw InvalidIndexException();

  CONST_ITERATOR i;
  i._sz = end - start;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = data() + start;
  return i;
}

CONST_ITERATOR begin() const
{
  CONST_ITERATOR i;
  i._sz = _len;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = data();
  return i;
}

CONST_ITERATOR end() const
{
  CONST_ITERATOR i;
  i._sz = _len;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = data();
  i._current_data  = data() + _len;
  return i;
}

ITERATOR begin()
{
  ITERATOR i;
  i._sz = _len;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = data();
  return i;
}

ITERATOR end()
{
  ITERATOR i;
  i._sz = _len;
  i._ld = _len;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = data();
  i._current_data  = data() + _len;
  return i;
}

/// Sets the vector to the zero vector
XVECTORN& set_zero()
{
  std::fill(begin(), end(), (REAL) 0.0);
  return *this;
}

/// Sets the vector to the zero vector
XVECTORN& set_one()
{
  std::fill(begin(), end(), (REAL) 1.0);
  return *this;
}

/// Returns true if all components of this vector are not infinite and not NaN, and false otherwise
bool is_finite() const
{
  CONST_ITERATOR i = begin();
  while (i != i.end())
  {
    if (std::isnan(*i) || std::isinf(*i))
      return false;
    i++;
  }

  return true;
}

/// Negates the matrix in place
XVECTORN& negate()
{
  ITERATOR i = begin();
  while (i != i.end())
  {
    *i = -*i;
    i++;
  }

  return *this;
}

/// Multiplies this vector in place by a scalar
XVECTORN& operator*=(REAL scalar)
{
  if (_len > 0)
    CBLAS::scal(_len, scalar, data(), inc());  

  return *this;
}

/// Sets a sub-vector of this vector
/**
 * \param start_idx the starting index (inclusive)
 * \param v a (end_idx - start_idx + 1)-dimensional vector
 */
template <class V>
XVECTORN& set_sub_vec(unsigned start, const V& v)
{
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // get size of v
  const unsigned SZ = v.size();

  // see whether index is valid
  if (start + SZ > _len)
    throw InvalidIndexException();

  // check whether to exit now
  if (SZ == 0)
    return *this;

  // copy using BLAS
  CBLAS::copy(SZ, v.data(), v.inc(), data(start), inc());

  return *this;
}

/// Adds a vector from this one in place
template <class V>
XVECTORN& operator+=(const V& v)
{  
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  if (v.size() != _len)
    throw MissizeException();

  // use the BLAS routine for subtraction
  if (_len > 0)
    CBLAS::axpy(_len, (REAL) 1.0, v.data(), v.inc(), data(), inc());

  return *this;
}

/// Subtracts a vector from this one in place
template <class V>
XVECTORN& operator-=(const V& v)
{  
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  if (v.size() != _len)
    throw MissizeException();

  // use the BLAS routine for subtraction
  if (_len > 0)
    CBLAS::axpy(_len, (REAL) -1.0, v.data(), v.inc(), data(), inc());

  return *this;
}

/// Sets a subvector; other components are unchanged
template <class ForwardIterator, class V>
XVECTORN& set(ForwardIterator idx_begin, ForwardIterator idx_end, const V& v)
{
  unsigned len = std::distance(idx_begin, idx_end);
  if (len != v.size() || len > _len)
    throw MissizeException();
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  // use iterator
  ITERATOR iter = begin();
  CONST_ITERATOR viter = v.begin();
  for (ForwardIterator i = idx_begin; i != idx_end; i++)
    iter[*i] = *viter++;;

  return *this;
}

/// Gets a sub-vector from this vector
/**
 * \param start_idx the starting index (inclusive)
 * \param end_idx the ending index (exclusive)
 */
template <class V>
V& get_sub_vec(unsigned start, unsigned end, V& v) const
{
  if (start > end || end > _len)
    throw InvalidIndexException();
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

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
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  CONST_ITERATOR data = begin();
  CONST_ITERATOR vdata = v.begin();

  // compare up to the shorter of the two
  unsigned shorter = std::min(size(), v.size());
  for (unsigned i=0; i< shorter; i++, data++, vdata++)
  {
    if (rel_equal(*data, *vdata))
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
  if (sizeof(data()) != sizeof(v.data()))
    throw DataMismatchException();

  if (_len != v.size())
    return false;

  // use iterator  
  CONST_ITERATOR data = begin();
  CONST_ITERATOR vdata = v.begin();
  for (unsigned i=0; i< _len; i++, data++, vdata++)
    if (!rel_equal(*data, *vdata))
      return false;
    
  return true;
}

/// Computes the dot-product between two vectors
template <class V1, class V2>
static REAL dot(const V1& v1, const V2& v2)
{    
  if (v1.size() != v2.size())
    throw MissizeException();
  if (sizeof(v1.data()) != sizeof(v2.data()))
    throw DataMismatchException();

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
  if (len > _len)
    throw MissizeException();
  v.resize(len);

  // use iterator
  CONST_ITERATOR iter = begin();
  ITERATOR viter = v.begin();
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
  if (len > _len)
    throw MissizeException();
  v.resize(len);

  // use iterator
  ITERATOR viter = v.begin();
  CONST_ITERATOR iter = begin();
  unsigned vidx = 0;
  for (ForwardIterator i = idx_begin; i != idx_end; i++)
    *viter++ = iter[*i];

  return v;
}

/// Computes the infinity norm of this vector
static REAL norm_inf(const XVECTORN& v) 
{
  REAL nrm = (REAL) 0.0;
  CONST_ITERATOR i = v.begin();
  while (i != i.end())
    nrm = std::max(nrm , std::fabs(*i++));  
  return nrm;
}

/// Computes the l1-norm of this vector
static REAL norm1(const XVECTORN& v) 
{
  REAL nrm = (REAL) 0.0;
  CONST_ITERATOR i = v.begin();
  while (i != i.end())
    nrm += std::fabs(*i++);  
  return nrm;
}

