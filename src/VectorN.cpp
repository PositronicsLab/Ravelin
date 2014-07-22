/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 *
 * This file contains code specific to VectorNf/VectorNd.
 ****************************************************************************/

/// Default constructor - constructs an empty vector
VECTORN::VECTORN()
{
}

/// Constructs an uninitialized N-dimensional vector
VECTORN::VECTORN(unsigned N)
{
  _data.resize(N);
}

/// Constructs a vector from a VECTOR2
VECTORN::VECTORN(const VECTOR2& v)
{
  const unsigned LEN = 2;
  _data.resize(LEN);
  for (unsigned i=0; i< LEN; i++)
    _data[i] = v[i];
}

/// Constructs a vector from a VECTOR3
VECTORN::VECTORN(const VECTOR3& v)
{
  const unsigned LEN = 3;
  _data.resize(LEN);
  for (unsigned i=0; i< LEN; i++)
    _data[i] = v[i];
}

/// Constructs a N-dimensional vector from the given array
/**
 * \param array a N-dimensional (or larger) array
 */
VECTORN::VECTORN(unsigned N, const REAL* array)
{
  _data.resize(N);
  CBLAS::copy(N,array,1,data(),inc());
}

/// Copy constructor
VECTORN::VECTORN(const VECTORN& source)
{
  _data.resize(source.size());
  operator=(source);
}

/// Copy constructor
VECTORN::VECTORN(const MATRIXN& source)
{
  operator=(source);
}

/// Copy constructor
VECTORN::VECTORN(const SHAREDMATRIXN& source)
{
  operator=(source);
}

/// Copy constructor
VECTORN::VECTORN(const CONST_SHAREDMATRIXN& source)
{
  operator=(source);
}

/// Copy constructor
VECTORN::VECTORN(const SHAREDVECTORN& source)
{
  operator=(source);
}

/// Copy constructor
VECTORN::VECTORN(const CONST_SHAREDVECTORN& source)
{
  operator=(source);
}

/// Constructs a N-dimension vector from the list of double values
/**
 * \note There is no means in C++ to check the types of a list of variable
 * arguments.  If the variable arguments are not of type double, then
 * unexpected values will result. Constructing a vector using the
 * statement VECTORN::construct_variable(3, 1.0, 1.0, 0) is incorrect because
 * the programmer has assumed that the integer 0 will be converted to a double
 * type; this is not the case.
 */
VECTORN VECTORN::construct_variable(unsigned N, ...)
{
  VECTORN v(N);
  std::va_list lst;
  va_start(lst, N);
  for (unsigned i=0; i< N; i++)
    v[i] = (REAL) va_arg(lst, double);
  va_end(lst);
  return v;
}

/// Returns a N-dimensional one vector
VECTORN VECTORN::one(unsigned N)
{
  VECTORN v(N);
  v.set_one();
  return v;
}

/// Returns a N-dimensional zero vector
VECTORN VECTORN::zero(unsigned N)
{
  VECTORN v(N);
  v.set_zero();
  return v;
}

/// Resizes the vector
VECTORN& VECTORN::resize(unsigned m, unsigned n, bool preserve) 
{ 
  if (n == 1)
    resize(m, preserve); 
  else if (m == 1)
    resize(n, preserve);
  else
    assert(false);

  return *this; 
}

/// Copies another vector
VECTORN& VECTORN::operator=(const VECTOR2& source)
{
  // resize the vector if necessary
  if (size() != source.size())
    resize(source.size());

  // don't even worry about BLAS for copying
  _data[0] = source[0];
  _data[1] = source[1];

  return *this;
}

/// Copies another vector
VECTORN& VECTORN::operator=(const VECTOR3& source)
{
  // resize the vector if necessary
  if (size() != source.size())
    resize(source.size());

  // don't even worry about BLAS for copying
  _data[0] = source[0];
  _data[1] = source[1];
  _data[2] = source[2];

  return *this;
}

/// Copies another vector 
VECTORN& VECTORN::operator=(const VECTORN& source)
{  
  // resize this vector if necessary
  if (size() != source.size())
    resize(source.size());

  // use the BLAS routine for copying
  if (size() > 0)
    CBLAS::copy(source.size(),source.data(),1,_data.get(),1);

  return *this;
}

/// Copies another vector 
VECTORN& VECTORN::operator=(const SHAREDVECTORN& source)
{  
  // resize this vector if necessary
  if (size() != source.size())
    resize(source.size());

  // use the BLAS routine for copying
  if (size() > 0)
    CBLAS::copy(source.size(),source.data(),source.inc(),_data.get(),1);

  return *this;
}

/// Copies another vector 
VECTORN& VECTORN::operator=(const CONST_SHAREDVECTORN& source)
{  
  // resize this vector if necessary
  if (size() != source.size())
    resize(source.size());

  // use the BLAS routine for copying
  if (size() > 0)
    CBLAS::copy(source.size(),source.data(),source.inc(),_data.get(),1);

  return *this;
}

/// Sets this vector to a matrix 
VECTORN& VECTORN::operator=(const MATRIXN& source)
{  
  // resize this vector if necessary
  resize(source.rows(), source.columns());

  // use the BLAS routine for copying
  if (size() > 0)
    CBLAS::copy(source.rows(),source.data(),1,_data.get(),1);

  return *this;
}

/// Sets this vector to a matrix 
VECTORN& VECTORN::operator=(const SHAREDMATRIXN& source)
{  
  // resize this vector if necessary
  resize(source.rows(), source.columns());

  // use the C++ routine for copying
  if (size() > 0)
    std::copy(source.column_iterator_begin(),source.column_iterator_end(),column_iterator_begin());

  return *this;
}

/// Sets this vector to a matrix 
VECTORN& VECTORN::operator=(const CONST_SHAREDMATRIXN& source)
{  
  // resize this vector if necessary
  resize(source.rows(), source.columns());

  // use the C++ routine for copying
  if (size() > 0)
    std::copy(source.column_iterator_begin(),source.column_iterator_end(),column_iterator_begin());

  return *this;
}

/// Gets a subvector of this vector as a shared vector
SHAREDVECTORN VECTORN::segment(unsigned start_idx, unsigned end_idx) 
{
  // verify that we cna get the subvector
  #ifndef NEXCEPT
  if (start_idx > end_idx || end_idx > size())
    throw InvalidIndexException();
  #endif

  SHAREDVECTORN x;
  x._data = _data;
  x._start = start_idx;
  x._len = end_idx - start_idx;
  x._inc = 1;
  return x; 
}

/// Gets a subvector of this vector as a shared vector
CONST_SHAREDVECTORN VECTORN::segment(unsigned start_idx, unsigned end_idx) const 
{
  // verify that we cna get the subvector
  #ifndef NEXCEPT
  if (start_idx > end_idx || end_idx > size())
    throw InvalidIndexException();
  #endif

  CONST_SHAREDVECTORN x;
  x._data = _data;
  x._start = start_idx;
  x._len = end_idx - start_idx;
  x._inc = 1;
  return x; 
}

/// Assigns this vector to a scalar
VECTORN& VECTORN::operator=(REAL scalar)
{
  resize(1);
  _data[0] = scalar;
  return *this;
}

/// Concatenates two vectors together and stores the result in a third
VECTORN& VECTORN::concat(const VECTORN& v1, const VECTORN& v2, VECTORN& result)
{
  result.resize(v1.size() + v2.size());
  result.set_sub_vec(0, v1);
  result.set_sub_vec(v1.size(), v2);
  return result;
}

/// Parses a string for a vector value
VECTORN& VECTORN::parse(const std::string& s, VECTORN& values)
{
  std::list<std::string> plist;

  // make a copy of the string
  std::string copy = s;

  while (true)
  {
    // get the portion of the string before the delimiter
    size_t space_idx = copy.find_first_of(" \t");
    size_t comma_idx = copy.find(',');
        
    // case 1: no delimiter found
    if (space_idx == std::string::npos && comma_idx == std::string::npos)
    {
      plist.push_back(copy);
      break;
    }
        
    // case 2: delimiter found
    if ((space_idx != std::string::npos && space_idx < comma_idx) || comma_idx == std::string::npos)
    {
      plist.push_back(copy.substr(0,space_idx));
      copy = copy.substr(space_idx);
    }
    else
    {
      plist.push_back(copy.substr(0,comma_idx));
      copy = copy.substr(comma_idx);
    }
        
    // get the new string
    size_t firstidx = copy.find_first_not_of(" ,");
    if (firstidx == std::string::npos)
      break;
    else
      copy = copy.substr(firstidx);        
  }

  // convert the list to a Vector
  values.resize(plist.size());
  unsigned idx = 0;
  for (std::list<std::string>::const_iterator i=plist.begin(); i != plist.end(); i++)
  {
    if (strcasecmp(i->c_str(), "inf") == 0)
      values[idx] = std::numeric_limits<REAL>::infinity();
    else if (strcasecmp(i->c_str(), "-inf") == 0)
      values[idx] = -std::numeric_limits<REAL>::infinity();
    else
      values[idx] = (REAL) atof(i->c_str());
    idx++;
  }
  
  return values;  
}

#define XVECTORN VECTORN
#include "XVectorN.cpp"
#undef XVECTORN

