/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 *
 * This file contains inline code specific to SharedVectorNf/SharedVectorNd.
 ****************************************************************************/

/// Returns the desired component of this vector
REAL& operator[](unsigned i)
{
  #ifndef NEXCEPT
  if (i > _len)
    throw InvalidIndexException();
  #endif
  return _data[i*inc()+_start];
}

/// Returns the desired component of this vector
const REAL& operator[](unsigned i) const
{
  #ifndef NEXCEPT
  if (i > _len)
    throw InvalidIndexException();
  #endif
  return _data[(i*inc())+_start];
}

/// Gets the appropriate data element
REAL* data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= _len)
    throw InvalidIndexException();
  #endif
  return &_data[i*inc()+_start];
}

/// Gets the appropriate data element
const REAL* data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= _len)
    throw InvalidIndexException();
  #endif
  return &_data[i*inc()+_start];
}


