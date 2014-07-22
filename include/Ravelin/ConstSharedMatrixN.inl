/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// This file consists of routines for constant shared dynamically allocated 
// matrices
//////////////////////////////////////////////////////////////////////////////

/// Gets a constant shared vector for a row
CONST_SHAREDVECTORN row(unsigned i) const
{
  // get the starting offset
  const unsigned OFFSET = data() - _data.get();

  CONST_SHAREDVECTORN v;
  v._data = _data;
  v._start = OFFSET + i;
  v._inc = leading_dim();
  v._len = columns();
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
  v._len = rows();
  return v;
}


