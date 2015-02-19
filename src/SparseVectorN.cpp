/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Creates a sparse vector from a dense vector
SPARSEVECTORN::SPARSEVECTORN(const VECTORN& x)
{
  // setup the vector size
  _size = x.size();

  // determine the nonzero elements
  _nelm = 0;
  shared_array<bool> nz_elms(new bool[x.size()]);
  for (unsigned i=0; i< x.size(); i++)
    if (std::fabs(x[i]) < EPS)
    {
      _nelm++;
      nz_elms[i] = true;
    }
    else
      nz_elms[i] = false;

  // declare memory
  _indices = shared_array<unsigned>(new unsigned[_nelm]);
  _data = shared_array<REAL>(new REAL[_nelm]);

  // setup data
  for (unsigned i=0, j=0; i< x.size(); i++)
    if (nz_elms[i])
    {
      _indices[j] = i;
      _data[j] = x[i];
      j++;
    }
}

SPARSEVECTORN::SPARSEVECTORN(unsigned n, const map<unsigned, REAL>& values)
{
  // declare memory
  _nelm = values.size();
  _size = n;
  _indices = shared_array<unsigned>(new unsigned[_nelm]);
  _data = shared_array<REAL>(new REAL[_nelm]);
  
  unsigned j=0;
  for (map<unsigned, REAL>::const_iterator i = values.begin(); i != values.end(); i++)
  {
    _indices[j] = i->first;
    _data[j] = i->second;
    j++;
  }
}

SPARSEVECTORN::SPARSEVECTORN(unsigned n, unsigned nelms, shared_array<unsigned> indices, shared_array<REAL> data)
{
  _size = n;
  _nelm = nelms;
  _indices = indices;
  _data = data;
}

/// Computes the dot product between a sparse vector and itself
REAL SPARSEVECTORN::square() const
{
  REAL result = 0;

  for (unsigned i=0; i< _nelm; i++)
    result += _data[i] * _data[i];

  return result;
}

/// Computes the dot product between a sparse vector and a dense vector
REAL SPARSEVECTORN::dot(const VECTORN& x) const
{
  REAL result = 0;

  #ifndef NEXCEPT
  if (x.size() != _size)
    throw MissizeException();
  #endif

  for (unsigned i=0; i< _nelm; i++)
    result += x[_indices[i]] * _data[i];

  return result;
}

/// Gets the dense version of the vector 
VECTORN& SPARSEVECTORN::to_dense(VECTORN& result) const
{
  result.set_zero(_size);
  for (unsigned i=0; i< _nelm; i++)
    result[_indices[i]] = _data[i];

  return result;
}

/// Multiplies this sparse vector by a scalar (in place)
SPARSEVECTORN& SPARSEVECTORN::operator*=(REAL scalar)
{
  for (unsigned i=0; i< _nelm; i++)
    _data[i] *= scalar;
  return *this;
}

/// Multiplies this sparse vector by a scalar and returns the result in a new vector
SPARSEVECTORN& SPARSEVECTORN::mult(REAL scalar, SPARSEVECTORN& result) const
{
  result._size = this->_size;
  result._nelm = this->_nelm;
  result._indices = shared_array<unsigned>(new unsigned[this->_nelm]);
  result._data = shared_array<REAL>(new REAL[this->_nelm]);
  for (unsigned i=0; i< this->_nelm; i++)
  {
    result._indices[i] = this->_indices[i];
    result._data[i] = this->_data[i] * scalar;
  }

  return result;
}

/// Negates this sparse vector in place 
SPARSEVECTORN& SPARSEVECTORN::negate()
{
  std::transform(_data.get(), _data.get()+_nelm, _data.get(), std::negate<REAL>());
  return *this;
}

