/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

SPARSEMATRIXN::SPARSEMATRIXN()
{
  _rows = _columns = 0;
  _nnz = 0;
  _nnz_capacity = 0;
  _ptr_capacity = 0;
  _stype = SPARSEMATRIXN::eCSR;
}

SPARSEMATRIXN::SPARSEMATRIXN(StorageType stype)
{
  _rows = _columns = 0;
  _nnz = 0;
  _nnz_capacity = 0;
  _ptr_capacity = 0;
  _stype = stype;
}

SPARSEMATRIXN::SPARSEMATRIXN(StorageType stype, unsigned m, unsigned n) 
{ 
  _rows = m; 
  _columns = n; 
  _stype = stype;
}

SPARSEMATRIXN::SPARSEMATRIXN(StorageType stype, unsigned m, unsigned n, const map<pair<unsigned, unsigned>, REAL>& values)
{
  _rows = _columns = 0;
  _nnz = 0;
  _nnz_capacity = 0;
  _ptr_capacity = 0;
  _stype = stype;
  set(m, n, values);
}

SPARSEMATRIXN::SPARSEMATRIXN(StorageType stype, unsigned m, unsigned n, shared_array<unsigned> ptr, shared_array<unsigned> indices, shared_array<REAL> data) 
{
  _rows = m;
  _columns = n;
  _stype = stype;
  _nnz = _ptr[_rows];
  _nnz_capacity = _nnz;
  _ptr_capacity = _rows+1;
  _data = data; 
  _ptr = ptr; 
  _indices = indices; 
}

/// Creates a sparse matrix from a dense matrix
SPARSEMATRIXN::SPARSEMATRIXN(const MATRIXN& m, REAL tol)
{
  _rows = m.rows();
  _columns = m.columns();
  _nnz = 0;
  _nnz_capacity = 0;
  _ptr_capacity = 0;
  _stype = eCSR;

  // determine number of non-zero values
  unsigned nv = 0;
  for (CONST_ROW_ITERATOR i=m.row_iterator_begin(); i!= i.end(); i++)
    if (std::fabs(*i) > tol)
      nv++;

  // setup arrays
  _nnz_capacity = nv;
  _ptr_capacity = m.rows()+1;
  _data = shared_array<REAL>(new REAL[nv]);
  _ptr = shared_array<unsigned>(new unsigned[m.rows()+1]);
  _indices = shared_array<unsigned>(new unsigned[nv]);

  // setup the arrays
  unsigned j = 0;
  unsigned k=0;
  _ptr[0] = j;
  CONST_COLUMN_ITERATOR i = m.column_iterator_begin();
  for (unsigned r=0; r< m.rows(); r++)
  {
    for (unsigned s=0; s< m.columns(); s++, i++)
      if (std::fabs(*i) > tol)
      {
        j++;
        _data[k] = *i;
        _indices[k] = s;
        k++;
      }
    _ptr[r+1] = j;
  }

  // set nnz
  _nnz = nv;
}

/// Creates a sparse matrix from a dense matrix
SPARSEMATRIXN::SPARSEMATRIXN(StorageType stype, const MATRIXN& m, REAL tol)
{
  _rows = m.rows();
  _columns = m.columns();
  _stype = stype;
  _nnz = 0;
  _nnz_capacity = 0;
  _ptr_capacity = 0;

  // determine number of non-zero values
  unsigned nv = 0;
  for (CONST_ROW_ITERATOR i=m.row_iterator_begin(); i!= i.end(); i++)
    if (std::fabs(*i) > tol)
      nv++;

  if (_stype == eCSR)
  {
    // setup arrays
    _nnz_capacity = nv;
    _ptr_capacity = m.rows()+1;
    _data = shared_array<REAL>(new REAL[nv]);
    _ptr = shared_array<unsigned>(new unsigned[m.rows()+1]);
    _indices = shared_array<unsigned>(new unsigned[nv]);

    // setup the arrays
    unsigned j=0;
    unsigned k=0;
    _ptr[0] = j;
    CONST_COLUMN_ITERATOR i = m.column_iterator_begin();
    for (unsigned r=0; r< m.rows(); r++)
    {
      for (unsigned s=0; s< m.columns(); s++, i++)
        if (std::fabs(*i) > tol)
        {
          j++;
          _data[k] = *i;
          _indices[k] = s;
          k++;
        }
      _ptr[r+1] = j;
    }

    // set nnz
    _nnz = nv;
  }
  else
  {
    assert(_stype == eCSC);

    // setup arrays
    _nnz_capacity = nv;
    _ptr_capacity = m.columns()+1;
    _data = shared_array<REAL>(new REAL[nv]);
    _ptr = shared_array<unsigned>(new unsigned[m.columns()+1]);
    _indices = shared_array<unsigned>(new unsigned[nv]);

    // setup ptr
    unsigned j = 0;
    unsigned k=0;
    _ptr[0] = j;
    CONST_ROW_ITERATOR i = m.row_iterator_begin();
    for (unsigned col=0; col< m.columns(); col++)
    {
      for (unsigned row=0; row< m.rows(); row++, i++)
        if (std::fabs(*i) > tol)
        {
          j++;
          _indices[k] = row;
          _data[k] = *i;
          k++;
        }
      _ptr[col+1] = j;
    }

    // set nnz
    _nnz = nv;
  }
}

/// Computes the infinity norm of this sparse matrix
REAL SPARSEMATRIXN::norm_inf() const
{
  if (_nnz == 0)
    return (REAL) 0.0;
  else
    return std::accumulate(_data.get(), _data.get()+_nnz, (REAL) 0.0);
}

/// Sets up an identity matrix from a sparse matrix
SPARSEMATRIXN SPARSEMATRIXN::identity(unsigned n)
{
  SPARSEMATRIXN m;

  // init matrix dimensions
  m._rows = m._columns = n;

  // init capacities
  m._nnz_capacity = n;
  m._ptr_capacity = n+1;

  // initialize arrays
  m._data = shared_array<REAL>(new REAL[n]);
  m._ptr = shared_array<unsigned>(new unsigned[n+1]);
  m._indices = shared_array<unsigned>(new unsigned[n]);

  // populate the matrix data, indices, and row pointers
  for (unsigned i=0; i< n; i++)
  {
    m._data[i] = (REAL) 1.0;
    m._indices[i] = i;
    m._ptr[i] = i;
  }
  m._ptr[n] = n;  

  // set nnz
  m._nnz = n;

  return m;
}

/// Sets up an identity matrix from a sparse matrix
SPARSEMATRIXN SPARSEMATRIXN::identity(StorageType stype, unsigned n)
{
  SPARSEMATRIXN m;

  // init matrix dimensions
  m._rows = m._columns = n;
  m._stype = stype;

  // init capacities
  m._nnz_capacity = n;
  m._ptr_capacity = n+1;

  // initialize arrays
  m._data = shared_array<REAL>(new REAL[n]);
  m._ptr = shared_array<unsigned>(new unsigned[n+1]);
  m._indices = shared_array<unsigned>(new unsigned[n]);

  // populate the matrix data, indices, and row pointers
  for (unsigned i=0; i< n; i++)
  {
    m._data[i] = (REAL) 1.0;
    m._indices[i] = i;
    m._ptr[i] = i;
  }
  m._ptr[n] = n;  

  // set nnz
  m._nnz = n;

  return m;
}

/// Gets all values from the matrix
void SPARSEMATRIXN::get_values(map<pair<unsigned, unsigned>, REAL>& values) const
{
  values.clear();

  // iterate through each row
  if (_stype == eCSR)
  {
    for (unsigned i=0; i< _rows; i++)
    {
      const REAL* data = _data.get() + _ptr[i];
      const unsigned* indices = _indices.get() + _ptr[i];
      const unsigned N = _ptr[i+1] - _ptr[i];
      for (unsigned j=0; j< N; j++)
        values[make_pair(i, indices[j])] = data[j];
    }
  }
  else
  {
    for (unsigned i=0; i< _columns; i++)
    {
      const REAL* data = _data.get() + _ptr[i];
      const unsigned* indices = _indices.get() + _ptr[i];
      const unsigned N = _ptr[i+1] - _ptr[i];
      for (unsigned j=0; j< N; j++)
        values[make_pair(indices[j], i)] = data[j];
    }
  }
}

/// Sets the column with the particular index
void SPARSEMATRIXN::set_column(unsigned col, const VECTORN& v)
{
  if (_stype == eCSR)
  {
    // get the values
    map<pair<unsigned, unsigned>, REAL> values; 
    get_values(values);

    // remove all values from column
    for (unsigned i=0; i< _rows; i++)
      values.erase(make_pair(i, col));

    // update the map
    for (unsigned i=0; i< v.size(); i++)
      if (v[i] > EPS || v[i] < -EPS)
        values[make_pair(i,col)] = v[i];

    // update the matrix
    set(_rows, _columns, values); 
  }
  else
  {
    assert(_stype == eCSC);

    // calculate the number of nonzeros in the vector
    unsigned nnz_v = std::count_if(v.column_iterator_begin(), v.column_iterator_end(), _1 > EPS || _1 < -EPS);

    // get the number of nonzeros in the column of the matrix
    unsigned nnz_col = _ptr[col+1] - _ptr[col];

    // three cases: number of zeros =, >, <
    if (nnz_v > nnz_col)
    {
      // determine how many more data entries we need
      unsigned nextra = nnz_v - nnz_col;

      // expand the capacity if necessary, preserving data
      if (_nnz_capacity < _nnz + nextra)
        set_capacities(_nnz + nextra, _ptr_capacity, true);
    
      // shift arrays
      REAL* data = _data.get();
      unsigned* indices = _indices.get();
      std::copy_backward(data+_ptr[col+1], data+_nnz, data+_ptr[col+1]+nextra);
      std::copy_backward(indices+_ptr[col+1], indices+_nnz, indices+_ptr[col+1]+nextra);
    
      // update ptr
      unsigned* ptr = _ptr.get();
      std::transform(ptr+col+1, ptr+_nnz+1, ptr+col+1, _1 + nextra);

      // update the number of nonzero entries
      _nnz += nextra;
    } 
    else if (nnz_v < nnz_col)
    {
      // determine how many fewer data entries there are
      unsigned nfewer = nnz_col - nnz_v;

      // shift arrays
      REAL* data = _data.get();
      unsigned* indices = _indices.get();
      std::copy(data+_ptr[col+1], data+_nnz, data+_ptr[col+1]-nfewer);
      std::copy(indices+_ptr[col+1], indices+_nnz, indices+_ptr[col+1]-nfewer);
    
      // update ptr
      unsigned* ptr = _ptr.get();
      std::transform(ptr+col+1, ptr+_nnz+1, ptr+col+1, _1 - nfewer);

      // update the number of nonzero entries
      _nnz -= nfewer;
    }

    // replace the non-zero values
    for (unsigned i=0, j=_ptr[col]; i< v.size(); i++)
      if (v[i] > EPS || v[i] < -EPS)
      {
        _data[j] = v[i];
        _indices[j] = i;
        j++;
      }
  }
}

/// Sets the row with the particular index
void SPARSEMATRIXN::set_row(unsigned row, const VECTORN& v)
{
  #ifndef NEXCEPT
  if (row >= _rows)
    throw InvalidIndexException();
  #endif 

  if (_stype == eCSR)
  {
    // calculate the number of nonzeros in the vector
    unsigned nnz_v = std::count_if(v.column_iterator_begin(), v.column_iterator_end(), _1 > EPS || _1 < -EPS);

    // get the number of nonzeros in the matrix
    unsigned nnz_row = _ptr[row+1] - _ptr[row];

    // three cases: number of zeros =, >, <
    if (nnz_v > nnz_row)
    {
      // determine how many more data entries we need
      unsigned nextra = nnz_v - nnz_row;

      // expand the capacity if necessary, preserving data
      if (_nnz_capacity < _nnz + nextra)
        set_capacities(_nnz + nextra, _ptr_capacity, true);
    
      // shift arrays
      REAL* data = _data.get();
      unsigned* indices = _indices.get();
      std::copy_backward(data+_ptr[row+1], data+_nnz, data+_ptr[row+1]+nextra);
      std::copy_backward(indices+_ptr[row+1], indices+_nnz, indices+_ptr[row+1]+nextra);
    
      // update ptr
      unsigned* ptr = _ptr.get();
      std::transform(ptr+row+1, ptr+_nnz+1, ptr+row+1, _1 + nextra);

      // update the number of nonzero entries
      _nnz += nextra;
    } 
    else if (nnz_v < nnz_row)
    {
      // determine how many fewer data entries there are
      unsigned nfewer = nnz_row - nnz_v;

      // shift arrays
      REAL* data = _data.get();
      unsigned* indices = _indices.get();
      std::copy(data+_ptr[row+1], data+_nnz, data+_ptr[row+1]-nfewer);
      std::copy(indices+_ptr[row+1], indices+_nnz, indices+_ptr[row+1]-nfewer);
    
      // update ptr
      unsigned* ptr = _ptr.get();
      std::transform(ptr+row+1, ptr+_nnz+1, ptr+row+1, _1 - nfewer);

      // update the number of nonzero entries
      _nnz -= nfewer;
    }

    // replace the non-zero values
    for (unsigned i=0, j=_ptr[row]; i< v.size(); i++)
      if (v[i] > EPS || v[i] < -EPS)
      {
        _data[j] = v[i];
        _indices[j] = i;
        j++;
      }
  }
  else
  {
    assert(_stype == eCSC);

    // get the values
    map<pair<unsigned, unsigned>, REAL> values; 
    get_values(values);

    // remove all values from row 
    for (unsigned i=0; i< _columns; i++)
      values.erase(make_pair(row, i));

    // update the map
    for (unsigned i=0; i< v.size(); i++)
      if (v[i] > EPS || v[i] < -EPS)
        values[make_pair(row, i)] = v[i];

    // update the matrix
    set(_rows, _columns, values); 
  }
}

/// Sets up a sparse matrix from a map 
void SPARSEMATRIXN::set(unsigned m, unsigned n, const map<pair<unsigned, unsigned>, REAL>& values)
{
  const unsigned nv = values.size();

  // setup rows and columns
  _rows = m;
  _columns = n;

  if (_stype == eCSR)
  {
    // setup arrays
    _nnz_capacity = nv;
    _ptr_capacity = m+1;
    _data = shared_array<REAL>(new REAL[nv]);
    _ptr = shared_array<unsigned>(new unsigned[m+1]);
    _indices = shared_array<unsigned>(new unsigned[nv]);

    // populate the matrix data
    map<pair<unsigned, unsigned>, REAL>::const_iterator i = values.begin();
    unsigned j;
    for (j=0; j< nv; j++, i++)
      _data[j] = i->second;

    // setup ptr
    j = 0;
    unsigned k=0;
    _ptr[0] = j;
    for (unsigned r=0; r< m; r++)
    {
      for (unsigned s=0; s< n; s++)
        if (values.find(make_pair(r, s)) != values.end())
        {
          j++;
          _indices[k++] = s;
        }
      _ptr[r+1] = j;
    }

    // set nnz
    _nnz = nv;
  }
  else
  {
    assert(_stype == eCSC);

    // setup a new values map
    map<pair<unsigned, unsigned>, REAL> values2;
    for (map<pair<unsigned, unsigned>, REAL>::const_iterator i = values.begin(); i != values.end(); i++)
      values2[make_pair(i->first.second, i->first.first)] = i->second;

    // setup arrays
    _nnz_capacity = nv;
    _ptr_capacity = n+1;
    _data = shared_array<REAL>(new REAL[nv]);
    _ptr = shared_array<unsigned>(new unsigned[n+1]);
    _indices = shared_array<unsigned>(new unsigned[nv]);

    // populate the matrix data
    map<pair<unsigned, unsigned>, REAL>::const_iterator i = values2.begin();
    unsigned j;
    for (j=0; j< nv; j++, i++)
      _data[j] = i->second;

    // setup ptr
    j = 0;
    unsigned k=0;
    _ptr[0] = j;
    for (unsigned col=0; col< n; col++)
    {
      for (unsigned row=0; row< m; row++)
        if (values.find(make_pair(row, col)) != values.end())
        {
          j++;
          _indices[k++] = row;
        }
      _ptr[col+1] = j;
    }

    // set nnz
    _nnz = nv;
  }
}

/// Gets a column of the sparse matrix as a sparse vector
SPARSEVECTORN& SPARSEMATRIXN::get_column(unsigned i, SPARSEVECTORN& result) const
{
  shared_array<unsigned> indices;
  shared_array<REAL> data;
  unsigned nelm;

  #ifndef NEXCEPT
  if (i >= _columns)
    throw InvalidIndexException();
  #endif

  // separate code depending on storage scheme
  if (_stype == eCSR)
  {
    // determine the number of elements
    nelm = 0;
    for (unsigned j=0; j< _ptr[_rows]; j++)
      if (_indices[j] == i)
        nelm++;

    // create new arrays
    indices = shared_array<unsigned>(new unsigned[nelm]);
    data = shared_array<REAL>(new REAL[nelm]);

    // get the data
    unsigned elm = 0;
    for (unsigned row=0; row< _rows; row++)
      for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
      {
        unsigned col = _indices[k];
  
        // look for early exit
        if (col > i)
          break; 

        // look for match
        if (col == i)
        {
          indices[elm] = row;
          data[elm] = _data[k];
          elm++;
        }
      }
  }
  else
  {
    assert(_stype == eCSC);

    // get the number of elements
    nelm = _ptr[i+1] - _ptr[i];
  
    // create new arrays
    indices = shared_array<unsigned>(new unsigned[nelm]);
    data = shared_array<REAL>(new REAL[nelm]);

    // setup the data
    for (unsigned j=_ptr[i], k=0; j< _ptr[i+1]; j++, k++)
    {
      indices[k] = _indices[j];
      data[k] = _data[j];
    }
  }

  // create the sparse vector
  result = SPARSEVECTORN(_rows, nelm, indices, data);
  return result;
}

/// Gets a column of the sparse matrix as a sparse vector
VECTORN& SPARSEMATRIXN::get_column(unsigned i, VECTORN& result) const
{
  // resize the resulting vector
  result.set_zero(_rows);

  // get the data
  REAL* rdata = result.data();

  #ifndef NEXCEPT
  if (i >= _columns)
    throw InvalidIndexException();
  #endif

  // separate code depending on storage scheme
  if (_stype == eCSR)
  {
    // get the data
    for (unsigned row=0; row< _rows; row++)
      for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
      {
        unsigned col = _indices[k];
  
        // look for early exit
        if (col > i)
          break; 

        // look for match
        if (col == i)
          rdata[row] = _data[k];
      }
  }
  else
  {
    assert(_stype == eCSC);

    // setup the data
    for (unsigned j=_ptr[i]; j< _ptr[i+1]; j++)
      rdata[_indices[j]] = _data[j];
  }

  return result;
}

/// Gets a row of the sparse matrix as a sparse vector
SPARSEVECTORN& SPARSEMATRIXN::get_row(unsigned i, SPARSEVECTORN& result) const
{
  shared_array<unsigned> indices;
  shared_array<REAL> data;
  unsigned nelm;

  #ifndef NEXCEPT
  if (i >= _rows)
    throw InvalidIndexException();
  #endif

  // different code depending on storage scheme
  if (_stype == eCSR)
  {
    // get the number of elements
    nelm = _ptr[i+1] - _ptr[i];
  
    // create new arrays
    indices = shared_array<unsigned>(new unsigned[nelm]);
    data = shared_array<REAL>(new REAL[nelm]);

    // setup the data
    for (unsigned j=_ptr[i], k=0; j< _ptr[i+1]; j++, k++)
    {
      indices[k] = _indices[j];
      data[k] = _data[j];
    }
  }
  else
  {
    assert(_stype == eCSC);

    // determine the number of elements
    nelm = 0;
    for (unsigned j=0; j< _ptr[_columns]; j++)
      if (_indices[j] == i)
        nelm++;

    // create new arrays
    indices = shared_array<unsigned>(new unsigned[nelm]);
    data = shared_array<REAL>(new REAL[nelm]);

    // get the data
    unsigned elm = 0;
    for (unsigned column=0; column < _columns; column++)
      for (unsigned k=_ptr[column]; k< _ptr[column+1]; k++)
      {
        unsigned row = _indices[k];
  
        // look for early exit
        if (row > i)
          break; 

        // look for match
        if (row == i)
        {
          indices[elm] = column;
          data[elm] = _data[k];
          elm++;
        }
      }
  }

  result = SPARSEVECTORN(_columns, nelm, indices, data);
  return result;
}

/// Gets a row of the sparse matrix as a sparse vector
VECTORN& SPARSEMATRIXN::get_row(unsigned i, VECTORN& result) const
{
  // resize the resulting vector
  result.set_zero(_columns);

  // get the data
  REAL* rdata = result.data();

  #ifndef NEXCEPT
  if (i >= _rows)
    throw InvalidIndexException();
  #endif

  // separate code depending on storage scheme
  if (_stype == eCSR)
  {
    // setup the data
    for (unsigned j=_ptr[i]; j< _ptr[i+1]; j++)
      rdata[_indices[j]] = _data[j];
  }
  else
  {
    assert(_stype == eCSC);

    // get the data
    for (unsigned col=0; col< _columns; col++)
      for (unsigned k=_ptr[col]; k< _ptr[col+1]; k++)
      {
        unsigned row = _indices[k];
  
        // look for early exit
        if (row > i)
          break; 

        // look for match
        if (row == i)
          rdata[col] = _data[k];
      }
  }

  return result;
}

/// Gets a submatrix of the sparse matrix
SPARSEMATRIXN SPARSEMATRIXN::get_sub_mat(unsigned rstart, unsigned rend, unsigned cstart, unsigned cend) const
{
  shared_array<REAL> data;
  shared_array<unsigned> ptr, indices;
  unsigned nv;

  // initialize the matrix
  SPARSEMATRIXN sub(_stype, rend - rstart, cend - cstart);

  #ifndef NEXCEPT
  if (rend < rstart || cend < cstart)
    throw InvalidIndexException();
  if (rend > _rows || cend > _columns)
    throw InvalidIndexException();
  #endif

  if (_stype == eCSR)
  {
    // determine how many values are in the sparse matrix
    nv = 0;
    for (unsigned row=rstart; row < rend; row++)
      for (unsigned col=_ptr[row]; col < _ptr[row+1]; col++)
        if (_indices[col] >= cstart && _indices[col] < cend)
          nv++;

    // setup arrays
    data = shared_array<REAL>(new REAL[nv]);
    ptr = shared_array<unsigned>(new unsigned[rend - rstart + 1]);
    indices = shared_array<unsigned>(new unsigned[nv]);

    // copy the data
    for (unsigned row=rstart, didx=0; row < rend; row++)
      for (unsigned col=_ptr[row]; col < _ptr[row+1]; col++)
        if (_indices[col] >= cstart && _indices[col] < cend)
          data[didx++] = _data[col];

    // setup ptr and indices
    ptr[0] = 0;
    unsigned j = 0;
    for (unsigned row=rstart; row < rend; row++)
    {
      ptr[row - rstart + 1] = ptr[row - rstart];
      for (unsigned col=_ptr[row]; col < _ptr[row+1]; col++)
        if (_indices[col] >= cstart && _indices[col] < cend)
        {
          ptr[row - rstart + 1]++;
          indices[j++] = _indices[col] - cstart;
        }
    }
  }
  else
  {
    assert(_stype == eCSC);

    // determine how many values are in the sparse matrix
    nv = 0;
    for (unsigned col=cstart; col < cend; col++)
      for (unsigned row=_ptr[col]; row < _ptr[col+1]; row++)
        if (_indices[row] >= rstart && _indices[row] < rend)
          nv++;

    // setup arrays
    data = shared_array<REAL>(new REAL[nv]);
    ptr = shared_array<unsigned>(new unsigned[cend - cstart + 1]);
    indices = shared_array<unsigned>(new unsigned[nv]);

    // copy the data
    for (unsigned col=cstart, didx=0; col < cend; col++)
      for (unsigned row=_ptr[col]; row < _ptr[col+1]; row++)
        if (_indices[row] >= rstart && _indices[row] < rend)
          data[didx++] = _data[row];

    // setup ptr and indices
    ptr[0] = 0;
    unsigned j = 0;
    for (unsigned col=cstart; col < cend; col++)
    {
      ptr[col - cstart + 1] = ptr[col - cstart];
      for (unsigned row=_ptr[col]; row < _ptr[col+1]; row++)
        if (_indices[row] >= rstart && _indices[row] < rend)
        {
          ptr[col - cstart + 1]++;
          indices[j++] = _indices[row] - rstart;
        }
    }
  }

  // setup pointers
  sub._ptr = ptr;
  sub._indices = indices;
  sub._data = data;
  sub._nnz = nv;

  return sub;
}

/// Sets the capacities of the arrays
void SPARSEMATRIXN::set_capacities(unsigned nnz_capacity, unsigned ptr_capacity, bool preserve = false)
{
  // increase ptr_capacity, just in case
  ptr_capacity++;

  // increase capacities, as necessary, if preserve is selected
  if (preserve)
  {
    ptr_capacity = std::max(ptr_capacity, _ptr_capacity);
    nnz_capacity = std::max(nnz_capacity, _nnz_capacity);
  }

  // create arrays
  shared_array<REAL> new_data(new REAL[nnz_capacity]);
  shared_array<unsigned> new_ptr(new unsigned[ptr_capacity]);
  shared_array<unsigned> new_indices(new unsigned[nnz_capacity]);

  // if there is no preservation, just setup the data
  if (!preserve)
  {
    _nnz_capacity = nnz_capacity;
    _ptr_capacity = ptr_capacity;
    _data = new_data;
    _ptr = new_ptr;
    _indices = new_indices;
    _nnz = 0;
  }
  else
  {
    unsigned len = (_stype == eCSR) ? _rows : _columns;
    std::copy(_data.get(), _data.get()+_nnz, new_data.get());
    std::copy(_ptr.get(), _ptr.get()+len+1, new_ptr.get());
    std::copy(_indices.get(), _indices.get()+_nnz, new_indices.get());
    _nnz_capacity = nnz_capacity;
    _ptr_capacity = ptr_capacity;
    _data = new_data;
    _ptr = new_ptr;
    _indices = new_indices;
  }
}

/// Multiplies this sparse matrix by a dense matrix
MATRIXN& SPARSEMATRIXN::mult(const MATRIXN& m, MATRIXN& result) const
{
  #ifndef NEXCEPT
  if (_columns != m.rows())
    throw MissizeException();
  #endif

  // setup the result matrix
  result.set_zero(_rows, m.columns());
  REAL* rdata = result.data();
  const REAL* mdata = m.data();

  if (_stype == eCSR)
  {
    // do the calculation
    for (unsigned col=0, idx=0, minc = 0; col < m.columns(); col++, minc += m.leading_dim())
      for (unsigned row=0; row < _rows; row++)
      {
        REAL dot = (REAL) 0.0;
        unsigned col_start = _ptr[row];
        unsigned col_end = _ptr[row+1];
        for (unsigned jj= col_start; jj< col_end; jj++)
          dot += _data[jj] * mdata[minc+_indices[jj]];
        rdata[idx++] += dot;
      }
  }
  else
  {
    assert(_stype == eCSC);

    for (unsigned k=0, minc=0, rinc=0; k< m.columns(); k++, minc += m.leading_dim(), rinc += result.leading_dim())
    for (unsigned col=0; col< _columns; col++)
      for (unsigned jj= _ptr[col]; jj < _ptr[col+1]; jj++)
        rdata[rinc + _indices[jj]] += _data[jj] * mdata[minc + col]; 
  }

  return result;
}

/// Multiplies this sparse matrix by a dense vector
VECTORN& SPARSEMATRIXN::mult(const VECTORN& x, VECTORN& result) const
{
  #ifndef NEXCEPT
  if (_columns != x.size())
    throw MissizeException();
  #endif

  // setup the result matrix
  result.set_zero(_rows);

  // get data vectors
  const REAL* xdata = x.data();
  REAL* rdata = result.data(); 

  if (_stype == eCSR)
  {
    for (unsigned row=0; row < _rows; row++)
    {
      REAL dot = (REAL) 0.0;
      unsigned row_start = _ptr[row];
      unsigned row_end = _ptr[row+1];
      for (unsigned jj= row_start; jj< row_end; jj++)
        dot += _data[jj] * xdata[_indices[jj]];

      rdata[row] += dot;
    }
  }
  else
  {
    assert(_stype == eCSC);

    for (unsigned col=0; col< _columns; col++)
      for (unsigned jj= _ptr[col]; jj < _ptr[col+1]; jj++)
        rdata[_indices[jj]] += _data[jj] * xdata[col]; 
  }

  return result;
}

/// Multiplies the transpose of this sparse matrix by a dense vector
VECTORN& SPARSEMATRIXN::transpose_mult(const VECTORN& x, VECTORN& result) const
{
  #ifndef NEXCEPT
  if (_rows != x.size())
    throw MissizeException();
  #endif

  // setup the result vector
  result.set_zero(_columns);

  // get data vectors
  const REAL* xdata = x.data();
  REAL* rdata = result.data(); 

  if (_stype == eCSR)
  {
    for (unsigned row=0; row< _rows; row++)
      for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
        rdata[_indices[k]] += _data[k] * xdata[row];
  }
  else
  {
    assert(_stype == eCSC);
    for (unsigned col=0; col< _columns; col++)
      for (unsigned k=_ptr[col]; k< _ptr[col+1]; k++)
        rdata[col] += _data[k] * xdata[_indices[k]];
  }

  return result;
}

/// Multiplies the transpose of this sparse matrix by a dense matrix
MATRIXN& SPARSEMATRIXN::transpose_mult(const MATRIXN& m, MATRIXN& result) const
{
  #ifndef NEXCEPT
  if (_rows != m.rows())
    throw MissizeException();
  #endif

  result.set_zero(_columns, m.columns());

  // get data vectors
  const REAL* mdata = m.data();
  REAL* rdata = result.data(); 

  if (_stype == eCSR)
  {
    for (unsigned col=0, idx=0, incr=0; col< m.columns(); col++, incr += result.leading_dim())
      for (unsigned row=0; row< _rows; row++, idx++)
        for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
          rdata[incr + _indices[k]] += _data[k] * mdata[idx];
  }
  else
  {
    assert(_stype == eCSC);
    for (unsigned idx=0, incr = 0, incm = 0; idx< m.columns(); idx++, incr += result.leading_dim(), incm += m.leading_dim())
      for (unsigned col=0; col< _columns; col++)
        for (unsigned k=_ptr[col]; k< _ptr[col+1]; k++)
          rdata[incr + col] += _data[k] * mdata[incm + _indices[k]];
  }

  return result;
}

/// Multiplies this matrix by the transpose of a dense matrix
MATRIXN& SPARSEMATRIXN::mult_transpose(const MATRIXN& m, MATRIXN& result) const
{
  #ifndef NEXCEPT
  if (_columns != m.columns())
    throw MissizeException();
  #endif

  // setup the result matrix
  result.set_zero(_rows, m.rows());

  // get data vectors
  const REAL* mdata = m.data();
  REAL* rdata = result.data(); 

  if (_stype == eCSR)
  {
    for (unsigned col=0, idx=0; col < m.rows(); col++)
      for (unsigned row=0; row < _rows; row++, idx++)
      {
        REAL dot = (REAL) 0.0;
        unsigned row_start = _ptr[row];
        unsigned row_end = _ptr[row+1];
        for (unsigned jj= row_start; jj< row_end; jj++)
          dot += _data[jj] * mdata[col+_indices[jj]*m.leading_dim()];
        rdata[idx] += dot;
      }
  }
  else
  {
    assert(_stype == eCSC);

    for (unsigned idx=0, rinc=0; idx< m.rows(); idx++, rinc+= result.leading_dim())
      for (unsigned col=0, minc=0; col< _columns; col++, minc += m.leading_dim())
        for (unsigned jj= _ptr[col]; jj < _ptr[col+1]; jj++)
          rdata[rinc+_indices[jj]] += _data[jj] * mdata[idx + minc]; 
  }

  return result;
}

/// Multiplies the transpose of this sparse matrix by the transpose of a dense matrix
MATRIXN& SPARSEMATRIXN::transpose_mult_transpose(const MATRIXN& m, MATRIXN& result) const
{
  #ifndef NEXCEPT
  if (_rows != m.columns())
    throw MissizeException();
  #endif

  // setup the result matrix
  result.set_zero(_columns, m.rows());

  // get data vectors
  const REAL* mdata = m.data();
  REAL* rdata = result.data(); 

  if (_stype == eCSR)
  {
    for (unsigned rowm=0, rinc=0; rowm< m.rows(); rowm++, rinc += result.leading_dim())
      for (unsigned row=0, minc=0; row< _rows; row++, minc += m.leading_dim())
        for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
        {
          unsigned col = _indices[k];
          rdata[col+rinc] += _data[k] * mdata[rowm+minc];
        }
  }
  else
  {
    assert(_stype == eCSC);

    for (unsigned idx=0, incr = 0; idx< m.rows(); idx++, incr += result.leading_dim())
      for (unsigned col=0; col< _columns; col++)
        for (unsigned k=_ptr[col]; k< _ptr[col+1]; k++)
          rdata[incr + col] += _data[k] * mdata[idx+_indices[k]*m.leading_dim()];
  }

  return result;
}

/// Gets a dense matrix from this sparse matrix
MATRIXN& SPARSEMATRIXN::to_dense(MATRIXN& m) const
{
  // resize m and make it zero
  m.set_zero(_rows, _columns);

  if (_stype == eCSR)
  {
    for (unsigned row=0; row< _rows; row++)
      for (unsigned k=_ptr[row]; k< _ptr[row+1]; k++)
        m(row, _indices[k]) = _data[k];
  }
  else
  {
    for (unsigned col=0; col< _columns; col++)
      for (unsigned k=_ptr[col]; k< _ptr[col+1]; k++)
        m(_indices[k], col) = _data[k];
  }

  return m;
}

/// Subtracts a sparse matrix from this one -- attempts to do it in place
SPARSEMATRIXN& SPARSEMATRIXN::operator-=(const SPARSEMATRIXN& m)
{
  // check rows/columns match up
  #ifndef NEXCEPT
  if (_rows != m._rows || _columns != m._columns)
    throw MissizeException();
  #endif

  typedef std::map<std::pair<unsigned, unsigned>, REAL> ValueMap; 
  ValueMap v1, v2;
  get_values(v1);
  m.get_values(v2);
  for (ValueMap::iterator i = v1.begin(); i != v1.end(); i++)
  {
    ValueMap::iterator j = v2.find(i->first);
    if (j == v2.end())
      continue;
    else
      i->second -= j->second;
  }

  // still must iterate through v2 in case keys in v2 aren't in v1 
  for (ValueMap::iterator i = v2.begin(); i != v2.end(); i++)
  {
    ValueMap::iterator j = v1.find(i->first);
    if (j != v1.end())
      continue;
    else
      v1[i->first] = -i->second;
  }

  // reset the values
  set(_rows, _columns, v1);

  return *this;
}

/// Adds a sparse matrix to this one -- attempts to do it in place
SPARSEMATRIXN& SPARSEMATRIXN::operator+=(const SPARSEMATRIXN& m)
{
  // check rows/columns match up
  #ifndef NEXCEPT
  if (_rows != m._rows || _columns != m._columns)
    throw MissizeException();
  #endif

  typedef std::map<std::pair<unsigned, unsigned>, REAL> ValueMap; 
  ValueMap v1, v2;
  get_values(v1);
  m.get_values(v2);
  for (ValueMap::iterator i = v1.begin(); i != v1.end(); i++)
  {
    ValueMap::iterator j = v2.find(i->first);
    if (j == v2.end())
      continue;
    else
      i->second += j->second;
  }

  // still must iterate through v2 in case keys in v2 aren't in v1 
  for (ValueMap::iterator i = v2.begin(); i != v2.end(); i++)
  {
    ValueMap::iterator j = v1.find(i->first);
    if (j != v1.end())
      continue;
    else
      v1[i->first] = i->second;
  }

  // reset the values
  set(_rows, _columns, v1);

  return *this;
}

/// Copies a sparse matrix to this
SPARSEMATRIXN& SPARSEMATRIXN::operator=(const SPARSEMATRIXN& m)
{
  // look for easiest exit
  if (m._rows == 0 || m._columns == 0)
  {
    _rows = m._rows;
    _columns = m._columns;
    _stype = m._stype;
    _nnz_capacity = 0;
    _ptr_capacity = 0;
    _data.reset();
    _ptr.reset();
    _indices.reset();
  }
  // see whether we can just copy the data without re-initing
  else
  {
    // make new nnz-dependent arrays if necessary
    if (_nnz_capacity < m._nnz)
    {
      _nnz_capacity = m._nnz;
      _data = shared_array<REAL>(new REAL[_nnz_capacity]);
      _indices = shared_array<unsigned>(new unsigned[_nnz_capacity]);
    }

    // make new ptr array, if necessary
    const unsigned PTR_SZ = (m._stype == eCSR) ? m._rows+1 : m._columns+1;
    if (_ptr_capacity < PTR_SZ)
    {
      _ptr_capacity = PTR_SZ;
      _ptr = shared_array<unsigned>(new unsigned[PTR_SZ]);
    }

    // copy everything
    std::copy(m._ptr.get(), m._ptr.get()+PTR_SZ, _ptr.get());
    std::copy(m._indices.get(), m._indices.get()+m._nnz, _indices.get());
    std::copy(m._data.get(), m._data.get()+m._nnz, _data.get());

    // update storage type, rows, and columns
    _nnz = m._nnz;
    _rows = m._rows;
    _columns = m._columns;
    _stype = m._stype;
  }

  return *this;
}

/// Multiplies a sparse matrix by a scalar
SPARSEMATRIXN& SPARSEMATRIXN::operator*=(REAL scalar)
{
  CBLAS::scal(_nnz, scalar, _data.get(), 1);
  return *this;
}

/// Negates this sparse matrix
SPARSEMATRIXN& SPARSEMATRIXN::negate()
{
  CBLAS::scal(_nnz, (REAL) -1.0, _data.get(), 1);
  return *this;
}

/// Calculates the outer product of a vector with itself and stores the result in a sparse matrix
SPARSEMATRIXN& SPARSEMATRIXN::outer_square(const SPARSEVECTORN& v, SPARSEMATRIXN& result)
{
  #ifdef REENTRANT
  FastThreadable<VECTORN> tmp;
  #else
  static FastThreadable<VECTORN> tmp;
  #endif

  // determine the size of the matrix
  unsigned n = v.size();

  // get the number of non-zero elements of v
  unsigned nz = v.num_elements();

  // get the indices of non-zero elements of v
  const unsigned* nz_indices = v.get_indices();
  const REAL* v_data = v.get_data();
  shared_array<bool> nz_elements(shared_array<bool>(new bool[n]));
  for (unsigned i=0; i< n; i++) nz_elements[i] = false;
  for (unsigned i=0; i< nz; i++) nz_elements[nz_indices[i]] = true;

  // setup ptr, indices, data
  shared_array<unsigned> ptr(new unsigned[n+1]);
  shared_array<unsigned> indices(new unsigned[nz*nz]);
  shared_array<REAL> data(new REAL[nz*nz]);

  // setup ptr
  ptr[0] = 0;
  for (unsigned i=0; i< n; i++)
  {
    ptr[i+1] = ptr[i];
    if (nz_elements[i])
      ptr[i+1] += nz;
  }

  // setup indices and data
  v.to_dense(tmp());
  for (unsigned i=0, k=0; i< n; i++)
  {
    // see whether to skip row
    if (!nz_elements[i])
      continue;

    for (unsigned j=0; j< n; j++)
    {
      // see whether to skip column
      if (!nz_elements[j])
        continue;

      // update indices
      assert(k < nz*nz);
      indices[k] = j;
      data[k] = tmp()[i]*tmp()[j];
      k++;
    }
  }

  // setup arrays
  result._rows = n;
  result._columns = n;
  result._ptr = ptr;
  result._indices = indices;
  result._data = data;
  result._nnz = nz*nz;

  return result;
}

/// Calculates the outer product of a vector with itself and stores the result in a sparse matrix
SPARSEMATRIXN& SPARSEMATRIXN::outer_square(const VECTORN& x, SPARSEMATRIXN& result)
{
  // determine the size of the matrix
  unsigned n = x.size();

  // determine the non-zero elements in x
  shared_array<bool> nz_elements(shared_array<bool>(new bool[n]));
  for (unsigned i=0; i< n; i++) nz_elements[i] = false;
  unsigned nz = 0;
  for (unsigned i=0; i< n; i++)
    if (std::fabs(x[i]) > EPS)
    {
      nz_elements[i] = true;
      nz++;
    }

  // setup ptr, indices, data
  shared_array<unsigned> ptr(new unsigned[n+1]);
  shared_array<unsigned> indices(new unsigned[nz*nz]);
  shared_array<REAL> data(new REAL[nz*nz]);

  // setup ptr
  ptr[0] = 0;
  for (unsigned i=0; i< n; i++)
  {
    ptr[i+1] = ptr[i];
    if (nz_elements[i])
      ptr[i+1] += nz;
  }

  // setup indices and data
  for (unsigned i=0, k=0; i< n; i++)
  {
    // see whether to skip row
    if (!nz_elements[i])
      continue;

    for (unsigned j=0; j< n; j++)
    {
      // see whether to skip column
      if (!nz_elements[j])
        continue;

      // update indices
      assert(k < nz*nz);
      indices[k] = j;
      data[k] = x[i]*x[j];
      k++;
    }
  }

  // setup arrays
  result._rows = n;
  result._columns = n;
  result._ptr = ptr;
  result._indices = indices;
  result._data = data;
  result._nnz = nz*nz;

  return result;
}

std::ostream& Ravelin::operator<<(std::ostream& out, const SPARSEMATRIXN& s)
{
  const unsigned* indices = s.get_indices();
  const unsigned* ptr = s.get_ptr();
  const REAL* data = s.get_data();
  vector<REAL> present(s.columns());

  out << "nnz: " << s.get_nnz() << std::endl;
  out << "ptr:";
  for (unsigned i=0; i<= s.rows(); i++)
    out << " " << ptr[i];
  out << std::endl;

  out << "indices:";
  for (unsigned i=0; i< ptr[s.rows()]; i++)
    out << " " << indices[i];
  out << std::endl;

  out << "data:";
  for (unsigned i=0; i< ptr[s.rows()]; i++)
    out << " " << data[i];
  out << std::endl;

  if (s.get_storage_type() == SPARSEMATRIXN::eCSR)
  {
    for (unsigned i=0; i< s.rows(); i++)
    {
      // mark all as not present
      for (unsigned j=0; j< present.size(); j++) present[j] = (REAL) 0.0;

      // mark ones that are present
      for (unsigned j=ptr[i]; j< ptr[i+1]; j++)
        present[indices[j]] = data[j];

      for (unsigned j=0; j< s.columns(); j++)
        out << present[j] << " ";
      out << std::endl;
    }
  }
  else
  {
    assert(s.get_storage_type() == SPARSEMATRIXN::eCSC);
    MATRIXN dense;
    out << s.to_dense(dense) << std::endl;
  }

  return out;
}

