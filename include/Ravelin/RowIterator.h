/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef ROW_ITERATOR 
#error This class is not to be included by the user directly. Use RowIteratorf.h or RowIteratord.h instead. 
#endif

/// A construct for iterating over a rectangular block of a matrix
class ROW_ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
{
  friend class MATRIXN;
  friend class VECTORN;
  friend class VECTOR3;
  friend class ORIGIN3;
  friend class VECTOR2;
  friend class ORIGIN2;
  friend class MATRIX3;
  friend class MATRIX2;
  friend class SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDVECTORN;
  friend class CONST_SHAREDMATRIXN;
  friend class SVECTOR6;
  friend class SFORCE;
  friend class SVELOCITY;
  friend class SACCEL;
  friend class CONST_ROW_ITERATOR;

  public:
    ROW_ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    ROW_ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    ROW_ITERATOR& operator+=(int n) 
    {
      // if there are no columns, verify that n = 0
      if (_columns == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count += n;

      // update the current data
      _current_data = _data_start + (_count / _columns) + (_count % _columns)*_ld; 

      return *this; 
    }

    ROW_ITERATOR& operator-=(int n) 
    { 
      // if there are no columns, verify that n = 0
      if (_columns == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count -= n;

      // update the current data
      _current_data = _data_start + (_count / _columns) + (_count % _columns)*_ld; 

      return *this; 
    }

    REAL& operator[](int i) const
    {
      assert(_columns > 0);
      int j = i + _count;
      if (j < 0 && j > _sz)
        throw std::runtime_error("Data outside of scope!");
      return _data_start[(j / _columns) + (j % _columns)*_ld];
    }

    ROW_ITERATOR operator+(int n) const
    {
      ROW_ITERATOR b = *this;
      b += n;
      return b;
    }

    ROW_ITERATOR operator-(int n) const
    {
      ROW_ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const ROW_ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const ROW_ITERATOR& j) const
    {
      return _count > j._count;
    }

    REAL& operator*() const 
    {
      if (_count < 0 || _count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const ROW_ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const ROW_ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _rows == b._rows &&
             _columns == b._columns);

      return (_count == b._count);
    }

    bool operator!=(const ROW_ITERATOR& j) const { return !operator==(j); }

    // prefix--
    ROW_ITERATOR& operator--() 
    { 
      assert(_columns > 0);
      if (--_count % _columns == 0)
        _current_data = _data_start + (_count / _columns);
      else
        _current_data -= _ld;

       return *this; 
    }

    // prefix++
    ROW_ITERATOR& operator++() 
    {
      assert(_columns > 0); 
      if (++_count % _columns == 0)
        _current_data = _data_start + (_count / _columns);
      else
        _current_data += _ld;

      return *this;
    }

    // postfix--
    ROW_ITERATOR operator--(int) 
    { 
      ROW_ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    ROW_ITERATOR operator++(int) 
    {
      ROW_ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    ROW_ITERATOR& operator=(const ROW_ITERATOR& i)
    {
      _current_data = i._current_data;
      _count = i._count;
      _sz = i._sz;
      _data_start = i._data_start;
      _ld = i._ld;
      _columns = i._columns;
      _rows = i._rows;
      return *this;
    }

  protected:
    int _count;            // number of iterator increments
    int _sz;               // size of the data block
    REAL* _data_start;     // pointer to the start of the data block
    REAL* _current_data;   // pointer to data corresponding to iterator state
    unsigned _ld;          // leading dimension of matrix
    unsigned _columns;     // columns of the matrix
    unsigned _rows;        // rows of the matrix
}; // end class

/// A construct for iterating over a rectangular block of a matrix
class CONST_ROW_ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
{
  friend class MATRIXN;
  friend class VECTORN;
  friend class VECTOR3;
  friend class ORIGIN3;
  friend class VECTOR2;
  friend class ORIGIN2;
  friend class MATRIX3;
  friend class MATRIX2;
  friend class SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDVECTORN;
  friend class CONST_SHAREDMATRIXN;
  friend class SVECTOR6;
  friend class SFORCE;
  friend class SVELOCITY;
  friend class SACCEL;

  public:
    CONST_ROW_ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    /// Converts a non-constant row iterator to a constant one
    CONST_ROW_ITERATOR(ROW_ITERATOR i)
    {
      _data_start = i._data_start;
      _current_data = i._current_data;
      _count = i._count;
      _sz = i._sz;
      _rows = i._rows;
      _ld = i._ld;
      _columns = i._columns;
    }

    /// Gets the iterator at the end of this block
    CONST_ROW_ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    CONST_ROW_ITERATOR& operator+=(int n) 
    {
      // if there are no columns, verify that n = 0
      if (_columns == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count += n;

      // update the current data
      _current_data = _data_start + (_count / _columns) + (_count % _columns)*_ld; 

      return *this; 
    }

    CONST_ROW_ITERATOR& operator-=(int n) 
    { 
      // if there are no columns, verify that n = 0
      if (_columns == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count -= n;

      // update the current data
      _current_data = _data_start + (_count / _columns) + (_count % _columns)*_ld; 

      return *this; 
    }

    const REAL& operator[](int i) const
    {
      int j = i + _count;
      if (j < 0 && j > _sz)
        throw std::runtime_error("Data outside of scope!");
      return _data_start[(j / _columns) + (j % _columns)*_ld];
    }

    CONST_ROW_ITERATOR operator+(int n) const
    {
      CONST_ROW_ITERATOR b = *this;
      b += n;
      return b;
    }

    CONST_ROW_ITERATOR operator-(int n) const
    {
      CONST_ROW_ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const CONST_ROW_ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const CONST_ROW_ITERATOR& j) const
    {
      return _count > j._count;
    }

    const REAL& operator*() const 
    {
      if (_count >= _sz || _count < 0)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const CONST_ROW_ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const CONST_ROW_ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _columns == b._columns &&
             _rows == b._rows);

      return (_count == b._count);
    }

    bool operator!=(const CONST_ROW_ITERATOR& j) const { return !operator==(j); }

    // prefix--
    CONST_ROW_ITERATOR& operator--() 
    { 
      assert(_columns > 0);

      if (--_count % _columns == 0)
        _current_data = _data_start + (_count / _columns);
      else
        _current_data -= _ld;

       return *this; 
    }

    // prefix++
    CONST_ROW_ITERATOR& operator++() 
    {
      assert(_columns > 0);
 
      if (++_count % _columns == 0)
        _current_data = _data_start + (_count / _columns);
      else
        _current_data += _ld;

      return *this;
    }

    // postfix--
    CONST_ROW_ITERATOR operator--(int) 
    { 
      CONST_ROW_ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    CONST_ROW_ITERATOR operator++(int) 
    {
      CONST_ROW_ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    CONST_ROW_ITERATOR& operator=(const CONST_ROW_ITERATOR& i)
    {
      _current_data = i._current_data;
      _count = i._count;
      _sz = i._sz;
      _data_start = i._data_start;
      _ld = i._ld;
      _columns = i._columns;
      _rows = i._rows;
      return *this;
    }

  protected:
    int _count;                  // number of iterator increments
    int _sz;                     // size of the data block
    const REAL* _data_start;     // pointer to the start of the data block
    const REAL* _current_data;   // pointer to data corresponding to iterator state
    unsigned _ld;                // leading dimension of matrix
    unsigned _columns;           // columns of the matrix
    unsigned _rows;              // rows of the matrix
}; // end class


