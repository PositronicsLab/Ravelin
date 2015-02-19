/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef COLUMN_ITERATOR 
#error This class is not to be included by the user directly. Use ColumnIteratorf.h or ColumnIteratord.h instead. 
#endif

/// A construct for iterating over a rectangular block of a matrix
class COLUMN_ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
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
  friend class CONST_COLUMN_ITERATOR;

  public:
    COLUMN_ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    COLUMN_ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    COLUMN_ITERATOR& operator+=(int n) 
    {
      // if there are no rows, verify that n = 0
      if (_rows == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count += n;

      // update the current data
      _current_data = _data_start + (_count / _rows)*_ld + (_count % _rows); 

      return *this; 
    }

    COLUMN_ITERATOR& operator-=(int n) 
    { 
      // if there are no rows, verify that n = 0
      if (_rows == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count -= n;

      // update the current data
      _current_data = _data_start + (_count / _rows)*_ld + (_count % _rows); 

      return *this; 
    }

    REAL& operator[](int i) const
    {
      // verify user not doing something wrong
      assert (_rows > 0);

      int j = i + _count;       
      if (j > _sz || j < 0)
        throw std::runtime_error("Data outside of scope!");
      return _data_start[(j / _rows)*_ld + (j % _rows)];
    }

    COLUMN_ITERATOR operator+(int n) const
    {
      COLUMN_ITERATOR b = *this;
      b += n;
      return b;
    }

    COLUMN_ITERATOR operator-(int n) const
    {
      COLUMN_ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const COLUMN_ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const COLUMN_ITERATOR& j) const
    {
      return _count > j._count;
    }

    REAL& operator*() const 
    {
      if (_count < 0 || _count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const COLUMN_ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const COLUMN_ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _rows == b._rows &&
             _columns == b._columns);

      return (_count == b._count);
    }

    bool operator!=(const COLUMN_ITERATOR& j) const { return !operator==(j); }

    // prefix--
    COLUMN_ITERATOR& operator--() 
    { 
      // verify user not doing something wrong
      assert(_rows > 0);

      _count--; 
      _current_data--;
      if (_count % _rows == 0)
        _current_data -= (_ld - _rows);

       return *this; 
    }

    // prefix++
    COLUMN_ITERATOR& operator++() 
    { 
      // verify user not doing something wrong
      assert(_rows > 0);

      _count++;
      _current_data++;
      if (_count % _rows == 0)
        _current_data += (_ld - _rows);

      return *this;
    }

    // postfix--
    COLUMN_ITERATOR operator--(int) 
    { 
      COLUMN_ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    COLUMN_ITERATOR operator++(int) 
    {
      COLUMN_ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    COLUMN_ITERATOR& operator=(const COLUMN_ITERATOR& i)
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
    int _count, _sz;
    REAL* _data_start;
    REAL* _current_data;
    unsigned _ld;
    unsigned _columns;
    unsigned _rows;
}; // end class

/// A construct for iterating over a rectangular block of a matrix
class CONST_COLUMN_ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
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
    CONST_COLUMN_ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    /// Converts a non-constant column iterator to a constant one
    CONST_COLUMN_ITERATOR(COLUMN_ITERATOR i)
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
    CONST_COLUMN_ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    CONST_COLUMN_ITERATOR& operator+=(int n) 
    {
      // if there are no rows, verify that n = 0
      if (_rows == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count += n;

      // update the current data
      _current_data = _data_start + (_count / _rows)*_ld + (_count % _rows); 

      return *this;    
    }

    CONST_COLUMN_ITERATOR& operator-=(int n) 
    { 
      // if there are no rows, verify that n = 0
      if (_rows == 0)
      {
        assert(n == 0);
        return *this;
      }

      // update the count
      _count -= n;

      // update the current data
      _current_data = _data_start + (_count / _rows)*_ld + (_count % _rows); 

      return *this;    
    }

    const REAL& operator[](int i) const
    {
      // verify user not doing something wrong
      assert(_rows > 0);

      int j = i + _count;       
      if (j > _sz || j < 0)
        throw std::runtime_error("Data outside of scope!");
      return _data_start[(j / _rows)*_ld + (j % _rows)];
    }

    CONST_COLUMN_ITERATOR operator+(int n) const
    {
      CONST_COLUMN_ITERATOR b = *this;
      b += n;
      return b;
    }

    CONST_COLUMN_ITERATOR operator-(int n) const
    {
      CONST_COLUMN_ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const CONST_COLUMN_ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const CONST_COLUMN_ITERATOR& j) const
    {
      return _count > j._count;
    }

    const REAL& operator*() const 
    {
      if (_count < 0 || _count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const CONST_COLUMN_ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const CONST_COLUMN_ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _columns == b._columns &&
             _rows == b._rows);

      return (_count == b._count);
    }

    bool operator!=(const CONST_COLUMN_ITERATOR& j) const { return !operator==(j); }

    // prefix--
    CONST_COLUMN_ITERATOR& operator--() 
    { 
      // verify user not doing something wrong
      assert(_rows > 0);

      _count--; 
      _current_data--;
      if (_count % _rows == 0)
        _current_data -= (_ld - _rows);

       return *this; 
    }

    // prefix++
    CONST_COLUMN_ITERATOR& operator++() 
    { 
      // verify user not doing something wrong
      assert(_rows > 0);

      _count++;
      _current_data++;
      if (_count % _rows == 0)
        _current_data += (_ld - _rows);

      return *this;
    }

    // postfix--
    CONST_COLUMN_ITERATOR operator--(int n) 
    { 
      // verify user not doing something wrong
      assert(_rows > 0);

      CONST_COLUMN_ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    CONST_COLUMN_ITERATOR operator++(int n) 
    {
      CONST_COLUMN_ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    CONST_COLUMN_ITERATOR& operator=(const CONST_COLUMN_ITERATOR& i)
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
    int _count, _sz;
    const REAL* _data_start;
    const REAL* _current_data;
    unsigned _ld;
    unsigned _columns;
    unsigned _rows;
}; // end class


