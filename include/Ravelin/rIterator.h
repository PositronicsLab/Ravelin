/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef ITERATOR 
#error This class is not to be included by the user directly. Use fIterator.h or dIterator.h instead. 
#endif

/// A construct for iterating over a rectangular block of a matrix
class ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
{
  friend class MATRIXN;
  friend class VECTORN;
  friend class VECTOR3;
  friend class VECTOR2;
  friend class MATRIX3;
  friend class MATRIX2;
  friend class SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDVECTORN;
  friend class CONST_SHAREDMATRIXN;
  friend class SVECTOR6;
  friend class WRENCH;
  friend class TWIST;

  public:
    ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    ITERATOR& operator+=(int n) 
    {
      assert(n >= 0);
      const unsigned NSKIP = _ld - _rows + 1;
      for (int i=0; i< n; i++)
      {
        if (++_count % _rows == 0)
          _current_data += NSKIP;
        else
          _current_data++;
      }

      return *this; 
    }

    ITERATOR& operator-=(int n) 
    { 
      assert(n >= 0);
      const unsigned NSKIP = _ld - _rows + 1;
      for (int i=0; i< n; i++)  
      {
        if (--_count % _rows == 0)
          _current_data -= NSKIP;
        else
          _current_data--;
      }

      return *this; 
    }

    REAL& operator[](unsigned i) const
    {
      if (i > _sz)
        throw std::runtime_error("Data outside of scope!");
      const unsigned NSKIP = _ld - _rows + 1;
      REAL* data = _current_data - _count;
      for (unsigned j=0; j< i; )
      {
        if (++j % _rows == 0)
          data += NSKIP;
        else
          data++;
      }

      return *data;
    }

    ITERATOR operator+(int n) const
    {
      ITERATOR b = *this;
      b += n;
      return b;
    }

    ITERATOR operator-(int n) const
    {
      ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const ITERATOR& j) const
    {
      return _count > j._count;
    }

    REAL& operator*() const 
    {
      if (_count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _rows == b._rows &&
             _columns == b._columns);

      return (_count == b._count);
    }

    bool operator!=(const ITERATOR& j) const { return !operator==(j); }

    // prefix--
    ITERATOR& operator--() 
    { 
      _count--; 
      _current_data--;
      if (_count % _rows == 0)
        _current_data -= (_ld - _rows);

       return *this; 
    }

    // prefix++
    ITERATOR& operator++() 
    { 
      _count++;
      _current_data++;
      if (_count % _rows == 0)
        _current_data += (_ld - _rows);

      return *this;
    }

    // postfix--
    ITERATOR operator--(int n) 
    { 
      ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    ITERATOR operator++(int n) 
    {
      ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    ITERATOR& operator=(const ITERATOR& i)
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
class CONST_ITERATOR : public std::iterator<std::random_access_iterator_tag, REAL>
{
  friend class MATRIXN;
  friend class VECTORN;
  friend class VECTOR3;
  friend class VECTOR2;
  friend class MATRIX3;
  friend class MATRIX2;
  friend class SHAREDVECTORN;
  friend class SHAREDMATRIXN;
  friend class CONST_SHAREDVECTORN;
  friend class CONST_SHAREDMATRIXN;
  friend class SVECTOR6;
  friend class WRENCH;
  friend class TWIST;

  public:
    CONST_ITERATOR()
    {
      _data_start = _current_data = NULL;
      _count = 0;
      _sz = 0;
      _rows = 0;
      _ld = 0;
      _columns = 0;
    }

    /// Gets the iterator at the end of this block
    CONST_ITERATOR end() const
    {
      if (_count <= _sz)
        return *this + (_sz - _count);
      else
        return *this - (_count - _sz);
    }

    CONST_ITERATOR& operator+=(int n) 
    {
      assert(n >= 0);
      const unsigned NSKIP = _ld - _rows + 1;
      for (int i=0; i< n; i++)
      {
        if (++_count % _rows == 0)
          _current_data += NSKIP;
        else
          _current_data++;
      }

      return *this; 
    }

    CONST_ITERATOR& operator-=(int n) 
    { 
      assert(n >= 0);
      const unsigned NSKIP = _ld - _rows + 1;
      for (int i=0; i< n; i++)  
      {
        if (--_count % _rows == 0)
          _current_data -= NSKIP;
        else
          _current_data--;
      }

      return *this; 
    }

    const REAL& operator[](unsigned i) const
    {
      if (i > _sz)
        throw std::runtime_error("Data outside of scope!");
      const unsigned NSKIP = _ld - _rows + 1;
      const REAL* data = _current_data - _count;
      for (unsigned j=0; j< i; )
      {
        if (++j % _rows == 0)
          data += NSKIP;
        else
          data++;
      }

      return *data;
    }

    CONST_ITERATOR operator+(int n) const
    {
      CONST_ITERATOR b = *this;
      b += n;
      return b;
    }

    CONST_ITERATOR operator-(int n) const
    {
      CONST_ITERATOR b = *this;
      b -= n;
      return b;
    }
    
    bool operator<(const CONST_ITERATOR& j) const
    {
      return _count < j._count;
    }

    bool operator>(const CONST_ITERATOR& j) const
    {
      return _count > j._count;
    }

    const REAL& operator*() const 
    {
      if (_count >= _sz)
        throw std::runtime_error("Iterator outside of range!");
      return *_current_data; 
    }

    int operator-(const CONST_ITERATOR& b) const
    {
      return _count - b._count;
    }
 
    bool operator==(const CONST_ITERATOR& b) const
    {
      // verify that we're not comparing two dissimilar iterators
      assert(_data_start == b._data_start &&
             _sz == b._sz &&
             _ld == b._ld &&
             _columns == b._columns &&
             _rows == b._rows);

      return (_count == b._count);
    }

    bool operator!=(const CONST_ITERATOR& j) const { return !operator==(j); }

    // prefix--
    CONST_ITERATOR& operator--() 
    { 
      _count--; 
      _current_data--;
      if (_count % _rows == 0)
        _current_data -= (_ld - _rows);

       return *this; 
    }

    // prefix++
    CONST_ITERATOR& operator++() 
    { 
      _count++;
      _current_data++;
      if (_count % _rows == 0)
        _current_data += (_ld - _rows);

      return *this;
    }

    // postfix--
    CONST_ITERATOR operator--(int n) 
    { 
      CONST_ITERATOR b = *this; 
      this->operator--(); 
      return b; 
    }

    // postfix++ 
    CONST_ITERATOR operator++(int n) 
    {
      CONST_ITERATOR b = *this; 
      this->operator++(); 
      return b; 
    }

    // assignment operator
    CONST_ITERATOR& operator=(const CONST_ITERATOR& i)
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


