/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Class for spatial arithmetic
class SPARITH
{
  public:

/// Converts a matrix (type X) to a vector of spatial axes
template <class X>
static void from_matrix(const X& m, std::vector<SVELOCITY>& w)
{
  const unsigned SPATIAL_DIM = 6;
  assert(m.rows() == SPATIAL_DIM);
  w.resize(m.columns());
  CONST_COLUMN_ITERATOR data = m.column_iterator_begin();
  for (unsigned k=0, i=0; i< w.size(); i++, k+= SPATIAL_DIM)
    w[i] = SVELOCITY(data[k+0], data[k+1], data[k+2], data[k+3], data[k+4], data[k+5]);
}

/// Converts an STL vector of axes to a matrix (type X)
template <class X>
static X& to_matrix(const std::vector<SVELOCITY>& w, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, w.size());
  COLUMN_ITERATOR data = m.column_iterator_begin();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    VECTOR3 f = w[i].get_angular();  
    VECTOR3 t = w[i].get_linear();
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
  }

  return m;
}

/// Converts an STL vector of forces to a matrix (type X)
template <class X>
static X& to_matrix(const std::vector<SFORCE>& w, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, w.size());
  COLUMN_ITERATOR data = m.column_iterator_begin();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    VECTOR3 f = w[i].get_force();  
    VECTOR3 t = w[i].get_torque();
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
  }

  return m;
}

/// Converts an STL vector of momenta to a matrix (type X)
template <class X>
static X& to_matrix(const std::vector<SMOMENTUM>& w, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(SPATIAL_DIM, w.size());
  COLUMN_ITERATOR data = m.column_iterator_begin();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    VECTOR3 f = w[i].get_linear();  
    VECTOR3 t = w[i].get_angular();
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
  }

  return m;
}

/// Converts an STL vector of momenta to a matrix (type X)
template <class X>
static X& spatial_transpose_to_matrix(const std::vector<SMOMENTUM>& w, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(w.size(), SPATIAL_DIM);
  ROW_ITERATOR data = m.row_iterator_begin();
  for (unsigned k=0, i=0; i< w.size(); i++)
  {
    VECTOR3 f = w[i].get_linear();  
    VECTOR3 t = w[i].get_angular();
    data[k++] = t[0];  data[k++] = t[1];  data[k++] = t[2];
    data[k++] = f[0];  data[k++] = f[1];  data[k++] = f[2];
  }

  return m;
}

/// Converts an STL vector of spatial velocities to a force matrix (type X)
template <class X>
static X& transpose_to_matrix(const std::vector<SVELOCITY>& t, X& m)
{
  const unsigned SPATIAL_DIM = 6;
  m.resize(t.size(), SPATIAL_DIM);
  for (unsigned i=0; i< t.size(); i++)
  {
    COLUMN_ITERATOR data = m.block_column_iterator_begin(i, i+1, 0, SPATIAL_DIM);
    VECTOR3 lin = t[i].get_linear();  
    VECTOR3 ang = t[i].get_angular();
    data[0] = lin[0]; 
    data[1] = lin[1]; 
    data[2] = lin[2]; 
    data[3] = ang[0]; 
    data[4] = ang[1]; 
    data[5] = ang[2];
  }

  return m;
}

/// Computes the "spatial dot product" between a vector of velocities and a vector of forces and returns the result in the matrix container (X)
template <class X>
static X& transpose_mult(const std::vector<SVELOCITY>& t, const std::vector<SFORCE>& w, X& result)
{
  result.resize(t.size(), w.size());
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< w.size(); j++)
      data[k++] = t[i].dot(w[j]);

  return result;
}

/// Computes the "spatial dot product" between a vector of velocities and a force and returns the result in the matrix container (X)
template <class X>
static X& transpose_mult(const std::vector<SVELOCITY>& t, const SFORCE& w, X& result)
{
  result.resize(t.size(), 1, false);
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0, k=0; i< t.size(); i++)
    data[k++] = t[i].dot(w);

  return result;
}

/// Computes the "spatial dot product" between a vector of velocities and a momentum and returns the result in the matrix container (X)
template <class X>
static X& transpose_mult(const std::vector<SVELOCITY>& t, const SMOMENTUM& w, X& result)
{
  result.resize(t.size(), 1, false);
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0, k=0; i< t.size(); i++)
    data[k++] = t[i].dot(w);

  return result;
}

/// Computes the "spatial dot product" between a vector of axes and a vector of momenta and returns the result in the matrix container (X)
template <class X>
static X& transpose_mult(const std::vector<SVELOCITY>& t, const std::vector<SMOMENTUM>& w, X& result)
{
  result.resize(t.size(), w.size());
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< w.size(); j++)
      data[k++] = t[i].dot(w[j]);

  return result;
}

/// Computes the "spatial dot product" between a vector of axes and a matrix or vector and returns the result in the matrix container (X)
template <class Y, class X>
static X& transpose_mult(const std::vector<SVELOCITY>& t, const Y& y, X& result)
{
  const unsigned SPATIAL_DIM = 6;
  result.resize(t.size(), y.columns(), false);
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0, k=0; i< t.size(); i++)
    for (unsigned j=0; j< y.columns(); j++)
      data[k++] = t[i].dot(y.column(j));

  return result;
}

/// Computes the "spatial dot product" between a vector of momenta and an axis and returns the result in the matrix container (X)
template <class X>
static X& transpose_mult(const std::vector<SMOMENTUM>& w, const SVELOCITY& t, X& result)
{
  result.resize(w.size());
  COLUMN_ITERATOR data = result.column_iterator_begin();
  for (unsigned i=0; i< w.size(); i++)
    data[i] = w[i].dot(t);

  return result;
}

  static MATRIXN& mult(const std::vector<SMOMENTUM>& Is, const MATRIXN& m, MATRIXN& result);
  static VECTORN& mult(const std::vector<SMOMENTUM>& Is, const VECTORN& v, VECTORN& result);
  static std::vector<SMOMENTUM>& mult(const SPATIAL_AB_INERTIA& I, const std::vector<SVELOCITY>& s, std::vector<SMOMENTUM>& result);
  static MATRIXN& mult(const SPATIAL_AB_INERTIA& I, const std::vector<SVELOCITY>& s, MATRIXN& result);
  static std::vector<SMOMENTUM>& mult(const SPATIAL_RB_INERTIA& I, const std::vector<SVELOCITY>& s, std::vector<SMOMENTUM>& result);
  static MATRIXN& mult(const SPATIAL_RB_INERTIA& I, const std::vector<SVELOCITY>& s, MATRIXN& result);
  static VECTORN& concat(const VECTORN& v, const SFORCE& w, VECTORN& result);
  static VECTORN& concat(const VECTORN& v, const SMOMENTUM& w, VECTORN& result);
  static SVELOCITY mult(const std::vector<SVELOCITY>& a, const VECTORN& v);
  static SACCEL transform_accel(boost::shared_ptr<const POSE3> target, const SACCEL& a);
  static std::vector<SACCEL>& transform_accel(boost::shared_ptr<const POSE3> target, const std::vector<SACCEL>& asrc, std::vector<SACCEL>& atgt);
  static void transform_accel(boost::shared_ptr<const POSE3> target, const SACCEL& w, const VECTOR3& r, const MATRIX3& E, SACCEL& result);

}; // end class 

