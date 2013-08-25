/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const XMATRIXN& m)
{
  const unsigned OUTPUT_PRECISION = 8;

  if (m.rows() == 0 || m.columns() == 0)
  {
    out << "(empty)" << std::endl;
    return out;
  }
  
  for (unsigned i=0; i< m.rows(); i++)
  {
    for (unsigned j=0; j< m.columns()-1; j++)
      out << std::setprecision(OUTPUT_PRECISION) << m(i,j) << " ";
    out << std::setprecision(OUTPUT_PRECISION) << m(i,m.columns()-1) << std::endl;
  }
   
  return out;
}

/// Reads a matrix from the stream
std::istream& Ravelin::operator>>(std::istream& in, XMATRIXN& m)
{
  unsigned rows, columns;

  // read in the size
  in >> rows;
  in >> columns;

  // create the matrix
  m.resize(rows, columns);
  for (unsigned i=0; i< rows; i++)
    for (unsigned j=0; j< columns; j++) 
      in >> m(i,j);

  return in;
} 

/// Computes the l-infinity norm of this matrix
REAL XMATRIXN::norm_inf() const
{
  CONST_COLUMN_ITERATOR i = column_iterator_begin();
  REAL nrm = (REAL) 0.0;
  while (i != i.end())
    nrm = std::max(nrm, std::fabs(*i++));

  return nrm;
}


