#ifndef SPARSEVECTORN
#error This class is not to be included by the user directly. Use SparseVectorNf.h or SparseVectorNd.h instead. 
#endif

/// A sparse vector represented in 'CSR' format
class SPARSEVECTORN
{
  public:
    SPARSEVECTORN() { _size = _nelm = 0; }
    SPARSEVECTORN(unsigned n, const std::map<unsigned, REAL>& values);
    SPARSEVECTORN(unsigned n, unsigned nnz, boost::shared_array<unsigned> indices, boost::shared_array<REAL> data);
    SPARSEVECTORN(const VECTORN& v);
    REAL dot(const VECTORN& x) const;
    REAL square() const;
    unsigned size() const { return _size; }
    unsigned num_elements() const { return _nelm; }
    unsigned* get_indices() { return _indices.get(); }
    REAL* get_data() { return _data.get(); }
    const unsigned* get_indices() const { return _indices.get(); }
    const REAL* get_data() const { return _data.get(); }
    VECTORN& to_dense(VECTORN& result) const;
    SPARSEVECTORN& negate();
    SPARSEVECTORN& operator*=(REAL scalar);
    SPARSEVECTORN& mult(REAL scalar, SPARSEVECTORN& result) const;

  protected:
    boost::shared_array<unsigned> _indices;   // indices of the data
    boost::shared_array<REAL> _data;          // 
    unsigned _size;
    unsigned _nelm;
}; // end class

std::ostream& operator<<(std::ostream& out, const SPARSEVECTORN& s);

