/// A sparse matrix represented in 'CSR' format
class SPARSEMATRIXN
{
  public:
    SPARSEMATRIXN();
    SPARSEMATRIXN(unsigned m, unsigned n, const std::map<std::pair<unsigned, unsigned>, REAL>& values);
    SPARSEMATRIXN(unsigned m, unsigned n, boost::shared_array<unsigned> ptr, boost::shared_array<unsigned> indices, boost::shared_array<REAL> data);
    SPARSEMATRIXN(const MATRIXN& m);
    static SPARSEMATRIXN identity(unsigned n);
    VECTORN& mult(const VECTORN& x, VECTORN& result) const;
    VECTORN& transpose_mult(const VECTORN& x, VECTORN& result) const;
    MATRIXN& mult(const MATRIXN& m, MATRIXN& result) const;
    MATRIXN& mult_transpose(const MATRIXN& m, MATRIXN& result) const;
    MATRIXN& transpose_mult(const MATRIXN& m, MATRIXN& result) const;
    MATRIXN& transpose_mult_transpose(const MATRIXN& m, MATRIXN& result) const;
    unsigned rows() const { return _rows; }
    unsigned columns() const { return _columns; }
    SPARSEMATRIXN get_sub_mat(unsigned rstart, unsigned rend, unsigned cstart, unsigned cend) const;
    SPARSEVECTORN& get_row(unsigned i, SPARSEVECTORN& row) const;
    SPARSEVECTORN& get_column(unsigned i, SPARSEVECTORN& column) const;
    unsigned* get_indices() { return _indices.get(); }
    unsigned* get_ptr() { return _ptr.get(); }
    REAL* get_data() { return _data.get(); }
    const unsigned* get_indices() const { return _indices.get(); }
    const unsigned* get_ptr() const { return _ptr.get(); }
    const REAL* get_data() const { return _data.get(); }
    SPARSEMATRIXN& operator-=(const SPARSEMATRIXN& m);
    SPARSEMATRIXN& operator+=(const SPARSEMATRIXN& m);
    SPARSEMATRIXN& operator*=(REAL scalar);
    SPARSEMATRIXN& negate();
    static SPARSEMATRIXN& outer_square(const VECTORN& g, SPARSEMATRIXN& result);
    static SPARSEMATRIXN& outer_square(const SPARSEVECTORN& v, SPARSEMATRIXN& result);
    MATRIXN& to_dense(MATRIXN& m) const;

  protected:
    SPARSEMATRIXN(unsigned m, unsigned n) { _rows = m; _columns = n; }
    boost::shared_array<unsigned> _indices;  // column indices
    boost::shared_array<unsigned> _ptr;      // starting indices for row i 
    boost::shared_array<REAL> _data;         // actual data
    unsigned _rows;
    unsigned _columns;

  private:
    void set(unsigned rows, unsigned columns, const std::map<std::pair<unsigned, unsigned>, REAL>& values);
}; // end class

std::ostream& operator<<(std::ostream& out, const SPARSEMATRIXN& s);

