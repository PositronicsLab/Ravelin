/// A sparse matrix
class SPARSEMATRIXN
{
  public:
    enum StorageType { eCSR, eCSC };

    SPARSEMATRIXN();
    SPARSEMATRIXN(StorageType s);
    SPARSEMATRIXN(StorageType s, unsigned m, unsigned n, const std::map<std::pair<unsigned, unsigned>, REAL>& values);
    SPARSEMATRIXN(StorageType s, unsigned m, unsigned n, boost::shared_array<unsigned> ptr, boost::shared_array<unsigned> indices, boost::shared_array<REAL> data);
    SPARSEMATRIXN(const MATRIXN& m, REAL tol=EPS);
    SPARSEMATRIXN(StorageType s, const MATRIXN& m, REAL tol=EPS);
    REAL norm_inf() const;
    static SPARSEMATRIXN identity(unsigned n);
    static SPARSEMATRIXN identity(StorageType stype, unsigned n);
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
    VECTORN& get_row(unsigned i, VECTORN& row) const;
    VECTORN& get_column(unsigned i, VECTORN& column) const;
    const unsigned* get_indices() const { return _indices.get(); }
    const unsigned* get_ptr() const { return _ptr.get(); }
    const REAL* get_data() const { return _data.get(); }
    void set_row(unsigned i, const VECTORN& v);
    void set_column(unsigned i, const VECTORN& v);
    SPARSEMATRIXN& operator=(const SPARSEMATRIXN& m);
    SPARSEMATRIXN& operator-=(const SPARSEMATRIXN& m);
    SPARSEMATRIXN& operator+=(const SPARSEMATRIXN& m);
    SPARSEMATRIXN& operator*=(REAL scalar);
    SPARSEMATRIXN& negate();
    static SPARSEMATRIXN& outer_square(const VECTORN& g, SPARSEMATRIXN& result);
    static SPARSEMATRIXN& outer_square(const SPARSEVECTORN& v, SPARSEMATRIXN& result);
    MATRIXN& to_dense(MATRIXN& m) const;
    void set_capacities(unsigned nnz_capacity, unsigned ptr_capacity, bool preserve);
    void get_values(std::map<std::pair<unsigned, unsigned>, REAL>& values) const;

    /// Gets the storage type
    StorageType get_storage_type() const { return _stype; }

    /// Gets the column (row, if CSC) indices of the nonzeros (sized get_nnz()
    unsigned* get_indices() { return _indices.get(); }

    /// Gets the row (column, if CSC) pointers
    /*
     * Element j denotes the location in the nonzero data and column (row, if
       CSC) indices that starts row (column, if CSC) j. This array has rows()+1 
       (columns()+1, if CSC) entries and the final element in the array is
       equal to get_nnz().
     */
    unsigned* get_ptr() { return _ptr.get(); }

    /// Gets the array of nonzeros (sized get_nnz())
    REAL* get_data() { return _data.get(); }

    /// Gets the number of nonzeros
    unsigned get_nnz() const { return _nnz; }

  protected:
    SPARSEMATRIXN(StorageType stype, unsigned m, unsigned n);

    boost::shared_array<unsigned> _indices;  // column (row) indices of nonzeros
    boost::shared_array<unsigned> _ptr;      // starting indices for row (col) i 
    boost::shared_array<REAL> _data;         // actual data (nonzeros)
    unsigned _nnz;                           // the number of nonzeros
    unsigned _rows;
    unsigned _columns;
    unsigned _nnz_capacity;                  // the nnz capacity
    unsigned _ptr_capacity;                  // the row capacity 
    StorageType _stype;                      // the storage capacity

  private:
    void set(unsigned rows, unsigned columns, const std::map<std::pair<unsigned, unsigned>, REAL>& values);
}; // end class

std::ostream& operator<<(std::ostream& out, const SPARSEMATRIXN& s);

