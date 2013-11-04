template <typename Elt> class Global_Transposition {
  Dist_Matrix<Elt>& m_1st;
  Dist_Matrix<Elt>& m_2nd;
  int* sendcounts; int* recvcounts; int* sdispls; int* rdispls;
public:
  Global_Transposition(Dist_Matrix<Elt>& _1st, Dist_Matrix<Elt>& _2nd):
    m_1st(_1st),
    m_2nd(_2nd){
    { // sanity checks
      assert(&m_1st != &m_2nd and "global transposition is out-of-place, cannot accept the same array for both input AND output");
      assert(m_1st.partitions() == m_2nd.partitions() and "arrays must be split to an equal number of parts");
      assert(m_1st.global_rows() == m_2nd.global_cols() and m_2nd.global_rows() == m_1st.global_cols() and "array dimensions must agree when interchanged"); }
    int nodes = m_1st.partitions();
    sendcounts = new int[nodes](); recvcounts = new int[nodes](); sdispls = new int[nodes](); rdispls = new int [nodes]();
    for(int rank = 0; rank < nodes; rank++){
      sendcounts[rank] = m_1st.local_rows() * m_2nd.get_divs_of_rank(rank);
      recvcounts[rank] = m_1st.get_divs_of_rank(rank) * m_2nd.local_rows(); }
    for(int rank = 1; rank < nodes; rank++){
      sdispls[rank] = sdispls[rank - 1] + sendcounts[rank - 1];
      rdispls[rank] = rdispls[rank - 1] + recvcounts[rank - 1]; }
  }
  ~Global_Transposition(){ delete[] sendcounts; delete[] recvcounts; delete[] sdispls; delete[] rdispls; }
  void exec(){
    Elt* __restrict__ arr_1st = m_1st.base_address();
    Elt* __restrict__ arr_2nd = m_2nd.base_address();
    
    { // prepare for mpi transposition, so that slowest dimension is indexed by rank
      int idx = 0; const int nodes = m_1st.partitions();
      const int local_rows = m_1st.local_rows(); const int local_cols = m_1st.local_cols();
      for(int rank = 0, width, c_offset = 0; rank < nodes; rank++, c_offset += width){
        width = m_2nd.get_divs_of_rank(rank);
        for(int r = 0; r < local_rows; r++){
          for(int c = 0; c < width; c++){
            arr_2nd[idx++] = arr_1st[r * local_cols + c_offset + c]; }}}}
    MPI::COMM_WORLD.Alltoallv(arr_2nd, sendcounts, sdispls, m_2nd.get_mpi_datatype(),
                     arr_1st, recvcounts, rdispls, m_1st.get_mpi_datatype());
    { // local transposition
      const int rows = m_2nd.global_cols(); const int cols = m_2nd.local_rows();
      for(int r = 0; r < rows; r++)
        for(int c = 0; c < cols; c++)
          arr_2nd[c * rows + r] = arr_1st[r * cols + c]; }
  }
};
