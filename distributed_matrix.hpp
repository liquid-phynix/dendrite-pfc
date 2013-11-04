#include <string>

template <typename T> struct RealWitness;
template <> struct RealWitness<float>{};
template <> struct RealWitness<double>{};

template <typename Elt> class Dist_Matrix {
private:
  Elt* m_array;
  int* m_rows_all;
  int m_rows, m_rows_offset;
  const int m_cols, m_rank, m_partitions, m_global_rows;
  const char* dtype_to_string();
  int m_read_offset;
public:
  typedef Elt Elt_t;
  std::pair<std::pair<int, int>, int> parse_dimensions(const char* fn){
    int buf[3] = {0, 0, 0};
    MASTER { 
      char typestring[10];
      char firstline[256];
      std::fstream fs(fn,  std::fstream::in);
      fs.getline(firstline, 256);
      buf[2] = fs.tellg();
      if(firstline[0] != '#') abort();
      std::sscanf(firstline, "# rows=%d cols=%d dtype=<%s", &buf[0], &buf[1], typestring);
      if(strcmp(typestring, "f8") != 0){
        std::cerr << "problem with input matrix format" << std::endl;
        abort(); }}
    MPI::COMM_WORLD.Bcast(buf, 3, MPI::INT, 0);
    //    std::cout << "rank " << MPI::COMM_WORLD.Get_rank() << " rows=" << buf[0] << ", cols=" << buf[1] << std::endl;
    return std::make_pair(std::make_pair(buf[0], buf[1]), buf[2]);
  }
  Dist_Matrix(const char* fn):
    Dist_Matrix(parse_dimensions(fn)){
    MPI::COMM_WORLD.Barrier();
    MPI::File fh = MPI::File::Open(MPI::COMM_WORLD, fn, MPI::MODE_RDONLY, MPI::INFO_NULL);
    fh.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    fh.Set_view(m_read_offset + m_rows_offset * local_cols() * sizeof(Elt), get_mpi_datatype(), get_mpi_datatype(), "native", MPI::INFO_NULL);
    
    //    std::cout << "rank " << MPI::COMM_WORLD.Get_rank() << " reads " << local_rows() * local_cols() << std::endl;
    fh.Read_all(m_array, local_rows() * local_cols(), get_mpi_datatype());
    fh.Close();
    MASTER std::cout << "matrix initialized from file '" << fn << "'" << std::endl;
    MPI::COMM_WORLD.Barrier();
  }
  Dist_Matrix(std::pair<std::pair<int, int>, int> dims):
    Dist_Matrix(dims.first.first, dims.first.second){ m_read_offset = dims.second; }
  Dist_Matrix(int g_rows, int g_cols):
    Dist_Matrix(g_rows, g_cols, MPI::COMM_WORLD.Get_size(), MPI::COMM_WORLD.Get_rank()){}
  template <typename Other_Elt> Dist_Matrix(MatchRealMatrix<Other_Elt> m):
    Dist_Matrix(m.m_rows, m.m_cols, MPI::COMM_WORLD.Get_size(), MPI::COMM_WORLD.Get_rank()){}
  Dist_Matrix(int g_rows, int g_cols, int partitions, int rank):
    m_array(0), m_rows_all(0), m_cols(g_cols), m_rank(rank), m_partitions(partitions), m_global_rows(g_rows){
    { // sanity checks
      assert(g_cols > 0 and g_rows > 0 and partitions > 0 and "number of {rows, columns, partitions} must be > 0");
      assert(g_rows >= partitions and "number of rows must be >= number of partitions");
      assert(rank >= 0 and rank < partitions and "rank must be >= 0 and < partitions"); }
    m_rows_all = new int[m_partitions]();
    int rows = m_global_rows / m_partitions; int rem = m_global_rows - rows * partitions;
    m_rows_offset = 0;
    for(int r = 0; r < m_partitions; r++){
      m_rows_all[r] = rows + (rem-- > 0 ? 1 : 0);
      m_rows_offset += (r > 0 and r <= m_rank) ? m_rows_all[r - 1] : 0; }
    m_rows = m_rows_all[m_rank];
    m_array = new Elt[std::max((g_cols / partitions + (g_cols % partitions != 0)) * g_rows,
                               (g_rows / partitions + (g_rows % partitions != 0)) * g_cols)]();
  }
  ~Dist_Matrix(){ delete[] m_array; delete[] m_rows_all; }
  void overwrite_with(Dist_Matrix<Elt>& m){
    assert(local_rows() == m.local_rows() and local_cols() == m.local_cols() and "matrix dimensions dont match");
    memcpy(m_array, m.m_array, local_rows() * local_cols() * sizeof(Elt)); }
  const MPI::Datatype& get_mpi_datatype();
  Elt& at(int r, int c){ return m_array[r * m_cols + c]; }
  int partitions(){ return m_partitions; }
  int local_rows(){ return m_rows; }
  int local_cols(){ return m_cols; }
  int global_rows(){ return m_global_rows; }
  int global_cols(){ return m_cols; }
  Elt* base_address(){ return m_array; }
  int get_divs_of_rank(int r){ return m_rows_all[r]; }
  int get_rows_start(){ return m_rows_offset; }
  void save_to_disk(const char* fmt, int it){
    char fn[128]; sprintf(fn, fmt, it);    
    MASTER { 
      std::fstream fs(fn,  std::fstream::out | std::fstream::trunc);
      fs << "# rows=" << global_rows()
         << " cols=" << global_cols()
         << " dtype=" << NUMPY_ENDIAN << dtype_to_string() << std::endl; }
    MPI::COMM_WORLD.Barrier();
    MPI::File fh = MPI::File::Open(MPI::COMM_WORLD, fn, MPI::MODE_WRONLY | MPI::MODE_APPEND, MPI::INFO_NULL);
    fh.Set_errhandler(MPI::ERRORS_ARE_FATAL);
    fh.Set_view(fh.Get_position() + m_rows_offset * local_cols() * sizeof(Elt), get_mpi_datatype(), get_mpi_datatype(), "native", MPI::INFO_NULL);
    fh.Write_all(m_array, local_rows() * local_cols(), get_mpi_datatype());
    fh.Close();
    if(m_rank == 0) std::cout << "matrix saved to '" << fn << "'" << std::endl;
  }
  Elt max(){
    typedef RealWitness<Elt> _;
    Elt local_max = at(0, 0);
    int max_i = local_rows() * local_cols();
    for(int i = 0; i < max_i; i++) local_max = local_max > m_array[i] ? local_max : m_array[i];
    Elt sendbuf[1] = { local_max };
    Elt recvbuf[1];
    MPI::COMM_WORLD.Allreduce(sendbuf, recvbuf, 1, get_mpi_datatype(), MPI::MAX);
    return recvbuf[0];
  }
  friend std::ostream& operator<<(std::ostream& o, Dist_Matrix<Elt>& m){
    o << "<rank=" << m.m_rank
      << ", grows=" << m.global_rows()
      << ", gcols=" << m.global_cols()
      << ", rows=" << m.local_rows()
      << ", cols=" << m.local_cols() << ">" << std::endl;
    return o;
  }
};

template <> const char* Dist_Matrix<float>::dtype_to_string(){ return "f4"; }
template <> const char* Dist_Matrix<double>::dtype_to_string(){ return "f8"; }
template <> const char* Dist_Matrix<std::complex<float> >::dtype_to_string(){ return "c8"; }
template <> const char* Dist_Matrix<std::complex<double> >::dtype_to_string(){ return "c16"; }

template <> const MPI::Datatype& Dist_Matrix<double>::get_mpi_datatype() { return MPI::DOUBLE; }
template <> const MPI::Datatype& Dist_Matrix<float>::get_mpi_datatype(){ return MPI::FLOAT; }
template <> const MPI::Datatype& Dist_Matrix<std::complex<double> >::get_mpi_datatype(){ return MPI::DOUBLE_COMPLEX; }
template <> const MPI::Datatype& Dist_Matrix<std::complex<float> >::get_mpi_datatype(){ return MPI::COMPLEX; }
