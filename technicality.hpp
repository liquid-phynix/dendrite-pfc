#if __BYTE_ORDER == __LITTLE_ENDIAN
#define NUMPY_ENDIAN "<"
#elif __BYTE_ORDER == __BIG_ENDIAN
#define NUMPY_ENDIAN ">"
#else
#error "__BYTE_ORDER needs to be defined"
#endif

template <typename T> class Dist_Matrix;

template <typename R_Elt> struct MatchRealMatrix {
  int m_rows, m_cols;
  MatchRealMatrix(int rows, int cols):
    m_rows(rows), m_cols(cols){}
};

template <typename R_Elt>
MatchRealMatrix<std::complex<R_Elt> >
match(Dist_Matrix<R_Elt>& dmat){
  return MatchRealMatrix<std::complex<R_Elt> >(dmat.global_rows(), dmat.global_cols() / 2 + 1);
}
template <typename R_Elt>
MatchRealMatrix<std::complex<R_Elt> >
match_transposed(Dist_Matrix<R_Elt>& dmat){
  return MatchRealMatrix<std::complex<R_Elt> >(dmat.global_cols() / 2 + 1, dmat.global_rows());
}
