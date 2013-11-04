#include <fftw3.h>

template <typename R_Elt>
class Dist_2D_FFT {
protected:
  Global_Transposition<std::complex<R_Elt> >* m_tr;
  fftw_plan m_plan_1, m_plan_2;
  void sanity(int Rrows, int Rcols, int Crows, int Ccols){
    assert(Rrows == Crows and Ccols == (Rcols / 2 + 1)
           and "array dimensions dont work out for distributed 2D R2C transform"); }
public:
  Dist_2D_FFT(): m_tr(0){};
  ~Dist_2D_FFT(){
    fftw_destroy_plan(m_plan_1);
    fftw_destroy_plan(m_plan_2);
    delete m_tr; }
  void exec(){ fftw_execute(m_plan_1); m_tr->exec(); fftw_execute(m_plan_2); }
};

template <typename R_Elt>
class Dist_2D_R2C_FFT: public Dist_2D_FFT<R_Elt> {
  R_Elt* m_kvec_row; R_Elt* m_kvec_col;
public:
  Dist_2D_R2C_FFT(Dist_Matrix<R_Elt>&                real,
                  Dist_Matrix<std::complex<R_Elt> >& comp,
                  Dist_Matrix<std::complex<R_Elt> >& real_like,
                  Dist_Matrix<std::complex<R_Elt> >& comp_like,
                  unsigned int fftw_flags = FFTW_ESTIMATE):
    m_kvec_row(0), m_kvec_col(0){
    //    std::cout << "in R2C\nreal: " << real << "\ncomp: " << comp << "\nreal_like: " << real_like << "\ncomp_like: " << comp_like << std::endl;
    this->sanity(real.global_rows(), real.global_cols(), comp.global_cols(),      comp.global_rows());
    this->sanity(real.global_rows(), real.global_cols(), comp_like.global_cols(), comp_like.global_rows());
    this->sanity(real.global_rows(), real.global_cols(), real_like.global_rows(), real_like.global_cols());

    const int length_r = real.local_cols();
    this->m_plan_1 = fftw_plan_many_dft_r2c(1, &length_r, real.local_rows(),
                                            reinterpret_cast<R_Elt*>(real.base_address()),             NULL, 1, length_r,
                                            reinterpret_cast<fftw_complex*>(real_like.base_address()), NULL, 1, length_r / 2 + 1,
                                            fftw_flags);
    this->m_tr = new Global_Transposition<std::complex<R_Elt> >(real_like, comp_like);
    const int length_c = comp.local_cols();
    this->m_plan_2 = fftw_plan_many_dft(1, &length_c, comp.local_rows(),
                                        reinterpret_cast<fftw_complex*>(comp_like.base_address()), NULL, 1, length_c,
                                        reinterpret_cast<fftw_complex*>(comp.base_address()),      NULL, 1, length_c,
                                        FFTW_FORWARD, fftw_flags);
    m_kvec_row = new R_Elt[comp.local_rows()]();
    m_kvec_col = new R_Elt[comp.local_cols()]();
    for(int r = 0; r < comp.local_rows(); r++) m_kvec_row[r] = comp.get_rows_start() + r;
    { int i = 0;
      for(int j = 0; j < comp.local_cols() / 2 + 1; j++, i++) m_kvec_col[i] = j;
      for(int j = - (comp.local_cols() - comp.local_cols() / 2 - 1); j < 0; j++, i++) m_kvec_col[i] = j; }
  }
  ~Dist_2D_R2C_FFT(){ delete m_kvec_row; delete m_kvec_col; }
  R_Elt k_at_row(int r){ return m_kvec_row[r]; }
  R_Elt k_at_col(int c){ return m_kvec_col[c]; }
};

template <typename R_Elt>
class Dist_2D_C2R_FFT: public Dist_2D_FFT<R_Elt> {
public:
  Dist_2D_C2R_FFT(Dist_Matrix<std::complex<R_Elt> >& comp,
                  Dist_Matrix<R_Elt>&                real,
                  Dist_Matrix<std::complex<R_Elt> >& comp_like,
                  Dist_Matrix<std::complex<R_Elt> >& real_like,
                  unsigned int fftw_flags = FFTW_ESTIMATE){
    this->sanity(real.global_rows(), real.global_cols(), comp.global_cols(),      comp.global_rows());
    this->sanity(real.global_rows(), real.global_cols(), comp_like.global_cols(), comp_like.global_rows());
    this->sanity(real.global_rows(), real.global_cols(), real_like.global_rows(), real_like.global_cols());
    
    const int length_c = comp.local_cols();
    this->m_plan_1 = fftw_plan_many_dft(1, &length_c, comp.local_rows(),
                                        reinterpret_cast<fftw_complex*>(comp.base_address()),      NULL, 1, length_c,
                                        reinterpret_cast<fftw_complex*>(comp_like.base_address()), NULL, 1, length_c,
                                        FFTW_BACKWARD, fftw_flags);
    this->m_tr = new Global_Transposition<std::complex<R_Elt> >(comp_like, real_like);
    const int length_r = real.local_cols();
    this->m_plan_2 = fftw_plan_many_dft_c2r(1, &length_r, real.local_rows(),
                                            reinterpret_cast<fftw_complex*>(real_like.base_address()), NULL, 1, length_r / 2 + 1,
                                            reinterpret_cast<R_Elt*>(real.base_address()),             NULL, 1, length_r,
                                            fftw_flags);
  }
};
