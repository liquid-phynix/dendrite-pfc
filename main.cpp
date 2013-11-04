#include <fstream>
#include <iostream>
#include <cassert>
#include <complex>
#include <mpi.h>
#include <cstdio>

#define MASTER if(MPI::COMM_WORLD.Get_rank() == 0)
#define LAST if(MPI::COMM_WORLD.Get_rank() == MPI::COMM_WORLD.Get_size() - 1)

#include "timing.hpp"
#include "technicality.hpp"
#include "distributed_matrix.hpp"
#include "global_transposition.hpp"
#include "distributed_fft.hpp"
#include "params.hpp"
#include "contour.hpp"

#ifndef FLOAT_TYPE
#define FLOAT_TYPE                 double
#endif

typedef FLOAT_TYPE                 Real;
typedef std::complex<Real>         Complex;
typedef Dist_Matrix<Real>          Matrix;
typedef Dist_Matrix<Complex>       CMatrix;
typedef Dist_2D_R2C_FFT<Real>      FFTR2C;
typedef Dist_2D_C2R_FFT<Real>      FFTC2R;
typedef Params_t<Real>             Par;

#include "solution.hpp"

TimeAcc timer;

int main(int argc, char* argv[]){
  MPI::Init(argc, argv);
  MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);

  // if(argc != 4){
  //   std::cout << "usage: " << argv[0] << " <rows> <cols> <iterations>" << std::endl;
  //   MPI::Finalize();
  //   return 0; }

  // const int rows = atoi(argv[1]);
  // const int cols = atoi(argv[2]);
  // const int iters = atoi(argv[3]);
  const int iters = 1;

  // for PFC
  //  Matrix r_nonlin(rows, cols);
  Matrix r_nonlin("dendrite.mat.bin");

  CMatrix c_fft_nonlin(match_transposed(r_nonlin));
  CMatrix c_fft_prev(match_transposed(r_nonlin));
  CMatrix c_fft_new(match_transposed(r_nonlin));
  CMatrix c_tmp(match(r_nonlin));
  //         real      comp          real_l comp_l
  FFTR2C r2c(r_nonlin, c_fft_nonlin, c_tmp, c_fft_new);
  //         comp       real      comp_l        real_l
  FFTC2R c2r(c_fft_new, r_nonlin, c_fft_nonlin, c_tmp);

  // for dendrite boundary computation
  Matrix r_square(r_nonlin.global_rows(), r_nonlin.global_cols());
  Matrix r_gauss(match_transposed(r_square));
  FFTR2C gauss_r2c(r_square, c_fft_nonlin, c_tmp, c_fft_new);
  FFTC2R gauss_c2r(c_fft_nonlin, r_square, c_fft_new, c_tmp);
  init_gauss_kernel(r_gauss, gauss_r2c);

  r_gauss.save_to_disk("gauss.mat.bin", 0);

  //  init_dummy(r_nonlin);
  //  r_nonlin.save_to_disk("start.mat.bin", 0);
  //  contour(r_nonlin, 5000);
  //  init_nucleus(r_nonlin);
  //  contour(r_nonlin, -0.4, "%s_%06d.txt", 0);

  r_nonlin.save_to_disk("start.mat.bin", 0);
  r2c.exec();
  c_fft_prev.overwrite_with(c_fft_nonlin);

  for(int it = 0; it <= iters; it++){ // read-only: c_fft_prev, r_nonlin
    // saving
//    if(it % Par::SAVEFREQ == 0) r_nonlin.save_to_disk("out_%06d.mat.bin", it);
    // computing dendrite boundary
    if(it % Par::SAVEFREQ == 0){
      c_fft_nonlin.overwrite_with(c_fft_prev);
      // (1 - G) * G spektrumban
      const int rows = r_gauss.local_rows(); const int cols = r_gauss.local_cols();
      for(int r = 0; r < rows; r++) for(int c = 0; c < cols; c++) c_fft_nonlin.at(r, c) *= 1 - r_gauss.at(r, c);
      gauss_c2r.exec();
      norm(r_square);
      nonlin_square(r_square);
      gauss_r2c.exec();
      for(int r = 0; r < rows; r++) for(int c = 0; c < cols; c++) c_fft_nonlin.at(r, c) *= r_gauss.at(r, c);
      gauss_c2r.exec();
      norm(r_square);

      r_square.save_to_disk("r_square.mat.bin", it);

      contour(r_square, r_square.max() / 2, "%s_%06d.txt", it);
    }

    timer.start();
    nonlin(r_nonlin);
    r2c.exec();
    solution(c_fft_prev, c_fft_nonlin, c_fft_new, r2c);
    c_fft_prev.overwrite_with(c_fft_new);
    c2r.exec();
    norm(r_nonlin);
    timer.stop();
  }

  r_nonlin.save_to_disk("end.mat.bin", 0);
  //  contour(r_nonlin, -0.508, "%s_%06d.txt", 0);


  MASTER timer.report("problem took %d ms per cycle\n");
  
  MPI::Finalize();
  return 0;
}
