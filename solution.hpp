
void norm(Matrix& m){
  const int rows = m.local_rows(); const int cols = m.local_cols();
  const double nfact = 1.0 / (m.global_rows() * m.global_cols());
  for(int r = 0; r < rows; r++) for(int c = 0; c < cols; c++) m.at(r, c) *= nfact;
}

void init_dummy(Matrix& m){
  const int rows = m.local_rows(), cols = m.local_cols();
  for(int r = m.get_rows_start(), rl = 0; rl < rows; r++, rl++)
    for(int c = 0; c < cols; c++)
      m.at(rl, c) = r * cols + c;
}

void init_rand(Matrix& m){
  const int rows = m.local_rows(), cols = m.local_cols();
  for(int r = m.get_rows_start(), rl = 0; rl < rows; r++, rl++)
    for(int c = 0; c < cols; c++)
      m.at(rl, c) = Par::PSI0 + Par::AMP * (rand() / (Matrix::Elt_t) RAND_MAX - 0.5);
}

void init_nucleus(Matrix& m){
  const Matrix::Elt_t sqrt3p2 = sqrt(3.0) / 2.0;
  const Matrix::Elt_t dx = Par::LX / m.global_cols();
  const Matrix::Elt_t dy = Par::LY / m.global_rows();
  const Matrix::Elt_t x0 = m.global_cols() / 2 * dx;
  const Matrix::Elt_t y0 = m.global_rows() / 2 * dy;
  const Matrix::Elt_t rad0 = Par::R0_IN_DX * dx;
  
  const int rows = m.local_rows(), cols = m.local_cols();
  for(int r = m.get_rows_start(), rl = 0; rl < rows; r++, rl++){
    const Matrix::Elt_t yrel = r * dy - y0;
    for(int c = 0; c < cols; c++){
      const Matrix::Elt_t xrel = c * dx - x0;
      const Matrix::Elt_t rr = sqrt(xrel * xrel + yrel * yrel);
      const Matrix::Elt_t envelope = 0.5 * (1.0 - tanh(0.1 * (rr - rad0)));
      m.at(rl, c) = Par::PSI0 + Par::AMP * envelope * (2.0 * cos(sqrt3p2 * xrel) * cos(0.5 * yrel) + cos(yrel)); }}
}

void nonlin(Matrix& m){
  const int rows = m.local_rows(); const int cols = m.local_cols();
  for(int r = 0; r < rows; r++)
    for(int c = 0; c < cols; c++) {
      const Matrix::Elt_t x = m.at(r, c);
      m.at(r, c) = x * x * x - x; }
}

void nonlin_square(Matrix& m){
  const int rows = m.local_rows(); const int cols = m.local_cols();
  for(int r = 0; r < rows; r++)
    for(int c = 0; c < cols; c++) {
      const Matrix::Elt_t x = m.at(r, c);
      m.at(r, c) = x * x; }
}

void solution(CMatrix& prev_step, CMatrix& nonlin, CMatrix& new_step, FFTR2C& fft){
  const int rows = prev_step.local_rows(), cols = prev_step.local_cols();
  for(int r = 0; r < rows; r++){
    const double kr = Par::DKX * fft.k_at_row(r);
    for(int c = 0; c < cols; c++){
      const double kc = Par::DKY * fft.k_at_col(c);

      const Matrix::Elt_t ik2 = - (kr * kr + kc * kc);
      const Matrix::Elt_t ik4 = ik2 * ik2;
      const Matrix::Elt_t ik6 = ik2 * ik4;
      const Matrix::Elt_t fk2 = (2.0 - Par::EPS) * ik2 + 2.0 * ik4 + ik6;
      new_step.at(r, c) =
        prev_step.at(r, c) / (1.0 - Par::DT * fk2)
        + ik2 * nonlin.at(r, c) / (1.0 / Par::DT - fk2); }}
}

void init_gauss_kernel(Matrix& gauss, FFTR2C& fft){
  const int rows = gauss.local_rows(), cols = gauss.local_cols();
  //  const double s = Par::SIGMA * Par::SIGMA;
  for(int r = 0; r < rows; r++){
    const double kr = Par::DKX * fft.k_at_row(r);
    for(int c = 0; c < cols; c++){
      const double kc = Par::DKY * fft.k_at_col(c);
      gauss.at(r, c) = exp(- 0.5 * 2.0 * Par::SIGMA * (kr * kr + kc * kc)); }}
}
