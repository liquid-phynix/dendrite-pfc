
template <typename Elt>
struct Params_t {
  static constexpr Elt LX = 4000.0;
  static constexpr Elt LY = 4000.0;
  static constexpr Elt DT = 1e-1;
  static constexpr Elt EPS = 0.75;
  static constexpr Elt PSI0 = -0.503; // -0.503355;
  static constexpr Elt AMP = 0.5;
  static constexpr Elt R0_IN_DX = 20;
  static constexpr int SAVEFREQ = 100;
  static constexpr Elt THS = 0.5;

  static constexpr Elt DKX = 2.0 * M_PI / LX;
  static constexpr Elt DKY = 2.0 * M_PI / LY;
  static constexpr Elt SIGMA = 4.0 * M_PI / sqrt(3.0);
};

/*
namespace Par {
  const double LX = 1000.0;
  const double LY = 1000.0;
  const double DT = 1e-2;
  const double EPS = 0.75;
  const double PSI0 = -0.503; // -0.503355;
  const double AMP = 0.5;
  const double R0_IN_DX = 50;
  const int SAVEFREQ = 1000;
  const double THS = 0.5;

  const double FKX = 2.0 * M_PI / LX;
  const double FKY = 2.0 * M_PI / LY;
};
*/
