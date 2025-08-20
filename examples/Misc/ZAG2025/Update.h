#ifndef RODIN_EXAMPLES_MISC_ZAG2025_UPDATE_H
#define RODIN_EXAMPLES_MISC_ZAG2025_UPDATE_H

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>

namespace Rodin::Misc::ZAG2025
{
  //––– Helper roots for Q₃ and Q₁ –––
  double theta_Q3(double alpha, double beta, double mu1, double mu2, double l);

  double theta_Q1(double alpha, double beta, double mu1, double mu2, double l);

  double theta_Q2(double alpha, double beta, double mu1, double mu2, double l);

  double updateThetaM(const std::array<double, 6>& data, double numTol = 1e-10);

  double updateThetaN(const std::array<double,6>& data, double numTol);

  // updateConductivityM: returns {θ, a11, a12, a22, opt}
  std::array<double, 5> updateConductivityM(const std::array<double, 6>& data);
}

#endif
