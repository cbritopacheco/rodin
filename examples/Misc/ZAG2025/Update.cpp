#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>

#include "Update.h"

namespace Rodin::Misc::ZAG2025
{
  //––– Helper roots for Q₃ and Q₁ –––
  double theta_Q1(double alpha, double beta, double mu1, double mu2, double l)
  {
    double res =
      -alpha - beta + (std::sqrt(mu1) + std::sqrt(mu2)) * std::sqrt(
          -beta * (beta - alpha) * (alpha + beta) / l);
    res /= beta - alpha;
    return res;
  }

  double theta_Q2(double alpha, double beta, double mu1, double mu2, double l)
  {
    double res =
      -2 * alpha + (std::sqrt(-mu1) + std::sqrt(-mu2)) * std::sqrt(
          alpha * (beta - alpha) * (alpha + beta) / l);
    res /= beta - alpha;
    return res;
  }

  double theta_Q3(double alpha, double beta, double mu1, double mu2, double l)
  {
    double res =
      -alpha + std::sqrt(
          -alpha * beta * (beta - alpha) * mu1 / (l + (beta - alpha) * mu2));
    res /= beta - alpha;
    return res;
  }

  //–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  // updateThetaM: exactly your “M” version, returns θ ∈ [0,1]
  //–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  double updateThetaM(const std::array<double, 6>& data, double numTol)
  {
    const double alpha = data[0];
    const double beta  = data[1];
    const double l     = data[2];
    const double M11   = data[3];
    const double M12   = data[4];
    const double M22   = data[5];

    Eigen::Matrix2d M;
    M << M11, M12,
         M12, M22;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M);
    if (es.info() != Eigen::Success)
      throw std::runtime_error("Eigen decomposition failed");

    double mu1 = es.eigenvalues()[0];
    double mu2 = es.eigenvalues()[1];

    if (std::hypot(mu1, mu2) < numTol * numTol)
      throw std::runtime_error("Null M matrix");

    // region tests
    if (mu1 > 0)
    {
      double theta1 = (alpha - beta * std::sqrt(mu1) / std::sqrt(mu2)) / (alpha - beta);
      if (theta1 <= 0)
      {
        // Q = Q3, strictly decreasing
        double Q03 = l - mu2 * (alpha - beta) - beta * mu1 * (alpha - beta) / alpha;
        double Q13 = l - mu2 * (alpha - beta) - alpha * mu1 * (alpha - beta) / beta;
        if (Q13 >= 0)
          return 0.0;
        else if (Q03 <= 0)
          return 1.0;
        else
          return theta_Q3(alpha, beta, mu1, mu2, l);
      }
      else if (theta1 >= 1)
      {
        // Q = Q1, strictly decreasing
        double Q01 =
          l - beta * (alpha - beta) * (std::sqrt(mu1) + std::sqrt(mu2)) * (
              std::sqrt(mu1) + std::sqrt(mu2)) / (beta + alpha);
        double Q11 =
          l - (std::sqrt(mu1) + std::sqrt(mu2)) * (
              std::sqrt(mu1) + std::sqrt(mu2)) * (alpha * alpha - beta * beta) / (4 * beta);
        if (Q11 >= 0)
          return 0.0;
        else if (Q01 <= 0)
          return 1.0;
        else
          return theta_Q1(alpha,beta,mu1,mu2,l);
      }
      else
      {
        // Mixed: choose between Q1 and Q3
        double Q01 =
          l - beta * (alpha - beta) * (std::sqrt(mu1) + std::sqrt(mu2)) * (
              std::sqrt(mu1) + std::sqrt(mu2)) / (beta + alpha);
        double Q13 = l - mu2 * (alpha - beta) - alpha * mu1 * (alpha - beta) / beta;
        double Qtheta = l - mu2 * (alpha * alpha - beta * beta) / beta;
        if (Q13 >= 0)
          return 0.0;
        else if (Q01 <= 0)
          return 1.0;
        else if (Qtheta >= 0)
          return theta_Q3(alpha,beta,mu1,mu2,l);
        else
          return theta_Q1(alpha,beta,mu1,mu2,l);
      }
    }
    else if (mu2 < 0)
    {
      double theta2 = alpha * (std::sqrt(-mu1) / std::sqrt(-mu2) - 1) / (beta - alpha);
      if (theta2 <= 0)
      {
        // Q = Q2, strictly increasing
        double Q02 = l - (beta * beta - alpha * alpha)
                         *(std::sqrt(-mu1)+std::sqrt(-mu2))
                         *(std::sqrt(-mu1)+std::sqrt(-mu2)) / (4*alpha);
        double Q12 = l + alpha*(alpha - beta)
                         *(std::sqrt(-mu1)+std::sqrt(-mu2))
                         *(std::sqrt(-mu1)+std::sqrt(-mu2)) / (alpha+beta);
        if (Q02 >= 0)      return 0.0;
        else if (Q12 <= 0) return 1.0;
        else               return theta_Q2(alpha,beta,mu1,mu2,l);
      }
      else if (theta2 >= 1)
      {
        // Q = Q3, strictly increasing
        double Q03 = l - mu2*(alpha - beta) - beta*mu1*(alpha - beta)/alpha;
        double Q13 = l - mu2*(alpha - beta) - alpha*mu1*(alpha - beta)/beta;
        if (Q03 >= 0)      return 0.0;
        else if (Q13 <= 0) return 1.0;
        else               return theta_Q3(alpha,beta,mu1,mu2,l);
      }
      else
      {
        // Mixed Q3/Q2
        double Q03    = l - mu2*(alpha - beta) - beta*mu1*(alpha - beta)/alpha;
        double Q12    = l + alpha*(alpha - beta)
                           *(std::sqrt(-mu1)+std::sqrt(-mu2))
                           *(std::sqrt(-mu1)+std::sqrt(-mu2)) / (alpha+beta);
        double Qtheta = (-mu2*alpha*alpha + l*alpha + mu2*beta*beta)/alpha;
        if (Q03 >= 0)      return 0.0;
        else if (Q12 <= 0) return 1.0;
        else if (Qtheta < 0) return theta_Q2(alpha,beta,mu1,mu2,l);
        else                 return theta_Q3(alpha,beta,mu1,mu2,l);
      }
    }
    else
    {
      // μ₁≤0≤μ₂ ⇒ pure Q3 (increasing)
      double Q03 = l - mu2*(alpha - beta) - beta*mu1*(alpha - beta)/alpha;
      double Q13 = l - mu2*(alpha - beta) - alpha*mu1*(alpha - beta)/beta;
      if (Q03 >= 0)      return 0.0;
      else if (Q13 <= 0) return 1.0;
      else               return theta_Q3(alpha,beta,mu1,mu2,l);
    }
  }

  //–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  // updateConductivityM: returns {θ, a11, a12, a22, opt}
  //–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  std::array<double,5> updateConductivityM(const std::array<double,6>& data)
  {
    const double alpha = data[0];
    const double beta  = data[1];

    // compute theta and eigen­stuff
    double theta = updateThetaM(data);
    Eigen::Matrix2d M;
    M << data[3], data[4], data[4], data[5];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M);
    if (es.info() != Eigen::Success)
      throw std::runtime_error("Eigen decomposition failed");
    double mu1 = es.eigenvalues()[0], mu2 = es.eigenvalues()[1];
    auto Q = es.eigenvectors();

    int opt = 1;

    std::array<double,5> R;
    R[0] = theta;
    if (opt == 1)
    {
      // pure phases
      R[1] = theta*alpha + (1 - theta) * beta;
      R[2] = 0.0;
      R[3] = R[1];
    }
    else
    {
      double lam1, lam2;
      if (opt == 11)
      {
        lam1 = 1.0/(theta/alpha - (theta-1)/beta);
        lam2 = alpha*theta - beta*(theta-1);
      }
      else if (opt == 12)
      {
        lam1 = alpha + (std::sqrt(-mu1)+std::sqrt(-mu2))/std::sqrt(-mu1)
               * (alpha*(beta-alpha)*(1-theta))/(theta*(beta-alpha)+2*alpha);
        lam2 = alpha + (std::sqrt(-mu1)+std::sqrt(-mu2))/std::sqrt(-mu2)
               * (alpha*(beta-alpha)*(1-theta))/(theta*(beta-alpha)+2*alpha);
      }
      else // opt == 22
      {
        lam1 = beta - (std::sqrt(mu1)+std::sqrt(mu2))/std::sqrt(mu1)
               * (beta*(beta-alpha)*theta)/(theta*(beta-alpha)+alpha+beta);
        lam2 = beta - (std::sqrt(mu1)+std::sqrt(mu2))/std::sqrt(mu2)
               * (beta*(beta-alpha)*theta)/(theta*(beta-alpha)+alpha+beta);
      }
      // A = Q * diag(lam1,lam2) * Qᵀ
      R[1] = lam1*Q(0,0)*Q(0,0) + lam2*Q(0,1)*Q(0,1);
      R[2] = lam1*Q(0,0)*Q(1,0) + lam2*Q(0,1)*Q(1,1);
      R[3] = lam1*Q(1,0)*Q(1,0) + lam2*Q(1,1)*Q(1,1);
    }
    R[4] = opt;
    return R;
  }
}

