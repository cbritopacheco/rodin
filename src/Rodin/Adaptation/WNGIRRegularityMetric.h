/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRREGULARITYMETRIC_H
#define RODIN_ADAPTATION_WNGIRREGULARITYMETRIC_H

#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>
#include <vector>
#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "CellDeformation.h"

namespace Rodin::Adaptation::Detail
{
  /// Frozen inverse deformation gradient, expressed as a form-language coefficient.
  template <class Displacement>
  class WNGIRCurrentInverse final
    : public Variational::MatrixFunctionBase<Real, WNGIRCurrentInverse<Displacement>>
  {
    public:
      WNGIRCurrentInverse(const Displacement& current, size_t dimension)
        : m_gradient(Variational::Jacobian(current)),
          m_dimension(dimension)
      {}
      template <class Point>
      Math::SpatialMatrix<Real> getValue(const Point& point) const
      {
        CellDeformation deformation(m_dimension);
        deformation.setDisplacementGradient(m_gradient.getValue(point));
        return Math::SpatialMatrix<Real>(deformation.getInverseTranspose().transpose());
      }
      size_t getRows() const
      {
        return m_dimension;
      }
      size_t getColumns() const
      {
        return m_dimension;
      }
      Optional<size_t> getOrder(const Geometry::Polytope&) const
      {
        return std::nullopt;
      }
      WNGIRCurrentInverse* copy() const noexcept override
      {
        return new WNGIRCurrentInverse(*this);
      }

    private:
      decltype(Variational::Jacobian(std::declval<const Displacement&>())) m_gradient;
      size_t m_dimension;
  };

  template <class Displacement>
  auto wngirCurrentVolumeWeight(const Displacement& current, size_t dimension)
  {
    return Variational::RealFunction(
      [gradient = Variational::Jacobian(current), dimension](const auto& point) -> Real {
        CellDeformation deformation(dimension);
        deformation.setDisplacementGradient(gradient.getValue(point));
        return deformation.getJacobian();
      });
  }

  /// Pullback of sym(grad_y v), with y=x+u_k(x).
  template <class Function, class Displacement>
  auto wngirCurrentStrain(
    const Function& function, const Displacement& current, size_t dimension)
  {
    const auto gradient =
      Variational::Jacobian(function) * WNGIRCurrentInverse(current, dimension);
    return Real(0.5) * (gradient + Variational::Transpose(gradient));
  }

  /// K_strain = integral j*eps(v):eps(z) - b(v)b(z)/(d*integral j).
  /// Only uniform current-configuration dilation is projected out.
  template <class TestFunction, class Displacement>
  std::vector<Math::Vector<Real>> wngirCurrentStrainCouplings(const TestFunction& test,
    const Displacement& current, size_t dimension, Real coefficient, size_t order)
  {
    if (coefficient == Real(0))
      return {};
    if (!(coefficient > Real(0)))
      Alert::Exception() << "Current strain regularity requires positive coefficient."
                         << Alert::Raise;
    const auto weight = wngirCurrentVolumeWeight(current, dimension);
    auto integral = Variational::Integral((Real(1) / std::sqrt(Real(dimension))) *
      weight * Variational::Trace(wngirCurrentStrain(test, current, dimension)));
    integral.setOrder(order);
    Variational::LinearForm form(test);
    form = integral;
    form.assemble();
    Variational::GridFunction dilation(test.getFiniteElementSpace());
    dilation = Variational::VectorFunction(dimension, [](const Geometry::Point& point) {
      return Math::SpatialVector<Real>(point.getCoordinates());
    });
    dilation += current;
    dilation *= Real(1) / std::sqrt(Real(dimension));
    const Real measure = form.getVector().dot(dilation.getData());
    if (!(measure > Real(0)) || !std::isfinite(measure))
      Alert::Exception()
        << "Current strain regularity requires positive finite current volume."
        << Alert::Raise;
    return {Math::Vector<Real>(std::sqrt(coefficient / measure) * form.getVector())};
  }

  /**
   * @brief Audit inertia of A + sum weights[k] modes[k] modes[k]^T, without repair.
   * Uses inertia(A)+inertia(-D^-1-U^T A^-1 U)-inertia(-D^-1).
   * Returns the negative eigenvalue count, or -1 if the unpivoted base
   * factorization or small Schur block cannot resolve inertia reliably.
   * No positive-definiteness claim follows from an unresolved check.
   */
  inline Integer wngirRegularityInertia(const Math::SparseMatrix<Real>& A,
    const std::vector<Math::Vector<Real>>& modes, const std::vector<Real>& weights)
  {
    Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> factor;
    factor.compute(A);
    if (factor.info() != Eigen::Success || !factor.vectorD().allFinite())
      return -1;
    const auto pivots = factor.vectorD();
    const Real margin =
      Real(64) * std::numeric_limits<Real>::epsilon() * pivots.cwiseAbs().maxCoeff();
    if (pivots.cwiseAbs().minCoeff() <= margin)
      return -1;
    Integer negative = (pivots.array() < Real(0)).count();
    if (weights.empty())
      return negative;
    const auto rank = static_cast<Eigen::Index>(weights.size());
    Math::Matrix<Real> U(A.rows(), rank), diagonal = Math::Matrix<Real>::Zero(rank, rank);
    Integer auxiliaryNegative = 0;
    for (Eigen::Index k = 0; k < rank; ++k)
    {
      U.col(k) = std::sqrt(std::abs(weights[k])) * modes[k];
      diagonal(k, k) = weights[k] > Real(0) ? Real(-1) : Real(1);
      auxiliaryNegative += weights[k] > Real(0);
    }
    const Math::Matrix<Real> inverseU = factor.solve(U);
    if (!inverseU.allFinite() || (A * inverseU - U).norm() > Real(1e-7) * U.norm())
      return -1;
    Math::Matrix<Real> schur = diagonal - U.transpose() * inverseU;
    schur = (Real(0.5) * (schur + schur.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(schur);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite() ||
      eigen.eigenvalues().cwiseAbs().minCoeff() <=
        Real(1e-9) * std::max(Real(1), eigen.eigenvalues().cwiseAbs().maxCoeff()))
      return -1;
    negative += (eigen.eigenvalues().array() < Real(0)).count() - auxiliaryNegative;
    return negative >= 0 ? negative : -1;
  }
}
#endif
