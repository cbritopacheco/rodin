/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRACTIVESETKKT_H
#define RODIN_ADAPTATION_WNGIRACTIVESETKKT_H

#include <Eigen/SparseCholesky>

#include "CellDeformation.h"
#include "WNGIRParameters.h"

namespace Rodin::Adaptation::Detail
{
  /**
   * @brief Eigen reference solver for the linearized admissibility QP.
   *
   * This class is a comparison instrument, not the production WNGIR path. It
   * solves the active equality-constrained systems through the Schur complement
   * of the Rodin-assembled SPD metric. Its cost grows with the active validation
   * points and therefore records the expense avoided by the metric models.
   */
  template <class TrialFunction, class Displacement>
  class WNGIRActiveSetKKT
  {
      struct Constraint
      {
          std::vector<std::pair<Index, Real>> entries;
          Real bound;
      };

    public:
      /// @brief Constructs the w n g i r active set k k t.
      WNGIRActiveSetKKT(const TrialFunction& trial, const Displacement& current,
        const WNGIRParameters& parameters)
        : m_trial(trial),
          m_current(current),
          m_parameters(parameters)
      {}

      template <class Mesh>
      /// @brief Solves the local constrained system.
      bool solve(const Mesh& mesh, const Math::SparseMatrix<Real>& metric,
        const Math::Vector<Real>& force, Math::Vector<Real>& direction,
        std::size_t& activeCount) const
      {
        Eigen::SimplicialLDLT<Math::SparseMatrix<Real>> factor;
        factor.compute(metric);
        if (factor.info() != Eigen::Success)
          return false;
        const Math::Vector<Real> unconstrained = factor.solve(force);
        if (factor.info() != Eigen::Success)
          return false;

        const auto constraints = getConstraints(mesh);
        direction = unconstrained;
        std::vector<std::size_t> active;
        std::vector<Math::Vector<Real>> responses;
        Math::Vector<Real> multipliers;
        constexpr Real tolerance = Real(1e-10);
        const std::size_t maxActive = std::min<std::size_t>(constraints.size(), 256);

        for (std::size_t iteration = 0; iteration < 2 * maxActive + 1; ++iteration)
        {
          std::size_t violated = constraints.size();
          Real maxViolation = tolerance;
          for (std::size_t row = 0; row < constraints.size(); ++row)
          {
            if (std::find(active.begin(), active.end(), row) != active.end())
              continue;
            const Real violation =
              dot(constraints[row], direction) - constraints[row].bound;
            if (violation > maxViolation)
            {
              maxViolation = violation;
              violated = row;
            }
          }

          if (violated != constraints.size())
          {
            if (active.size() == maxActive)
              return false;
            active.push_back(violated);
          }
          else if (active.empty() || multipliers.minCoeff() >= -tolerance)
          {
            activeCount = active.size();
            return true;
          }
          else
          {
            Eigen::Index remove = 0;
            multipliers.minCoeff(&remove);
            active.erase(active.begin() + remove);
          }

          const Eigen::Index count = static_cast<Eigen::Index>(active.size());
          responses.clear();
          responses.reserve(active.size());
          for (const auto row : active)
          {
            Math::Vector<Real> rhs = Math::Vector<Real>::Zero(force.size());
            for (const auto& [column, value] : constraints[row].entries)
              rhs(static_cast<Eigen::Index>(column)) = value;
            responses.push_back(factor.solve(rhs));
            if (factor.info() != Eigen::Success)
              return false;
          }

          Math::Matrix<Real> schur(count, count);
          Math::Vector<Real> rhs(count);
          for (Eigen::Index i = 0; i < count; ++i)
          {
            rhs(i) =
              dot(constraints[active[static_cast<std::size_t>(i)]], unconstrained) -
              constraints[active[static_cast<std::size_t>(i)]].bound;
            for (Eigen::Index j = 0; j < count; ++j)
              schur(i, j) = dot(constraints[active[static_cast<std::size_t>(i)]],
                responses[static_cast<std::size_t>(j)]);
            schur(i, i) += Real(1e-12);
          }
          multipliers = schur.ldlt().solve(rhs);
          if (!multipliers.allFinite())
            return false;
          direction = unconstrained;
          for (Eigen::Index i = 0; i < count; ++i)
            direction -= multipliers(i) * responses[static_cast<std::size_t>(i)];
        }
        return false;
      }

    private:
      static Real dot(const Constraint& row, const Math::Vector<Real>& vector)
      {
        Real value = 0;
        for (const auto& [column, coefficient] : row.entries)
          value += coefficient * vector(static_cast<Eigen::Index>(column));
        return value;
      }

      template <class Mesh>
      std::vector<Constraint> getConstraints(const Mesh& mesh) const
      {
        const auto& fes = m_trial.get().getFiniteElementSpace();
        const auto& params = m_parameters.get();
        const std::size_t dim = mesh.getDimension();
        std::vector<Constraint> constraints;
        auto basisJacobian = Variational::Jacobian(m_trial.get());
        CellDeformation deformation(dim);

        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const Index index = cell->getIndex();
          const auto& fe = fes.getFiniteElement(dim, index);
          const std::size_t order = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(2, 2 * fe.getOrder());
          const auto& qf = QF::PolytopeQuadratureFormula::get(order, cell->getGeometry());
          const auto& quad = cell->getQuadrature(qf);
          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const Variational::IntegrationPoint ip(quad.getPoint(qp), &qf, qp);
            deformation.setDisplacementGradient(
              Variational::Jacobian(m_current.get()).getValue(ip));
            if (!deformation.isAdmissible())
              continue;
            basisJacobian.setIntegrationPoint(ip);
            Constraint jacobian;
            Constraint quality;
            jacobian.bound = deformation.getJacobian() - params.jSafe;
            quality.bound = params.qMax - deformation.getRelativeDistortion();
            for (std::size_t local = 0; local < fe.getCount(); ++local)
            {
              const Index global = fes.getGlobalIndex({dim, index}, local);
              const auto gradient = basisJacobian.getBasis(local);
              jacobian.entries.emplace_back(
                global, -deformation.getJacobianAction(gradient));
              quality.entries.emplace_back(
                global, deformation.getRelativeDistortionAction(gradient));
            }
            constraints.push_back(std::move(jacobian));
            constraints.push_back(std::move(quality));
          }
        }
        return constraints;
      }

      std::reference_wrapper<const TrialFunction> m_trial;
      std::reference_wrapper<const Displacement> m_current;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
  };
}

#endif
