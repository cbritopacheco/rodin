/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGE_MODEL_H
#define RODIN_ADAPTATION_SWIFT_HINGE_MODEL_H

#include <functional>
#include <memory>
#include <optional>
#include <variant>
#include <utility>
#include <vector>
#include "Rodin/Variational/Jacobian.h"
#include "State.h"
#include "Weights.h"

namespace Rodin::Adaptation::SWIFT::Hinge
{
  /// @brief Cell-local hinge evaluation with its weights in the covering's own order.
  class Cell
  {
    public:
      /// @returns The same geometric witnesses used for unweighted admissibility.
      const Witnesses& getWitnesses() const
      {
        return m_witnesses;
      }

      /// @returns Frozen measures for these witnesses.
      const Weights& getWeights() const
      {
        if (const auto* local = std::get_if<Weights>(&m_weights))
          return *local;
        return std::get<std::reference_wrapper<const Weights>>(m_weights).get();
      }

      /**
       * @brief Visits each point with the corresponding two penalty measures.
       * @param evaluate Callable accepting a mapped point and its two weights.
       */
      template <class Evaluate>
      void forEach(Evaluate&& evaluate) const
      {
        const auto& weights = getWeights();
        for (size_t q = 0; q < m_witnesses.getSize(); ++q)
        {
          const auto measure = weights[q];
          evaluate(m_witnesses.getPoint(q), measure.jacobian, measure.distortion);
        }
      }

    private:
      template <class>
      friend class Model;
      using Storage = std::variant<std::reference_wrapper<const Weights>, Weights>;

      Cell(const Witnesses& witnesses, Storage weights)
        : m_witnesses(witnesses),
          m_weights(std::move(weights))
      {}

      Witnesses m_witnesses;
      Storage m_weights;
  };

  /**
   * @brief Affine hinge model at one outer state.
   *
   * The displacement, scaled predictor and parameters are borrowed and must
   * remain fixed during the inner solve. Form copies share the immutable
   * table; the inner increment is not part of this model. Cell views borrow
   * the mesh geometry and, for precomputed measures, the model's table.
   */
  template <class Displacement>
  class Model
  {
    public:
      /**
       * @brief Binds the outer fields, computing weights on demand.
       * @param current Frozen outer displacement.
       * @param predictor Scaled predictor increment.
       * @param parameters Quality and sampling policy for this model.
       * @param coefficient Predicted decrease per unit background measure.
       */
      Model(const Displacement& current, const Displacement& predictor,
        const Parameters& parameters, Real coefficient)
        : m_current(current),
          m_predictor(predictor),
          m_parameters(parameters),
          m_coefficient(coefficient)
      {}

      /**
       * @brief Precomputes the weights shared by assembly and merit evaluation.
       * @param current Frozen outer displacement.
       * @param predictor Scaled predictor increment.
       * @param parameters Quality and sampling policy for this model.
       * @param coefficient Predicted decrease per unit background measure.
       * @param cells Distinct indices of the cells on which to freeze weights.
       */
      Model(const Displacement& current, const Displacement& predictor,
        const Parameters& parameters, Real coefficient, const std::vector<Index>& cells)
        : Model(current, predictor, parameters, coefficient)
      {
        const auto& fes = current.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        auto table = std::make_shared<Table>(mesh.getCellCount());
#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          const auto currentJacobian = Variational::Jacobian(current);
          const auto predictorJacobian = Variational::Jacobian(predictor);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(cells.size()); ++i)
          {
            const Index index = cells[i];
            const auto cell = mesh.getCell(index);
            const Witnesses witnesses(*cell, getDegree(*cell), parameters.sampling);
            (*table)[index].emplace(witnesses,
              Weights(witnesses, currentJacobian, predictorJacobian, parameters.model));
          }
        }
        m_table = std::move(table);
      }

      /// @returns The frozen displacement.
      const Displacement& getDisplacement() const
      {
        return m_current;
      }

      /// @returns The frozen model and sampling parameters.
      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      /// @returns The predicted decrease per unit background measure.
      Real getCoefficient() const
      {
        return m_coefficient;
      }

      /**
       * @param cell Background cell.
       * @returns A view on the displacement space's covering for this cell.
       */
      Cell on(const Geometry::Polytope& cell) const
      {
        return on(cell, getDegree(cell));
      }

      /**
       * @param cell Background cell.
       * @param degree Order requested by the form.
       * @returns A view for a form's degree; reuse requires the identical mapped covering.
       */
      Cell on(const Geometry::Polytope& cell, size_t degree) const
      {
        const Witnesses witnesses(cell, degree, m_parameters.sampling);
        const auto index = cell.getIndex();
        if (m_table && index < m_table->size())
        {
          const auto& entry = (*m_table)[index];
          if (entry && &entry->first.getQuadrature() == &witnesses.getQuadrature())
            return Cell(witnesses, std::cref(entry->second));
        }
        return Cell(witnesses,
          Weights(witnesses, Variational::Jacobian(m_current),
            Variational::Jacobian(m_predictor), m_parameters.model));
      }

      /**
       * @param increment Proposed inner displacement increment.
       * @param cells Background cells included in the penalty.
       * @returns The penalty contribution to the inner merit at the proposed increment.
       */
      Real getEnergy(const Displacement& increment, const std::vector<Index>& cells) const
      {
        const auto& mesh = m_current.getFiniteElementSpace().getMesh();
        Real energy = 0;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(+ : energy)
#endif
        {
          auto currentJacobian = Variational::Jacobian(m_current);
          auto innerJacobian = Variational::Jacobian(increment);
          CellDeformation deformation(mesh.getDimension());
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(cells.size()); ++i)
          {
            const auto cell = mesh.getCell(cells[i]);
            on(*cell).forEach(
              [&](const Variational::IntegrationPoint& ip, Real wj, Real wq) {
                deformation.setDisplacementGradient(currentJacobian.getValue(ip));
                const State state(
                  deformation, innerJacobian.getValue(ip), m_parameters, m_coefficient);
                energy += state.getEnergy(m_parameters, m_coefficient, wj, wq);
              });
          }
        }
        return energy;
      }

    private:
      using Table = std::vector<std::optional<std::pair<Witnesses, Weights>>>;

      size_t getDegree(const Geometry::Polytope& cell) const
      {
        const auto& fe = m_current.getFiniteElementSpace().getFiniteElement(
          cell.getDimension(), cell.getIndex());
        return fe.getDegree().value_or(2 * fe.getOrder());
      }

      const Displacement& m_current;
      const Displacement& m_predictor;
      const Parameters& m_parameters;
      Real m_coefficient;
      std::shared_ptr<const Table> m_table;
  };
}

#endif
