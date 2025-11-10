/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PeriodicBC.h
 * @brief Periodic boundary condition implementation.
 *
 * This file defines classes for imposing periodic boundary conditions in
 * finite element problems. Periodic conditions enforce that the solution
 * values match across opposite boundaries of the domain.
 *
 * ## Mathematical Foundation
 * A periodic boundary condition specifies:
 * @f[
 *   u(x + L) = u(x) \quad \text{for} \quad x \in \Gamma_-
 * @f]
 * where @f$ \Gamma_- @f$ and @f$ \Gamma_+ @f$ are matching boundary pairs
 * separated by periodicity vector @f$ L @f$.
 *
 * ## Applications
 * - Crystal lattice simulations
 * - Unit cell problems in homogenization
 * - Fluid flow in periodic geometries
 * - Wave propagation in periodic media
 *
 * ## Implementation
 * Periodic conditions are enforced by:
 * 1. Identifying matching DOF pairs on opposite boundaries
 * 2. Constraining one DOF to equal its periodic partner
 * 3. Eliminating dependent DOFs from the system
 *
 * ## Usage Example
 * ```cpp
 * // Periodic BC between boundaries 1 and 2
 * auto pbc = PeriodicBC(u).from(1).to(2);
 * ```
 */
#ifndef RODIN_VARIATIONAL_PERIODICBC_H
#define RODIN_VARIATIONAL_PERIODICBC_H

#include <set>
#include <variant>

#include "Rodin/Utility.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup PeridodicBCSpecializations PeriodicBC Template Specializations
   * @brief Template specializations of the PeriodicBC class.
   * @see PeriodicBC
   */

  /**
   * @brief Abstract base class for a periodic boundary condition.
   *
   * @see PeriodicBC
   */
  template <class Scalar>
  class PeriodicBCBase : public FormLanguage::Base
  {
    public:
      using ScalarType = Scalar;

      using DOFs = IndexMap<std::pair<IndexArray, Math::Vector<ScalarType>>>;

      /**
       * @brief Assembles the periodic boundary condition.
       *
       * Computes the DOF constraints enforcing periodicity across specified
       * boundaries. Creates relationships between DOFs on opposite periodic
       * boundaries such that they share the same value.
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the map of periodic DOF constraints.
       * @return Map containing DOF relationships for periodicity
       *
       * Returns pairs of DOF indices and constraint coefficients defining
       * the periodic relationships between degrees of freedom.
       */
      virtual const DOFs& getDOFs() const = 0;

      /**
       * @brief Gets the operand (trial function) of the periodic BC.
       * @return Reference to the trial function being constrained
       */
      virtual const FormLanguage::Base& getOperand() const = 0;

      /**
       * @brief Checks if this is a component-wise periodic BC.
       * @return True if BC applies to a single component
       */
      virtual Boolean isComponent() const = 0;

      /**
       * @brief Creates a polymorphic copy of this periodic BC.
       * @return Pointer to a new copy
       */
      virtual PeriodicBCBase* copy() const noexcept override = 0;
  };

  /// Alias for a list of peridodic boundary conditions
  template <class Scalar>
  using PeriodicBoundary = FormLanguage::List<PeriodicBCBase<Scalar>>;

  /**
   * @ingroup PeridodicBCSpecializations
   * @brief Represents a Peridodic boundary condition on a ShapeFunction
   * object.
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Type of value
   *
   */
  template <class Solution, class FES>
  class PeriodicBC<TrialFunction<Solution, FES>, IndexMap<IndexSet>> final
    : public PeriodicBCBase<typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType = FES;

      using SolutionType = Solution;

      /// Operand type
      using OperandType = TrialFunction<SolutionType, FESType>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DOFs = IndexMap<std::pair<IndexArray, Math::Vector<ScalarType>>>;

      /// Parent class
      using Parent = PeriodicBCBase<ScalarType>;

      PeriodicBC(const OperandType& u, const IndexMap<IndexSet>& adjacency)
        : m_u(u),
          m_adjacency(adjacency)
      {}

      PeriodicBC(const OperandType& u, IndexMap<IndexSet>&& adjacency)
        : m_u(u),
          m_adjacency(std::move(adjacency))
      {}

      /**
       * @brief Copy constructor
       */
      PeriodicBC(const PeriodicBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_adjacency(other.m_adjacency)
      {}

      /**
       * @brief Move constructor
       */
      PeriodicBC(PeriodicBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_adjacency(std::move(other.m_adjacency))
      {}

      /**
       * @brief Computes the indices of the degrees of freedoms associated to
       * the boundary region.
       */
      void assemble() override
      {
        const auto& adjacency = getAdjacency();
        m_dofs.reserve(adjacency.size());
        for (const auto& [k, v] : adjacency)
        {
          IndexArray dofs(v.size());
          size_t i = 0;
          for (const auto& child : v)
            dofs.coeffRef(i++) = child;
          Math::Vector<ScalarType> coeffs(v.size());
          coeffs.array() /= v.size();
          m_dofs.emplace(k, std::pair{ std::move(dofs), std::move(coeffs) });
        }
      }

      const IndexMap<IndexSet>& getAdjacency() const
      {
        return m_adjacency;
      }

      const OperandType& getOperand() const override
      {
        return m_u;
      }

      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      Boolean isComponent() const override
      {
        return false;
      }

      PeriodicBC* copy() const noexcept override
      {
        return new PeriodicBC(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      IndexMap<IndexSet> m_adjacency;
      DOFs m_dofs;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for PeriodicBC
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Derived type of FunctionBase
   */
  template <class Solution, class FES>
  PeriodicBC(const TrialFunction<Solution, FES>&, const IndexMap<IndexSet>&)
    -> PeriodicBC<TrialFunction<Solution, FES>, IndexMap<IndexSet>>;
}

#endif

