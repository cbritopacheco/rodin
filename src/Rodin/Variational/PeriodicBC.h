/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
       * @brief Assembles the Peridodic boundary condition.
       *
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the global degree of freedom map.
       */
      virtual const DOFs& getDOFs() const = 0;

      /**
       * @brief Gets the associated operand.
       */
      virtual const FormLanguage::Base& getOperand() const = 0;

      virtual Boolean isComponent() const = 0;

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

