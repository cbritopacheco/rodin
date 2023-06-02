/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <set>
#include <variant>

#include "Rodin/Utility.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DirichletBCSpecializations DirichletBC Template Specializations
   * @brief Template specializations of the DirichletBC class.
   * @see DirichletBC
   */

  class DirichletBCBase : public FormLanguage::Base
  {
    public:
      virtual void dofs() = 0;
      virtual void project() const = 0;
      virtual const IndexSet& getDOFs() const = 0;
      virtual const FormLanguage::Base& getOperand() const = 0;
      virtual DirichletBCBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Represents a Dirichlet boundary condition on a ShapeFunction
   * object.
   *
   * When utilized in a Problem construction, it will impose the Dirichlet
   * condition
   * @f[
   *   u = g \quad \text{ on } \quad \Gamma_D
   * @f]
   * on the subset of the boundary @f$ \Gamma_D \subset \mathcal{B}_h @f$.
   */
  template <class FES, class ValueDerived>
  class DirichletBC<TrialFunction<FES>, FunctionBase<ValueDerived>> final
    : public DirichletBCBase
  {
    public:
      /// Operand type
      using Operand = TrialFunction<FES>;

      /// Value type
      using Value = FunctionBase<ValueDerived>;

      /// Parent class
      using Parent = DirichletBCBase;

      /**
       * @brief Constructs the object given the Operand and Value.
       * @param[in] u ShapeFunction object
       * @param[in] v Value object
       */
      DirichletBC(Operand& u, const Value& v)
        : m_u(u), m_value(v.copy())
      {}

      /**
       * @brief Copy constructor
       */
      DirichletBC(const DirichletBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_value(other.m_value->copy()),
          m_essBdr(other.m_essBdr),
          m_dofs(other.m_dofs)
      {}

      /**
       * @brief Move constructor
       */
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_value(std::move(other.m_value)),
          m_essBdr(std::move(other.m_essBdr)),
          m_dofs(std::move(other.m_dofs))
      {}

      /**
       * @brief Specifies the region of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAttr Attribute associated to the boundary region
       */
      inline
      constexpr
      DirichletBC& on(Geometry::Attribute bdrAtr)
      {
        return on(std::set<Geometry::Attribute>{bdrAtr});
      }

      /**
       * @brief Specifies the regions of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAttrs Attributes associated to the boundary regions
       */
      inline
      constexpr
      DirichletBC& on(const std::set<Geometry::Attribute>& bdrAttrs)
      {
        assert(bdrAttrs.size() > 0);
        m_essBdr = bdrAttrs;
        return *this;
      }

      /**
       * @returns Returns reference to the value of the boundary condition.
       */
      inline
      constexpr
      const Value& getValue() const
      {
        assert(m_value);
        return *m_value;
      }

      /**
       * @returns Attributes over which the boundary condition is imposed.
       */
      inline
      constexpr
      const std::set<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

      /**
       * @brief Computes the indices of the degrees of freedoms associated to
       * the boundary region.
       */
      inline
      void dofs() override
      {
        const auto& fes = m_u.get().getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        for (auto it = mesh.getBoundary(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (m_essBdr.size() == 0 || m_essBdr.count(polytope.getAttribute()))
          {
            const size_t d = polytope.getDimension();
            const size_t i = polytope.getIndex();
            const auto& fe = fes.getFiniteElement(d, i);
            for (Index local = 0; local < fe.getCount(); local++)
              m_dofs.insert(fes.getGlobalIndex({ d, i }, local));
          }
        }
      }

      /**
       * @brief Projects the value over the associated boundary region.
       */
      inline
      void project() const override
      {
        m_u.get().getSolution().projectOnBoundary(getValue(), m_essBdr);
      }

      inline
      const Operand& getOperand() const override
      {
        return m_u;
      }

      inline
      const IndexSet& getDOFs() const override
      {
        return m_dofs;
      }

      inline
      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<Operand> m_u;
      std::unique_ptr<Value> m_value;
      std::set<Geometry::Attribute> m_essBdr;
      IndexSet m_dofs;
  };

  template <class FES, class ValueDerived>
  DirichletBC(TrialFunction<FES>&, const FunctionBase<ValueDerived>&)
    -> DirichletBC<TrialFunction<FES>, FunctionBase<ValueDerived>>;
}

#endif
