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
#include "Rodin/FormLanguage/List.h"

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

  /**
   * @brief Abstract base class for a Dirichlet boundary condition.
   *
   * Used as a base class to represent the Dirichlet boundary condition:
   * @f[
   *   \mathrm{Operand} = \mathrm{Value} \ \text{ on } \ \Gamma_D
   * @f]
   * on some subset of the boundary @f$ \Gamma_D \subset \mathcal{B}_h @f$.
   *
   * @see DirichletBC
   */
  class DirichletBCBase : public FormLanguage::Base
  {
    public:
      using DOFs = IndexMap<Scalar>;

      /**
       * @brief Assembles the Dirichlet boundary condition.
       *
       * This method computes the global degree of freedom map associated to
       * the Dirichlet boundary. In other words, it computes the IndexMap which
       * has a keys the global indices of the DOFs, and as values @f$ \ell_i
       * @f$ the
       * @f[
       *  \ell_i(\mathrm{Value}), \quad i = 1, \ldots, n
       * @f]
       * where @f$ \ell_i @f$ is the i-th linear form on the associated finite
       * element space.
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

      /**
       * @brief Gets the associated value.
       */
      virtual const FormLanguage::Base& getValue() const = 0;

      virtual DirichletBCBase* copy() const noexcept override = 0;
  };

  /// Alias for a list of Dirichlet boundary conditions
  using EssentialBoundary = FormLanguage::List<DirichletBCBase>;

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Represents a Dirichlet boundary condition on a ShapeFunction
   * object.
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Type of value
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
      DirichletBC(const Operand& u, const Value& v)
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
        return on(FlatSet<Geometry::Attribute>{bdrAtr});
      }

      template <class A1, class A2, class ... As>
      inline
      constexpr
      DirichletBC& on(A1 a1, A2 a2, As... as)
      {
        return on(FlatSet<Geometry::Attribute>{ a1, a2, as... });
      }

      /**
       * @brief Specifies the regions of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAttrs Attributes associated to the boundary regions
       */
      inline
      constexpr
      DirichletBC& on(const FlatSet<Geometry::Attribute>& bdrAttrs)
      {
        assert(bdrAttrs.size() > 0);
        m_essBdr = bdrAttrs;
        return *this;
      }

      /**
       * @returns Attributes over which the boundary condition is imposed.
       */
      inline
      constexpr
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

      /**
       * @brief Computes the indices of the degrees of freedoms associated to
       * the boundary region.
       *
       * This will compute the degrees of freedom over the incidence set
       * @f[
       *  D - 1 \longrightarrow 0 ~.
       * @f]
       *
       * If the set of specified attributes is empty, this will
       * compute the degrees of freedom over the boundary, in which case the
       * incidence set
       * @f[
       *  D - 1 \longrightarrow D
       * @f]
       * is also required.
       */
      inline
      void assemble() override
      {
        const auto& fes = m_u.get().getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        Geometry::FaceIterator it;
        if (m_essBdr.size() > 0)
          it = mesh.getFace();
        else
          it = mesh.getBoundary();
        for (; !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (m_essBdr.size() == 0 || m_essBdr.count(polytope.getAttribute()))
          {
            const size_t d = polytope.getDimension();
            const size_t i = polytope.getIndex();
            const auto& fe = fes.getFiniteElement(d, i);
            const auto& mapping = fes.getMapping({ d, i }, getValue());
            for (Index local = 0; local < fe.getCount(); local++)
            {
              const Index global = fes.getGlobalIndex({ d, i }, local);
              auto find = m_dofs.find(global);
              if (find == m_dofs.end())
              {
                const auto& lf = fe.getLinearForm(local);
                const Scalar s = lf(mapping);
                m_dofs.insert(find, std::pair{ global, s });
              }
            }
          }
        }
      }

      inline
      const Operand& getOperand() const override
      {
        return m_u;
      }

      inline
      const Value& getValue() const override
      {
        assert(m_value);
        return *m_value;
      }

      inline
      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      inline
      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
      std::unique_ptr<Value> m_value;
      FlatSet<Geometry::Attribute> m_essBdr;
      IndexMap<Scalar> m_dofs;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for DirichletBC
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Derived type of FunctionBase
   */
  template <class FES, class FunctionDerived>
  DirichletBC(const TrialFunction<FES>&, const FunctionBase<FunctionDerived>&)
    -> DirichletBC<TrialFunction<FES>, FunctionBase<FunctionDerived>>;
}

#endif
