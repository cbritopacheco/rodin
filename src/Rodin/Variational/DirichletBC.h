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
      virtual void project() const = 0;
      virtual const mfem::Array<int>& getDOFs() const = 0;
      virtual const FormLanguage::Base& getOperand() const = 0;
      virtual DirichletBCBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Represents a Dirichlet boundary condition on an H1 TrialFunction
   * object.
   *
   * When utilized in a Problem construction, it will impose the Dirichlet
   * condition
   * @f[
   *   u = g \text{ on } \Gamma_D
   * @f]
   * on the segment of the boundary @f$ \Gamma_D \subset \partial \Omega @f$
   * specified by the boundary attribute.
   */
  template <class OperandDerived, class ValueDerived, ShapeFunctionSpaceType Space, class ... Ts>
  class DirichletBC<ShapeFunction<OperandDerived, H1<Ts...>, Space>, FunctionBase<ValueDerived>> final
    : public DirichletBCBase
  {
    public:
      using Operand = ShapeFunction<OperandDerived, H1<Ts...>, Space>;
      using Value = FunctionBase<ValueDerived>;
      using Parent = DirichletBCBase;

      constexpr
      DirichletBC(Operand& u, const Value& v)
        : m_u(u), m_value(v.copy())
      {}

      constexpr
      DirichletBC(const DirichletBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_value(other.m_value->copy()),
          m_essBdr(other.m_essBdr)
      {}

      constexpr
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_value(std::move(other.m_value)),
          m_essBdr(std::move(other.m_essBdr))
      {}

      inline
      constexpr
      DirichletBC& on(Geometry::Attribute bdrAtr)
      {
        return on(std::set<Geometry::Attribute>{bdrAtr});
      }

      inline
      constexpr
      DirichletBC& on(const std::set<Geometry::Attribute>& bdrAtr)
      {
        m_essBdr = bdrAtr;
        return *this;
      }

      /**
       * @returns Returns reference to the value of the boundary condition
       * at the boundary
       */
      inline
      constexpr
      const Value& getValue() const
      {
        assert(m_value);
        return *m_value;
      }

      inline
      constexpr
      const std::set<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

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
      const mfem::Array<int>& getDOFs() const override
      {
        if (!m_dofs.has_value())
          m_dofs.emplace(m_u.get().getFiniteElementSpace().getEssentialTrueDOFs(m_essBdr));
        assert(m_dofs.has_value());
        return m_dofs.value();
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
      mutable std::optional<const mfem::Array<int>> m_dofs;
  };

  template <class OperandDerived, class ValueDerived, ShapeFunctionSpaceType Space, class ... Ts>
  DirichletBC(ShapeFunction<OperandDerived, H1<Ts...>, Space>&, const FunctionBase<ValueDerived>&)
    -> DirichletBC<ShapeFunction<OperandDerived, H1<Ts...>, Space>, FunctionBase<ValueDerived>>;

  // /**
  //  * @ingroup DirichletBCSpecializations
  //  */
  // template <class OperandDerived, class ValueDerived, ShapeFunctionSpaceType Space, class ... Ts>
  // class DirichletBC<Component<ShapeFunction<OperandDerived, H1<Ts...>, Space>>, FunctionBase<ValueDerived>>
  //   : public DirichletBCBase
  // {
  //   public:
  //     using Operand = Component<ShapeFunction<OperandDerived, H1<Ts...>, Space>>;
  //     using Value = FunctionBase<ValueDerived>;
  //     using Parent = DirichletBCBase;

  //     constexpr
  //     DirichletBC(const Operand& ux, const Value& v)
  //       : m_ux(ux), m_value(v)
  //     {}

  //     constexpr
  //     DirichletBC(const DirichletBC& other)
  //       : Parent(other),
  //         m_ux(other.m_ux),
  //         m_value(other.m_value),
  //         m_essBdr(other.m_essBdr)
  //     {}

  //     constexpr
  //     DirichletBC(DirichletBC&& other)
  //       : Parent(std::move(other)),
  //         m_ux(std::move(other.m_ux)),
  //         m_value(std::move(other.m_value)),
  //         m_essBdr(std::move(other.m_essBdr))
  //     {}

  //     inline
  //     constexpr
  //     DirichletBC& on(int bdrAtr)
  //     {
  //       return on(std::set<int>{bdrAtr});
  //     }

  //     inline
  //     constexpr
  //     DirichletBC& on(const std::set<Geometry::Attribute>& bdrAtr)
  //     {
  //       m_essBdr = bdrAtr;
  //       return *this;
  //     }

  //     /**
  //      * @returns Returns reference to the value of the boundary condition
  //      * at the boundary
  //      */
  //     inline
  //     constexpr
  //     Value& getValue() const
  //     {
  //       return m_value;
  //     }

  //     inline
  //     void project() const override
  //     {
  //       assert(false);
  //     }

  //     inline
  //     const Operand& getOperand() const override
  //     {
  //       return m_ux;
  //     }

  //     inline
  //     const std::set<size_t>& getDOFs() const override
  //     {
  //       assert(false);
  //       return m_dofs;
  //     }

  //     inline
  //     constexpr
  //     const std::set<Geometry::Attribute>& getAttributes() const
  //     {
  //       return m_essBdr;
  //     }

  //     inline
  //     DirichletBC* copy() const noexcept
  //     override
  //     {
  //       return new DirichletBC(*this);
  //     }
  //   private:
  //     Operand m_ux;
  //     Value m_value;
  //     std::set<Geometry::Attribute> m_essBdr;
  //     mutable std::set<size_t> m_dofs;
  // };

  // template <class OperandDerived, class ValueDerived, ShapeFunctionSpaceType Space, class ... Ts>
  // DirichletBC(const Component<ShapeFunction<OperandDerived, H1<Ts...>, Space>>&, const FunctionBase<ValueDerived>&)
  //   -> DirichletBC<Component<ShapeFunction<OperandDerived, H1<Ts...>, Space>>, FunctionBase<ValueDerived>>;
}

#endif
