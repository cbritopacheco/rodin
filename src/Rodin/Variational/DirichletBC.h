/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DirichletBC.h
 * @brief Dirichlet boundary condition implementation.
 *
 * This file defines classes for imposing Dirichlet (essential) boundary
 * conditions in finite element problems. Dirichlet conditions prescribe the
 * solution values on specified boundaries.
 *
 * ## Mathematical Foundation
 * A Dirichlet boundary condition specifies:
 * @f[
 *   u = g \quad \text{on} \quad \Gamma_D
 * @f]
 * where @f$ u @f$ is the solution, @f$ g @f$ is the prescribed boundary value,
 * and @f$ \Gamma_D @f$ is a portion of the domain boundary.
 *
 * ## Implementation
 * Dirichlet conditions are typically enforced by:
 * 1. Identifying boundary degrees of freedom
 * 2. Setting their values to @f$ g @f$
 * 3. Modifying the system matrix and right-hand side
 *
 * ## Usage Example
 * ```cpp
 * // Homogeneous Dirichlet BC (u = 0 on boundary)
 * auto bc = DirichletBC(u, Zero());
 * 
 * // Inhomogeneous Dirichlet BC
 * auto g = [](const Point& p) { return sin(p.x()); };
 * auto bc = DirichletBC(u, g).on(1);  // On boundary attribute 1
 * ```
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

#include <set>
#include <variant>

#include "Rodin/Utility.h"
#include "Rodin/FormLanguage/List.h"

#include "Rodin/Assembly/ForwardDecls.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"
#include "VectorFunction.h"

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
  template <class Scalar>
  class DirichletBCBase : public FormLanguage::Base
  {
    public:
      using ScalarType = Scalar;

      using DOFs = IndexMap<ScalarType>;

      /**
       * @brief Assembles the Dirichlet boundary condition.
       *
       * Computes the global DOF map for the Dirichlet boundary by evaluating
       * the prescribed value at each boundary DOF. The result is a map
       * @f$ \{(i, g(x_i))\} @f$ where @f$ i @f$ is the global DOF index and
       * @f$ g(x_i) @f$ is the prescribed value at that DOF.
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the map of constrained DOFs and their values.
       * @return Map from global DOF index to prescribed value
       *
       * Returns the assembled DOF map containing pairs @f$ (i, g_i) @f$ where
       * @f$ i @f$ is the global DOF index and @f$ g_i @f$ is the prescribed
       * boundary value.
       */
      virtual const DOFs& getDOFs() const = 0;

      /**
       * @brief Checks if this is a component-wise boundary condition.
       * @return True if BC applies to a single component of a vector function
       */
      virtual bool isComponent() const = 0;

      /**
       * @brief Gets the operand (trial function) of the BC.
       * @return Reference to the trial function being constrained
       */
      virtual const FormLanguage::Base& getOperand() const = 0;

      /**
       * @brief Gets the prescribed boundary value.
       * @return Reference to the function defining the BC value
       */
      virtual const FormLanguage::Base& getValue() const = 0;

      /**
       * @brief Creates a polymorphic copy of this BC.
       * @return Pointer to a new copy
       */
      virtual DirichletBCBase* copy() const noexcept override = 0;
  };

  /// Alias for a list of Dirichlet boundary conditions
  template <class Scalar>
  using EssentialBoundary = FormLanguage::List<DirichletBCBase<Scalar>>;

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
  template <class Solution, class FES, class ValueDerived>
  class DirichletBC<TrialFunction<Solution, FES>, FunctionBase<ValueDerived>> final
    : public DirichletBCBase<typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType =
        FES;

      /// Operand type
      using OperandType =
        TrialFunction<Solution, FESType>;

      /// Scalar type
      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      using DOFs =
        IndexMap<ScalarType>;

      /// Value type
      using ValueType =
        FunctionBase<ValueDerived>;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using DefaultAssemblyType =
        typename Assembly::Default<FESMeshContextType>::template Type<DOFs, DirichletBC>;

      using AssemblyType =
        DefaultAssemblyType;

      /// Parent class
      using Parent = DirichletBCBase<ScalarType>;

      /**
       * @brief Constructs the object given the Operand and Value.
       * @param[in] u ShapeFunction object
       * @param[in] v Value object
       */
      DirichletBC(const OperandType& u, const ValueType& v)
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
          m_dofs(other.m_dofs),
          m_assembly(other.m_assembly)
      {}

      /**
       * @brief Move constructor
       */
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_value(std::move(other.m_value)),
          m_essBdr(std::move(other.m_essBdr)),
          m_dofs(std::move(other.m_dofs)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Specifies the region of the boundary over which the condition
       * will be imposed.
       * @param[in] bdrAttr Attribute associated to the boundary region
       */
      constexpr
      DirichletBC& on(Geometry::Attribute bdrAtr)
      {
        return on(FlatSet<Geometry::Attribute>{bdrAtr});
      }

      template <class A1, class A2, class ... As>
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
      void assemble() override
      {
        m_assembly.execute(m_dofs, { m_u.get(), *m_value, m_essBdr });
      }

      bool isComponent() const override
      {
        return false;
      }

      const OperandType& getOperand() const override
      {
        return m_u.get();
      }

      const ValueType& getValue() const override
      {
        assert(m_value);
        return *m_value;
      }

      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      const Assembly::AssemblyBase<IndexMap<ScalarType>, DirichletBC>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      std::unique_ptr<ValueType> m_value;
      FlatSet<Geometry::Attribute> m_essBdr;
      IndexMap<ScalarType> m_dofs;
      AssemblyType m_assembly;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for DirichletBC
   * @tparam FES Type of finite element space
   * @tparam ValueDerived Derived type of FunctionBase
   */
  template <class Solution, class FES, class FunctionDerived>
  DirichletBC(const TrialFunction<Solution, FES>&, const FunctionBase<FunctionDerived>&)
    -> DirichletBC<TrialFunction<Solution, FES>, FunctionBase<FunctionDerived>>;
}

#endif
