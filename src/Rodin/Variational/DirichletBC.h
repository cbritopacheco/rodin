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
#include "Rodin/Math/Vector.h"
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

      /**
       * @brief DOF map for the value-based Dirichlet condition `u = g`.
       *
       * Maps each constrained global FES DOF index to its prescribed scalar
       * value @f$ g(x_i) @f$.
       */
      using ValueDOFs = IndexMap<ScalarType>;

      /**
       * @brief DOF map for the identification-based Dirichlet condition
       *        `u = A(v)`.
       *
       * Maps each slave DOF (in @f$ u @f$'s FES) to a pair of arrays:
       *  - the master DOF indices in @f$ v @f$'s FES (FES-global, pre-offset),
       *  - the corresponding scalar coefficients @f$ C_{sj} = [A(\varphi_j^v)](x_s) @f$.
       *
       * The constraint encoded is @f$ u_s = \sum_j C_{sj}\, v_j @f$.
       */
      using IdentifiedDOFs =
        IndexMap<std::pair<IndexArray, Math::Vector<ScalarType>>>;

      /**
       * @brief Variant DOF type covering both Dirichlet semantics.
       *
       * The active alternative depends on which subclass produced it: a
       * value-prescribing BC writes the @ref ValueDOFs alternative, an
       * identification BC writes the @ref IdentifiedDOFs alternative.
       */
      using DOFs = std::variant<ValueDOFs, IdentifiedDOFs>;

      /**
       * @brief Assembles the Dirichlet boundary condition.
       *
       * Either fills the @ref ValueDOFs alternative (for value-prescribing
       * BCs) or the @ref IdentifiedDOFs alternative (for identification BCs).
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the map of constrained DOFs.
       * @return Variant holding either the ValueDOFs or the IdentifiedDOFs map.
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
       * @brief Returns the UUID of @f$ v @f$'s leaf trial function for an
       *        identification BC, or @c nullopt for a value-prescribing BC.
       *
       * Consumer assemblies use this to locate @f$ v @f$'s trial block and
       * apply the correct global offset when assembling identified-DOF
       * constraints.
       */
      virtual Optional<Identifiable::UUID> getValueUUID() const { return {}; }

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

      /// Parent class
      using Parent = DirichletBCBase<ScalarType>;

      /// Value DOF map type (the alternative populated by this specialization)
      using ValueDOFs = typename Parent::ValueDOFs;

      /// Variant DOFs type, inherited from Parent
      using DOFs = typename Parent::DOFs;

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
        typename Assembly::Default<FESMeshContextType>::template Type<ValueDOFs, DirichletBC>;

      using AssemblyType =
        DefaultAssemblyType;

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
       * @param[in] bdrAtr Attribute associated to the boundary region
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
        // Ensure variant holds the ValueDOFs alternative.
        if (!std::holds_alternative<ValueDOFs>(m_dofs))
          m_dofs = ValueDOFs{};
        m_assembly.execute(
            std::get<ValueDOFs>(m_dofs), { m_u.get(), *m_value, m_essBdr });
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

      const Assembly::AssemblyBase<ValueDOFs, DirichletBC>& getAssembly() const
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
      DOFs m_dofs{ValueDOFs{}};
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

  /**
   * @ingroup DirichletBCSpecializations
   * @brief Identification Dirichlet boundary condition `u = A(v)`.
   *
   * Constrains the DOFs of @f$ u @f$ on a boundary subset to equal a linear
   * expression @f$ A(v) @f$, where @f$ v @f$ is itself a (trial) shape
   * function and @f$ A @f$ is any expression that yields a
   * @ref ShapeFunctionBase (e.g. @f$ v @f$ itself, a component
   * @f$ v_x @f$, a product @f$ f \cdot v @f$, a matrix product
   * @f$ R \cdot v @f$, etc.).
   *
   * The constraint encoded for each slave DOF @f$ s @f$ is
   * @f[
   *   u_s \;=\; \sum_j C_{sj}\, v_j \,, \qquad
   *   C_{sj} \;=\; [A(\varphi_j^v)](x_s)
   * @f]
   * with all DOF pairings determined exactly by the FES connectivity
   * (`getDOFs(face)`) — no geometric search and no tolerance.
   *
   * @tparam Solution Solution type of the trial @f$ u @f$
   * @tparam FES1 Finite element space of @f$ u @f$
   * @tparam Derived2 CRTP-derived type of @f$ A(v) @f$
   * @tparam FES2 Finite element space of @f$ v @f$
   * @tparam Sp Shape function space type (Trial or Test) of @f$ A(v) @f$
   */
  template <class Solution, class FES1,
            class Derived2, class FES2, ShapeFunctionSpaceType Sp>
  class DirichletBC<TrialFunction<Solution, FES1>,
                    ShapeFunctionBase<Derived2, FES2, Sp>> final
    : public DirichletBCBase<typename FormLanguage::Traits<FES1>::ScalarType>
  {
    public:
      using FESType = FES1;

      /// Operand type (the slave trial function)
      using OperandType = TrialFunction<Solution, FESType>;

      /// Value type (the right-hand-side shape function expression)
      using ValueType = ShapeFunctionBase<Derived2, FES2, Sp>;

      /// Scalar type
      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      /// Parent class
      using Parent = DirichletBCBase<ScalarType>;

      /// Identified-DOFs alternative populated by this specialization
      using IdentifiedDOFs = typename Parent::IdentifiedDOFs;

      /// Variant DOFs type
      using DOFs = typename Parent::DOFs;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESRangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using DefaultAssemblyType =
        typename Assembly::Default<FESMeshContextType>::template Type<IdentifiedDOFs, DirichletBC>;

      using AssemblyType =
        DefaultAssemblyType;

      /**
       * @brief Constructs the identification Dirichlet BC.
       * @param[in] u Slave trial function
       * @param[in] v Right-hand-side shape function expression @f$ A(v) @f$
       */
      DirichletBC(const OperandType& u, const ValueType& v)
        : m_u(u), m_v(v.copy())
      {}

      /// Copy constructor
      DirichletBC(const DirichletBC& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v->copy()),
          m_essBdr(other.m_essBdr),
          m_dofs(other.m_dofs),
          m_assembly(other.m_assembly)
      {}

      /// Move constructor
      DirichletBC(DirichletBC&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_essBdr(std::move(other.m_essBdr)),
          m_dofs(std::move(other.m_dofs)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Specifies the boundary attribute over which the BC applies.
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
       * @brief Specifies the set of boundary attributes for the BC.
       */
      constexpr
      DirichletBC& on(const FlatSet<Geometry::Attribute>& bdrAttrs)
      {
        assert(bdrAttrs.size() > 0);
        m_essBdr = bdrAttrs;
        return *this;
      }

      constexpr
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_essBdr;
      }

      void assemble() override
      {
        if (!std::holds_alternative<IdentifiedDOFs>(m_dofs))
          m_dofs = IdentifiedDOFs{};
        m_assembly.execute(
            std::get<IdentifiedDOFs>(m_dofs), { m_u.get(), *m_v, m_essBdr });
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
        assert(m_v);
        return *m_v;
      }

      const DOFs& getDOFs() const override
      {
        return m_dofs;
      }

      Optional<Identifiable::UUID> getValueUUID() const override
      {
        assert(m_v);
        return m_v->getLeaf().getUUID();
      }

      const Assembly::AssemblyBase<IdentifiedDOFs, DirichletBC>& getAssembly() const
      {
        return m_assembly;
      }

      DirichletBC* copy() const noexcept override
      {
        return new DirichletBC(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      std::unique_ptr<ValueType> m_v;
      FlatSet<Geometry::Attribute> m_essBdr;
      DOFs m_dofs{IdentifiedDOFs{}};
      AssemblyType m_assembly;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for the identification DirichletBC.
   */
  template <class Solution, class FES1,
            class Derived2, class FES2, ShapeFunctionSpaceType Sp>
  DirichletBC(const TrialFunction<Solution, FES1>&,
              const ShapeFunctionBase<Derived2, FES2, Sp>&)
    -> DirichletBC<TrialFunction<Solution, FES1>,
                   ShapeFunctionBase<Derived2, FES2, Sp>>;
}

#endif
