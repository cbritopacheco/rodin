/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_H
#define RODIN_VARIATIONAL_DIRICHLETBC_H

/**
 * @file
 * @brief Dirichlet (essential) boundary conditions for finite element problems.
 *
 * This file provides facilities for imposing Dirichlet boundary conditions,
 * also known as essential boundary conditions, which prescribe solution values
 * on portions of the domain boundary.
 *
 * # Mathematical Foundation
 *
 * A Dirichlet boundary condition specifies the value of the solution on a
 * boundary subset @f$ \Gamma_D \subset \partial\Omega @f$:
 * @f[
 *   u = g \quad \text{on } \Gamma_D
 * @f]
 * where @f$ g @f$ is a given function representing the prescribed boundary value.
 *
 * ## Implementation in FEM
 * In the discrete setting, Dirichlet conditions are enforced by:
 * 1. Identifying DOFs associated with boundary vertices/edges/faces in @f$ \Gamma_D @f$
 * 2. Setting these DOFs to match the prescribed values @f$ g(x_i) @f$
 * 3. Modifying the linear system to preserve these constraints
 *
 * ## Types of Boundary Conditions
 * - **Homogeneous**: @f$ u = 0 @f$ on @f$ \Gamma_D @f$
 * - **Inhomogeneous**: @f$ u = g(x) @f$ on @f$ \Gamma_D @f$ for non-zero @f$ g @f$
 * - **Component-wise**: Different values for different components (vector problems)
 *
 * # Usage Examples
 *
 * ## Homogeneous Dirichlet BC (Zero Boundary)
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * P1 Vh(mesh);
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Poisson with u = 0 on entire boundary
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - Integral(f * v)
 *         + DirichletBC(u, 0.0);  // Homogeneous BC
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Inhomogeneous Dirichlet BC (Non-Zero Function)
 * @code{.cpp}
 * // Boundary value function
 * auto g = ScalarFunction([](const Math::SpatialPoint& x) {
 *   return std::sin(x(0)) * std::cos(x(1));
 * });
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - Integral(f * v)
 *         + DirichletBC(u, g);  // Inhomogeneous BC: u = g on boundary
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Selective Boundary Conditions (Specific Attributes)
 * @code{.cpp}
 * // Apply BC only on boundary regions with specific attributes
 * std::set<int> boundaryAttrs = {1, 2, 3};  // Bottom, left, right walls
 *
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - Integral(f * v)
 *         + DirichletBC(u, 0.0).on(boundaryAttrs);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Vector-Valued Dirichlet BC (Elasticity)
 * @code{.cpp}
 * // 2D elasticity with vector displacement
 * P1 Vh(mesh, 2);  // 2-component space
 *
 * TrialFunction u(Vh);  // Displacement vector
 * TestFunction  v(Vh);
 *
 * // Prescribed displacement on boundary
 * auto g = VectorFunction([](const Math::SpatialPoint& x) {
 *   return Math::Vector2D(0.0, -0.1);  // Downward displacement
 * });
 *
 * Problem problem(u, v);
 * problem = /* elasticity bilinear form */
 *         + DirichletBC(u, g);  // Vector BC
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Component-Wise BC (Mixed Conditions)
 * @code{.cpp}
 * // Fix only x-component, leave y-component free
 * Problem problem(u, v);
 * problem = /* bilinear form */
 *         + DirichletBC(Component(u, 0), 0.0);  // u_x = 0
 *         // u_y is free (natural/Neumann BC)
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Time-Dependent BC
 * @code{.cpp}
 * // Time-varying boundary condition
 * double t = 0.0;  // Current time
 *
 * auto g = ScalarFunction([&t](const Math::SpatialPoint& x) {
 *   return std::sin(t) * x(0);  // Captures time variable
 * });
 *
 * // In time stepping loop:
 * for (t = 0; t < T; t += dt) {
 *   Problem problem(u, v);
 *   problem = /* form */
 *           + DirichletBC(u, g);  // BC depends on current t
 *
 *   problem.solve(solver);
 * }
 * @endcode
 *
 * @see PeriodicBC, Problem, TrialFunction
 * @see ScalarFunction, VectorFunction, Component
 */


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

      virtual bool isComponent() const = 0;

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
