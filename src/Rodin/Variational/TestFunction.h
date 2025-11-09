/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

/**
 * @file
 * @brief Test function for variational formulations.
 */

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<Variational::TestFunction<FES>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TestSpace;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Represents a test function in variational formulations.
   *
   * A TestFunction represents the test function @f$ v @f$ in a variational problem.
   * In finite element analysis, test functions are elements of the test space
   * @f$ V_h @f$ used to convert a PDE into its weak form. They appear as the
   * second argument in bilinear forms @f$ a(u,v) @f$ and as the argument in
   * linear forms @f$ l(v) @f$.
   *
   * ## Mathematical Foundation
   *
   * In the weak formulation: Find @f$ u \in V_h @f$ such that
   * @f[
   *   a(u,v) = l(v) \quad \forall v \in V_h
   * @f]
   * the TestFunction represents @f$ v @f$, while TrialFunction represents @f$ u @f$.
   *
   * Test functions are used to build the discrete system @f$ Au = b @f$ where:
   * - @f$ A_{ij} = a(\phi_j, \psi_i) @f$ (bilinear form)
   * - @f$ b_i = l(\psi_i) @f$ (linear form)
   * 
   * Here @f$ \psi_i @f$ are the test basis functions.
   *
   * ## Usage Example
   * ```cpp
   * P1 Vh(mesh);
   * TrialFunction u(Vh);  // Unknown solution
   * TestFunction v(Vh);   // Test function
   * 
   * // Poisson equation: -Δu = f
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v))  // a(u,v)
   *         - Integral(f, v)               // l(v)
   *         + DirichletBC(u, Zero());
   * problem.solve(solver);
   * ```
   *
   * ## Galerkin Method
   * In the Galerkin method, the test and trial spaces are the same (@f$ V_h @f$),
   * which is the most common choice in finite element methods. This leads to
   * symmetric matrices for self-adjoint operators.
   *
   * @tparam FES Type of finite element space
   *
   * @see TrialFunction, BilinearForm, LinearForm, Problem
   */
  template <class FES>
  class TestFunction
    : public ShapeFunction<TestFunction<FES>, FES, TestSpace>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;
      
      /// @brief Test space identifier
      static constexpr ShapeFunctionSpaceType Space = TestSpace;

      /// @brief Parent class type
      using Parent = ShapeFunction<TestFunction<FESType>, FESType, TestSpace>;

      /**
       * @brief Constructs a test function on a finite element space.
       * @param[in] fes Finite element space on which to define the test function
       */
      constexpr
      TestFunction(const FES& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      TestFunction(const TestFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr
      TestFunction(TestFunction&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Copy assignment is deleted.
       */
      void operator=(const TestFunction&) = delete;

      /**
       * @brief Move assignment is deleted.
       */
      void operator=(TestFunction&&) = delete;

      /**
       * @brief Gets the x-component of a vector test function.
       * @returns Component expression for x-coordinate
       *
       * Requires that the finite element space has vector dimension @f$ \geq 1 @f$.
       */
      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      /**
       * @brief Gets the y-component of a vector test function.
       * @returns Component expression for y-coordinate
       *
       * Requires that the finite element space has vector dimension @f$ \geq 2 @f$.
       */
      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      /**
       * @brief Gets the z-component of a vector test function.
       * @returns Component expression for z-coordinate
       *
       * Requires that the finite element space has vector dimension @f$ \geq 3 @f$.
       */
      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      /**
       * @brief Gets the underlying test function (returns itself).
       * @returns Reference to this test function
       */
      constexpr
      const TestFunction& getLeaf() const
      {
        return *this;
      }

      /**
       * @brief Creates a copy of this test function.
       * @returns Pointer to the copy
       */
      TestFunction* copy() const noexcept override
      {
        return new TestFunction(*this);
      }
  };

  template <class FES>
  TestFunction(const FES&) -> TestFunction<FES>;
}

#endif
