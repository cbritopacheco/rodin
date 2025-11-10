/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TestFunction.h
 * @brief Test functions for variational formulations.
 *
 * This file defines the TestFunction class, which represents the test functions
 * (also known as weighting functions) in the weak formulation of partial 
 * differential equations using the finite element method.
 */
#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

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
   * @brief Test function in a finite element space.
   *
   * In the weak formulation of a partial differential equation, test functions
   * (also called weighting functions) are used to convert the PDE into a 
   * variational problem. Given a finite element space @f$ V_h @f$, the test
   * function @f$ v \in V_h @f$ is chosen from this space.
   *
   * ## Mathematical Foundation
   * The weak formulation seeks @f$ u \in U_h @f$ such that:
   * @f[
   *   a(u, v) = l(v) \quad \forall v \in V_h
   * @f]
   * where:
   * - @f$ v @f$ is the test function
   * - @f$ a(\cdot, \cdot) @f$ is a bilinear form
   * - @f$ l(\cdot) @f$ is a linear functional
   * - @f$ V_h @f$ is the test space (typically equal to the trial space @f$ U_h @f$)
   *
   * ## Usage Example
   * ```cpp
   * P1 Vh(mesh);                    // Define finite element space
   * TestFunction v(Vh);              // Define test function
   * 
   * // Use in variational formulation
   * auto bilinearForm = Integral(Grad(u), Grad(v));
   * auto linearForm = Integral(f, v);
   * ```
   *
   * @tparam FES Finite element space type
   *
   * @see TrialFunction, ShapeFunction
   */
  template <class FES>
  class TestFunction
    : public ShapeFunction<TestFunction<FES>, FES, TestSpace>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;
      
      /// @brief Space type (test space)
      static constexpr ShapeFunctionSpaceType Space = TestSpace;

      /// @brief Parent class type
      using Parent = ShapeFunction<TestFunction<FESType>, FESType, TestSpace>;

      /**
       * @brief Constructs a test function in the given finite element space.
       * @param[in] fes Finite element space
       */
      constexpr
      TestFunction(const FES& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Test function to copy
       */
      constexpr
      TestFunction(const TestFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Test function to move
       */
      constexpr
      TestFunction(TestFunction&& other)
        : Parent(std::move(other))
      {}

      /// @brief Copy assignment is deleted
      void operator=(const TestFunction&) = delete;

      /// @brief Move assignment is deleted
      void operator=(TestFunction&&) = delete;

      /**
       * @brief Gets the x-component of a vector-valued test function.
       * @returns Component expression representing @f$ v_x @f$
       */
      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      /**
       * @brief Gets the y-component of a vector-valued test function.
       * @returns Component expression representing @f$ v_y @f$
       */
      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      /**
       * @brief Gets the z-component of a vector-valued test function.
       * @returns Component expression representing @f$ v_z @f$
       */
      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      /**
       * @brief Gets the leaf (innermost) expression in the shape function tree.
       * @returns Reference to this test function
       */
      constexpr
      const TestFunction& getLeaf() const
      {
        return *this;
      }

      /**
       * @brief Creates a copy of this test function.
       * @returns Pointer to newly allocated copy
       */
      TestFunction* copy() const noexcept override
      {
        return new TestFunction(*this);
      }
  };

  /**
   * @brief Deduction guide for TestFunction.
   *
   * Allows automatic type deduction when constructing a TestFunction:
   * ```cpp
   * P1 Vh(mesh);
   * TestFunction v(Vh);  // Type deduced as TestFunction<P1>
   * ```
   */
  template <class FES>
  TestFunction(const FES&) -> TestFunction<FES>;
}

#endif