/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TrialFunction.h
 * @brief Trial functions (solution functions) for variational formulations.
 *
 * This file defines the TrialFunction class, which represents the trial or
 * solution functions in the weak formulation of partial differential equations.
 * Trial functions are the unknowns being solved for in the finite element method.
 */
#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include <functional>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class Solution, class FES>
  struct Traits<Variational::TrialFunction<Solution, FES>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TrialSpace;

    using SolutionType = Solution;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Reference wrapper for trial functions.
   *
   * TrialFunctionReference provides a lightweight reference to a TrialFunction
   * without owning the associated solution data. This is useful for creating
   * expressions involving trial functions without unnecessary copies.
   *
   * @tparam Solution Type of the solution grid function
   * @tparam FES Finite element space type
   *
   * @see TrialFunction
   */
  template <class Solution, class FES>
  class TrialFunctionReference
    : public ShapeFunction<TrialFunctionReference<Solution, FES>, FES, ShapeFunctionSpaceType::Trial>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;

      /// @brief Space type (trial space)
      static constexpr ShapeFunctionSpaceType SpaceType = ShapeFunctionSpaceType::Trial;

      /// @brief Parent class type
      using Parent =
        ShapeFunction<TrialFunctionReference<Solution, FESType>, FESType, SpaceType>;

      /**
       * @brief Constructs a reference to a trial function.
       * @param[in] ref Trial function to reference
       * @param[in] fes Finite element space
       */
      explicit
      TrialFunctionReference(const TrialFunction<Solution, FESType>& ref, const FESType& fes)
        : Parent(fes),
          m_ref(std::cref(ref))
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Reference to copy
       */
      TrialFunctionReference(const TrialFunctionReference& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Reference to move
       */
      TrialFunctionReference(TrialFunctionReference&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      /**
       * @brief Sets the evaluation point for the shape function.
       * @param[in] p Point at which to evaluate
       * @returns Reference to this object
       */
      TrialFunctionReference& setIntegrationPoint(const IntegrationPoint& p)
      {
        m_ref.get().setIntegrationPoint(p);
        return *this;
      }

      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_ref.get().getIntegrationPoint();
      }

      /**
       * @brief Gets the x-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_x @f$
       */
      constexpr
      auto x() const
      {
        return m_ref.get().x();
      }

      /**
       * @brief Gets the y-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_y @f$
       */
      constexpr
      auto y() const
      {
        return m_ref.get().y();
      }

      /**
       * @brief Gets the z-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_z @f$
       */
      constexpr
      auto z() const
      {
        return m_ref.get().z();
      }

      /**
       * @brief Gets the leaf (innermost) expression in the shape function tree.
       * @returns Reference to the underlying trial function
       */
      constexpr
      const auto& getLeaf() const
      {
        return m_ref.get().getLeaf();
      }

      /**
       * @brief Gets the solution grid function (mutable).
       * @returns Reference to solution
       */
      constexpr
      auto& getSolution()
      {
        return m_ref.get().getSolution();
      }

      /**
       * @brief Gets the basis function at the specified local index.
       * @param[in] local Local index of the basis function
       * @returns Basis function expression
       */
      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return m_ref.get().getBasis(local);
      }

      /**
       * @brief Gets the solution grid function (const).
       * @returns Const reference to solution
       */
      constexpr
      const auto& getSolution() const
      {
        return m_ref.get().getSolution();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_ref.get().getOrder(poly);
      }

      TrialFunctionReference& setName(const std::string& name)
      {
        m_ref.get().setName(name);
        return *this;
      }

      Optional<StringView> getName() const override
      {
        return m_ref.get().getName();
      }

      /**
       * @brief Creates a copy of this reference.
       * @returns Pointer to newly allocated copy
       */
      TrialFunctionReference* copy() const noexcept final override
      {
        return new TrialFunctionReference(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<Solution, FESType>> m_ref;
  };

  /**
   * @ingroup RodinVariational
   * @brief Trial function (unknown solution) in a finite element space.
   *
   * In the weak formulation of a partial differential equation, the trial function
   * represents the unknown solution being sought. Given a finite element space 
   * @f$ U_h @f$, the trial function @f$ u \in U_h @f$ is the function to be 
   * determined by solving the variational problem.
   *
   * ## Mathematical Foundation
   * The weak formulation seeks @f$ u \in U_h @f$ such that:
   * @f[
   *   a(u, v) = l(v) \quad \forall v \in V_h
   * @f]
   * where:
   * - @f$ u @f$ is the trial function (unknown)
   * - @f$ v @f$ is a test function
   * - @f$ a(\cdot, \cdot) @f$ is a bilinear form
   * - @f$ l(\cdot) @f$ is a linear functional
   * - @f$ U_h @f$ is the trial space
   *
   * The TrialFunction internally manages a GridFunction that stores the degrees
   * of freedom of the solution after the problem is solved.
   *
   * ## Usage Example
   * ```cpp
   * P1 Uh(mesh);                    // Define finite element space
   * TrialFunction u(Uh);             // Define trial function
   * TestFunction v(Uh);              // Define test function
   * 
   * // Define problem
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
   * 
   * // Solve
   * problem.solve(solver);
   * 
   * // Access solution
   * auto& solution = u.getSolution();
   * ```
   *
   * @tparam Solution Type of the solution grid function
   * @tparam FES Finite element space type
   *
   * @see TestFunction, GridFunction, ShapeFunction
   */
  template <class Solution, class FES>
  class TrialFunction : public TrialFunctionReference<Solution, FES>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;

      /// @brief Space type (trial space)
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      /// @brief Solution grid function type
      using SolutionType = Solution;

      /// @brief Parent class type
      using Parent = TrialFunctionReference<SolutionType, FESType>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>);

      /**
       * @brief Constructs a trial function in the given finite element space.
       * @param[in] fes Finite element space
       *
       * Creates a trial function with an associated solution grid function
       * initialized to zero.
       */
      constexpr
      TrialFunction(const FES& fes)
        : Parent(*this, fes),
          m_gf(fes)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Trial function to copy
       */
      constexpr
      TrialFunction(const TrialFunction& other)
        : Parent(other),
          m_gf(other.m_gf)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Trial function to move
       */
      constexpr
      TrialFunction(TrialFunction&& other)
        : Parent(std::move(other)),
          m_gf(std::move(other.m_gf))
      {}

      /// @brief Copy assignment is deleted
      void operator=(const TrialFunction&) = delete;

      /// @brief Move assignment is deleted
      void operator=(TrialFunction&&) = delete;

      /**
       * @brief Gets the x-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_x @f$
       */
      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      /**
       * @brief Gets the y-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_y @f$
       */
      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      /**
       * @brief Gets the z-component of a vector-valued trial function.
       * @returns Component expression representing @f$ u_z @f$
       */
      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      /**
       * @brief Gets the leaf (innermost) expression in the shape function tree.
       * @returns Reference to this trial function
       */
      constexpr
      const TrialFunction& getLeaf() const
      {
        return *this;
      }

      /**
       * @brief Gets the solution grid function (mutable).
       * @returns Reference to the solution
       *
       * The solution grid function stores the degrees of freedom of the 
       * computed solution after solving the variational problem.
       */
      constexpr
      SolutionType& getSolution()
      {
        return m_gf;
      }

      /**
       * @brief Gets the solution grid function (const).
       * @returns Const reference to the solution
       */
      constexpr
      const SolutionType& getSolution() const
      {
        return m_gf;
      }

      TrialFunction& setName(const std::string& name)
      {
        m_name = name;
        return *this;
      }

      Optional<StringView> getName() const override
      {
        return m_name;
      }

    private:
      std::string m_name; ///< Optional name for the trial function
      SolutionType m_gf;
  };

  /**
   * @brief Deduction guide for TrialFunction.
   *
   * Allows automatic type deduction when constructing a TrialFunction:
   * ```cpp
   * P1 Uh(mesh);
   * TrialFunction u(Uh);  // Type deduced automatically
   * ```
   */
  template <class FES>
  TrialFunction(const FES&)
    -> TrialFunction<
        GridFunction<FES, Math::Vector<
          typename FormLanguage::Traits<FES>::ScalarType>>, FES>;
}
#endif

