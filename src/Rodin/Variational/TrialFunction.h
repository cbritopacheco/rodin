/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

/**
 * @file
 * @brief Trial function for variational formulations.
 */

#include <functional>

#include "Rodin/Variational/ForwardDecls.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Traits specialization for TrialFunction.
   */
  template <class Solution, class FES>
  struct Traits<Variational::TrialFunction<Solution, FES>>
  {
    /// @brief Finite element space type
    using FESType = FES;
    /// @brief Space type identifier (Trial space)
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TrialSpace;

    /// @brief Solution storage type
    using SolutionType = Solution;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Reference wrapper for trial functions in variational formulations.
   *
   * TrialFunctionReference provides a lightweight reference wrapper around
   * TrialFunction objects, enabling efficient pass-by-reference semantics
   * in variational form expressions without unnecessary copying.
   *
   * @tparam Solution Type of the solution (typically GridFunction)
   * @tparam FES Type of finite element space
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
      
      /// @brief Trial space identifier
      static constexpr ShapeFunctionSpaceType SpaceType = ShapeFunctionSpaceType::Trial;

      /// @brief Parent class type
      using Parent =
        ShapeFunction<TrialFunctionReference<Solution, FESType>, FESType, SpaceType>;

      /**
       * @brief Constructs reference to a trial function.
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
       */
      TrialFunctionReference(const TrialFunctionReference& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      /**
       * @brief Move constructor.
       */
      TrialFunctionReference(TrialFunctionReference&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      /**
       * @brief Sets the evaluation point.
       * @param[in] p Point at which to evaluate
       * @returns Reference to this instance
       */
      TrialFunctionReference& setPoint(const Geometry::Point& p)
      {
        m_ref.get().setPoint(p);
        return *this;
      }

      /**
       * @brief Gets the x-component.
       * @returns Component expression for x-coordinate
       */
      constexpr
      auto x() const
      {
        return m_ref.get().x();
      }

      /**
       * @brief Gets the y-component.
       * @returns Component expression for y-coordinate
       */
      constexpr
      auto y() const
      {
        return m_ref.get().y();
      }

      /**
       * @brief Gets the z-component.
       * @returns Component expression for z-coordinate
       */
      constexpr
      auto z() const
      {
        return m_ref.get().z();
      }

      /**
       * @brief Gets the underlying trial function.
       * @returns Reference to the trial function
       */
      constexpr
      const auto& getLeaf() const
      {
        return m_ref.get().getLeaf();
      }

      /**
       * @brief Gets the solution object.
       * @returns Reference to the solution
       */
      constexpr
      auto& getSolution()
      {
        return m_ref.get().getSolution();
      }

      /**
       * @brief Gets the basis function value at a local DOF.
       * @param[in] local Local degree of freedom index
       * @returns Basis function value
       */
      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return m_ref.get().getBasis(local);
      }

      /**
       * @brief Gets the solution object (const version).
       * @returns Const reference to the solution
       */
      constexpr
      const auto& getSolution() const
      {
        return m_ref.get().getSolution();
      }

      /**
       * @brief Creates a copy of this reference.
       * @returns Pointer to the copy
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
   * @brief Represents a trial function in variational formulations.
   *
   * A TrialFunction represents the unknown solution @f$ u @f$ in a variational
   * problem. In finite element analysis, trial functions are elements of the
   * approximation space @f$ V_h @f$ where the solution is sought. They appear
   * as the first argument in bilinear forms @f$ a(u,v) @f$ and provide the
   * degrees of freedom that will be determined by solving the discrete system.
   *
   * ## Mathematical Foundation
   *
   * In the weak formulation: Find @f$ u \in V_h @f$ such that
   * @f[
   *   a(u,v) = l(v) \quad \forall v \in V_h
   * @f]
   * the TrialFunction represents @f$ u @f$, while TestFunction represents @f$ v @f$.
   *
   * The discrete representation is:
   * @f[
   *   u_h(x) = \sum_{j=1}^N u_j \phi_j(x)
   * @f]
   * where @f$ \phi_j @f$ are the trial basis functions and @f$ u_j @f$ are the
   * unknown coefficients stored in the solution vector.
   *
   * ## Usage Example
   * ```cpp
   * P1 Vh(mesh);
   * TrialFunction u(Vh);  // Unknown solution
   * TestFunction v(Vh);   // Test function
   * 
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
   * problem.solve(solver);
   * 
   * // Access solution
   * const auto& solution = u.getSolution();
   * ```
   *
   * @tparam Solution Type of solution storage (typically GridFunction)
   * @tparam FES Type of finite element space
   *
   * @see TestFunction, GridFunction, Problem
   */
  template <class Solution, class FES>
  class TrialFunction : public TrialFunctionReference<Solution, FES>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;
      
      /// @brief Trial space identifier
      static constexpr ShapeFunctionSpaceType Space = TrialSpace;

      /// @brief Solution storage type
      using SolutionType = Solution;

      /// @brief Parent class type
      using Parent = TrialFunctionReference<SolutionType, FESType>;

      static_assert(std::is_base_of_v<FiniteElementSpaceBase, FES>);

      /**
       * @brief Constructs a trial function on a finite element space.
       * @param[in] fes Finite element space on which to define the trial function
       */
      constexpr
      TrialFunction(const FES& fes)
        : Parent(*this, fes),
          m_gf(fes)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      TrialFunction(const TrialFunction& other)
        : Parent(other),
          m_gf(other.m_gf)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr
      TrialFunction(TrialFunction&& other)
        : Parent(std::move(other)),
          m_gf(std::move(other.m_gf))
      {}

      /**
       * @brief Copy assignment is deleted.
       */
      void operator=(const TrialFunction&) = delete;

      /**
       * @brief Move assignment is deleted.
       */
      void operator=(TrialFunction&&) = delete;

      /**
       * @brief Gets the x-component of a vector trial function.
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
       * @brief Gets the y-component of a vector trial function.
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
       * @brief Gets the z-component of a vector trial function.
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
       * @brief Gets the underlying trial function (returns itself).
       * @returns Reference to this trial function
       */
      constexpr
      const TrialFunction& getLeaf() const
      {
        return *this;
      }

      /**
       * @brief Gets the solution object containing DOF values.
       * @returns Reference to the solution storage
       *
       * After solving, this contains the computed coefficients @f$ u_j @f$.
       */
      constexpr
      SolutionType& getSolution()
      {
        return m_gf;
      }

      /**
       * @brief Gets the solution object (const version).
       * @returns Const reference to the solution storage
       */
      constexpr
      const SolutionType& getSolution() const
      {
        return m_gf;
      }

    private:
      SolutionType m_gf;
  };

  template <class FES>
  TrialFunction(const FES&)
    -> TrialFunction<
        GridFunction<FES, Math::Vector<
          typename FormLanguage::Traits<FES>::ScalarType>>, FES>;
}
#endif

