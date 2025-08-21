/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_H1A_H
#define RODIN_MODELS_DISTANCE_H1A_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Models::Hilbert
{
  /**
   * @brief Hilbertian @f$ H^1 @f$ extension-regularization procedure with
   * a regularization length-scale parameter.
   *
   * This class implements a Hilbert space extension operator that solves 
   * regularized problems of the form:
   * @f[
   * a(u, v) := \alpha^2 \int_{\mathbb{R}^d} \nabla u : \nabla v \, dx +
   * \int_{\mathbb{R}^d} u \cdot v \, dx
   * @f]
   * where @f$ \alpha @f$ is the regularization parameter controlling the 
   * length scale of the extension.
   *
   * @tparam Solution Solution type
   * @tparam FES Finite element space type
   */
  template <class Solution, class FES>
  class H1a
  {
    /// Scalar type from the finite element space
    using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;

    /// Range type from the finite element space  
    using FESRange = typename FormLanguage::Traits<FES>::RangeType;

    /// Linear system type for the regularized problem
    using LinearSystemType =
      Math::LinearSystem<Math::SparseMatrix<FESRange>, Math::Vector<FESRange>>;

    /// Solver type (Conjugate Gradient)
    using SolverType = Solver::CG<LinearSystemType>;

    static_assert(std::is_same_v<FESRange, Real> ||
        std::is_same_v<FESRange, Math::Vector<Real>>);

    public:
      /**
       * @brief Constructs an H1a extension operator.
       * @param fes Finite element space for the problem
       */
      H1a(const FES& fes)
        : m_fes(fes),
          m_trial(fes), m_test(fes),
          m_pb(m_trial, m_test),
          m_alpha(1)
      {}

      /**
       * @brief Sets the regularization parameter @f$ \alpha @f$.
       * @param alpha Regularization length-scale parameter
       * @return Reference to this object for method chaining
       */
      H1a& setAlpha(Real alpha)
      {
        m_alpha = alpha;
        return *this;
      }

      /**
       * @brief Applies the extension operator to a differential form.
       * @tparam Differential Type of the differential form
       * @param lf Linear functional to extend
       * @return Solution of the regularized extension problem
       */
      template <class Differential>
      auto operator()(const Differential& lf) const
      {
        const auto& fes = getFiniteElementSpace();
        if constexpr (std::is_same_v<FESRange, Real>)
        {
          Variational::TrialFunction g(fes);
          Variational::TestFunction  w(fes);
          Variational::Problem hilbert(g, w);
          hilbert = Integral(m_alpha * m_alpha * Grad(g), Grad(w))
                  + Integral(g, w)
                  - lf;
          hilbert.solve(m_solver);
          return g.getSolution();
        }
        else if constexpr (std::is_same_v<FESRange, Math::Vector<Real>>)
        {
          Variational::TrialFunction g(fes);
          Variational::TestFunction  w(fes);
          Variational::Problem hilbert(g, w);
          hilbert = Integral(m_alpha * m_alpha * Jacobian(g), Jacobian(w))
                  + Integral(g, w)
                  - lf;
          hilbert.solve(m_solver);
          return g.getSolution();
        }
        else
        {
          assert(false);
          return void();
        }
      }

      /**
       * @brief Adds a Dirichlet boundary condition to the problem.
       * @param dbc Dirichlet boundary condition to add
       * @return Reference to this object for method chaining
       */
      H1a& operator+=(const Variational::DirichletBCBase<ScalarType>& dbc)
      {
        m_pb += dbc;
        return *this;
      }

      /**
       * @brief Gets the underlying variational problem.
       * @return Reference to the variational problem
       */
      auto& getProblem()
      {
        return m_pb;
      }

      /**
       * @brief Gets the underlying variational problem.
       * @return Const reference to the variational problem
       */
      const auto& getProblem() const
      {
        return m_pb;
      }

      /**
       * @brief Gets the trial function used in the formulation.
       * @return Const reference to the trial function
       */
      const auto& getTrialFunction() const
      {
        return m_trial;
      }

      /**
       * @brief Gets the test function used in the formulation.
       * @return Const reference to the test function
       */
      const auto& getTestFunction() const
      {
        return m_test;
      }

      /**
       * @brief Gets the finite element space.
       * @return Const reference to the finite element space
       */
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

    private:
      std::reference_wrapper<const FES>           m_fes;    ///< Reference to finite element space
      Variational::TrialFunction<Solution, FES>   m_trial;  ///< Trial function for the formulation
      Variational::TestFunction<FES>              m_test;   ///< Test function for the formulation
      Variational::Problem<
        FES, FES, Context::Local, Math::SparseMatrix<Real>, Math::Vector<Real>> m_pb; ///< Underlying variational problem
      Real m_alpha;        ///< Regularization parameter @f$ \alpha @f$
      SolverType m_solver; ///< Conjugate gradient solver
  };
}

#endif




