/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file H1a.h
 * @brief Regularized Hilbert-space extension operator.
 *
 * Provides @ref Rodin::Hilbert::H1a for smoothing linear functionals and
 * extending boundary shape derivatives into volume fields.
 */
#ifndef RODIN_HILBERT_H1A_H
#define RODIN_HILBERT_H1A_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Hilbert
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
      Math::LinearSystem<Math::SparseMatrix<ScalarType>, Math::Vector<ScalarType>>;

    /// Solver type (Conjugate Gradient)
    using SolverType = Solver::CG<LinearSystemType>;

    static_assert(std::is_same_v<FESRange, Real> ||
      std::is_same_v<FESRange, Math::SpatialVector<Real>>);

    public:
      /**
       * @brief Constructs an H1a extension operator.
       * @param fes Finite element space for the problem
       */
      H1a(const FES& fes)
        : m_fes(fes),
          m_trial(fes),
          m_test(fes),
          m_pb(m_trial, m_test),
          m_alpha(1)
      {}

      /**
       * @brief Copies parameters, boundary conditions and the owned solution.
       * @param other Extension to copy.
       *
       * The working problem is bound to this object's trial and test functions.
       * Its formulation and linear system are rebuilt on the next application.
       */
      H1a(const H1a& other)
        : m_fes(other.m_fes),
          m_trial(other.m_trial),
          m_test(other.m_test),
          m_pb(m_trial, m_test),
          m_alpha(other.m_alpha)
      {
        typename decltype(m_pb)::ProblemBodyType body;
        body.getDBCs() = other.m_pb.getBody().getDBCs();
        for (auto& condition : body.getDBCs())
          condition.setOperand(m_trial);
        m_pb = body;
      }

      /**
       * @brief Moves an extension while binding its working problem to its new owner.
       * @param other Extension to move.
       *
       * Boundary conditions are retained; the formulation and linear system are
       * rebuilt on the next application, as in copy construction.
       */
      H1a(H1a&& other)
        : m_fes(other.m_fes),
          m_trial(std::move(other.m_trial)),
          m_test(std::move(other.m_test)),
          m_pb(m_trial, m_test),
          m_alpha(other.m_alpha)
      {
        typename decltype(m_pb)::ProblemBodyType body;
        body.getDBCs() = other.m_pb.getBody().getDBCs();
        for (auto& condition : body.getDBCs())
          condition.setOperand(m_trial);
        m_pb = body;
      }

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
       *
       * Replaces the retained formulation and solves into the owned trial
       * function. This operation updates the extension's working state.
       * @tparam Differential Type of the differential form
       * @param lf Linear functional to extend
       * @return Solution of the regularized extension problem
       */
      template <class Differential>
      auto operator()(const Differential& lf)
      {
        auto body = [&]() {
          if constexpr (std::is_same_v<FESRange, Real>)
          {
            return Variational::Integral(m_alpha * m_alpha * Variational::Grad(m_trial),
                     Variational::Grad(m_test)) +
              Variational::Integral(m_trial, m_test) - lf;
          }
          else
          {
            return Variational::Integral(
                     m_alpha * m_alpha * Variational::Jacobian(m_trial),
                     Variational::Jacobian(m_test)) +
              Variational::Integral(m_trial, m_test) - lf;
          }
        }();
        body.getDBCs() = m_pb.getBody().getDBCs();
        m_pb = body;
        SolverType solver(m_pb);
        m_pb.solve(solver);
        return m_trial.getSolution();
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
      std::reference_wrapper<const FES> m_fes;    ///< Reference to finite element space
      Variational::TrialFunction<Solution, FES>
        m_trial;  ///< Trial function for the formulation
      Variational::TestFunction<FES> m_test;   ///< Test function for the formulation
      Variational::Problem<LinearSystemType, Variational::TrialFunction<Solution, FES>,
        Variational::TestFunction<FES>>
        m_pb; ///< Underlying variational problem
      Real m_alpha;        ///< Regularization parameter @f$ \alpha @f$
  };
}

#endif
