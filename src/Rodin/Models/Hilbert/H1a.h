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
   * @f[
   * a(u, v) := \alpha^2 \int_{\mathbb{R}^d} \nabla u : \nabla v \ dx +
   * \int_{\mathbb{R}^d} u \cdot v \ dx
   * @f]
   */
  template <class Solution, class FES>
  class H1a
  {
    using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;

    using FESRange = typename FormLanguage::Traits<FES>::RangeType;

    using LinearSystemType =
      Math::LinearSystem<Math::SparseMatrix<FESRange>, Math::Vector<FESRange>>;

    using SolverType = Solver::CG<LinearSystemType>;

    static_assert(std::is_same_v<FESRange, Real> ||
        std::is_same_v<FESRange, Math::Vector<Real>>);

    public:
      H1a(const FES& fes)
        : m_fes(fes),
          m_trial(fes), m_test(fes),
          m_pb(m_trial, m_test),
          m_alpha(1)
      {}

      H1a& setAlpha(Real alpha)
      {
        m_alpha = alpha;
        return *this;
      }

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

      H1a& operator+=(const Variational::DirichletBCBase<ScalarType>& dbc)
      {
        m_pb += dbc;
        return *this;
      }

      auto& getProblem()
      {
        return m_pb;
      }

      const auto& getProblem() const
      {
        return m_pb;
      }

      const auto& getTrialFunction() const
      {
        return m_trial;
      }

      const auto& getTestFunction() const
      {
        return m_test;
      }

      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

    private:
      std::reference_wrapper<const FES>           m_fes;
      Variational::TrialFunction<Solution, FES>   m_trial;
      Variational::TestFunction<FES>              m_test;
      Variational::Problem<
        FES, FES, Context::Local, Math::SparseMatrix<Real>, Math::Vector<Real>> m_pb;
      Real m_alpha;
      SolverType m_solver;
  };
}

#endif




