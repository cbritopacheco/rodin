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
  template <class FES>
  class H1a
  {
    using FESRange = typename FormLanguage::Traits<FES>::RangeType;
    static_assert(std::is_same_v<FESRange, Scalar> || std::is_same_v<FESRange, Math::Vector<Scalar>>);

    public:
      H1a(const FES& fes)
        : m_fes(fes),
          m_trial(fes), m_test(fes),
          m_pb(m_trial, m_test),
          m_alpha(1)
      {}

      H1a& setAlpha(Scalar alpha)
      {
        m_alpha = alpha;
        return *this;
      }

      template <class Differential>
      auto operator()(const Differential& lf) const
      {
        const auto& fes = getFiniteElementSpace();
        if constexpr (std::is_same_v<FESRange, Scalar>)
        {
          Variational::TrialFunction g(fes);
          Variational::TestFunction  w(fes);
          Variational::Problem hilbert(g, w);
          hilbert = Integral(m_alpha * m_alpha * Grad(g), Grad(w))
                  + Integral(g, w)
                  - lf;
          hilbert.solve(m_cg);
          return g.getSolution();
        }
        else if constexpr (std::is_same_v<FESRange, Math::Vector<Scalar>>)
        {
          Variational::TrialFunction g(fes);
          Variational::TestFunction  w(fes);
          Variational::Problem hilbert(g, w);
          hilbert = Integral(m_alpha * m_alpha * Jacobian(g), Jacobian(w))
                  + Integral(g, w)
                  - lf;
          hilbert.solve(m_cg);
          return g.getSolution();
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline
      H1a& operator+=(const Variational::DirichletBCBase& dbc)
      {
        m_pb += dbc;
        return *this;
      }

      inline
      auto& getProblem()
      {
        return m_pb;
      }

      inline
      const auto& getProblem() const
      {
        return m_pb;
      }

      inline
      const auto& getTrialFunction() const
      {
        return m_trial;
      }

      inline
      const auto& getTestFunction() const
      {
        return m_test;
      }

      inline
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

    private:
      std::reference_wrapper<const FES>   m_fes;
      Variational::TrialFunction<FES>     m_trial;
      Variational::TestFunction<FES>      m_test;
      Variational::Problem<FES, FES, Context::Sequential, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> m_pb;
      Scalar m_alpha;
      Solver::CG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> m_cg;
  };
}

#endif




