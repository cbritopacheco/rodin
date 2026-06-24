/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalVirtualWork.h
 * @brief Façade for the internal virtual work nonlinear form in hyperelastic formulations.
 *
 * Provides the unified API for constructing the internal virtual work
 * contributions to a nonlinear solid mechanics problem:
 *
 * @code
 * auto ivw = Solid::InternalVirtualWork(law, displacement);
 *
 * // Residual contribution (linear form):
 * LinearForm F(fes);
 * F += ivw.Residual(v);
 *
 * // Tangent contribution (bilinear form):
 * BilinearForm K(trialFES, testFES);
 * K += ivw.Tangent(u, v);
 *
 * // Both at once for problem definition:
 * Problem newton(u, v);
 * newton = ivw(u, v) - Integral(f, v) + DirichletBC(u, g).on(bc);
 * @endcode
 *
 * @par Sign convention with Problem and NewtonSolver
 * @c Residual(v) returns the internal virtual work with its physical sign,
 * @f[
 *   R(\mathbf{u}; \mathbf{v})
 *     = \int_{\Omega_0} \mathbf{P}(\mathbf{u}) : \nabla_0\mathbf{v}\,dX.
 * @f]
 * @c Tangent(du, v) returns its positive Gateaux derivative @f$K = DR@f$.
 * When a @c Problem assembles a @c ProblemBody into a linear system, linear
 * form contributions are moved to the right-hand side with a minus sign.
 * Therefore
 * @code
 * Problem newton(du, v);
 * newton = ivw(du, v) + DirichletBC(du, zero);
 * @endcode
 * assembles the Newton system
 * @f[
 *   K(\mathbf{u}^k)\,\Delta\mathbf{u} = -R(\mathbf{u}^k),
 * @f]
 * and @c NewtonSolver then applies @f$\mathbf{u}^{k+1}
 * = \mathbf{u}^k + \Delta\mathbf{u}@f$. Do not write
 * @c Tangent(du, v) - Residual(v) for this Newton form: the @c Problem
 * assembly will negate the already-negated residual and produce
 * @f$K\Delta\mathbf{u}=+R@f$.
 *
 * All three forms share the same constitutive law, displacement state,
 * quadrature order, and input hook — ensuring consistency of the Newton
 * linearization.
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORK_H
#define RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORK_H

#include <functional>
#include <type_traits>

#include "Rodin/Variational/NonLinearFormIntegrator.h"
#include "Rodin/Solid/Local/Input.h"

#include "InternalVirtualWorkResidual.h"
#include "InternalVirtualWorkTangent.h"

namespace Rodin::Solid
{
  /**
   * @brief Façade for the internal virtual work nonlinear form.
   *
   * Owns the constitutive law, a reference to the current displacement, and
   * optional configuration (quadrature order, auxiliary input). Exposes
   * `.Residual(v)`, `.Tangent(u, v)`, and `operator()(u, v)` (which combines
   * both into a ProblemBody for direct use in a Problem assignment).
   *
   * @note The combined @c operator()(u, v) intentionally returns
   * @c Tangent(u, v) + Residual(v). @c Problem assembly places linear forms on
   * the right-hand side with a minus sign, so this is the correct convention
   * for @c NewtonSolver.
   *
   * @tparam Law  The hyperelastic constitutive law type
   * @tparam State The displacement grid-function type
   */
  template <class Law, class State>
  class InternalVirtualWork
    : public Variational::NonLinearFormIntegrator<InternalVirtualWork<Law, State>>
  {
    public:
      using LawType   = Law;
      using StateType = State;

      /**
       * @brief Constructs the internal virtual work object.
       * @param law The constitutive law (stored by value)
       * @param displacement The current displacement (stored by reference)
       */
      InternalVirtualWork(const Law& law, const State& displacement)
        : m_law(law),
          m_displacement(displacement),
          m_quadOrder(0)
      {}

      /// @brief Rebinds the current displacement state.
      InternalVirtualWork& setDisplacement(const State& displacement)
      {
        m_displacement = std::cref(displacement);
        return *this;
      }

      /**
       * @brief Sets the quadrature order for both residual and tangent.
       * @param order Polynomial order (0 = auto: 2 * max(fe order))
       */
      InternalVirtualWork& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /**
       * @brief Sets the auxiliary input for both residual and tangent.
       *
       * The input is propagated to every quadrature point in both
       * @c Residual() and @c Tangent() integrators, ensuring they see the
       * same fiber directions, activation parameters, etc.
       */
      InternalVirtualWork& setInput(InputFunction input)
      {
        m_input = std::move(input);
        return *this;
      }

      /// @brief Gets the constitutive law.
      const Law& getLaw() const { return m_law; }

      /**
       * @brief Returns the residual integrator (linear form contribution).
       *
       * Produces:
       * @f[
       *   \delta W^{\text{int}}(\mathbf{v})
       *     = \int_{\Omega_0} \mathbf{P}(\mathbf{u}) : \nabla_0 \mathbf{v} \, dX
       * @f]
       */
      template <class Test>
      auto Residual(const Test& v) const
      {
        InternalVirtualWorkResidual integrator(m_law, v, m_displacement.get());
        if (m_quadOrder > 0)
          integrator.setQuadratureOrder(m_quadOrder);
        if (m_input)
          integrator.setInput(m_input);
        return integrator;
      }

      /**
       * @brief Returns the tangent integrator (bilinear form contribution).
       *
       * Produces:
       * @f[
       *   D(\delta W^{\text{int}})[\Delta\mathbf{u}, \mathbf{v}]
       *     = \int_{\Omega_0} \mathbf{A}(\mathbf{F}) : \nabla_0 \Delta\mathbf{u}
       *       : \nabla_0 \mathbf{v} \, dX
       * @f]
       * where @f$ \mathbf{A} = \partial\mathbf{P}/\partial\mathbf{F} @f$.
       */
      template <class Trial, class Test>
      auto Tangent(const Trial& u, const Test& v) const
      {
        InternalVirtualWorkTangent integrator(m_law, u, v, m_displacement.get());
        if (m_quadOrder > 0)
          integrator.setQuadratureOrder(m_quadOrder);
        if (m_input)
          integrator.setInput(m_input);
        return integrator;
      }

      // operator()(u, v) → ProblemBody{Tangent, Residual} is provided by the
      // NonLinearFormIntegrator CRTP base.

    private:
      Law m_law;
      std::reference_wrapper<const State> m_displacement;
      size_t m_quadOrder;
      InputFunction m_input;
  };

  /// CTAD deduction guide for InternalVirtualWork
  template <class Law, class State>
  InternalVirtualWork(const Law&, const State&)
    -> InternalVirtualWork<std::decay_t<Law>, std::decay_t<State>>;
}

#endif
