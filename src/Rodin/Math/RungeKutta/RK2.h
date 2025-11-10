/**
 * @file RK2.h
 * @brief Second-order Runge-Kutta method (midpoint method).
 *
 * This file implements the classical second-order Runge-Kutta method,
 * also known as the midpoint method or modified Euler method.
 */
#ifndef RODIN_MATH_RUNGEKUTTA_RK2_H
#define RODIN_MATH_RUNGEKUTTA_RK2_H

#include "Rodin/Types.h"

namespace Rodin::Math::RungeKutta
{
  /**
   * @brief Classical second-order Runge-Kutta method with midpoint rule.
   *
   * This method provides second-order accurate numerical integration for
   * ordinary differential equations (ODEs) of the form:
   * @f[
   *   \frac{dq}{dt} = f(q, t)
   * @f]
   *
   * ## Algorithm
   * Given a time step @f$ \Delta t @f$ and current state @f$ q_n @f$, the
   * method computes @f$ q_{n+1} @f$ using:
   * @f[
   *   \begin{align}
   *   k_1 &= f(q_n) \\
   *   k_2 &= f(q_n + 0.5 \Delta t \, k_1) \\
   *   q_{n+1} &= q_n + \Delta t \, k_2
   *   \end{align}
   * @f]
   *
   * ## Properties
   * - **Order of accuracy**: 2 (error @f$ \sim O(\Delta t^3) @f$ per step)
   * - **Stages**: 2 function evaluations per step
   * - **Stability**: A-stable for linear problems
   *
   * ## Typical Use Cases
   * - Time-dependent PDEs (heat equation, wave equation)
   * - Particle dynamics
   * - Systems requiring better accuracy than forward Euler
   *
   * ## Example Usage
   * ```cpp
   * RK2 integrator;
   * Vector<Real> q = initial_state;
   * Real dt = 0.01;
   * auto f = [](const auto& state) { return compute_derivative(state); };
   * integrator.step(q, dt, q, f);
   * ```
   */
  struct RK2
  {
    /**
     * @brief Performs one step of the RK2 method.
     *
     * Advances the state @f$ q @f$ by one time step using the RK2 scheme.
     *
     * @tparam T Type of the state vector
     * @tparam G Type of the current state (typically same as T)
     * @tparam F Type of the function object computing @f$ f(q) @f$
     * @param[in,out] q State vector to update (output)
     * @param[in] dt Time step size @f$ \Delta t @f$
     * @param[in] p Current state (input)
     * @param[in] f Function computing the derivative @f$ f(q) @f$
     */
    template <class T, class G, class F>
    void step(T& q, Real dt, const G& p, const F& f) const
    {
      static thread_local T s_k1, s_k2;
      s_k1 = f(p);
      s_k2 = f(p + Real(0.5) * dt * s_k1);
      q += dt * s_k2;
    }
  };
}

#endif
