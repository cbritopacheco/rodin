/**
 * @file RK4.h
 * @brief Fourth-order Runge-Kutta method (classical RK4).
 *
 * This file implements the classical fourth-order Runge-Kutta method,
 * one of the most widely used numerical integration methods for ODEs.
 */
#ifndef RODIN_MATH_RUNGEKUTTA_RK4_H
#define RODIN_MATH_RUNGEKUTTA_RK4_H

#include "Rodin/Types.h"

namespace Rodin::Math::RungeKutta
{
  /**
   * @brief Classical fourth-order Runge-Kutta method (RK4).
   *
   * This method provides fourth-order accurate numerical integration for
   * ordinary differential equations (ODEs) of the form:
   * @f[
   *   \frac{dq}{dt} = f(q, t)
   * @f]
   *
   * ## Algorithm
   * Given a time step @f$ \Delta t @f$ and current state @f$ q_n @f$, the
   * method computes @f$ q_{n+1} @f$ using four stages:
   * @f[
   *   \begin{aligned}
   *   k_1 &= f(q_n) \\
   *   k_2 &= f(q_n + 0.5 \Delta t \, k_1) \\
   *   k_3 &= f(q_n + 0.5 \Delta t \, k_2) \\
   *   k_4 &= f(q_n + \Delta t \, k_3) \\
   *   q_{n+1} &= q_n + \frac{\Delta t}{6} (k_1 + 2k_2 + 2k_3 + k_4)
   *   \end{aligned}
   * @f]
   *
   * ## Properties
   * - **Order of accuracy**: 4 (error @f$ \sim O(\Delta t^5) @f$ per step)
   * - **Stages**: 4 function evaluations per step
   * - **Stability**: Good stability properties for moderate time steps
   *
   * ## Advantages
   * - High accuracy with relatively few function evaluations
   * - Well-balanced between accuracy and computational cost
   * - Suitable for smooth solutions
   *
   * ## Typical Use Cases
   * - Celestial mechanics and orbital dynamics
   * - Chemical kinetics
   * - Time-dependent PDEs requiring high accuracy
   * - General-purpose ODE integration
   *
   * ## Example Usage
   * ```cpp
   * RK4 integrator;
   * Vector<Real> q = initial_state;
   * Real dt = 0.01;
   * auto f = [](const auto& state) { return compute_derivative(state); };
   * integrator.step(q, dt, q, f);
   * ```
   *
   * @note Thread-safe when used with different state vectors. Internal
   * stage vectors are thread-local.
   */
  struct RK4
  {
    /**
     * @brief Performs one step of the RK4 method.
     *
     * Advances the state by one time step using the classical four-stage
     * Runge-Kutta method.
     *
     * @tparam T Type of the state vector
     * @tparam G Type of the current state (typically same as T)
     * @tparam F Type of the function object computing @f$ f(q) @f$
     * @param[out] q State vector to update (output)
     * @param[in] dt Time step size @f$ \Delta t @f$
     * @param[in] p Current state (input)
     * @param[in] f Function computing the derivative @f$ f(q) @f$
     */
    template <class T, class G, class F>
    void step(T& q, Real dt, const G& p, const F& f) const
    {
      static thread_local T s_k1, s_k2, s_k3, s_k4;
      s_k1 = f(p);
      s_k2 = f(p + 0.5 * dt * s_k1);
      s_k3 = f(p + 0.5 * dt * s_k2);
      s_k4 = f(p + dt * s_k3);
      q = p + (dt / 6) * (s_k1 + 2 * s_k2 + 2 * s_k3 + s_k4);
    }
  };
}

#endif
