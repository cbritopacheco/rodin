#ifndef RODIN_MATH_RUNGEKUTTA_RK4_H
#define RODIN_MATH_RUNGEKUTTA_RK4_H

#include "Rodin/Types.h"

namespace Rodin::Math::RungeKutta
{
  struct RK4
  {
    template <class T, class G, class F>
    void step(T& q, Real dt, const G& p, const F& f) const
    {
      static thread_local T s_k1, s_k2, s_k3, s_k4;
      s_k1 = f(p);
      s_k2 = f(p + 0.5 * dt * s_k1);
      s_k3 = f(p + 0.5 * dt * s_k2);
      s_k4 = f(p + dt * s_k3);
      q += (dt / 6) * (s_k1 + 2 * s_k2 + 2 * s_k3 + s_k4);
    }

    // -------- Dense output --------
    template <class T, class F>
    struct DenseOut
    {
      T x0, x1;     // endpoints in state space
      T v0, v1;     // endpoint velocities f(x0), f(x1)
      Real h;       // step size
      F f;          // velocity field (captured by value)

      // Cubic Hermite interpolant X(t) for t in [0, h]
      T X(Real t) const
      {
        const Real s  = t / h;
        const Real s2 = s * s;
        const Real s3 = s2 * s;

        // Hermite basis
        const Real H00 = (2 * s3 - 3 * s2 + 1);
        const Real H10 = (    s3 - 2 * s2 + s);
        const Real H01 = (-2 * s3 + 3 * s2);
        const Real H11 = (    s3 -     s2);

        return H00 * x0 + H01 * x1 + h * (H10 * v0 + H11 * v1);
      }

      // Physical velocity along the interpolant
      T V(Real t) const
      {
        const T xt = X(t);
        return f(xt);
      }
    };

    // Build dense output on [0, dt]
    template <class T, class F>
    DenseOut<T, F> dense(T&& p, Real dt, F&& f) const
    {
      // RK4 stages (local, non-static to keep this re-entrant)
      T k1 = f(p);
      T k2 = f(p + 0.5 * dt * k1);
      T k3 = f(p + 0.5 * dt * k2);
      T k4 = f(p + dt * k3);

      // End state (RK4)
      const T x1 = p + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

      // Endpoint velocities
      const T v0 = k1;
      const T v1 = f(x1);

      return DenseOut<T, F>{ p, x1, v0, v1, dt, std::forward<F>(f) };
    }
  };
}

#endif
