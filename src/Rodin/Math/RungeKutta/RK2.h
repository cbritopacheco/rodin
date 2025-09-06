#ifndef RODIN_MATH_RUNGEKUTTA_RK2_H
#define RODIN_MATH_RUNGEKUTTA_RK2_H

#include "Rodin/Types.h"

namespace Rodin::Math::RungeKutta
{
  // Classical RK2 (midpoint) with cubic-Hermite dense output.
  struct RK2
  {
    // Step: q ← q + dt * f(p + 0.5*dt*f(p))
    template <class T, class G, class F>
    void step(T& q, Real dt, const G& p, const F& f) const
    {
      static thread_local T s_k1, s_k2;
      s_k1 = f(p);
      s_k2 = f(p + Real(0.5) * dt * s_k1);
      q += dt * s_k2;
    }

    // Dense output object for evaluating X(t), V(t) with t in (0, h).
    template <class T, class F>
    struct DenseOut
    {
      T    x0, x1;   // endpoints
      T    v0, v1;   // endpoint velocities f(x0), f(x1)
      Real h;        // step size
      F    f;        // vector field

      // Cubic Hermite interpolation
      T X(Real t) const
      {
        const Real s  = t / h;
        const Real s2 = s * s;
        const Real s3 = s2 * s;

        const Real H00 = (2*s3 - 3*s2 + 1);
        const Real H10 = (   s3 - 2*s2 + s);
        const Real H01 = (-2*s3 + 3*s2);
        const Real H11 = (   s3 -    s2);

        return H00 * x0 + H01 * x1 + h * (H10 * v0 + H11 * v1);
      }

      // Velocity along the interpolant
      T V(Real t) const
      {
        const T xt = X(t);
        return f(xt);
      }
    };

    // Build dense output on [0, dt]
    template <class T, class F>
    DenseOut<T, F> dense(const T& p, Real dt, F&& f) const
    {
      T k1 = f(p);
      T k2 = f(p + Real(0.5) * dt * k1);
      const T x1 = p + dt * k2;

      const T v0 = k1;
      const T v1 = f(x1);

      return DenseOut<T, F>{ p, x1, v0, v1, dt, std::forward<F>(f) };
    }
  };
}

#endif
