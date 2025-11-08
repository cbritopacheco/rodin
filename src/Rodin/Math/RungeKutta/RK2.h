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
  };
}

#endif
