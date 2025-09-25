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
      q = p + (dt / 6) * (s_k1 + 2 * s_k2 + 2 * s_k3 + s_k4);
    }
  };
}

#endif
