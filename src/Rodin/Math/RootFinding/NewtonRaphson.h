#ifndef RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H
#define RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H

#include <cmath>
#include <optional>

#include <iostream>

namespace Rodin::Math::RootFinding
{
  template <class Scalar>
  class NewtonRaphson
  {
    public:
      NewtonRaphson(
          Scalar abs_t_tol = 1e-12,
          Scalar rel_t_tol = 1e-9,
          Scalar abs_g_tol = 1e-12,
          std::size_t max_iter = 25)
        : m_abs_t_tol(abs_t_tol),
          m_rel_t_tol(rel_t_tol),
          m_abs_g_tol(abs_g_tol),
          m_max_iter(max_iter)
      {}

      template <class Event>
      std::optional<Scalar>
      solve(const Event& event, Scalar t0, Scalar lo, Scalar hi) const
      {
        // Require an open interval (lo, hi) and try to keep it
        Scalar a = lo;
        Scalar b = hi;

        // Ensure initial guess is inside the open interval
        if (!(t0 > a && t0 < b))
          t0 = static_cast<Scalar>(0.5) * (a + b);

        // Evaluate at endpoints to establish a bracket (caller already checked sign change)
        Scalar ta = a; auto pa = event(ta); Scalar ga = pa.first;
        Scalar tb = b; auto pb = event(tb); Scalar gb = pb.first;

        Scalar t = t0;
        Scalar t_prev = t;
        Scalar g_prev = std::numeric_limits<Scalar>::quiet_NaN();
        bool have_prev = false;

        for (std::size_t it = 0; it < m_max_iter; ++it)
        {
          Scalar te = t; auto pg = event(te);
          Scalar g  = pg.first;
          Scalar gp = pg.second;

          // Convergence checks (strict)
          if (std::fabs(g) < m_abs_g_tol)
            return t;
          const Scalar at = m_abs_t_tol + m_rel_t_tol * std::fabs(t);
          if ((b - a) < at)
            return static_cast<Scalar>(0.5) * (a + b);

          // Keep the bracket with current sample
          if (g * ga < static_cast<Scalar>(0))
          {
            b = t; gb = g;
          }
          else
          {
            a = t; ga = g;
          }

          // Proposed Newton step
          Scalar dt;
          if (std::fabs(gp) > std::numeric_limits<Scalar>::epsilon())
          {
            dt = -g / gp;
          }
          else
          {
            // Derivative too small → fallback to secant/bisection
            if (have_prev && (g != g_prev)) // uses !=? forbidden. Avoid equality.
            {
              // Use secant only if denominator is not tiny (strict check via ratio)
              Scalar denom = (g - g_prev);
              if (std::fabs(denom) > std::numeric_limits<Scalar>::epsilon())
                dt = -g * (t - t_prev) / denom;
              else
                dt = static_cast<Scalar>(0); // will be rejected below, replaced by bisection
            }
            else
            {
              dt = static_cast<Scalar>(0); // will be rejected below, replaced by bisection
            }
          }

          Scalar t_new = t + dt;

          // Safeguard: keep inside (a, b). If outside or too aggressive, bisect.
          const Scalar mid = static_cast<Scalar>(0.5) * (a + b);
          const Scalar half = static_cast<Scalar>(0.5) * (b - a);
          bool outside = (t_new < a) || (t_new > b);
          bool too_large = std::fabs(t_new - t) > half;
          if (outside || too_large)
            t_new = mid;

          // Stopping on step size (strict)
          if (std::fabs(t_new - t) < at)
            return t_new;

          // Prepare next iteration
          t_prev  = t;
          g_prev  = g;
          have_prev = true;
          t = t_new;
          std::cout << "  it=" << it << " t=" << t << " g=" << g << " gp=" << gp << "\n";
        }

        // No convergence within max_iter
        return std::nullopt;
      }

    private:
      Scalar m_abs_t_tol;
      Scalar m_rel_t_tol;
      Scalar m_abs_g_tol;
      std::size_t m_max_iter;
  };
}

#endif
