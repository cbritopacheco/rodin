#ifndef RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H
#define RODIN_MATH_ROOTFINDING_NEWTONRAPHSON_H

#include <cmath>
#include <optional>

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

      template <class F>
      std::optional<Scalar> solve(const F& f, Scalar t0, Scalar lo, Scalar hi) const
      {
        // basic interval sanity
        assert(std::isfinite(lo) && std::isfinite(hi));
        assert(lo < hi);

        Scalar a = lo, b = hi;

        // Initial values and bracket check
        auto fa = f(a);
        Scalar fav = fa.first;

        auto fb = f(b);
        Scalar fbv = fb.first;

        // endpoints must be finite
        assert(std::isfinite(fav) && std::isfinite(fbv));

        // either already bracketed or caller misused the API
        assert(fav * fbv <= 0); // Root not bracketed at [lo, hi]

        if (!(fav * fbv <= 0))
          return {};

        if (!(t0 > a && t0 < b))
          t0 = 0.5 * (a + b);

        Scalar t = t0, t_prev = t;
        Scalar g_prev = std::numeric_limits<Scalar>::quiet_NaN();
        bool have_prev = false;

        for (std::size_t it = 0; it < m_max_iter; ++it)
        {
          const auto ft = f(t);
          const Scalar fv = ft.first;
          const Scalar fg = ft.second;

          // function and derivative must be finite
          assert(std::isfinite(fv));
          assert(std::isfinite(fg) || std::isnan(fg)); // allow NaN to trigger safeguards

          if (std::fabs(fv) < m_abs_g_tol) return t;

          const Scalar at = m_abs_t_tol + m_rel_t_tol * std::fabs(t);
          if ((b - a) < at)
            return 0.5 * (a + b);

          // Maintain bracket using *current* endpoint values
          if (fv * fav < 0)
          {
            b = t;
            fbv = fv;
          }
          else
          {
            a = t;
            fav = fv;
          }

          if (std::fabs(fav) < m_abs_g_tol)
            return a;

          if (std::fabs(fbv) < m_abs_g_tol)
            return b;

          if (fav * fbv > 0)
            return {};

          // Newton step with scale-aware guard on gp
          Scalar dt;
          const Scalar gp_tol = std::numeric_limits<Scalar>::epsilon() * (static_cast<Scalar>(1) + std::fabs(fg));
          if (std::fabs(fg) > gp_tol)
          {
            dt = -fv / fg;
          }
          else if (have_prev)
          {
            const Scalar denom = (fv - g_prev);
            if (std::fabs(denom) > gp_tol)
              dt = -fv * (t - t_prev) / denom; // secant
            else
              dt = 0;
          }
          else
          {
            dt = 0;
          }

          Scalar t_new = t + dt;

          // Safeguard: keep inside (a,b). If outside or too big, bisect.
          const Scalar mid = 0.5 * (a + b);
          const Scalar half = 0.5 * (b - a);

          if (t_new <= a || t_new >= b || std::fabs(t_new - t) > half)
            t_new = mid;

          if (std::fabs(t_new - t) < at)
            return t_new;

          t_prev = t;
          g_prev = fv;
          have_prev = true;
          t = t_new;
        }
        return {};
      }

    private:
      Scalar m_abs_t_tol;
      Scalar m_rel_t_tol;
      Scalar m_abs_g_tol;
      std::size_t m_max_iter;
  };
}

#endif
