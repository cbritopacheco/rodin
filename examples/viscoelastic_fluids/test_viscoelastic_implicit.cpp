/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file test_viscoelastic_implicit.cpp
 * @brief Convergence tests of the fully-implicit log-conformation Oldroyd-B
 *        solver of Heart/LeftAtrium2D_viscous_logarithmic_implicit.
 *
 * The flow solver below is the parent one, term by term (psi form, Newton,
 * P1/P1/P1, BDF1, lagged split-OSS: S1 convection, pressure and div sigma, S2
 * grad-div, S3 u-psi compatibility and psi advection), with the PTT factor
 * switched off (epsilon = 0): Oldroyd-B. Only the domain, the data and the
 * error measurement are new. Two studies, both in the channel
 * (0, L) x (-H, H), inlet x = 0 (tag 1), outlet x = L (tag 2), walls y = +-H
 * (tag 3), meshes built in the code:
 *
 * 1. space: Pacheco & Castillo, IJNME 124 (2023) 1908-1927, Sec. 5.1.4,
 *    "ramped up channel flow". nu_0 = 1, beta = 0.5, lambda = 0.1, U = 1,
 *    H = 1, L = 3, t* = 2, T = 10; inflow 1.5 U (1 - y^2/H^2) f(t) with the
 *    smooth step (29); h and dt refined together from dt = 1/24 and twelve
 *    triangles; relative errors of u (L2, H1), p (L2) and sigma (L2) at t = T
 *    against the stationary solution (their Fig. 5). The paper uses the
 *    half-channel with a symmetry line; here the full channel is meshed, which
 *    is the same problem. The stationary solution is re-derived: as printed
 *    it has p = 3 U nu0 (L - 2 x1)/H^2, sigma11 = 72 (...) and a lambda in
 *    sigma12, none of which satisfies (1)-(6); the Oldroyd-B solution is
 *      u1 = 1.5 U (1 - y^2/H^2),  gdot = -3 U y/H^2,
 *      sigma12 = eta_p gdot,  sigma11 = 2 lambda eta_p gdot^2,  sigma22 = 0,
 *      p = 3 U eta_0 (L - x)/H^2.
 *
 * 2. time: the paper's temporal test (Sec. 5.1.1) and its Taylor-Green test
 *    (5.1.2) use the closed-form family (27)-(28), for which the conformation
 *    I + (lambda/eta_p) sigma is IDENTICALLY ZERO: psi = log(c) does not exist
 *    and no log-conformation solver can run them. They are replaced by an
 *    exact, fully time-dependent Oldroyd-B solution with positive-definite
 *    conformation: pulsatile channel flow at q(t) = 1 + A sin(omega t).
 *    For unidirectional flow the upper-convected terms leave
 *      sigma12 + lambda d_t sigma12 = eta_p gdot     (Maxwell, linear),
 *      sigma11 + lambda (d_t sigma11 - 2 gdot sigma12) = 0,
 *    so u and sigma12 are Womersley modes with eta*(w) = eta_s + eta_p/(1 +
 *    i w lambda), and sigma11 is the finite Fourier sum (modes 0, +-1, +-2) of
 *    the product 2 lambda gdot sigma12 filtered by 1/(1 + i m w lambda).
 *    The run starts from the exact state at t = 0; dt is refined on a fixed
 *    mesh and each solution is compared both with the exact one (total error,
 *    which stalls at the spatial error) and with a reference computed with
 *    dt_min/4 on the same mesh (temporal error alone).
 *
 * Boundary data are exact: u and psi = log(I + (lambda0/eta_p) sigma) on the
 * inlet (the psi equation is hyperbolic), u = 0 on the walls, and the exact
 * traction (2 eta_s eps(u) + sigma - p I) n on the outlet. The backflow term
 * of the parent is inconsistent with an exact solution and is off by default.
 *
 * Run (sequential is enough; any number of ranks works):
 *   ./examples/viscoelastic_fluids/test_viscoelastic_implicit                 # both studies
 *   ./examples/viscoelastic_fluids/test_viscoelastic_implicit -vt_study space -vt_levels 5
 *   ./examples/viscoelastic_fluids/test_viscoelastic_implicit -vt_study time  -vt_n 16
 *   python3 ../examples/viscoelastic_fluids/plot_test_viscoelastic.py        # -> PNG
 *
 * Options (prefix -vt_): study space|time|all, levels, n (time-study mesh),
 * n0 (space-study coarsest mesh, cells across 2H), dt0 (space study),
 * steps0 (time study, steps per period at the coarsest level),
 * lambda, beta, amplitude, wo (time study), conformation_its,
 * conformation_tol, newton_max_step, backflow, vms, and the stabilisation
 * scales vms_scale, graddiv_scale, pressure_scale, stress_div_scale,
 * stress_scale.
 */
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "OrthogonalProjection.h"

namespace Rodin::Examples::ViscoelasticFluids
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  // The block below is copied verbatim from
  // Heart/LeftAtrium2D_viscous_logarithmic_implicit.cpp (its anonymous
  // namespace), so that the test assembles the same operators.
  namespace
  {
    constexpr int RootRank = 0;

    void setPrefixedDefault(
      const std::string& prefix, const char* suffix, const char* value)
    {
      const std::string name = "-" + prefix + suffix;
      PetscBool set = PETSC_FALSE;
      PetscErrorCode ierr =
        PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name.c_str(), &set);
      assert(ierr == PETSC_SUCCESS);
      if (!set)
      {
        ierr = PetscOptionsSetValue(PETSC_NULLPTR, name.c_str(), value);
        assert(ierr == PETSC_SUCCESS);
      }
      (void)ierr;
    }

    // ---- sigma in Voigt order (xx, xy, yy), 2D ---------------------------------
    // tensor() turns the stored vector into the full symmetric tensor, so every
    // pairing below is sigma:chi. The rest are the operators of the
    // constitutive equation; they take trial/test functions and finite element
    // fields alike, so a term and its projection share one definition.

    /// @brief E0 = e_x e_x, E1 = e_x e_y + e_y e_x, E2 = e_y e_y.
    MatrixFunction<Math::Matrix<Real>> voigtBasis(size_t k)
    {
      Math::Matrix<Real> E = Math::Matrix<Real>::Zero(2, 2);
      if (k == 0)
        E(0, 0) = 1;
      else if (k == 1)
        E(0, 1) = E(1, 0) = 1;
      else
        E(1, 1) = 1;
      return MatrixFunction(E);
    }

    /// @brief e_j. VectorFunction keeps a reference, hence the static storage.
    VectorFunction<Math::Vector<Real>> unit(size_t j)
    {
      static const Math::Vector<Real> e[2] = {
        Math::Vector<Real>::Unit(2, 0), Math::Vector<Real>::Unit(2, 1) };
      return VectorFunction(e[j]);
    }

    /// @brief sigma = E0 s_xx + E1 s_xy + E2 s_yy.
    template <class S>
    auto tensor(const S& s)
    {
      return voigtBasis(0) * Component(s, 0) + voigtBasis(1) * Component(s, 1) +
        voigtBasis(2) * Component(s, 2);
    }

    /// @brief The Voigt vector of a symmetric tensor function, for projections.
    template <class M>
    auto voigt(const M& m)
    {
      return VectorFunction{ Component(m, 0, 0), Component(m, 0, 1), Component(m, 1, 1) };
    }

    /// @brief (u . grad) sigma.
    template <class S, class U>
    auto advection(const S& s, const U& u)
    {
      return tensor(Mult(Jacobian(s), u));
    }

    /// @brief (div sigma)_i = d sigma_ij / d x_j.
    template <class S>
    auto divergence(const S& s)
    {
      const auto dx = Mult(Jacobian(s), unit(0));
      const auto dy = Mult(Jacobian(s), unit(1));
      return unit(0) * (Component(dx, 0) + Component(dy, 1)) +
        unit(1) * (Component(dx, 1) + Component(dy, 2));
    }

    // ---- log-conformation ---------------------------------------------------

    /// @brief Eigen-decomposition of a symmetric 2x2 psi (Voigt), with exp,
    ///        Dexp and Dexp^{-1} by Daleckii-Krein: in the eigenbasis R,
    ///        R^T Dexp[H] R = F o (R^T H R), F_ij = (e^li - e^lj)/(li - lj),
    ///        F_ii = e^li, and Dexp^{-1} divides by F instead. F > 0 always, so
    ///        the inverse needs no special case for equal eigenvalues.
    struct LogEig
    {
        Eigen::Matrix2d R;
        Eigen::Vector2d l;
        Eigen::Matrix2d F;

        Eigen::Matrix2d exp() const
        {
          return R * l.array().exp().matrix().asDiagonal() * R.transpose();
        }
        Eigen::Matrix2d expMinus() const
        {
          return R * (-l).array().exp().matrix().asDiagonal() * R.transpose();
        }
        Eigen::Matrix2d dexp(const Eigen::Matrix2d& H) const
        {
          return R * F.cwiseProduct(R.transpose() * H * R) * R.transpose();
        }
        Eigen::Matrix2d dexpInverse(const Eigen::Matrix2d& H) const
        {
          return R * (R.transpose() * H * R).cwiseQuotient(F) * R.transpose();
        }
    };

    LogEig logEig(const Eigen::Matrix2d& psi)
    {
      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(psi);
      LogEig out;
      out.R = eig.eigenvectors();
      out.l = eig.eigenvalues();
      const Real e0 = std::exp(out.l(0)), e1 = std::exp(out.l(1)), d = out.l(1) - out.l(0);
      const Real f = std::abs(d) > 1e-12 ? e0 * std::expm1(d) / d : e0 * (1.0 + 0.5 * d);
      out.F << e0, f, f, e1;
      return out;
    }

    template <class Voigt>
    Eigen::Matrix2d unpack(const Voigt& v)
    {
      Eigen::Matrix2d m;
      m << v(0), v(1), v(1), v(2);
      return m;
    }

    Math::SpatialVector<Real> pack(const Eigen::Matrix2d& m)
    {
      return Math::SpatialVector<Real>{{ m(0, 0), m(0, 1), m(1, 1) }};
    }

    const std::array<Eigen::Matrix2d, 3>& voigtBasisEigen()
    {
      static const std::array<Eigen::Matrix2d, 3> B = {
        (Eigen::Matrix2d() << 1, 0, 0, 0).finished(),
        (Eigen::Matrix2d() << 0, 1, 1, 0).finished(),
        (Eigen::Matrix2d() << 0, 0, 0, 1).finished() };
      return B;
    }

    /// @brief The constants of the psi equation.
    struct LogModel
    {
        Real lambda = 0.0;
        Real lambda0 = 0.0;
        Real epsilon = 0.0;   ///< sPTT
        Real c = 0.0;         ///< 2 (1 - lambda0/lambda), on eps(u)
    };

    /// @brief Everything the forms need at one quadrature point, from
    ///        (psi^k, L^k = grad u^k): exp, Dexp[B_j] (momentum, S1), and the
    ///        Newton data of the psi equation, N = -Dexp^{-1}[L T + T L^T -
    ///        c eps] + (f/lambda)(I - exp(-psi)), its FD Jacobian JN_j =
    ///        dN/dpsi_j, and Q_ab = -Dexp^{-1}[E_ab T + T E_ab^T - c sym(E_ab)],
    ///        the exact (linear) dependence on L_ab.
    struct LogPoint
    {
        Eigen::Matrix2d exp;
        std::array<Eigen::Matrix2d, 3> D;
        Real f = 1.0;
        Eigen::Matrix2d N;
        std::array<Eigen::Matrix2d, 3> JN;
        std::array<std::array<Eigen::Matrix2d, 2>, 2> Q;
    };

    Eigen::Matrix2d nonlinearMap(
      const Eigen::Matrix2d& psi, const Eigen::Matrix2d& L, const LogModel& m, Real* fOut = nullptr)
    {
      const LogEig e = logEig(psi);
      const Eigen::Matrix2d T = e.exp();
      const Eigen::Matrix2d H = L * T + T * L.transpose() - 0.5 * m.c * (L + L.transpose());
      const Real f = 1.0 + m.epsilon * (m.lambda / m.lambda0) * (T.trace() - 2.0);
      if (fOut)
        *fOut = f;
      return -e.dexpInverse(H) +
        (f / m.lambda) * (Eigen::Matrix2d::Identity() - e.expMinus());
    }

    LogPoint logPoint(const Eigen::Matrix2d& psi, const Eigen::Matrix2d& L, const LogModel& m)
    {
      LogPoint out;
      const LogEig e = logEig(psi);
      out.exp = e.exp();
      const auto& B = voigtBasisEigen();
      for (size_t j = 0; j < 3; ++j)
        out.D[j] = e.dexp(B[j]);
      out.N = nonlinearMap(psi, L, m, &out.f);
      // Central differences in the three Voigt directions; the map is smooth
      // and the eigenproblem is 2x2, so this costs six decompositions.
      const Real h = 1.0e-6 * std::max<Real>(1.0, psi.cwiseAbs().maxCoeff());
      for (size_t j = 0; j < 3; ++j)
        out.JN[j] = (nonlinearMap(psi + h * B[j], L, m) - nonlinearMap(psi - h * B[j], L, m)) / (2.0 * h);
      for (size_t a = 0; a < 2; ++a)
        for (size_t b = 0; b < 2; ++b)
        {
          Eigen::Matrix2d E = Eigen::Matrix2d::Zero();
          E(a, b) = 1.0;
          out.Q[a][b] = -e.dexpInverse(E * out.exp + out.exp * E.transpose() -
            0.5 * m.c * (E + E.transpose()));
        }
      return out;
    }

    /// @brief A 2x2 tensor field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseTensor final
      : public MatrixFunctionBase<Real, PointwiseTensor<F>>
    {
      public:
        using Parent = MatrixFunctionBase<Real, PointwiseTensor<F>>;

        explicit PointwiseTensor(F f) : m_f(std::move(f)) {}
        PointwiseTensor(const PointwiseTensor& other) : Parent(other), m_f(other.m_f) {}
        PointwiseTensor(PointwiseTensor&& other) : Parent(std::move(other)), m_f(std::move(other.m_f)) {}

        Math::Matrix<Real> getValue(const Point& p) const { return Math::Matrix<Real>(m_f(p)); }
        size_t getRows() const { return 2; }
        size_t getColumns() const { return 2; }
        /// Not a polynomial; quadratic is enough against P1 x P1.
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept { return 2; }
        PointwiseTensor* copy() const noexcept override { return new PointwiseTensor(*this); }

      private:
        F m_f;
    };

    /// @brief logPoint(psi(x), grad u(x)), cached per quadrature point: the
    ///        form language evaluates a coefficient once per basis function,
    ///        and without the cache the eigenproblems dominate assembly. Keyed
    ///        by the psi field, its revision (bumped on every change of the
    ///        pair), cell and reference coordinates.
    template <class GF, class UF>
    const LogPoint& logPointAt(
      const GF& psi, const UF& u, const LogModel& m, std::uint64_t revision, const Point& p)
    {
      struct Entry
      {
        const void* field = nullptr;
        std::uint64_t revision = 0;
        Index cell = 0;
        Real r0 = 0, r1 = 0;
        LogPoint value;
      };
      thread_local std::array<Entry, 16> cache;
      thread_local size_t next = 0;

      const Index cell = p.getPolytope().getIndex();
      const auto& rc = p.getReferenceCoordinates();
      for (const auto& e : cache)
        if (e.field == &psi && e.revision == revision && e.cell == cell &&
            e.r0 == rc(0) && e.r1 == rc(1))
          return e.value;

      auto& e = cache[next];
      next = (next + 1) % cache.size();
      const auto J = Jacobian(u).getValue(p);   // Math::SpatialMatrix, not Eigen
      Eigen::Matrix2d L;
      L << J(0, 0), J(0, 1), J(1, 0), J(1, 1);
      e = Entry{ &psi, revision, cell, rc(0), rc(1),
        logPoint(unpack(psi.getValue(p)), L, m) };
      return e.value;
    }

    /// @brief A pointwise matrix drawn from the cached LogPoint by `pick`.
    template <class GF, class UF, class Pick>
    auto fromPoint(const GF& psi, const UF& u, const LogModel& m,
      const std::uint64_t& revision, Pick pick)
    {
      return PointwiseTensor([&psi, &u, m, &revision, pick](const Point& p) -> Eigen::Matrix2d {
        return pick(logPointAt(psi, u, m, revision, p)); });
    }

    /// @brief Dexp_{psi(x)}[w] = sum_j w_j Dexp[B_j] for a Voigt w: psi itself
    ///        or a derivative of it (chain rule); trial, test or field.
    template <class GF, class UF, class W>
    auto dexp(const GF& psi, const UF& u, const LogModel& m, const std::uint64_t& revision, const W& w)
    {
      const auto D = [&](size_t j) {
        return fromPoint(psi, u, m, revision, [j](const LogPoint& q) { return q.D[j]; });
      };
      return D(0) * Component(w, 0) + D(1) * Component(w, 1) + D(2) * Component(w, 2);
    }

    /// @brief JN_{psi(x)}[w] = sum_j w_j dN/dpsi_j.
    template <class GF, class UF, class W>
    auto jacobianN(const GF& psi, const UF& u, const LogModel& m, const std::uint64_t& revision, const W& w)
    {
      const auto J = [&](size_t j) {
        return fromPoint(psi, u, m, revision, [j](const LogPoint& q) { return q.JN[j]; });
      };
      return J(0) * Component(w, 0) + J(1) * Component(w, 1) + J(2) * Component(w, 2);
    }

    /// @brief -Dexp^{-1}[L(v) T + T L(v)^T - c eps(v)] = sum_ab L_ab(v) Q_ab, for
    ///        a vector v (trial or field): the exact velocity dependence of the
    ///        psi equation at psi^k.
    template <class GF, class UF, class LogM, class V>
    auto velocityMap(const GF& psi, const UF& u, const LogM& m, const std::uint64_t& revision, const V& v)
    {
      const auto Q = [&](size_t a, size_t b) {
        return fromPoint(psi, u, m, revision, [a, b](const LogPoint& q) { return q.Q[a][b]; });
      };
      // L_ab = d v_a / d x_b = (J(v) e_b)_a
      const auto col0 = Mult(Jacobian(v), unit(0));
      const auto col1 = Mult(Jacobian(v), unit(1));
      return Q(0, 0) * Component(col0, 0) + Q(1, 0) * Component(col0, 1) +
        Q(0, 1) * Component(col1, 0) + Q(1, 1) * Component(col1, 1);
    }

    /// @brief The projections invert mass matrices: CG + Jacobi, never the
    ///        global direct solver used by the coupled flow system.
    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  // ==========================================================================
  // Exact solutions of the Oldroyd-B channel problem
  // ==========================================================================
  /// @brief Velocity, velocity gradient (grad u)_ij = du_i/dx_j, pressure and
  ///        polymeric stress sigma, as functions of (x, y, t).
  class ExactSolution
  {
    public:
      virtual ~ExactSolution() = default;
      virtual Eigen::Vector2d velocity(Real x, Real y, Real t) const = 0;
      virtual Eigen::Matrix2d velocityGradient(Real x, Real y, Real t) const = 0;
      virtual Real pressure(Real x, Real y, Real t) const = 0;
      virtual Eigen::Matrix2d stress(Real x, Real y, Real t) const = 0;
  };

  /**
   * @brief Pacheco & Castillo (2023), Sec. 5.1.4: Poiseuille flow switched on
   *        by the smooth step f(t) of their Eq. (29).
   *
   * @details Exact (stationary) for t >= t* once the start-up has relaxed; for
   *          t < t* the same expressions scaled by f(t) are used as inlet and
   *          outlet data, as the paper does for the velocity. The error is
   *          only measured at T = 10 = 100 lambda, long after both.
   */
  class RampedChannel final : public ExactSolution
  {
    public:
      RampedChannel(Real U, Real H, Real L, Real etaS, Real etaP, Real lambda, Real tStar)
        : m_U(U), m_H(H), m_L(L), m_etaS(etaS), m_etaP(etaP), m_lambda(lambda), m_tStar(tStar)
      {}

      /// @brief f(t) = 6 s^5 - 15 s^4 + 10 s^3, s = t/t*, and 1 for t >= t*.
      Real ramp(Real t) const
      {
        if (t <= 0.0)
          return 0.0;
        if (t >= m_tStar)
          return 1.0;
        const Real s = t / m_tStar;
        return s * s * s * (10.0 - 15.0 * s + 6.0 * s * s);
      }

      Eigen::Vector2d velocity(Real, Real y, Real t) const override
      {
        return Eigen::Vector2d(1.5 * m_U * (1.0 - y * y / (m_H * m_H)) * ramp(t), 0.0);
      }

      Eigen::Matrix2d velocityGradient(Real, Real y, Real t) const override
      {
        Eigen::Matrix2d g = Eigen::Matrix2d::Zero();
        g(0, 1) = shearRate(y, t);
        return g;
      }

      Real pressure(Real x, Real, Real t) const override
      {
        return 3.0 * m_U * (m_etaS + m_etaP) * ramp(t) * (m_L - x) / (m_H * m_H);
      }

      Eigen::Matrix2d stress(Real, Real y, Real t) const override
      {
        const Real gd = shearRate(y, t);
        Eigen::Matrix2d s;
        s << 2.0 * m_lambda * m_etaP * gd * gd, m_etaP * gd, m_etaP * gd, 0.0;
        return s;
      }

    private:
      Real shearRate(Real y, Real t) const
      {
        return -3.0 * m_U * y * ramp(t) / (m_H * m_H);
      }

      Real m_U, m_H, m_L, m_etaS, m_etaP, m_lambda, m_tStar;
  };

  /**
   * @brief Fully developed pulsatile Oldroyd-B channel flow at the flow rate
   *        2 h U (1 + A sin(omega t)): an exact solution of the full problem.
   *
   * @details With X(t) = sum_m X_m e^{i m omega t}, X_{-m} = conj(X_m):
   *          u_0 = (3/2) U (1 - y^2/h^2), u_1 = a_1 (1 - cosh(ky)/cosh(kh)),
   *          k^2 = i omega rho/eta*, eta* = eta_s + eta_p/(1 + i omega lambda),
   *          a_1 = U c_1/(1 - tanh(kh)/(kh)), c_1 = -iA/2;
   *          sigma12_m = eta_p/(1 + i m omega lambda) du_m/dy;
   *          sigma11_m = [2 lambda sum_{a+b=m} sigma12_a du_b/dy]/(1 + i m omega lambda);
   *          sigma22 = 0; dp/dx = G_0 = -3 eta_0 U/h^2 and G_1 = -i omega rho a_1,
   *          p = G(t) (x - L). Every hyperbolic function is written with
   *          decaying exponentials only.
   */
  class MaxwellWomersley final : public ExactSolution
  {
    public:
      using Complex = std::complex<Real>;

      MaxwellWomersley(Real U, Real A, Real h, Real L, Real omega, Real rho,
        Real etaS, Real etaP, Real lambda)
        : m_U(U), m_h(h), m_L(L), m_omega(omega), m_etaP(etaP), m_lambda(lambda)
      {
        const Complex I(0.0, 1.0);
        const Complex etaStar = etaS + etaP / (1.0 + I * omega * lambda);
        m_k = std::sqrt(I * omega * rho / etaStar);
        m_e2 = std::exp(-2.0 * m_k * h);
        const Complex tanhKh = (1.0 - m_e2) / (1.0 + m_e2);
        m_a1 = U * (-0.5 * I * A) / (1.0 - tanhKh / (m_k * h));
        m_mu1 = etaP / (1.0 + I * omega * lambda);
        m_G0 = -3.0 * (etaS + etaP) * U / (h * h);
        m_G1 = -I * omega * rho * m_a1;
      }

      Eigen::Vector2d velocity(Real, Real y, Real t) const override
      {
        const Real u = 1.5 * m_U * (1.0 - y * y / (m_h * m_h)) +
          2.0 * (m_a1 * (1.0 - coshRatio(y)) * phase(1, t)).real();
        return Eigen::Vector2d(u, 0.0);
      }

      Eigen::Matrix2d velocityGradient(Real, Real y, Real t) const override
      {
        Eigen::Matrix2d g = Eigen::Matrix2d::Zero();
        g(0, 1) = shearRate0(y) + 2.0 * (shearRate1(y) * phase(1, t)).real();
        return g;
      }

      Real pressure(Real x, Real, Real t) const override
      {
        return (m_G0 + 2.0 * (m_G1 * phase(1, t)).real()) * (x - m_L);
      }

      Eigen::Matrix2d stress(Real, Real y, Real t) const override
      {
        const Complex I(0.0, 1.0);
        const Real g0 = shearRate0(y);
        const Complex g1 = shearRate1(y);
        const Real s0 = m_etaP * g0;
        const Complex s1 = m_mu1 * g1;
        const Real sigma12 = s0 + 2.0 * (s1 * phase(1, t)).real();
        // Modes 0, 1, 2 of 2 lambda sigma12 gdot, filtered by 1/(1 + i m w lambda).
        const Real p0 = 2.0 * m_lambda * (s0 * g0 + 2.0 * (s1 * std::conj(g1)).real());
        const Complex p1 = 2.0 * m_lambda * (s0 * g1 + s1 * g0);
        const Complex p2 = 2.0 * m_lambda * s1 * g1;
        const Real sigma11 = p0 +
          2.0 * (p1 / (1.0 + I * m_omega * m_lambda) * phase(1, t)).real() +
          2.0 * (p2 / (1.0 + 2.0 * I * m_omega * m_lambda) * phase(2, t)).real();
        Eigen::Matrix2d s;
        s << sigma11, sigma12, sigma12, 0.0;
        return s;
      }

      Real getPeriod() const
      {
        return 2.0 * M_PI / m_omega;
      }

    private:
      Complex phase(int m, Real t) const
      {
        const Real theta = m * m_omega * t;
        return Complex(std::cos(theta), std::sin(theta));
      }

      /// @brief cosh(ky)/cosh(kh), |y| <= h.
      Complex coshRatio(Real y) const
      {
        return (std::exp(m_k * (y - m_h)) + std::exp(-m_k * (y + m_h))) / (1.0 + m_e2);
      }

      /// @brief sinh(ky)/cosh(kh), |y| <= h.
      Complex sinhRatio(Real y) const
      {
        return (std::exp(m_k * (y - m_h)) - std::exp(-m_k * (y + m_h))) / (1.0 + m_e2);
      }

      Real shearRate0(Real y) const
      {
        return -3.0 * m_U * y / (m_h * m_h);
      }

      Complex shearRate1(Real y) const
      {
        return -m_a1 * m_k * sinhRatio(y);
      }

      Real m_U, m_h, m_L, m_omega, m_etaP, m_lambda;
      Complex m_k, m_e2, m_a1, m_mu1, m_G1;
      Real m_G0 = 0.0;
  };

  /// @brief psi = log(I + (lambda0/eta_p) sigma), the log-conformation of an
  ///        exact stress (Castillo's scaling sigma = (eta_p/lambda0)(exp(psi) - I)).
  Eigen::Matrix2d logConformation(const Eigen::Matrix2d& sigma, Real lambda0, Real etaP)
  {
    const Eigen::Matrix2d tau = Eigen::Matrix2d::Identity() + (lambda0 / etaP) * sigma;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(tau);
    if (!(eig.eigenvalues().minCoeff() > 0.0))
      throw std::runtime_error("The exact conformation is not positive definite.");
    return eig.eigenvectors() * eig.eigenvalues().array().log().matrix().asDiagonal() *
      eig.eigenvectors().transpose();
  }

  // ==========================================================================
  // The solver under test
  // ==========================================================================
  class ViscoelasticChannelTest
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;
      using AttributeSet = FlatSet<Attribute>;

      using VectorFESType = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFESType = H1<1, Real, MeshType>;
      using VectorGridFunctionType = PETSc::Variational::GridFunction<VectorFESType>;
      using ScalarGridFunctionType = PETSc::Variational::GridFunction<ScalarFESType>;
      using VectorTrialFunctionType =
        PETSc::Variational::TrialFunction<VectorGridFunctionType, VectorFESType>;
      using ScalarTrialFunctionType =
        PETSc::Variational::TrialFunction<ScalarGridFunctionType, ScalarFESType>;
      using VectorTestFunctionType = PETSc::Variational::TestFunction<VectorFESType>;
      using ScalarTestFunctionType = PETSc::Variational::TestFunction<ScalarFESType>;
      using LinearSystemType = PETSc::Math::LinearSystem;
      using FlowProblemType = Problem<LinearSystemType, VectorTrialFunctionType,
        ScalarTrialFunctionType, VectorTrialFunctionType, VectorTestFunctionType,
        ScalarTestFunctionType, VectorTestFunctionType>;
      using CellFESType = P0<Real, MeshType>;
      using CellGridFunctionType = PETSc::Variational::GridFunction<CellFESType>;

      static constexpr Attribute Inlet = 1;
      static constexpr Attribute Outlet = 2;
      static constexpr Attribute Wall = 3;

      /// @brief Oldroyd-B: the parent's sPTT with epsilon = 0.
      struct OldroydB
      {
          Real etaS = 0.5;
          Real etaP = 0.5;
          Real lambda = 0.1;
          Real lambda0Factor = 1.0;
          Real lambda0Min = 0.0;
          Real pttEpsilon = 0.0;
      };

      struct Config
      {
          Real length = 3.0;       ///< L
          Real halfHeight = 1.0;   ///< H
          size_t nx = 3;
          size_t ny = 2;
          Real rho = 1.0;
          OldroydB oldroydB;
          Real pressurePenalty = 1.0e-12;
          /// @brief The parent's directional do-nothing term. Not consistent
          ///        with an exact solution, hence 0.
          Real outletBackflowStabilization = 0.0;
          Real vmsScale = 1.0;
          Real gradDivScale = 1.0;
          Real pressureScale = 1.0;
          Real stressDivScale = 1.0;
          Real stressScale = 1.0;
          bool useVMS = true;
          int conformationIterations = 10;
          Real conformationTolerance = 1.0e-10;
          Real newtonMaxStep = 2.0;
      };

      /// @brief Squared L2 norms of the errors and of the reference field.
      struct Errors
      {
          Real u = 0, gradU = 0, p = 0, sigma = 0;
          Real uRef = 0, gradURef = 0, pRef = 0, sigmaRef = 0;

          Real relU() const { return std::sqrt(u / uRef); }
          Real relH1() const { return std::sqrt((u + gradU) / (uRef + gradURef)); }
          Real relP() const { return std::sqrt(p / pRef); }
          Real relSigma() const { return std::sqrt(sigma / sigmaRef); }
      };

      /// @brief A stored discrete state (u, p, psi) on this solver's spaces.
      struct Snapshot
      {
          VectorGridFunctionType u;
          ScalarGridFunctionType p;
          VectorGridFunctionType psi;
      };

      ViscoelasticChannelTest(const Context::MPI& context, const Config& cfg,
        std::shared_ptr<const ExactSolution> exact)
        : m_cfg(cfg),
          m_exact(std::move(exact)),
          m_mesh(makeMesh(context, m_cfg)),
          m_wallSet{ Wall },
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_tauh(std::integral_constant<size_t, 1>{}, m_mesh,
            m_mesh.getSpaceDimension() * (m_mesh.getSpaceDimension() + 1) / 2),
          m_ch(m_mesh),
          m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh),
          m_uOld(m_vh), m_uIt(m_vh),
          m_psi(m_tauh), m_chi(m_tauh), m_psiOld(m_tauh), m_psiIt(m_tauh), m_sigma(m_tauh),
          m_tauK(m_ch), m_alpha1(m_ch), m_alpha2(m_ch), m_alpha3(m_ch), m_alphaPsi(m_ch),
          m_piConv(m_vh, "vt_proj_"), m_sub(m_vh, "vt_proj_"), m_subOld(m_vh),
          m_piDiv(m_sh, "vt_proj_"), m_piGradP(m_vh, "vt_proj_"),
          m_piDivSigma(m_vh, "vt_proj_"), m_piEps(m_tauh, "vt_proj_"),
          m_piAdvPsi(m_tauh, "vt_proj_"),
          m_flow(m_u, m_p, m_psi, m_v, m_q, m_chi),
          m_flowKSP(m_flow)
      {
        m_h = std::hypot(m_cfg.length / m_cfg.nx, 2.0 * m_cfg.halfHeight / m_cfg.ny);
      }

      Real getMeshSize() const { return m_h; }

      /// @brief From the exact state at t = 0 to t = T with steps dt.
      bool run(Real dt, Real T)
      {
        m_dt = dt;
        m_t = 0.0;
        setInitialState();
        setupFlow();
        m_flowFieldSplitsSet = false;
        const int steps = static_cast<int>(std::lround(T / dt));
        for (int n = 0; n < steps; ++n)
        {
          m_t += m_dt;
          if (!solveFlow())
          {
            if (isRoot())
              Alert::Warning() << "[test] the flow solve failed at t = " << m_t
                               << " (max|u| = " << m_speed << ")" << Alert::Raise;
            return false;
          }
        }
        return true;
      }

      Snapshot snapshot()
      {
        Snapshot s{ VectorGridFunctionType(m_vh), ScalarGridFunctionType(m_sh),
          VectorGridFunctionType(m_tauh) };
        s.u.setData(m_u.getSolution().getData());
        s.p.setData(m_p.getSolution().getData());
        s.psi.setData(m_psi.getSolution().getData());
        return s;
      }

      /// @brief Errors of the current state against the exact solution at t.
      Errors errorToExact(Real t)
      {
        const Real s = m_cfg.oldroydB.etaP / lambda0();
        const auto& uSol = m_u.getSolution();
        const auto& pSol = m_p.getSolution();
        const auto& psiSol = m_psi.getSolution();
        const auto J = Jacobian(uSol);
        return integrate([&](const Point& pt, Errors& e, Real w) {
          const auto& x = pt.getPhysicalCoordinates();
          const Eigen::Vector2d uh = toEigen(uSol.getValue(pt));
          const auto Jh = J.getValue(pt);
          Eigen::Matrix2d gh;
          gh << Jh(0, 0), Jh(0, 1), Jh(1, 0), Jh(1, 1);
          const Real ph = pSol.getValue(pt);
          const Eigen::Matrix2d sh =
            s * (logEig(unpack(psiSol.getValue(pt))).exp() - Eigen::Matrix2d::Identity());

          const Eigen::Vector2d ue = m_exact->velocity(x(0), x(1), t);
          const Eigen::Matrix2d ge = m_exact->velocityGradient(x(0), x(1), t);
          const Real pe = m_exact->pressure(x(0), x(1), t);
          const Eigen::Matrix2d se = m_exact->stress(x(0), x(1), t);
          accumulate(e, w, uh - ue, gh - ge, ph - pe, sh - se, ue, ge, pe, se);
        });
      }

      /// @brief Errors of the current state against a stored one.
      Errors errorTo(const Snapshot& ref)
      {
        const Real s = m_cfg.oldroydB.etaP / lambda0();
        const auto& uSol = m_u.getSolution();
        const auto& pSol = m_p.getSolution();
        const auto& psiSol = m_psi.getSolution();
        const auto J = Jacobian(uSol);
        const auto JRef = Jacobian(ref.u);
        return integrate([&](const Point& pt, Errors& e, Real w) {
          const auto toMatrix = [](const auto& M) {
            Eigen::Matrix2d g;
            g << M(0, 0), M(0, 1), M(1, 0), M(1, 1);
            return g;
          };
          const Eigen::Vector2d uh = toEigen(uSol.getValue(pt));
          const Eigen::Vector2d ur = toEigen(ref.u.getValue(pt));
          const Eigen::Matrix2d gh = toMatrix(J.getValue(pt));
          const Eigen::Matrix2d gr = toMatrix(JRef.getValue(pt));
          const Real ph = pSol.getValue(pt);
          const Real pr = ref.p.getValue(pt);
          const Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
          const Eigen::Matrix2d sh = s * (logEig(unpack(psiSol.getValue(pt))).exp() - I);
          const Eigen::Matrix2d sr = s * (logEig(unpack(ref.psi.getValue(pt))).exp() - I);
          accumulate(e, w, uh - ur, gh - gr, ph - pr, sh - sr, ur, gr, pr, sr);
        });
      }

    private:
      bool isRoot() const
      {
        return m_mesh.getContext().getCommunicator().rank() == RootRank;
      }

      static Eigen::Vector2d toEigen(const Math::SpatialVector<Real>& v)
      {
        return Eigen::Vector2d(v(0), v(1));
      }

      static void accumulate(Errors& e, Real w,
        const Eigen::Vector2d& du, const Eigen::Matrix2d& dg, Real dp, const Eigen::Matrix2d& ds,
        const Eigen::Vector2d& u, const Eigen::Matrix2d& g, Real p, const Eigen::Matrix2d& s)
      {
        e.u += w * du.squaredNorm();
        e.gradU += w * dg.squaredNorm();
        e.p += w * dp * dp;
        e.sigma += w * ds.squaredNorm();
        e.uRef += w * u.squaredNorm();
        e.gradURef += w * g.squaredNorm();
        e.pRef += w * p * p;
        e.sigmaRef += w * s.squaredNorm();
      }

      /// @brief Sum of f over the owned cells with the degree-5 Dunavant rule,
      ///        reduced over the communicator.
      template <class F>
      Errors integrate(F&& f) const
      {
        static const std::array<std::array<Real, 2>, 7> xq = { {
          { 1.0 / 3.0, 1.0 / 3.0 },
          { 0.470142064105115, 0.470142064105115 }, { 0.059715871789770, 0.470142064105115 },
          { 0.470142064105115, 0.059715871789770 },
          { 0.101286507323456, 0.101286507323456 }, { 0.797426985353087, 0.101286507323456 },
          { 0.101286507323456, 0.797426985353087 } } };
        static const std::array<Real, 7> wq = { 0.225,
          0.132394152788506, 0.132394152788506, 0.132394152788506,
          0.125939180544827, 0.125939180544827, 0.125939180544827 };

        const size_t D = m_mesh.getDimension();
        const auto& shard = m_mesh.getShard();
        Errors local;
        for (auto it = m_mesh.getCell(); it; ++it)
        {
          if (!shard.isOwned(D, it->getIndex()))
            continue;
          const Real area = it->getMeasure();
          for (size_t k = 0; k < wq.size(); ++k)
          {
            Math::SpatialPoint rc(2);
            rc(0) = xq[k][0];
            rc(1) = xq[k][1];
            const Point pt(*it, rc);
            f(pt, local, wq[k] * area);
          }
        }

        const auto& comm = m_mesh.getContext().getCommunicator();
        const auto sum = [&](Real v) { return boost::mpi::all_reduce(comm, v, std::plus<Real>()); };
        Errors out;
        out.u = sum(local.u);
        out.gradU = sum(local.gradU);
        out.p = sum(local.p);
        out.sigma = sum(local.sigma);
        out.uRef = sum(local.uRef);
        out.gradURef = sum(local.gradURef);
        out.pRef = sum(local.pRef);
        out.sigmaRef = sum(local.sigmaRef);
        return out;
      }

      /// @brief The channel (0, L) x (-H, H), nx x ny squares cut in two.
      static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        MPI::Sharder sharder(context);
        if (comm.rank() == RootRank)
        {
          Geometry::Mesh<Context::Local> mesh;
          mesh = mesh.UniformGrid(Polytope::Type::Triangle, { cfg.nx + 1, cfg.ny + 1 });
          for (auto it = mesh.getVertex(); it; ++it)
          {
            Math::SpatialPoint x = mesh.getVertexCoordinates(it->getIndex());
            x(0) = x(0) * cfg.length / static_cast<Real>(cfg.nx);
            x(1) = -cfg.halfHeight + x(1) * 2.0 * cfg.halfHeight / static_cast<Real>(cfg.ny);
            mesh.setVertexCoordinates(it->getIndex(), x);
          }
          mesh.flush();

          const size_t D = mesh.getDimension();
          mesh.getConnectivity().compute(D, D);
          mesh.getConnectivity().compute(D, 0);
          mesh.getConnectivity().compute(D, D - 1);
          mesh.getConnectivity().compute(D - 1, D);
          mesh.getConnectivity().compute(D - 1, 0);

          const Real eps = 1.0e-9 * cfg.length;
          std::vector<std::pair<Index, Attribute>> tags;
          for (auto it = mesh.getBoundary(); it; ++it)
          {
            Math::SpatialPoint c(2);
            c.setZero();
            for (const auto& v : it->getVertices())
              c += mesh.getVertexCoordinates(v);
            c /= static_cast<Real>(it->getVertices().size());
            Attribute a = Wall;
            if (c(0) < eps)
              a = Inlet;
            else if (c(0) > cfg.length - eps)
              a = Outlet;
            tags.emplace_back(it->getIndex(), a);
          }
          for (const auto& [f, a] : tags)
            mesh.setAttribute({ D - 1, f }, a);

#ifdef RODIN_USE_SCOTCH
          Scotch::Partitioner partitioner(mesh);
#else
          Geometry::BalancedCompactPartitioner partitioner(mesh);
#endif
          partitioner.partition(static_cast<size_t>(comm.size()));
          sharder.shard(partitioner);
          sharder.scatter(RootRank);
        }

        MeshType mesh = sharder.gather(RootRank);
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, 0);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        mesh.reconcile(1);
        return mesh;
      }

      void setInitialState()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real l0 = lambda0();
        const Real etaP = m_cfg.oldroydB.etaP;
        m_uOld.project(VectorFunction(dim, [this](const Point& p) -> Math::SpatialVector<Real> {
          const auto& x = p.getPhysicalCoordinates();
          const Eigen::Vector2d u = m_exact->velocity(x(0), x(1), 0.0);
          return Math::SpatialVector<Real>{{ u(0), u(1) }};
        }));
        m_psiOld.project(VectorFunction(size_t(3),
          [this, l0, etaP](const Point& p) -> Math::SpatialVector<Real> {
            const auto& x = p.getPhysicalCoordinates();
            return pack(logConformation(m_exact->stress(x(0), x(1), 0.0), l0, etaP));
          }));
        m_p.getSolution().project(RealFunction([this](const Point& p) {
          const auto& x = p.getPhysicalCoordinates();
          return m_exact->pressure(x(0), x(1), 0.0);
        }));
        m_uIt.setData(m_uOld.getData());
        m_u.getSolution().setData(m_uOld.getData());
        m_psiIt.setData(m_psiOld.getData());
        m_psi.getSolution().setData(m_psiOld.getData());
        ++m_psiRevision;
        m_subOld = Math::SpatialVector<Real>{{0.0, 0.0}};
        m_sub.get() = Math::SpatialVector<Real>{{0.0, 0.0}};
        updateConformation();
      }

      // ---- The parent's operators, unchanged ---------------------------------
      static Real cellSize(const Point& p)
      {
        return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      }

      Real tau1At(const Point& p) const
      {
        const auto uc = m_uOld.getValue(p);
        const Real h = cellSize(p);
        const Real nu = m_cfg.oldroydB.etaS / m_cfg.rho;
        return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
      }

      Real pttFactor(Real traceExp) const
      {
        const auto& ob = m_cfg.oldroydB;
        if (ob.pttEpsilon == 0.0)
          return 1.0;
        return 1.0 + ob.pttEpsilon * (ob.lambda / lambda0()) * (traceExp - 2.0);
      }

      Real alpha3At(const Point& p) const
      {
        const auto& ob = m_cfg.oldroydB;
        const auto uc = m_uOld.getValue(p);
        const Real h = cellSize(p);
        const Real speed = std::sqrt(Math::dot(uc, uc));
        const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
        const Real f = pttFactor(logEig(unpack(m_psiOld.getValue(p))).exp().trace());
        return 1.0 / (4.0 * f / (2.0 * ob.etaP) +
          0.25 * (ob.lambda * speed / (2.0 * ob.etaP * h) + ob.lambda * gradNorm / ob.etaP));
      }

      Real lambda0() const
      {
        const auto& ob = m_cfg.oldroydB;
        const Real l0 = std::max(ob.lambda0Factor * ob.lambda, ob.lambda0Min);
        if (!(l0 > 0.0))
          throw std::runtime_error("lambda0 = max(k lambda, lambda0_min) must be positive.");
        return l0;
      }

      void updateConformation()
      {
        const Real s = m_cfg.oldroydB.etaP / lambda0();
        const auto field = [](auto f) { return VectorFunction(size_t(3), f); };
        m_sigma.project(field([&](const Point& p) -> Math::SpatialVector<Real> {
          return pack(s * (logEig(unpack(m_psiIt.getValue(p))).exp() - Eigen::Matrix2d::Identity()));
        }));
      }

      static void axpy(Real a, const ::Vec& x, ::Vec& y)
      {
        PetscErrorCode ierr = VecAXPY(y, a, x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      void setupFlow()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const auto normal = BoundaryNormal(m_mesh);
        const Real rho = m_cfg.rho;
        const Real dt = m_dt;
        const auto& ob = m_cfg.oldroydB;
        const Real etaS = ob.etaS;
        const Real l0 = lambda0();
        const Real s = ob.etaP / l0;
        const Real chiScale = 2.0;
        const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };

        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto convU = Mult(Jacobian(m_u), m_uOld);
        const auto temam = Div(m_uOld) * Dot(m_u, m_v);
        const auto outletBackflow =
          0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);

        // Exact data at the current time m_t.
        VectorFunction inletVelocity(dim, [this](const Point& p) -> Math::SpatialVector<Real> {
          const auto& x = p.getPhysicalCoordinates();
          const Eigen::Vector2d u = m_exact->velocity(x(0), x(1), m_t);
          return Math::SpatialVector<Real>{{ u(0), u(1) }};
        });
        VectorFunction inletPsi(size_t(3), [this, l0](const Point& p) -> Math::SpatialVector<Real> {
          const auto& x = p.getPhysicalCoordinates();
          return pack(logConformation(m_exact->stress(x(0), x(1), m_t), l0, m_cfg.oldroydB.etaP));
        });
        // t = (2 eta_s eps(u) + sigma - p I) n, n = e_x on the outlet x = L.
        VectorFunction outletTraction(dim, [this, etaS](const Point& p) -> Math::SpatialVector<Real> {
          const auto& x = p.getPhysicalCoordinates();
          const Eigen::Matrix2d g = m_exact->velocityGradient(x(0), x(1), m_t);
          const Eigen::Matrix2d T = etaS * (g + g.transpose()) + m_exact->stress(x(0), x(1), m_t) -
            m_exact->pressure(x(0), x(1), m_t) * Eigen::Matrix2d::Identity();
          return Math::SpatialVector<Real>{{ T(0, 0), T(1, 0) }};
        });

        const auto streamV = Mult(Jacobian(m_v), m_uOld);
        const auto X = tensor(m_chi);
        const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(2, 2));

        const auto at = [this, model](auto pick) {
          return fromPoint(m_psiIt, m_uIt, model, m_psiRevision, pick); };
        const auto Dexp = [this, model](const auto& w) {
          return dexp(m_psiIt, m_uIt, model, m_psiRevision, w); };
        const auto JN = [this, model](const auto& w) {
          return jacobianN(m_psiIt, m_uIt, model, m_psiRevision, w); };
        const auto Gu = [this, model](const auto& v) {
          return velocityMap(m_psiIt, m_uIt, model, m_psiRevision, v); };
        const auto divExp = [&](const auto& psi) {
          return Mult(Dexp(Mult(Jacobian(psi), unit(0))), unit(0)) +
            Mult(Dexp(Mult(Jacobian(psi), unit(1))), unit(1));
        };

        const auto T = Dexp(m_psi);
        const auto Tn = at([](const LogPoint& q) { return q.exp; });
        const auto C = Tn - Dexp(m_psiIt);
        const auto N = at([](const LogPoint& q) { return q.N; });
        const auto known = N - JN(m_psiIt) - Gu(m_uIt) - (1.0 / dt) * tensor(m_psiOld)
          - advection(m_psiIt, m_uIt);

        m_flow =
          // ---- Momentum, sigma = s (exp(psi) - I) ----------------------------
            (rho / dt) * Integral(m_u, m_v) - (rho / dt) * Integral(m_uOld, m_v)
          + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)
          + 2.0 * etaS * Integral(symU, symV)
          + s * Integral(T, symV) + s * Integral(C, symV) - s * Integral(I, symV)
          - Integral(m_p, Div(m_v))
          // ---- Continuity -----------------------------------------------------
          + Integral(Div(m_u), m_q) + m_cfg.pressurePenalty * Integral(m_p, m_q)
          // ---- Constitutive, psi form, BDF1, Newton ------------------------------
          + Integral((1.0 / dt) * tensor(m_psi) + JN(m_psi) + advection(m_psi, m_uIt), X)
          + Integral(advection(m_psiIt, m_u) + Gu(m_u), X)
          + Integral(known, X)
          // ---- S1 ---------------------------------------------------------------
          + rho * rho * Integral(m_tauK * convU, streamV)
          - rho * rho * Integral(m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get()), streamV)
          + m_cfg.pressureScale * Integral(m_alpha1 * Grad(m_p), Grad(m_q))
          - m_cfg.pressureScale * Integral(m_alpha1 * m_piGradP.get(), Grad(m_q))
          + chiScale * m_cfg.stressDivScale * s * Integral(m_alpha1 * divExp(m_psi), divergence(m_chi))
          - chiScale * m_cfg.stressDivScale * Integral(m_alpha1 * m_piDivSigma.get(), divergence(m_chi))
          // ---- S2 ---------------------------------------------------------------
          + Integral(m_alpha2 * Div(m_u), Div(m_v))
          - Integral(m_alpha2 * m_piDiv.get(), Div(m_v))
          // ---- S3 ---------------------------------------------------------------
          + Integral(m_alpha3 * symU, symV)
          - Integral(m_alpha3 * tensor(m_piEps.get()), symV)
          + Integral(m_alphaPsi * advection(m_psi, m_uOld), advection(m_chi, m_uOld))
          - Integral(m_alphaPsi * tensor(m_piAdvPsi.get()), advection(m_chi, m_uOld))
          // ---- Boundary: exact traction on the outlet ------------------------------
          - BoundaryIntegral(Dot(outletTraction, m_v)).over(Outlet)
          + BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(Outlet)
          + DirichletBC(m_u, inletVelocity).on(Inlet)
          + DirichletBC(m_u, Zero(dim)).on(m_wallSet)
          + DirichletBC(m_psi, inletPsi).on(Inlet);
      }

      void updateStabilization()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real rho = m_cfg.rho;
        const Real dt = m_dt;
        const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
        const Real gradDiv = m_cfg.useVMS ? m_cfg.gradDivScale : 0.0;

        m_tauK.project(RealFunction([this, rho, dt, vms](const Point& p) {
          return vms / (rho / dt + rho / tau1At(p)); }));
        m_alpha1.project(RealFunction([this, rho](const Point& p) { return tau1At(p) / rho; }));
        m_alpha2.project(RealFunction([this, rho, gradDiv](const Point& p) {
          const Real h = cellSize(p);
          return gradDiv * rho * h * h / (4.0 * tau1At(p)); }));
        m_alpha3.project(RealFunction([this](const Point& p) {
          return m_cfg.stressScale * alpha3At(p); }));
        {
          const auto& ob = m_cfg.oldroydB;
          const Real k = ob.lambda / (2.0 * ob.etaP);
          const Real s = ob.etaP / lambda0();
          m_alphaPsi.project(RealFunction([this, k, s](const Point& p) {
            return 2.0 * k * k * s * m_cfg.stressScale * alpha3At(p); }));
        }

        m_piConv.project(Mult(Jacobian(m_uOld), m_uOld));
        m_sub.project(VectorFunction(dim,
          [this, rho, dt, dim](const Point& p) -> Math::SpatialVector<Real> {
            const auto conv = Mult(Jacobian(m_uOld), m_uOld).getValue(p);
            const auto proj = m_piConv.get().getValue(p);
            const auto old = m_subOld.getValue(p);
            const Real tau = m_tauK.getValue(p);
            Math::SpatialVector<Real> out(dim);
            for (Index c = 0; c < static_cast<Index>(dim); ++c)
              out(c) = tau * rho * (old(c) / dt - (conv(c) - proj(c)));
            return out;
          }));
        m_piDiv.project(Div(m_uOld));
        m_piGradP.project(Grad(m_p.getSolution()));
        const auto& ob = m_cfg.oldroydB;
        const Real l0 = lambda0();
        const Real s = ob.etaP / l0;
        const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };
        const auto Dexp = [this, model](const auto& w) {
          return dexp(m_psiOld, m_uOld, model, m_psiRevision, w); };
        const auto twoEpsN = Jacobian(m_uOld) + Transpose(Jacobian(m_uOld));
        m_piDivSigma.project(s * (Mult(Dexp(Mult(Jacobian(m_psiOld), unit(0))), unit(0)) +
          Mult(Dexp(Mult(Jacobian(m_psiOld), unit(1))), unit(1))));
        m_piEps.project(voigt(0.5 * twoEpsN));
        m_piAdvPsi.project(Mult(Jacobian(m_psiOld), m_uOld));
      }

      bool solveFlow()
      {
        updateStabilization();

        ::KSPConvergedReason reason = KSP_CONVERGED_ITS;
        PetscErrorCode ierr = PETSC_SUCCESS;
        ::Vec increment = PETSC_NULLPTR;
        ierr = VecDuplicate(m_psiIt.getData(), &increment);
        assert(ierr == PETSC_SUCCESS);
        ::Vec uIncrement = PETSC_NULLPTR;
        ierr = VecDuplicate(m_uIt.getData(), &uIncrement);
        assert(ierr == PETSC_SUCCESS);
        m_conformationIts = 0;
        m_uIt.setData(m_uOld.getData());
        for (int k = 0; k < std::max(1, m_cfg.conformationIterations); ++k)
        {
          m_flow.assemble();
          if (!m_flowFieldSplitsSet)
          {
            m_flow.setFieldSplits();
            m_flowFieldSplitsSet = true;
          }
          m_flow.solve(m_flowKSP);

          ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
          assert(ierr == PETSC_SUCCESS);
          if (reason <= 0)
            break;

          const ::Vec& solution = m_psi.getSolution().getData();
          ierr = VecWAXPY(increment, -1.0, m_psiIt.getData(), solution);
          assert(ierr == PETSC_SUCCESS);
          Real stepSize = 0.0;
          ierr = VecNorm(increment, NORM_INFINITY, &stepSize);
          assert(ierr == PETSC_SUCCESS);
          Real psiScale = 0.0;
          ierr = VecNorm(m_psiIt.getData(), NORM_INFINITY, &psiScale);
          assert(ierr == PETSC_SUCCESS);
          const Real omega = (std::isfinite(stepSize) && stepSize > m_cfg.newtonMaxStep)
            ? m_cfg.newtonMaxStep / stepSize : 1.0;
          m_psiIncrement = stepSize / std::max<Real>(1.0, psiScale);
          axpy(omega, increment, m_psiIt.getData());
          ierr = VecWAXPY(uIncrement, -1.0, m_uIt.getData(), m_u.getSolution().getData());
          assert(ierr == PETSC_SUCCESS);
          axpy(omega, uIncrement, m_uIt.getData());
          ++m_psiRevision;
          ++m_conformationIts;
          if (!std::isfinite(m_psiIncrement) ||
              (omega == 1.0 && m_psiIncrement < m_cfg.conformationTolerance))
            break;
        }
        ierr = VecDestroy(&increment);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDestroy(&uIncrement);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;

        m_uOld.setData(m_uIt.getData());
        m_u.getSolution().setData(m_uIt.getData());
        m_psiOld.setData(m_psiIt.getData());
        m_psi.getSolution().setData(m_psiIt.getData());
        ++m_psiRevision;
        updateConformation();
        if (m_cfg.useVMS)
          m_subOld.setData(m_sub.get().getData());

        m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
        return reason > 0 && std::isfinite(m_speed) && m_speed < 1.0e3;
      }

      Config m_cfg;
      std::shared_ptr<const ExactSolution> m_exact;
      MeshType m_mesh;
      AttributeSet m_wallSet;
      Real m_h = 0.0;

      VectorFESType m_vh;
      ScalarFESType m_sh;
      VectorFESType m_tauh;
      CellFESType m_ch;

      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld;
      VectorGridFunctionType m_uIt;
      VectorTrialFunctionType m_psi;
      VectorTestFunctionType m_chi;
      VectorGridFunctionType m_psiOld;
      VectorGridFunctionType m_psiIt;
      std::uint64_t m_psiRevision = 0;
      VectorGridFunctionType m_sigma;

      CellGridFunctionType m_tauK;
      CellGridFunctionType m_alpha1;
      CellGridFunctionType m_alpha2;
      CellGridFunctionType m_alpha3;
      CellGridFunctionType m_alphaPsi;
      Heart::OrthogonalProjection<VectorFESType> m_piConv;
      Heart::OrthogonalProjection<VectorFESType> m_sub;
      VectorGridFunctionType m_subOld;
      Heart::OrthogonalProjection<ScalarFESType> m_piDiv;
      Heart::OrthogonalProjection<VectorFESType> m_piGradP;
      Heart::OrthogonalProjection<VectorFESType> m_piDivSigma;
      Heart::OrthogonalProjection<VectorFESType> m_piEps;
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsi;

      FlowProblemType m_flow;
      Solver::KSP m_flowKSP;

      Real m_t = 0.0;
      Real m_dt = 0.0;
      Real m_speed = 0.0;
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;
      bool m_flowFieldSplitsSet = false;
  };

  // ==========================================================================
  // Studies
  // ==========================================================================
  struct Options
  {
      std::string study = "all";
      int levels = 4;
      int n0 = 2;            ///< space study: cells across 2H on the coarsest mesh
      Real dt0 = 1.0 / 24.0; ///< space study: dt on the coarsest mesh
      int n = 16;            ///< time study: cells across 2H
      int steps0 = 20;       ///< time study: steps per period on the coarsest level
      Real lambda = -1.0;    ///< < 0: the default of each study
      Real beta = -1.0;
      Real amplitude = 0.5;
      Real womersley = std::sqrt(2.0 * M_PI);
      ViscoelasticChannelTest::Config solver;
  };

  Real order(Real coarse, Real fine, Real ratio)
  {
    if (!(coarse > 0.0) || !(fine > 0.0))
      return std::numeric_limits<Real>::quiet_NaN();
    return std::log(coarse / fine) / std::log(ratio);
  }

  std::string ord(Real v)
  {
    std::ostringstream os;
    if (std::isnan(v))
      os << "   - ";
    else
      os << std::fixed << std::setprecision(2) << std::setw(5) << v;
    return os.str();
  }

  /// @brief Pacheco & Castillo, Sec. 5.1.4: h and dt halved together.
  int spaceStudy(const Context::MPI& context, const Options& opt, bool root)
  {
    const Real nu0 = 1.0, U = 1.0, H = 1.0, L = 3.0, tStar = 2.0, T = 10.0;
    const Real beta = opt.beta > 0.0 ? opt.beta : 0.5;
    const Real lambda = opt.lambda > 0.0 ? opt.lambda : 0.1;

    auto cfg = opt.solver;
    cfg.length = L;
    cfg.halfHeight = H;
    cfg.rho = 1.0;
    cfg.oldroydB.etaS = beta * nu0;
    cfg.oldroydB.etaP = (1.0 - beta) * nu0;
    cfg.oldroydB.lambda = lambda;
    auto exact = std::make_shared<RampedChannel>(
      U, H, L, cfg.oldroydB.etaS, cfg.oldroydB.etaP, lambda, tStar);

    if (root)
      std::cout << "\n[space] ramped channel (Pacheco & Castillo 2023, Sec. 5.1.4): nu0=" << nu0
                << " beta=" << beta << " lambda=" << lambda << " U=" << U << " H=" << H
                << " L=" << L << " t*=" << tStar << " T=" << T << "  We=lambda U/H="
                << lambda * U / H << "\n";

    std::vector<std::array<Real, 7>> rows;   // n, h, dt, eU, eH1, eP, eS
    for (int l = 0; l < opt.levels; ++l)
    {
      const int n = opt.n0 << l;
      cfg.ny = static_cast<size_t>(n);
      cfg.nx = static_cast<size_t>(std::lround(n * L / (2.0 * H)));
      const Real dt = opt.dt0 / static_cast<Real>(1 << l);
      ViscoelasticChannelTest test(context, cfg, exact);
      if (root)
        std::cout << "  level " << l << ": " << cfg.nx << "x" << cfg.ny << " squares, h="
                  << test.getMeshSize() << ", dt=" << dt << ", " << std::lround(T / dt)
                  << " steps ..." << std::endl;
      if (!test.run(dt, T))
        return 1;
      const auto e = test.errorToExact(T);
      rows.push_back({ Real(n), test.getMeshSize(), dt, e.relU(), e.relH1(), e.relP(), e.relSigma() });
    }

    if (root)
    {
      std::ofstream csv("vt_space.csv");
      csv << "n,h,dt,e_u_L2,e_u_H1,e_p_L2,e_sigma_L2\n";
      std::cout << "\n   n      h         dt        u L2      ord   u H1      ord   p L2      ord"
                   "   sigma L2  ord   (relative, t = T)\n";
      for (size_t i = 0; i < rows.size(); ++i)
      {
        const auto& r = rows[i];
        const auto o = [&](size_t k) {
          return i == 0 ? ord(std::numeric_limits<Real>::quiet_NaN())
                        : ord(order(rows[i - 1][k], r[k], rows[i - 1][1] / r[1]));
        };
        std::cout << std::setw(4) << int(r[0]) << "  " << std::scientific << std::setprecision(2)
                  << r[1] << "  " << r[2] << "  " << r[3] << ' ' << o(3) << "  " << r[4] << ' '
                  << o(4) << "  " << r[5] << ' ' << o(5) << "  " << r[6] << ' ' << o(6) << '\n';
        csv << r[0] << ',' << r[1] << ',' << r[2] << ',' << r[3] << ',' << r[4] << ',' << r[5]
            << ',' << r[6] << '\n';
      }
      std::cout << "Written vt_space.csv" << std::endl;
    }
    return 0;
  }

  /// @brief Exact pulsatile Oldroyd-B channel flow: dt halved on a fixed mesh.
  int timeStudy(const Context::MPI& context, const Options& opt, bool root)
  {
    const Real rho = 1.0, eta0 = 1.0, U = 1.0, h = 1.0, L = 2.0;
    const Real beta = opt.beta > 0.0 ? opt.beta : 0.5;
    const Real lambda = opt.lambda > 0.0 ? opt.lambda : 0.5;
    // Wo = h sqrt(omega rho/eta_0)
    const Real omega = opt.womersley * opt.womersley * eta0 / (rho * h * h);

    auto cfg = opt.solver;
    cfg.length = L;
    cfg.halfHeight = h;
    cfg.rho = rho;
    cfg.oldroydB.etaS = beta * eta0;
    cfg.oldroydB.etaP = (1.0 - beta) * eta0;
    cfg.oldroydB.lambda = lambda;
    cfg.ny = static_cast<size_t>(opt.n);
    cfg.nx = static_cast<size_t>(std::lround(opt.n * L / (2.0 * h)));
    auto exact = std::make_shared<MaxwellWomersley>(U, opt.amplitude, h, L, omega, rho,
      cfg.oldroydB.etaS, cfg.oldroydB.etaP, lambda);
    const Real T = exact->getPeriod();

    ViscoelasticChannelTest test(context, cfg, exact);
    if (root)
      std::cout << "\n[time] pulsatile Oldroyd-B channel (exact): beta=" << beta
                << " lambda=" << lambda << " Wo=" << opt.womersley << " A=" << opt.amplitude
                << "  period T=" << T << "  De=lambda/T=" << lambda / T
                << "  mesh " << cfg.nx << "x" << cfg.ny << " (h=" << test.getMeshSize()
                << "), one period from the exact state\n";

    const int finest = opt.steps0 << (opt.levels - 1);
    const Real dtRef = T / static_cast<Real>(4 * finest);
    if (root)
      std::cout << "  reference: dt = " << dtRef << " ..." << std::endl;
    if (!test.run(dtRef, T))
      return 1;
    const auto reference = test.snapshot();
    const auto refExact = test.errorToExact(T);

    // dt, steps, exact (u, H1, p, sigma), reference (u, H1, p, sigma)
    std::vector<std::array<Real, 10>> rows;
    for (int l = 0; l < opt.levels; ++l)
    {
      const int steps = opt.steps0 << l;
      const Real dt = T / static_cast<Real>(steps);
      if (root)
        std::cout << "  level " << l << ": dt=" << dt << " (" << steps << " steps) ..." << std::endl;
      if (!test.run(dt, T))
        return 1;
      const auto ex = test.errorToExact(T);
      const auto rf = test.errorTo(reference);
      rows.push_back({ dt, Real(steps), ex.relU(), ex.relH1(), ex.relP(), ex.relSigma(),
        rf.relU(), rf.relH1(), rf.relP(), rf.relSigma() });
    }

    if (root)
    {
      std::ofstream csv("vt_time.csv");
      csv << "dt,steps,ex_u_L2,ex_u_H1,ex_p_L2,ex_sigma_L2,ref_u_L2,ref_u_H1,ref_p_L2,ref_sigma_L2\n";
      std::cout << "\n  errors against the reference (temporal error alone), relative, t = T\n"
                << "   dt        u L2      ord   u H1      ord   p L2      ord   sigma L2  ord\n";
      for (size_t i = 0; i < rows.size(); ++i)
      {
        const auto& r = rows[i];
        const auto o = [&](size_t k) {
          return i == 0 ? ord(std::numeric_limits<Real>::quiet_NaN())
                        : ord(order(rows[i - 1][k], r[k], rows[i - 1][0] / r[0]));
        };
        std::cout << "  " << std::scientific << std::setprecision(2) << r[0] << "  " << r[6]
                  << ' ' << o(6) << "  " << r[7] << ' ' << o(7) << "  " << r[8] << ' ' << o(8)
                  << "  " << r[9] << ' ' << o(9) << '\n';
      }
      std::cout << "\n  errors against the exact solution (temporal + spatial), relative, t = T\n"
                << "   dt        u L2      u H1      p L2      sigma L2\n";
      for (const auto& r : rows)
        std::cout << "  " << std::scientific << std::setprecision(2) << r[0] << "  " << r[2]
                  << "  " << r[3] << "  " << r[4] << "  " << r[5] << '\n';
      std::cout << "  spatial floor (reference vs exact): " << refExact.relU() << "  "
                << refExact.relH1() << "  " << refExact.relP() << "  " << refExact.relSigma()
                << '\n';
      for (const auto& r : rows)
      {
        for (size_t k = 0; k < r.size(); ++k)
          csv << (k ? "," : "") << r[k];
        csv << '\n';
      }
      std::cout << "Written vt_time.csv" << std::endl;
    }
    return 0;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setPETScDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };
  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");
  setPETScDefault("-mat_mumps_icntl_20", "0");
  setPETScDefault("-mat_mumps_icntl_21", "0");
  // The projections are mass solves: tight, so they do not pollute the rates.
  setPETScDefault("-vt_proj_ksp_rtol", "1e-12");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);
  const bool root = world.rank() == 0;

  int status = 0;
  try
  {
    using namespace Rodin::Examples::ViscoelasticFluids;
    Options opt;

    const auto getReal = [](const char* name, Rodin::Real& out) {
      PetscBool set = PETSC_FALSE;
      PetscReal value = 0.0;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &set);
      if (set)
        out = value;
    };
    const auto getInt = [](const char* name, int& out) {
      PetscBool set = PETSC_FALSE;
      PetscInt value = 0;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &set);
      if (set)
        out = static_cast<int>(value);
    };
    const auto getBool = [](const char* name, bool& out) {
      PetscBool set = PETSC_FALSE;
      PetscBool value = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &set);
      if (set)
        out = (value == PETSC_TRUE);
    };
    {
      char buffer[64];
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-vt_study", buffer, sizeof(buffer), &set);
      if (set)
        opt.study = buffer;
    }
    getInt("-vt_levels", opt.levels);
    getInt("-vt_n0", opt.n0);
    getReal("-vt_dt0", opt.dt0);
    getInt("-vt_n", opt.n);
    getInt("-vt_steps0", opt.steps0);
    getReal("-vt_lambda", opt.lambda);
    getReal("-vt_beta", opt.beta);
    getReal("-vt_amplitude", opt.amplitude);
    getReal("-vt_wo", opt.womersley);
    getInt("-vt_conformation_its", opt.solver.conformationIterations);
    getReal("-vt_conformation_tol", opt.solver.conformationTolerance);
    getReal("-vt_newton_max_step", opt.solver.newtonMaxStep);
    getReal("-vt_backflow", opt.solver.outletBackflowStabilization);
    getBool("-vt_vms", opt.solver.useVMS);
    getReal("-vt_vms_scale", opt.solver.vmsScale);
    getReal("-vt_graddiv_scale", opt.solver.gradDivScale);
    getReal("-vt_pressure_scale", opt.solver.pressureScale);
    getReal("-vt_stress_div_scale", opt.solver.stressDivScale);
    getReal("-vt_stress_scale", opt.solver.stressScale);

    if (root)
      std::cout << "Log-conformation Oldroyd-B (psi form, Newton, P1/P1/P1, BDF1, split OSS): "
                   "convergence tests" << std::endl;

    if (opt.study == "space" || opt.study == "all")
      status = spaceStudy(context, opt, root);
    if (status == 0 && (opt.study == "time" || opt.study == "all"))
      status = timeStudy(context, opt, root);
    if (opt.study != "space" && opt.study != "time" && opt.study != "all")
      throw std::runtime_error("-vt_study must be space, time or all.");
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_viscoelastic_implicit failed: " << e.what() << "\n";
    status = 1;
  }

  PetscFinalize();
  return status;
}
