/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file test_kolmogorov_viscoelastic.cpp
 * @brief Two-dimensional viscoelastic Kolmogorov flow (Berti & Boffetta, Phys.
 *        Rev. E 82, 036314, 2010) with the fully-implicit log-conformation
 *        Oldroyd-B solver of Heart/LeftAtrium2D_viscous_logarithmic_implicit.
 *
 * The flow solver is the parent one, term by term (psi form, Newton,
 * P1/P1/P1, BDF1, lagged split-OSS), with epsilon = 0 (Oldroyd-B). What is new:
 * the doubly periodic box [0, 2 pi]^2, the Kolmogorov forcing f = (F cos(y/L),
 * 0), and the diagnostics of the paper.
 *
 * Periodicity. Rodin's PeriodicBC is not consumed by any assembler (the PETSc
 * LinearSystem::merge throws "Unimplemented"), but every assembler, PETSc
 * included, applies identification DirichletBCs through the same
 * ConstraintMap (A <- T^T A T plus reconstruction rows). PeriodicIdentification
 * below is such a DirichletBCBase: it ties each DOF of a vertex on x = 2 pi or
 * y = 2 pi to the same DOF of its image on x = 0 or y = 0 (corners to the
 * origin), for u, p and psi. The OSS projections carry the same constraint.
 * The pressure is fixed at the vertex (pi, pi). With several MPI ranks the
 * master of a slave vertex usually lives on another rank, so every rank
 * gathers (all_gather) the global DOFs of the boundary vertices, indexed by
 * their grid position, and registers all slave -> master pairs: the
 * assembler needs them for ghost slaves too. The pairs are O(n) per field.
 *
 * Paper -> Rodin. Their conformation sigma relaxes as -2(sigma - 1)/tau and
 * the polymer stress is (2 eta nu/tau)(sigma - 1), i.e. Oldroyd-B with
 *   lambda = tau/2, eta_s = nu, eta_p = eta nu, rho = 1,
 * and c = exp(psi) is their sigma (lambda0 = lambda). U0 = 4 (their Table I)
 * and L = 1/4: the base flow has four periods in the box (their Figs. 5, 7)
 * and only then does their ordering omega_A <~ omega* < omega_tau << omega_Gamma
 * hold. nu_t = nu (1 + eta) = U0 L/Re0, tau = Wi0 L/U0, F = nu_t U0/L^2.
 * The laminar fixed point, their Eq. (3), is
 *   u = U0 cos(y/L) e_x,  c11 = 1 + tau^2 U0^2 sin^2(y/L)/(2 L^2),
 *   c12 = -tau U0 sin(y/L)/(2L),  c22 = 1,  p = const,
 * with <K> = U0^2/4, <tr c> = 2 + Wi^2/4 and <c11> - <c22> = Wi^2/4.
 * The paper does not give eta; 0.2 is used here (-kf_eta).
 *
 * Cases (-kf_case):
 *   steady      Wi0 = 8 < Wi_c ~ 10, from rest, marched until
 *               max|u^{n+1} - u^n|/(dt U0) < steady_tol; the relative L2
 *               errors of u and c against the laminar fixed point are printed.
 *               n = 48, dt = 0.1, T <= 15.
 *   transition  Wi0 = 16 > Wi_c, from the laminar fixed point plus a small
 *               periodic, divergence-free perturbation, as in the paper; the
 *               time series of K and Sigma (their Figs. 1, 3, 4) are written.
 *               n = 64, dt = 0.02, T = 100 (5000 steps).
 *
 * Cost. Each step is a Newton loop (4-6 iterations) on the parent's forms;
 * about 87% of the time is their quadrature, 10% the MUMPS factorisations.
 * On one core: ~15 s per step at n = 48 and ~30 s at n = 64, so the steady
 * case takes well under an hour and the transition case more than a day;
 * both parts scale with the number of MPI ranks (1.8x on 2 ranks at n = 32).
 * Memory at n = 64 (25k unknowns, six coupled fields) is ~300 MB on one
 * rank, dominated by the monolithic LU of the flow (~12M factor entries).
 * The mesh must resolve the forcing period 2 pi L: with fewer than about 8
 * cells per period (n < 32) the P1 interpolant of cos(y/L) degenerates
 * (n = 8 puts it at the Nyquist mode, where D(u) is L2-orthogonal to P1 and
 * psi stays exactly zero). A line is printed every -kf_sample_every steps.
 *
 * Output: kf_<case>_Wi<Wi0>.csv with t, K = <|u|^2>/2, Sigma = <tr c>,
 * N1 = <c11> - <c22>, <c12>, U = 2 <u_x cos(y/L)> (the measured mean-flow
 * amplitude, from which the a-posteriori Wi and Re follow), rms(u_y), and
 * XDMF snapshots (velocity, logConformation, trace, vorticity) every
 * -kf_xdmf_every steps. plot_kolmogorov.py compares them with the paper.
 *
 * Run (any number of ranks; results agree across ranks to round-off):
 *   mpirun -n 4 ./examples/viscoelastic_fluids/test_kolmogorov_viscoelastic -kf_case steady
 *   mpirun -n 4 ./examples/viscoelastic_fluids/test_kolmogorov_viscoelastic -kf_case transition
 *   python3 ../examples/viscoelastic_fluids/plot_kolmogorov.py kf_*.csv
 *
 * Options (prefix -kf_): case, wi (Wi0), re (Re0), eta, L, n (cells per side),
 * dt, T, steady_tol, perturbation, seed, sample_every, xdmf_every,
 * conformation_its, conformation_tol, newton_max_step, vms, vms_scale,
 * graddiv_scale, pressure_scale, stress_div_scale, stress_scale.
 */
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

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

    /// @brief The parts of a LogPoint, each computed only when a form asks
    ///        for it: the forms are assembled one integrator at a time over
    ///        the whole mesh, so each pass recomputes the point and should pay
    ///        only for what it reads (the FD Jacobian JN alone costs six
    ///        eigenproblems).
    enum LogPart : std::uint8_t
    {
      LogExp = 1,  ///< exp(psi)
      LogD = 2,    ///< Dexp[B_j]
      LogN = 4,    ///< N, and f
      LogJN = 8,   ///< dN/dpsi_j
      LogQ = 16    ///< Q_ab
    };

    /// @brief logPoint(psi(x), grad u(x)), cached per quadrature point: the
    ///        form language evaluates a coefficient once per basis function,
    ///        and without the cache the eigenproblems dominate assembly. Keyed
    ///        by the psi field, its revision (bumped on every change of the
    ///        pair), cell and reference coordinates. Only the parts in `need`
    ///        are guaranteed.
    template <class GF, class UF>
    const LogPoint& logPointAt(const GF& psi, const UF& u, const LogModel& m,
      std::uint64_t revision, const Point& p, std::uint8_t need)
    {
      struct Entry
      {
        const void* field = nullptr;
        std::uint64_t revision = 0;
        Index cell = 0;
        Real r0 = 0, r1 = 0;
        std::uint8_t have = 0;
        bool hasEig = false;
        bool hasL = false;
        Eigen::Matrix2d psi;
        Eigen::Matrix2d L;
        LogEig eig;
        LogPoint value;
      };
      thread_local std::array<Entry, 16> cache;
      thread_local size_t next = 0;

      const Index cell = p.getPolytope().getIndex();
      const auto& rc = p.getReferenceCoordinates();
      Entry* hit = nullptr;
      for (auto& e : cache)
        if (e.field == &psi && e.revision == revision && e.cell == cell &&
            e.r0 == rc(0) && e.r1 == rc(1))
        {
          hit = &e;
          break;
        }

      if (!hit)
      {
        hit = &cache[next];
        next = (next + 1) % cache.size();
        hit->field = &psi;
        hit->revision = revision;
        hit->cell = cell;
        hit->r0 = rc(0);
        hit->r1 = rc(1);
        hit->have = 0;
        hit->hasEig = false;
        hit->hasL = false;
        hit->psi = unpack(psi.getValue(p));
      }

      Entry& e = *hit;
      const std::uint8_t missing = need & ~e.have;
      if (!missing)
        return e.value;

      LogPoint& out = e.value;
      if ((missing & (LogN | LogJN)) && !e.hasL)
      {
        const auto J = Jacobian(u).getValue(p);   // Math::SpatialMatrix, not Eigen
        e.L << J(0, 0), J(0, 1), J(1, 0), J(1, 1);
        e.hasL = true;
      }
      if ((missing & (LogExp | LogD | LogQ)) && !e.hasEig)
      {
        e.eig = logEig(e.psi);
        e.hasEig = true;
      }
      if ((missing & (LogExp | LogQ)) && !(e.have & LogExp))
      {
        out.exp = e.eig.exp();
        e.have |= LogExp;
      }
      if (missing & LogD)
      {
        const auto& B = voigtBasisEigen();
        for (size_t j = 0; j < 3; ++j)
          out.D[j] = e.eig.dexp(B[j]);
      }
      if (missing & LogN)
        out.N = nonlinearMap(e.psi, e.L, m, &out.f);
      if (missing & LogJN)
      {
        // Central differences in the three Voigt directions; the map is smooth
        // and the eigenproblem is 2x2, so this costs six decompositions.
        const auto& B = voigtBasisEigen();
        const Real h = 1.0e-6 * std::max<Real>(1.0, e.psi.cwiseAbs().maxCoeff());
        for (size_t j = 0; j < 3; ++j)
          out.JN[j] = (nonlinearMap(e.psi + h * B[j], e.L, m) -
            nonlinearMap(e.psi - h * B[j], e.L, m)) / (2.0 * h);
      }
      if (missing & LogQ)
      {
        for (size_t a = 0; a < 2; ++a)
          for (size_t b = 0; b < 2; ++b)
          {
            Eigen::Matrix2d E = Eigen::Matrix2d::Zero();
            E(a, b) = 1.0;
            out.Q[a][b] = -e.eig.dexpInverse(E * out.exp + out.exp * E.transpose() -
              0.5 * m.c * (E + E.transpose()));
          }
      }
      e.have |= need;
      return out;
    }

    /// @brief A pointwise matrix drawn from the cached LogPoint by `pick`,
    ///        which reads only the parts in `need`.
    template <class GF, class UF, class Pick>
    auto fromPoint(const GF& psi, const UF& u, const LogModel& m,
      const std::uint64_t& revision, std::uint8_t need, Pick pick)
    {
      return PointwiseTensor([&psi, &u, m, &revision, need, pick](const Point& p) -> Eigen::Matrix2d {
        return pick(logPointAt(psi, u, m, revision, p, need)); });
    }

    /// @brief Dexp_{psi(x)}[w] = sum_j w_j Dexp[B_j] for a Voigt w: psi itself
    ///        or a derivative of it (chain rule); trial, test or field.
    template <class GF, class UF, class W>
    auto dexp(const GF& psi, const UF& u, const LogModel& m, const std::uint64_t& revision, const W& w)
    {
      const auto D = [&](size_t j) {
        return fromPoint(psi, u, m, revision, LogD, [j](const LogPoint& q) { return q.D[j]; });
      };
      return D(0) * Component(w, 0) + D(1) * Component(w, 1) + D(2) * Component(w, 2);
    }

    /// @brief JN_{psi(x)}[w] = sum_j w_j dN/dpsi_j.
    template <class GF, class UF, class W>
    auto jacobianN(const GF& psi, const UF& u, const LogModel& m, const std::uint64_t& revision, const W& w)
    {
      const auto J = [&](size_t j) {
        return fromPoint(psi, u, m, revision, LogJN, [j](const LogPoint& q) { return q.JN[j]; });
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
        return fromPoint(psi, u, m, revision, LogQ, [a, b](const LogPoint& q) { return q.Q[a][b]; });
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
  // Periodicity as identification constraints
  // ==========================================================================
  /**
   * @brief u_s = u_m for every (slave, master) pair of global DOF indices.
   *
   * @details A DirichletBCBase in identification mode, so that every
   *          assembler distributes the slave rows into the master rows
   *          (A <- T^T A T) and writes u_s - u_m = 0 as the slave row. The
   *          pairs are global (KolmogorovFlow::makeDOFPairs) and identical on
   *          every rank; the assembler writes each slave row on its owner.
   */
  template <class Trial>
  class PeriodicIdentification final : public DirichletBCBase<Real>
  {
    public:
      using Parent = DirichletBCBase<Real>;
      using Pairs = std::vector<std::pair<Index, Index>>;

      PeriodicIdentification(Trial& u, std::shared_ptr<const Pairs> pairs)
        : m_u(u), m_pairs(std::move(pairs))
      {}

      PeriodicIdentification(const PeriodicIdentification& other)
        : Parent(other), m_u(other.m_u), m_pairs(other.m_pairs), m_dofs(other.m_dofs)
      {}

      void assemble() override
      {
        IdentifiedDOFs dofs;
        for (const auto& [slave, master] : *m_pairs)
        {
          IndexArray masters(1);
          masters(0) = master;
          Math::Vector<Real> coeffs(1);
          coeffs(0) = 1.0;
          dofs[slave] = std::make_pair(masters, coeffs);
        }
        m_dofs = std::move(dofs);
      }

      const DOFs& getDOFs() const override { return m_dofs; }
      bool isComponent() const override { return false; }
      const FormLanguage::Base& getOperand() const override { return m_u.get(); }
      const FormLanguage::Base& getValue() const override { return m_u.get(); }
      Optional<Identifiable::UUID> getValueUUID() const override { return m_u.get().getUUID(); }
      PeriodicIdentification* copy() const noexcept override { return new PeriodicIdentification(*this); }

    private:
      std::reference_wrapper<Trial> m_u;
      std::shared_ptr<const Pairs> m_pairs;
      DOFs m_dofs{ IdentifiedDOFs{} };
  };

  /// @brief u = 0 at one global DOF (removes the pressure constant).
  template <class Trial>
  class PinnedDOF final : public DirichletBCBase<Real>
  {
    public:
      using Parent = DirichletBCBase<Real>;

      PinnedDOF(Trial& u, Index dof) : m_u(u), m_dof(dof) {}
      PinnedDOF(const PinnedDOF& other)
        : Parent(other), m_u(other.m_u), m_dof(other.m_dof), m_dofs(other.m_dofs)
      {}

      void assemble() override
      {
        ValueDOFs dofs;
        dofs[m_dof] = 0.0;
        m_dofs = std::move(dofs);
      }

      const DOFs& getDOFs() const override { return m_dofs; }
      bool isComponent() const override { return false; }
      const FormLanguage::Base& getOperand() const override { return m_u.get(); }
      const FormLanguage::Base& getValue() const override { return m_u.get(); }
      PinnedDOF* copy() const noexcept override { return new PinnedDOF(*this); }

    private:
      std::reference_wrapper<Trial> m_u;
      Index m_dof;
      DOFs m_dofs{ ValueDOFs{} };
  };

  /**
   * @brief Heart/OrthogonalProjection.h with the periodic constraint: the
   *        L2 projection onto the periodic P1 space.
   *
   * @details The identification rows make the mass system nonsymmetric;
   *          the default solve is GMRES + Jacobi (options prefix kf_proj_).
   */
  template <class FES>
  class PeriodicProjection
  {
    public:
      using GridFunctionType = PETSc::Variational::GridFunction<FES>;
      using TrialFunctionType = PETSc::Variational::TrialFunction<GridFunctionType, FES>;
      using TestFunctionType = PETSc::Variational::TestFunction<FES>;
      using ProblemType = Problem<PETSc::Math::LinearSystem, TrialFunctionType, TestFunctionType>;
      using Pairs = typename PeriodicIdentification<TrialFunctionType>::Pairs;

      PeriodicProjection(const FES& fes, std::shared_ptr<const Pairs> pairs)
        : m_trial(fes), m_test(fes), m_problem(m_trial, m_test), m_ksp(m_problem),
          m_projection(fes), m_pairs(std::move(pairs))
      {
        m_ksp.setPrefix("kf_proj_");
      }

      PeriodicProjection(const PeriodicProjection&) = delete;
      PeriodicProjection& operator=(const PeriodicProjection&) = delete;

      template <class Expression>
      PeriodicProjection& project(const Expression& f)
      {
        m_problem = Integral(m_trial, m_test) - Integral(f, m_test)
          + PeriodicIdentification<TrialFunctionType>(m_trial, m_pairs);
        m_problem.solve(m_ksp);
        m_projection.setData(m_trial.getSolution().getData());
        return *this;
      }

      const GridFunctionType& get() const { return m_projection; }
      GridFunctionType& get() { return m_projection; }

    private:
      TrialFunctionType m_trial;
      TestFunctionType m_test;
      ProblemType m_problem;
      Solver::KSP m_ksp;
      GridFunctionType m_projection;
      std::shared_ptr<const Pairs> m_pairs;
  };

  // ==========================================================================
  // The Kolmogorov flow
  // ==========================================================================
  class KolmogorovFlow
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;
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
      using Pairs = std::vector<std::pair<Index, Index>>;   ///< (slave, master) global DOFs

      struct Config
      {
          std::string name = "steady";
          Real U0 = 4.0;        ///< laminar amplitude (their Table I)
          Real L = 0.25;        ///< forcing scale, f ~ cos(y/L): four periods in 2 pi
          Real wi0 = 8.0;       ///< Wi0 = tau U0/L
          Real re0 = 1.0;       ///< Re0 = U0 L/nu_t
          Real eta = 0.2;       ///< nu_p/nu
          size_t n = 64;        ///< cells per side of [0, 2 pi]
          Real dt = 0.01;
          Real T = 30.0;
          bool fromRest = true;
          Real perturbation = 1.0e-3;  ///< relative amplitude of the initial perturbation
          unsigned seed = 1;
          Real steadyTol = 0.0;        ///< > 0: stop once max|du|/(dt U0) is below it
          int sampleEvery = 5;
          int xdmfEvery = 0;
          // The parent's solver parameters.
          Real pressurePenalty = 1.0e-12;
          Real vmsScale = 1.0, gradDivScale = 1.0, pressureScale = 1.0;
          Real stressDivScale = 1.0, stressScale = 1.0;
          bool useVMS = true;
          int conformationIterations = 10;
          Real conformationTolerance = 1.0e-8;
          Real newtonMaxStep = 2.0;
      };

      KolmogorovFlow(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_comm(context.getCommunicator()),
          m_mesh(makeMesh(context, cfg)),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_tauh(std::integral_constant<size_t, 1>{}, m_mesh,
            m_mesh.getSpaceDimension() * (m_mesh.getSpaceDimension() + 1) / 2),
          m_ch(m_mesh),
          m_vPairs(std::make_shared<Pairs>(makeDOFPairs(m_vh))),
          m_sPairs(std::make_shared<Pairs>(makeDOFPairs(m_sh))),
          m_tPairs(std::make_shared<Pairs>(makeDOFPairs(m_tauh))),
          m_pinnedPressure(globalVertexDOF(m_sh, m_cfg.n / 2, m_cfg.n / 2)),
          m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh),
          m_uOld(m_vh), m_uIt(m_vh), m_uPrev(m_vh),
          m_psi(m_tauh), m_chi(m_tauh), m_psiOld(m_tauh), m_psiIt(m_tauh),
          m_trace(m_sh), m_vorticity(m_sh),
          m_tauK(m_ch), m_alpha1(m_ch), m_alpha2(m_ch), m_alpha3(m_ch), m_alphaPsi(m_ch),
          m_piConv(m_vh, m_vPairs), m_sub(m_vh, m_vPairs), m_subOld(m_vh),
          m_piDiv(m_sh, m_sPairs), m_piGradP(m_vh, m_vPairs),
          m_piDivSigma(m_vh, m_vPairs), m_piEps(m_tauh, m_tPairs), m_piAdvPsi(m_tauh, m_tPairs),
          m_flow(m_u, m_p, m_psi, m_v, m_q, m_chi),
          m_flowKSP(m_flow),
          m_xdmf(context.getCommunicator(), "kf_" + cfg.name)
      {
        // Paper -> Rodin (file header).
        const Real nuT = m_cfg.U0 * m_cfg.L / m_cfg.re0;
        m_nu = nuT / (1.0 + m_cfg.eta);
        m_etaP = m_cfg.eta * m_nu;
        m_tau = m_cfg.wi0 * m_cfg.L / m_cfg.U0;
        m_lambda = 0.5 * m_tau;
        m_force = nuT * m_cfg.U0 / (m_cfg.L * m_cfg.L);

        if (m_comm.rank() == RootRank)
          std::cout << "[kolmogorov] " << m_cfg.name << ": Wi0=" << m_cfg.wi0 << " Re0=" << m_cfg.re0
                    << " eta=" << m_cfg.eta << " U0=" << m_cfg.U0 << " L=" << m_cfg.L
                    << "  ->  nu=" << m_nu << " nu_p=" << m_etaP << " tau=" << m_tau
                    << " (lambda=" << m_lambda << ") F=" << m_force << "\n"
                    << "  mesh " << m_cfg.n << "x" << m_cfg.n << " on [0, 2pi]^2, "
                    << m_sPairs->size() << " periodic vertex pairs, " << m_comm.size()
                    << " rank(s), dt=" << m_cfg.dt
                    << ", T=" << m_cfg.T << ", laminar <K>=" << 0.25 * m_cfg.U0 * m_cfg.U0
                    << " <tr c>=" << 2.0 + 0.25 * m_cfg.wi0 * m_cfg.wi0 << std::endl;
      }

      int run()
      {
        setInitialState();
        setupFlow();

        m_trace.setName("trace");
        m_vorticity.setName("vorticity");
        m_u.setName("velocity");
        m_psi.setName("logConformation");
        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("logConformation", m_psi.getSolution());
        m_xdmf.add("trace", m_trace);
        m_xdmf.add("vorticity", m_vorticity);

        std::ostringstream csvName;
        csvName << "kf_" << m_cfg.name << "_Wi" << m_cfg.wi0 << ".csv";
        const bool root = m_comm.rank() == RootRank;
        std::ofstream csv;
        if (root)
        {
          csv.open(csvName.str());
          csv << "t,K,Sigma,N1,c12,U,uy_rms,newtonIts\n";
        }

        const int steps = static_cast<int>(std::lround(m_cfg.T / m_cfg.dt));
        Real change = 0.0;
        for (int n = 0; n <= steps; ++n)
        {
          if (n > 0)
          {
            m_t += m_cfg.dt;
            m_uPrev.setData(m_uOld.getData());
            if (!solveFlow())
            {
              if (root)
                std::cerr << "[kolmogorov] the flow solve failed at t = " << m_t << "\n";
              return 1;
            }
            change = maxChange() / (m_cfg.dt * m_cfg.U0);
          }

          const bool last = (n == steps) ||
            (m_cfg.steadyTol > 0.0 && n > 0 && change < m_cfg.steadyTol);
          if (n % std::max(1, m_cfg.sampleEvery) == 0 || last)
          {
            const auto d = diagnostics();
            if (root)
            {
              csv << m_t << ',' << d.K << ',' << d.Sigma << ',' << d.N1 << ',' << d.c12 << ','
                  << d.U << ',' << d.uyRms << ',' << m_conformationIts << '\n';
              csv.flush();
              std::cout << "  t=" << std::fixed << std::setprecision(3) << m_t
                        << std::scientific << std::setprecision(4) << "  K=" << d.K
                        << "  Sigma=" << d.Sigma << "  N1=" << d.N1 << "  U=" << d.U
                        << "  rms(uy)=" << d.uyRms << "  |du|/(dt U0)=" << change
                        << "  newton=" << m_conformationIts << std::endl;
            }
          }
          if (m_cfg.xdmfEvery > 0 && (n % m_cfg.xdmfEvery == 0 || last))
          {
            updateOutputFields();
            m_xdmf.write(m_t).flush();
          }
          if (last)
          {
            if (root && m_cfg.steadyTol > 0.0 && n < steps)
              std::cout << "  steady state reached at t = " << m_t << std::endl;
            break;
          }
        }
        m_xdmf.close();

        const auto e = errorToLaminar();
        const auto d = diagnostics();
        const Real Wi = m_tau * d.U / m_cfg.L;
        if (root)
          std::cout << "[kolmogorov] final: <K>=" << d.K << " (laminar " << 0.25 * m_cfg.U0 * m_cfg.U0
                    << ")  <tr c>=" << d.Sigma << " (laminar " << 2.0 + 0.25 * m_cfg.wi0 * m_cfg.wi0
                    << ")  N1=" << d.N1 << " (laminar " << 0.25 * m_cfg.wi0 * m_cfg.wi0 << ")\n"
                    << "  a posteriori: U=" << d.U << "  Wi=tau U/L=" << Wi
                    << "  Re=U L/nu_t=" << d.U * m_cfg.L / (m_nu * (1.0 + m_cfg.eta))
                    << "\n  relative L2 distance to the laminar fixed point: u " << e.first
                    << "  c " << e.second << "\nWritten " << csvName.str() << std::endl;
        return 0;
      }

    private:
      struct Diagnostics
      {
          Real K = 0, Sigma = 0, N1 = 0, c12 = 0, U = 0, uyRms = 0;
      };

      // ---- Mesh and periodic pairs ------------------------------------------
      static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        MPI::Sharder sharder(context);
        if (comm.rank() == RootRank)
        {
          Geometry::Mesh<Context::Local> mesh;
          mesh = mesh.UniformGrid(Polytope::Type::Triangle, { cfg.n + 1, cfg.n + 1 });
          const Real h = 2.0 * M_PI / static_cast<Real>(cfg.n);
          for (auto it = mesh.getVertex(); it; ++it)
          {
            Math::SpatialPoint x = mesh.getVertexCoordinates(it->getIndex());
            x(0) *= h;
            x(1) *= h;
            mesh.setVertexCoordinates(it->getIndex(), x);
          }
          mesh.flush();
          const size_t D = mesh.getDimension();
          mesh.getConnectivity().compute(D, D);
          mesh.getConnectivity().compute(D, 0);
          mesh.getConnectivity().compute(D, D - 1);
          mesh.getConnectivity().compute(D - 1, D);
          mesh.getConnectivity().compute(D - 1, 0);
          Geometry::BalancedCompactPartitioner partitioner(mesh);
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

      /// @brief Grid position (i, j) of a vertex of the n x n grid.
      std::pair<long, long> gridIndex(Index vertex) const
      {
        const Real h = 2.0 * M_PI / static_cast<Real>(m_cfg.n);
        const auto& x = m_mesh.getVertexCoordinates(vertex);
        return { std::lround(x(0) / h), std::lround(x(1) / h) };
      }

      /**
       * @brief Global (slave, master) DOF pairs of @p fes: every DOF of a
       *        vertex on x = 2 pi or y = 2 pi tied to the same DOF of its
       *        image on x = 0 or y = 0 (corners to the origin).
       *
       * Each rank contributes the global DOFs of its boundary vertices, owned
       * or ghost, keyed by grid position; the table is gathered on every rank,
       * so each rank knows all pairs, including those whose master lives on
       * another rank.
       */
      template <class FES>
      Pairs makeDOFPairs(const FES& fes) const
      {
        const long n = static_cast<long>(m_cfg.n);
        const size_t k = fes.getVectorDimension();
        std::vector<long> local;   // records [grid, dof_0, ..., dof_{k-1}]
        for (auto it = m_mesh.getVertex(); it; ++it)
        {
          const auto [i, j] = gridIndex(it->getIndex());
          if (i != 0 && j != 0 && i != n && j != n)
            continue;
          const IndexArray dofs = fes.getDOFs(0, it->getIndex());   // copy: thread-local buffer
          assert(static_cast<size_t>(dofs.size()) == k);
          local.push_back(j * (n + 1) + i);
          for (size_t c = 0; c < k; ++c)
            local.push_back(static_cast<long>(dofs(static_cast<Index>(c))));
        }
        std::vector<std::vector<long>> all;
        boost::mpi::all_gather(m_comm, local, all);
        std::vector<long> table(static_cast<size_t>((n + 1) * (n + 1)) * k, -1);
        for (const auto& rec : all)
          for (size_t r = 0; r < rec.size(); r += k + 1)
            for (size_t c = 0; c < k; ++c)
              table[static_cast<size_t>(rec[r]) * k + c] = rec[r + 1 + c];

        Pairs pairs;
        for (long j = 0; j <= n; ++j)
          for (long i = 0; i <= n; ++i)
          {
            if (i < n && j < n)
              continue;
            const size_t gs = static_cast<size_t>(j * (n + 1) + i);
            const size_t gm = static_cast<size_t>((j % n) * (n + 1) + (i % n));
            for (size_t c = 0; c < k; ++c)
            {
              const long ds = table[gs * k + c], dm = table[gm * k + c];
              if (ds < 0 || dm < 0)
                throw std::runtime_error("makeDOFPairs: a boundary vertex is missing.");
              pairs.emplace_back(static_cast<Index>(ds), static_cast<Index>(dm));
            }
          }
        return pairs;
      }

      /// @brief Global DOF of the scalar space @p fes at grid vertex (i, j).
      template <class FES>
      Index globalVertexDOF(const FES& fes, size_t i, size_t j) const
      {
        long dof = -1;
        for (auto it = m_mesh.getVertex(); it; ++it)
        {
          const auto [gi, gj] = gridIndex(it->getIndex());
          if (gi == static_cast<long>(i) && gj == static_cast<long>(j))
          {
            dof = static_cast<long>(fes.getDOFs(0, it->getIndex())(0));
            break;
          }
        }
        long out = -1;
        boost::mpi::all_reduce(m_comm, dof, out, boost::mpi::maximum<long>());
        if (out < 0)
          throw std::runtime_error("globalVertexDOF: vertex not found.");
        return static_cast<Index>(out);
      }

      // ---- Laminar fixed point, their Eq. (3), in Rodin variables ---------------
      Eigen::Vector2d laminarVelocity(Real, Real y) const
      {
        return Eigen::Vector2d(m_cfg.U0 * std::cos(y / m_cfg.L), 0.0);
      }

      Eigen::Matrix2d laminarConformation(Real, Real y) const
      {
        const Real a = m_tau * m_cfg.U0 * std::sin(y / m_cfg.L) / m_cfg.L;
        Eigen::Matrix2d c;
        c << 1.0 + 0.5 * a * a, -0.5 * a, -0.5 * a, 1.0;
        return c;
      }

      static Eigen::Matrix2d logSPD(const Eigen::Matrix2d& c)
      {
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(c);
        return eig.eigenvectors() * eig.eigenvalues().array().log().matrix().asDiagonal() *
          eig.eigenvectors().transpose();
      }

      void setInitialState()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        if (m_cfg.fromRest)
        {
          m_uOld = Math::SpatialVector<Real>{{0.0, 0.0}};
          m_psiOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
        }
        else
        {
          // Laminar fixed point plus a divergence-free, periodic perturbation
          // u' = (d_y phi, -d_x phi), phi = sum of random low modes.
          std::mt19937 rng(m_cfg.seed);
          std::uniform_real_distribution<Real> phase(0.0, 2.0 * M_PI), amp(-1.0, 1.0);
          struct Mode { int kx, ky; Real a, ph; };
          std::vector<Mode> modes;
          for (int kx = 1; kx <= 4; ++kx)
            for (int ky = -4; ky <= 4; ++ky)
              modes.push_back({ kx, ky, amp(rng), phase(rng) });
          const Real eps = m_cfg.perturbation * m_cfg.U0;
          m_uOld.project(VectorFunction(dim, [this, modes, eps](const Point& p) -> Math::SpatialVector<Real> {
            const auto& x = p.getPhysicalCoordinates();
            Eigen::Vector2d u = laminarVelocity(x(0), x(1));
            for (const auto& m : modes)
            {
              const Real arg = m.kx * x(0) + m.ky * x(1) + m.ph;
              const Real k = std::hypot(Real(m.kx), Real(m.ky));
              u(0) += eps * m.a * m.ky * std::cos(arg) / k;
              u(1) -= eps * m.a * m.kx * std::cos(arg) / k;
            }
            return Math::SpatialVector<Real>{{ u(0), u(1) }};
          }));
          m_psiOld.project(VectorFunction(size_t(3), [this](const Point& p) -> Math::SpatialVector<Real> {
            const auto& x = p.getPhysicalCoordinates();
            return pack(logSPD(laminarConformation(x(0), x(1))));
          }));
        }
        m_p.getSolution() = Real(0);
        m_uIt.setData(m_uOld.getData());
        m_u.getSolution().setData(m_uOld.getData());
        m_psiIt.setData(m_psiOld.getData());
        m_psi.getSolution().setData(m_psiOld.getData());
        ++m_psiRevision;
        m_subOld = Math::SpatialVector<Real>{{0.0, 0.0}};
        m_sub.get() = Math::SpatialVector<Real>{{0.0, 0.0}};
      }

      // ---- The parent's operators, unchanged ----------------------------------
      static Real cellSize(const Point& p)
      {
        return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      }

      Real tau1At(const Point& p) const
      {
        const auto uc = m_uOld.getValue(p);
        const Real h = cellSize(p);
        return 1.0 / (4.0 * m_nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
      }

      Real alpha3At(const Point& p) const
      {
        const auto uc = m_uOld.getValue(p);
        const Real h = cellSize(p);
        const Real speed = std::sqrt(Math::dot(uc, uc));
        const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
        const Real f = 1.0;  // Oldroyd-B: no PTT factor
        return 1.0 / (4.0 * f / (2.0 * m_etaP) +
          0.25 * (m_lambda * speed / (2.0 * m_etaP * h) + m_lambda * gradNorm / m_etaP));
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
        const Real rho = 1.0;
        const Real dt = m_cfg.dt;
        const Real etaS = m_nu;
        const Real l0 = m_lambda;             // lambda0 = lambda: c = exp(psi)
        const Real s = m_etaP / l0;
        const Real chiScale = 2.0;
        const LogModel model{ m_lambda, l0, 0.0, 2.0 * (1.0 - l0 / m_lambda) };

        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto convU = Mult(Jacobian(m_u), m_uOld);
        const auto temam = Div(m_uOld) * Dot(m_u, m_v);

        const Real F = m_force, L = m_cfg.L;
        VectorFunction forcing(dim, [F, L](const Point& p) -> Math::SpatialVector<Real> {
          return Math::SpatialVector<Real>{{ F * std::cos(p.getPhysicalCoordinates()(1) / L), 0.0 }};
        });

        const auto streamV = Mult(Jacobian(m_v), m_uOld);
        const auto X = tensor(m_chi);
        const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(2, 2));

        const auto at = [this, model](std::uint8_t need, auto pick) {
          return fromPoint(m_psiIt, m_uIt, model, m_psiRevision, need, pick); };
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
        const auto Tn = at(LogExp, [](const LogPoint& q) { return q.exp; });
        const auto C = Tn - Dexp(m_psiIt);
        const auto N = at(LogN, [](const LogPoint& q) { return q.N; });
        const auto known = N - JN(m_psiIt) - Gu(m_uIt) - (1.0 / dt) * tensor(m_psiOld)
          - advection(m_psiIt, m_uIt);

        m_flow =
          // ---- Momentum, sigma = s (exp(psi) - I), Kolmogorov forcing -----------
            (rho / dt) * Integral(m_u, m_v) - (rho / dt) * Integral(m_uOld, m_v)
          + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)
          + 2.0 * etaS * Integral(symU, symV)
          + s * Integral(T, symV) + s * Integral(C, symV) - s * Integral(I, symV)
          - Integral(m_p, Div(m_v))
          - Integral(Dot(forcing, m_v))
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
          // ---- Doubly periodic box, pressure fixed at the centre ------------------
          + PeriodicIdentification<VectorTrialFunctionType>(m_u, m_vPairs)
          + PeriodicIdentification<ScalarTrialFunctionType>(m_p, m_sPairs)
          + PeriodicIdentification<VectorTrialFunctionType>(m_psi, m_tPairs)
          + PinnedDOF<ScalarTrialFunctionType>(m_p, m_pinnedPressure);
      }

      void updateStabilization()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real rho = 1.0;
        const Real dt = m_cfg.dt;
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
          const Real k = m_lambda / (2.0 * m_etaP);
          const Real s = m_etaP / m_lambda;
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
        const Real s = m_etaP / m_lambda;
        const LogModel model{ m_lambda, m_lambda, 0.0, 0.0 };
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
        if (m_cfg.useVMS)
          m_subOld.setData(m_sub.get().getData());

        const Real speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
        return reason > 0 && std::isfinite(speed) && speed < 100.0 * m_cfg.U0;
      }

      // ---- Diagnostics -----------------------------------------------------------
      Real maxChange()
      {
        ::Vec d = PETSC_NULLPTR;
        PetscErrorCode ierr = VecDuplicate(m_uOld.getData(), &d);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecWAXPY(d, -1.0, m_uPrev.getData(), m_uOld.getData());
        assert(ierr == PETSC_SUCCESS);
        Real out = 0.0;
        ierr = VecNorm(d, NORM_INFINITY, &out);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDestroy(&d);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        return out;
      }

      /// @brief Sum of f over this rank's owned cells, degree-5 Dunavant
      ///        rule (the callers reduce over ranks).
      template <class F>
      void quadrature(F&& f)
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
            f(pt, wq[k] * area);
          }
        }
      }

      Diagnostics diagnostics()
      {
        const auto& uSol = m_u.getSolution();
        const auto& psiSol = m_psi.getSolution();
        Diagnostics d;
        Real uy2 = 0.0;
        quadrature([&](const Point& pt, Real w) {
          const auto u = uSol.getValue(pt);
          const Eigen::Matrix2d c = logEig(unpack(psiSol.getValue(pt))).exp();
          const Real y = pt.getPhysicalCoordinates()(1);
          d.K += w * 0.5 * (u(0) * u(0) + u(1) * u(1));
          d.Sigma += w * c.trace();
          d.N1 += w * (c(0, 0) - c(1, 1));
          d.c12 += w * c(0, 1);
          d.U += w * 2.0 * u(0) * std::cos(y / m_cfg.L);
          uy2 += w * u(1) * u(1);
        });
        {
          const std::vector<Real> local = { d.K, d.Sigma, d.N1, d.c12, d.U, uy2 };
          std::vector<Real> sum(local.size(), 0.0);
          boost::mpi::all_reduce(m_comm, local.data(), static_cast<int>(local.size()),
            sum.data(), std::plus<Real>());
          d.K = sum[0]; d.Sigma = sum[1]; d.N1 = sum[2]; d.c12 = sum[3]; d.U = sum[4];
          uy2 = sum[5];
        }
        const Real area = 4.0 * M_PI * M_PI;
        d.K /= area;
        d.Sigma /= area;
        d.N1 /= area;
        d.c12 /= area;
        d.U /= area;
        d.uyRms = std::sqrt(uy2 / area);
        return d;
      }

      /// @brief Relative L2 distance of (u, c) to the laminar fixed point.
      std::pair<Real, Real> errorToLaminar()
      {
        const auto& uSol = m_u.getSolution();
        const auto& psiSol = m_psi.getSolution();
        Real eu = 0, nu = 0, ec = 0, nc = 0;
        quadrature([&](const Point& pt, Real w) {
          const auto& x = pt.getPhysicalCoordinates();
          const auto u = uSol.getValue(pt);
          const Eigen::Vector2d ue = laminarVelocity(x(0), x(1));
          const Eigen::Matrix2d c = logEig(unpack(psiSol.getValue(pt))).exp();
          const Eigen::Matrix2d ce = laminarConformation(x(0), x(1));
          eu += w * (Eigen::Vector2d(u(0), u(1)) - ue).squaredNorm();
          nu += w * ue.squaredNorm();
          ec += w * (c - ce).squaredNorm();
          nc += w * ce.squaredNorm();
        });
        const std::vector<Real> local = { eu, nu, ec, nc };
        std::vector<Real> sum(local.size(), 0.0);
        boost::mpi::all_reduce(m_comm, local.data(), static_cast<int>(local.size()),
          sum.data(), std::plus<Real>());
        eu = sum[0]; nu = sum[1]; ec = sum[2]; nc = sum[3];
        return { std::sqrt(eu / nu), std::sqrt(ec / nc) };
      }

      void updateOutputFields()
      {
        m_trace.project(RealFunction([this](const Point& p) {
          return logEig(unpack(m_psiIt.getValue(p))).exp().trace(); }));
        const auto J = Jacobian(m_u.getSolution());
        m_vorticity.project(RealFunction([J](const Point& p) {
          const auto g = J.getValue(p);
          return g(1, 0) - g(0, 1); }));
      }

      Config m_cfg;
      boost::mpi::communicator m_comm;
      MeshType m_mesh;
      Real m_nu = 0, m_etaP = 0, m_tau = 0, m_lambda = 0, m_force = 0;

      VectorFESType m_vh;
      ScalarFESType m_sh;
      VectorFESType m_tauh;
      CellFESType m_ch;
      std::shared_ptr<const Pairs> m_vPairs, m_sPairs, m_tPairs;   ///< global DOF pairs
      Index m_pinnedPressure;

      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      VectorGridFunctionType m_uOld, m_uIt, m_uPrev;
      VectorTrialFunctionType m_psi;
      VectorTestFunctionType m_chi;
      VectorGridFunctionType m_psiOld, m_psiIt;
      std::uint64_t m_psiRevision = 0;
      ScalarGridFunctionType m_trace, m_vorticity;

      CellGridFunctionType m_tauK, m_alpha1, m_alpha2, m_alpha3, m_alphaPsi;
      PeriodicProjection<VectorFESType> m_piConv;
      PeriodicProjection<VectorFESType> m_sub;
      VectorGridFunctionType m_subOld;
      PeriodicProjection<ScalarFESType> m_piDiv;
      PeriodicProjection<VectorFESType> m_piGradP;
      PeriodicProjection<VectorFESType> m_piDivSigma;
      PeriodicProjection<VectorFESType> m_piEps;
      PeriodicProjection<VectorFESType> m_piAdvPsi;

      FlowProblemType m_flow;
      Solver::KSP m_flowKSP;
      IO::XDMF m_xdmf;

      Real m_t = 0.0;
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;
      bool m_flowFieldSplitsSet = false;
  };
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
  // Workspace margin (%) over MUMPS's analysis estimate. The identification
  // rows are nonsymmetric, so pivoting adds fill beyond the estimate; with
  // several ranks the default margin can fail with INFOG(1) = -9.
  setPETScDefault("-mat_mumps_icntl_14", "100");
  // Periodic projections: mass matrices plus identification rows, which are
  // not symmetric but well conditioned. GMRES + Jacobi converges in < 20
  // iterations and avoids keeping seven LU factors in memory.
  setPETScDefault("-kf_proj_ksp_type", "gmres");
  setPETScDefault("-kf_proj_pc_type", "jacobi");
  setPETScDefault("-kf_proj_ksp_rtol", "1e-10");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Flow = Rodin::Examples::ViscoelasticFluids::KolmogorovFlow;
    Flow::Config cfg;

    std::string which = "steady";
    {
      char buffer[64];
      PetscBool set = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-kf_case", buffer, sizeof(buffer), &set);
      if (set)
        which = buffer;
    }
    if (which == "steady")
    {
      cfg.name = "steady";
      cfg.wi0 = 8.0;
      cfg.fromRest = true;
      cfg.n = 48;        // 12 cells per forcing period
      cfg.T = 15.0;
      cfg.dt = 0.1;      // only the fixed point matters: BDF1 with a large step
      cfg.steadyTol = 1.0e-6;
    }
    else if (which == "transition")
    {
      cfg.name = "transition";
      cfg.wi0 = 16.0;
      cfg.fromRest = false;
      cfg.n = 64;        // 16 cells per forcing period
      cfg.T = 100.0;
      cfg.dt = 0.02;
      cfg.steadyTol = 0.0;
      cfg.xdmfEvery = 250;
    }
    else
      throw std::runtime_error("-kf_case must be steady or transition.");

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
    int n = static_cast<int>(cfg.n), seed = static_cast<int>(cfg.seed);
    getReal("-kf_wi", cfg.wi0);
    getReal("-kf_re", cfg.re0);
    getReal("-kf_eta", cfg.eta);
    getReal("-kf_L", cfg.L);
    getInt("-kf_n", n);
    getReal("-kf_dt", cfg.dt);
    getReal("-kf_T", cfg.T);
    getReal("-kf_steady_tol", cfg.steadyTol);
    getReal("-kf_perturbation", cfg.perturbation);
    getInt("-kf_seed", seed);
    getInt("-kf_sample_every", cfg.sampleEvery);
    getInt("-kf_xdmf_every", cfg.xdmfEvery);
    getInt("-kf_conformation_its", cfg.conformationIterations);
    getReal("-kf_conformation_tol", cfg.conformationTolerance);
    getReal("-kf_newton_max_step", cfg.newtonMaxStep);
    getBool("-kf_vms", cfg.useVMS);
    getReal("-kf_vms_scale", cfg.vmsScale);
    getReal("-kf_graddiv_scale", cfg.gradDivScale);
    getReal("-kf_pressure_scale", cfg.pressureScale);
    getReal("-kf_stress_div_scale", cfg.stressDivScale);
    getReal("-kf_stress_scale", cfg.stressScale);
    cfg.n = static_cast<size_t>(n);
    cfg.seed = static_cast<unsigned>(seed);

    Flow flow(context, cfg);
    status = flow.run();
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_kolmogorov_viscoelastic failed: " << e.what() << "\n";
    status = 1;
  }

  PetscFinalize();
  return status;
}
