/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LeftAtrium2D_viscous_logarithmic_implicit.cpp
 * @brief Driver for the 2D left atrium with sPTT/Oldroyd-B blood in the
 *        fully-implicit log-conformation form; see the header.
 *
 * Run (from the build directory, mesh flattened as for LeftAtrium2D):
 *   mpirun -n 7 ./examples/Heart/LeftAtrium2D_viscous_logarithmic_implicit \
 *     -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -la2d_dt 1e-3
 *
 * Options: those of LeftAtrium2D_viscous_logarithmic (-la2d_ptt_epsilon,
 *          -la2d_lambda0_factor, -la2d_lambda0_min, -la2d_conformation_its,
 *          -la2d_conformation_tol) plus -la2d_newton_max_step and
 *          -la2d_inlet_conformation.
 */
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "LeftAtrium2D_viscous_logarithmic_implicit.h"
#include "CoronaryArtery/CoronaryArteryAlerts.h"
#include "CoronaryArtery/CoronaryArteryTiming.h"

namespace Rodin::Examples::Heart
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

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
  // PressureWaveform
  // ==========================================================================
  void PressureWaveform::load(const std::string& path)
  {
    std::ifstream file(path);
    if (!file)
      throw std::runtime_error("Failed to open the pressure file " + path);

    m_path = path;
    m_t.clear();
    m_p.clear();

    std::string line;
    while (std::getline(file, line))
    {
      const auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || line[first] == '#')
        continue;

      std::istringstream row(line.substr(first));
      Real t = 0.0;
      Real p = 0.0;
      if (!(row >> t >> p))
        throw std::runtime_error("Malformed line in " + path + ": " + line);

      if (!m_t.empty() && !(t > m_t.back()))
        throw std::runtime_error("Non-increasing time column in " + path);

      m_t.push_back(t);
      m_p.push_back(p);
    }

    if (m_t.size() < 2)
      throw std::runtime_error("Fewer than two samples in " + path);

    m_period = m_t.back() - m_t.front();
    if (!(m_period > 0.0))
      throw std::runtime_error("Zero period in " + path);

    m_min = *std::min_element(m_p.begin(), m_p.end());
    m_max = *std::max_element(m_p.begin(), m_p.end());

    // Trapezoidal mean over one period, so the reported figure is the mean of
    // the signal and not the mean of the samples: the file is uniformly
    // sampled, but nothing here requires it to be.
    Real integral = 0.0;
    for (size_t i = 0; i + 1 < m_t.size(); ++i)
      integral += 0.5 * (m_p[i] + m_p[i + 1]) * (m_t[i + 1] - m_t[i]);
    m_mean = integral / m_period;
  }

  PressureWaveform::Real PressureWaveform::operator()(Real t) const
  {
    assert(m_t.size() >= 2);

    const Real t0 = m_t.front();
    Real tau = t - t0;
    tau -= m_period * std::floor(tau / m_period);
    const Real x = t0 + tau;

    // upper_bound, then step back: the samples are strictly increasing, so the
    // bracketing interval is [it - 1, it) and the index arithmetic cannot run
    // off either end.
    const auto it = std::upper_bound(m_t.begin(), m_t.end(), x);
    if (it == m_t.begin())
      return m_p.front();
    if (it == m_t.end())
      return m_p.back();

    const size_t hi = static_cast<size_t>(it - m_t.begin());
    const size_t lo = hi - 1;
    const Real dt = m_t[hi] - m_t[lo];
    const Real s = (x - m_t[lo]) / dt;
    return (1.0 - s) * m_p[lo] + s * m_p[hi];
  }

  // ==========================================================================
  // LeftAtrium2DViscousLogImplicit
  // ==========================================================================
  LeftAtrium2DViscousLogImplicit::AttributeSet LeftAtrium2DViscousLogImplicit::makeInletSet(const Config& cfg)
  {
    return AttributeSet(cfg.labels.inlets.begin(), cfg.labels.inlets.end());
  }

  LeftAtrium2DViscousLogImplicit::AttributeSet LeftAtrium2DViscousLogImplicit::makeWallSet(const Config& cfg)
  {
    AttributeSet out(cfg.labels.wall.begin(), cfg.labels.wall.end());
    out.insert(cfg.labels.appendage.begin(), cfg.labels.appendage.end());
    return out;
  }

  LeftAtrium2DViscousLogImplicit::LeftAtrium2DViscousLogImplicit(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_inletSet(makeInletSet(m_cfg)),
      m_wallSet(makeWallSet(m_cfg)),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_tauh(std::integral_constant<size_t, 1>{}, m_mesh,
        m_mesh.getSpaceDimension() * (m_mesh.getSpaceDimension() + 1) / 2),
      m_ch(m_mesh),
      m_u(m_vh),
      m_p(m_sh),
      m_v(m_vh),
      m_q(m_sh),
      m_uOld(m_vh),
      m_uIt(m_vh),
      m_psi(m_tauh),
      m_chi(m_tauh),
      m_psiOld(m_tauh),
      m_psiIt(m_tauh),
      m_sigma(m_tauh),
      m_wTrial(m_vh),
      m_wTest(m_vh),
      m_tauK(m_ch),
      m_alpha1(m_ch),
      m_alpha2(m_ch),
      m_alpha3(m_ch),
      m_alphaPsi(m_ch),
      m_piConv(m_vh, "la2d_proj_"),
      m_sub(m_vh, "la2d_proj_"),
      m_subOld(m_vh),
      m_piDiv(m_sh, "la2d_proj_"),
      m_piGradP(m_vh, "la2d_proj_"),
      m_piDivSigma(m_vh, "la2d_proj_"),
      m_piEps(m_tauh, "la2d_proj_"),
      m_piAdvPsi(m_tauh, "la2d_proj_"),
      m_th(m_sh),
      m_fg(m_sh),
      m_fn(m_sh),
      m_vth(m_sh),
      m_vfg(m_sh),
      m_vfn(m_sh),
      m_thCur(m_sh),
      m_fgCur(m_sh),
      m_fnCur(m_sh),
      m_thPrev(m_sh),
      m_fgPrev(m_sh),
      m_fnPrev(m_sh),
      m_wss(m_vh),
      m_symRec0(m_vh),
      m_symRec1(m_vh),
      m_netShear(m_vh),
      m_absShear(m_sh),
      m_shearMagnitude(m_sh),
      m_tawss(m_sh),
      m_osi(m_sh),
      m_activation(m_sh),
      m_qFlux(m_sh),
      m_one(m_sh),
      m_flux(m_qFlux),
      m_flow(m_u, m_p, m_psi, m_v, m_q, m_chi),
      m_flowKSP(m_flow),
      m_species(m_th, m_fg, m_fn, m_vth, m_vfg, m_vfn),
      m_speciesKSP(m_species),
      m_vectorProjection(m_wTrial, m_wTest),
      m_vectorProjectionKSP(m_vectorProjection),
      m_wssTrial(m_vh),
      m_wssTest(m_vh),
      m_wssProjection(m_wssTrial, m_wssTest),
      m_wssKSP(m_wssProjection)
  {
    m_inletWave.load(m_cfg.inletPressurePath);
    m_outletWave.load(m_cfg.outletPressurePath);

    // Queried on every rank, printed on one: nothing that may reduce over the
    // communicator belongs inside a root-only branch.
    const auto cellCount = m_mesh.getCellCount();
    const auto vertexCount = m_mesh.getVertexCount();
    const auto velocityDOFs = m_vh.getSize();
    const auto pressureDOFs = m_sh.getSize();

    if (isRoot())
    {
      Alert::Info() << "[mesh] cells=" << cellCount << " vertices=" << vertexCount
                    << " velocity DOFs=" << velocityDOFs
                    << " pressure DOFs=" << pressureDOFs << Alert::Raise;

      const auto describe = [](const char* what, const PressureWaveform& w) {
        Alert::Info() << "[" << what << "] " << w.getPath()
                      << "  samples=" << w.getSampleCount() << "  T=" << w.getPeriod()
                      << " s  mean=" << w.getMean() << " Pa  range=[" << w.getMinimum()
                      << ", " << w.getMaximum() << "] Pa" << Alert::Raise;
      };
      describe("p_pv", m_inletWave);
      describe("p_mv", m_outletWave);

      const auto& ob = m_cfg.oldroydB;
      Alert::Info() << "[model] psi-form log-conformation, Newton, BDF1; "
                    << (ob.pttEpsilon > 0.0 ? "linear sPTT (xi = 0)" : "Oldroyd-B")
                    << "  eta_s=" << ob.etaS << "  eta_p=" << ob.etaP
                    << "  lambda=" << ob.lambda << " s  lambda0=" << lambda0()
                    << " s  epsilon=" << ob.pttEpsilon
                    << "  1/(2 lambda)=" << (0.5 / ob.lambda)
                    << " 1/s (Oldroyd-B extensional threshold)" << Alert::Raise;
    }

    // The two files must share a period, and the run's cycle length must be
    // that period: the wall-shear indices are accumulated over one cycle, and
    // a cycle that is not a period of the forcing averages two different
    // phases of the flow into the same TAWSS.
    const Real tolerance = 1.0e-9;
    if (std::abs(m_inletWave.getPeriod() - m_outletWave.getPeriod()) > tolerance)
      throw std::runtime_error("The inlet and outlet waveforms have different periods.");

    if (std::abs(m_inletWave.getPeriod() - m_cfg.period) > 1.0e-6)
    {
      if (isRoot())
        Alert::Warning() << "[waveform] Config::period = " << m_cfg.period
                         << " s does not match the file period "
                         << m_inletWave.getPeriod() << " s; taking the file's."
                         << Alert::Raise;
      m_cfg.period = m_inletWave.getPeriod();
    }
  }

  LeftAtrium2DViscousLogImplicit::~LeftAtrium2DViscousLogImplicit() = default;

  bool LeftAtrium2DViscousLogImplicit::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  LeftAtrium2DViscousLogImplicit::MeshType LeftAtrium2DViscousLogImplicit::makeMesh(
    const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 2 || mesh.getDimension() != 2)
        throw std::runtime_error(
          "LeftAtrium2D expects a planar triangular mesh written as MEDIT "
          "\"Dimension 2\". RectLAA.mesh is stored as \"Dimension 3\" "
          "with z = 0; flatten it first with "
          "examples/Heart/LeftAtrium2D/make_la2d_mesh.py.");

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);

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
    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.reconcile(1);

    return mesh;
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = m_cfg.oldroydB.etaS / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::pttFactor(Real traceExp) const
  {
    const auto& ob = m_cfg.oldroydB;
    if (ob.pttEpsilon == 0.0)
      return 1.0;
    return 1.0 + ob.pttEpsilon * (ob.lambda / lambda0()) * (traceExp - 2.0);
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::alpha3At(const Point& p) const
  {
    const auto& ob = m_cfg.oldroydB;
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real speed = std::sqrt(Math::dot(uc, uc));
    const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
    // f at psi^n, from the nodal exp (P1 interpolation is enough for a P0
    // parameter).
    const Real f = pttFactor(logEig(unpack(m_psiOld.getValue(p))).exp().trace());
    return 1.0 / (4.0 * f / (2.0 * ob.etaP) +
      0.25 * (ob.lambda * speed / (2.0 * ob.etaP * h) + ob.lambda * gradNorm / ob.etaP));
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::lambda0() const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = std::max(ob.lambda0Factor * ob.lambda, ob.lambda0Min);
    if (!(l0 > 0.0))
      throw std::runtime_error("lambda0 = max(k lambda, lambda0_min) must be positive.");
    return l0;
  }

  void LeftAtrium2DViscousLogImplicit::updateConformation()
  {
    // Nodal (P1 interpolation), from psi^k; output and wall traction only.
    const Real s = m_cfg.oldroydB.etaP / lambda0();
    const auto field = [](auto f) { return VectorFunction(size_t(3), f); };

    m_sigma.project(field([&](const Point& p) -> Math::SpatialVector<Real> {
      return pack(s * (logEig(unpack(m_psiIt.getValue(p))).exp() - Eigen::Matrix2d::Identity()));
    }));
  }

  void LeftAtrium2DViscousLogImplicit::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  LeftAtrium2DViscousLogImplicit& LeftAtrium2DViscousLogImplicit::initialize()
  {
    setupSpaces();
    setupFlow();
    setupSpecies();
    setupWallShear();

    if (isRoot())
    {
      m_csv.open(m_cfg.csvPath);
      if (!m_csv)
        throw std::runtime_error("Failed to open " + m_cfg.csvPath);
      writeCSVHeader();
    }

    m_initialized = true;
    return *this;
  }

  void LeftAtrium2DViscousLogImplicit::setupSpaces()
  {
    const auto zeroVector = Math::SpatialVector<Real>{{0.0, 0.0}};
    const auto zeroTensor = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_uOld = zeroVector;
    m_uIt = zeroVector;
    m_subOld = zeroVector;
    m_sub.get() = zeroVector;
    m_psiOld = zeroTensor;
    m_psiIt = zeroTensor;
    ++m_psiRevision;
    m_wss = zeroVector;
    m_netShear = zeroVector;
    m_symRec0 = zeroVector;
    m_symRec1 = zeroVector;

    m_absShear = Real(0);
    m_shearMagnitude = Real(0);
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);
    m_one = Real(1);

    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;
    m_thCur = Real(m_cfg.inletThrombin);
    m_thPrev = Real(m_cfg.inletThrombin);
    m_fnCur = Real(m_cfg.inletFibrin);
    m_fnPrev = Real(m_cfg.inletFibrin);
    m_fgCur = fg0;
    m_fgPrev = fg0;

    configureMassSolver(m_vectorProjectionKSP, "la2d_vproj_");
    configureMassSolver(m_wssKSP, "la2d_wss_");

    updateConformation();  // psi = 0: tau = I, sigma = 0

    m_u.setName("velocity");
    m_p.setName("pressure");
    m_psi.setName("logConformation");
    m_sigma.setName("elasticStress");
    m_thCur.setName("thrombin");
    m_fgCur.setName("fibrinogen");
    m_fnCur.setName("fibrin");
    m_tawss.setName("TAWSS");
    m_osi.setName("OSI");
    m_activation.setName("activation");
    m_wss.setName("shearStress");

    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("logConformation", m_psi.getSolution());
    m_xdmf.add("elasticStress", m_sigma);
    m_xdmf.add("thrombin", m_thCur);
    m_xdmf.add("fibrinogen", m_fgCur);
    m_xdmf.add("fibrin", m_fnCur);
    m_xdmf.add("TAWSS", m_tawss);
    m_xdmf.add("OSI", m_osi);
    m_xdmf.add("activation", m_activation);
    m_xdmf.add("shearStress", m_wss);

    // In 2D these are lengths, not areas: the "flux" through a boundary is a
    // volumetric flow per unit depth, m^2/s.
    m_outletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.outlet});
    m_inletMeasure = boundaryMeasure(m_inletSet);

    m_pIn = m_inletWave(0.0) + m_cfg.pressureOffset;
    m_pOut = m_outletWave(0.0) + m_cfg.pressureOffset;
    m_outletPressure = m_pOut;

    // The velocity scale the forcing can account for. With rigid walls, no body
    // force and both ends on prescribed pressure, sqrt(2 max|dp| / rho) is the
    // whole budget: a peak far above it has no source in the data, and says the
    // scheme is making energy rather than that the atrium is doing something
    // interesting. Printed next to the divergence guard so the two can be read
    // against each other.
    {
      Real maxDp = 0.0;
      const int samples = 2000;
      for (int i = 0; i <= samples; ++i)
      {
        const Real t = m_cfg.period * static_cast<Real>(i) / samples;
        maxDp = std::max<Real>(maxDp, std::abs(m_inletWave(t) - m_outletWave(t)));
      }
      m_velocityScale = std::sqrt(2.0 * maxDp / m_cfg.rho);

      if (isRoot())
      {
        Alert::Info() << "[scale] max|dp| = " << maxDp
                      << " Pa  ->  sqrt(2 dp/rho) = " << m_velocityScale
                      << " m/s; mean inlet velocity at that head ~ "
                      << (m_velocityScale * m_outletMeasure / m_inletMeasure)
                      << " m/s. Divergence guard at " << m_cfg.maxVelocity
                      << " m/s = " << (m_cfg.maxVelocity / m_velocityScale)
                      << "x the scale." << Alert::Raise;

        if (m_cfg.maxVelocity > 5.0 * m_velocityScale)
          Alert::Warning() << "[scale] the guard sits more than five times "
                              "above the velocity the forcing can account for, "
                              "so a blow-up will run a long way before it trips."
                           << Alert::Raise;
      }
    }

    if (isRoot())
      Alert::Info() << "[boundary] MV length = " << m_outletMeasure
                    << " m  PV length = " << m_inletMeasure
                    << " m  |  p_pv(0) = " << m_pIn << " Pa  p_mv(0) = " << m_pOut
                    << " Pa" << Alert::Raise;
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void LeftAtrium2DViscousLogImplicit::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;
    const auto& ob = m_cfg.oldroydB;
    const Real etaS = ob.etaS;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;                  // sigma = s (exp(psi) - I)
    /// The psi equation is 2x the sigma equation of the exp form (header).
    const Real chiScale = 2.0;
    const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = Div(m_uOld) * Dot(m_u, m_v);

    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);

    // Time-dependent scalars enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };

    const auto streamV = Mult(Jacobian(m_v), m_uOld);  // (grad v) u^n
    const auto X = tensor(m_chi);
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(2, 2));

    // Everything pointwise is at (psi^k, u^k), the Newton linearisation point.
    const auto at = [this, model](std::uint8_t need, auto pick) {
      return fromPoint(m_psiIt, m_uIt, model, m_psiRevision, need, pick); };
    const auto Dexp = [this, model](const auto& w) {
      return dexp(m_psiIt, m_uIt, model, m_psiRevision, w); };
    const auto JN = [this, model](const auto& w) {
      return jacobianN(m_psiIt, m_uIt, model, m_psiRevision, w); };
    const auto Gu = [this, model](const auto& v) {
      return velocityMap(m_psiIt, m_uIt, model, m_psiRevision, v); };
    const auto divExp = [&](const auto& psi) {   // div exp(psi), chain rule
      return Mult(Dexp(Mult(Jacobian(psi), unit(0))), unit(0)) +
        Mult(Dexp(Mult(Jacobian(psi), unit(1))), unit(1));
    };

    // Momentum: exp(psi) ~ exp(psi^k) + Dexp[psi - psi^k] = T + C, all at psi^k.
    const auto T = Dexp(m_psi);
    const auto Tn = at(LogExp, [](const LogPoint& q) { return q.exp; });
    const auto C = Tn - Dexp(m_psiIt);

    // psi equation, Newton about (psi^k, u^k):
    //   (psi - psi^n)/dt + u^k.grad psi + u.grad psi^k - u^k.grad psi^k
    //   + N + JN[psi - psi^k] + Gu(u) - Gu(u^k) = 0,
    // N = -Dexp^{-1}[L^k T + T L^kT - c eps^k] + (f/lambda)(I - exp(-psi^k)),
    // Gu(v) = -Dexp^{-1}[L(v) T + T L(v)^T - c eps(v)].
    const auto N = at(LogN, [](const LogPoint& q) { return q.N; });
    const auto known = N - JN(m_psiIt) - Gu(m_uIt) - (1.0 / dt) * tensor(m_psiOld)
      - advection(m_psiIt, m_uIt);

    const auto body =
      // ---- Momentum, Eq. (19), sigma = s (exp(psi) - I) ----------------------
        (rho / dt) * Integral(m_u, m_v) - (rho / dt) * Integral(m_uOld, m_v)
      + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)
      + 2.0 * etaS * Integral(symU, symV)
      + s * Integral(T, symV) + s * Integral(C, symV) - s * Integral(I, symV)
      - Integral(m_p, Div(m_v))

      // ---- Continuity, Eq. (20) -------------------------------------------
      + Integral(Div(m_u), m_q) + m_cfg.pressurePenalty * Integral(m_p, m_q)

      // ---- Constitutive, psi form, BDF1, Newton --------------------------------
      + Integral((1.0 / dt) * tensor(m_psi) + JN(m_psi)             // psi
          + advection(m_psi, m_uIt), X)                            //   u^k . grad psi
      + Integral(advection(m_psiIt, m_u) + Gu(m_u), X)             // u
      + Integral(known, X)                                         // known

      // ---- S1, momentum (Eqs. 52, 74) -------------------------------------
      + rho * rho * Integral(m_tauK * convU, streamV)
      - rho * rho * Integral(m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get()), streamV)
      + m_cfg.pressureScale * Integral(m_alpha1 * Grad(m_p), Grad(m_q))
      - m_cfg.pressureScale * Integral(m_alpha1 * m_piGradP.get(), Grad(m_q))
      + chiScale * m_cfg.stressDivScale * s * Integral(m_alpha1 * divExp(m_psi), divergence(m_chi))
      - chiScale * m_cfg.stressDivScale * Integral(m_alpha1 * m_piDivSigma.get(), divergence(m_chi))

      // ---- S2, continuity (Eq. 53) ----------------------------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v))
      - Integral(m_alpha2 * m_piDiv.get(), Div(m_v))

      // ---- S3: u-psi compatibility, and the psi advection (header) ----------
      + Integral(m_alpha3 * symU, symV)
      - Integral(m_alpha3 * tensor(m_piEps.get()), symV)
      + Integral(m_alphaPsi * advection(m_psi, m_uOld), advection(m_chi, m_uOld))
      - Integral(m_alphaPsi * tensor(m_piAdvPsi.get()), advection(m_chi, m_uOld))

      // ---- Boundary ----------------------------------------------------------
      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)

           // Source impedance of the pressure inlets: positive semidefinite
           // and assembled implicitly. Default 0, see Config.
      + m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet)

      + m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)

      + BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // The psi equation is hyperbolic in psi: relaxed blood enters at the veins.
    if (m_cfg.inletConformation)
      m_flow = body + DirichletBC(m_psi, Zero(size_t(3))).on(m_inletSet);
    else
      m_flow = body;
  }

  LeftAtrium2DViscousLogImplicit::Real LeftAtrium2DViscousLogImplicit::crosswind(const Point& p, Real cur, Real prev,
    const Math::SpatialVector<Real>& gradient, Real diffusivity, Real reaction) const
  {
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = m_uOld.getValue(p);
    const Real residual = (cur - prev) / m_cfg.dt + Math::dot(uc, gradient) + reaction;
    return std::max<Real>(0.0,
      m_cfg.crosswindC * std::abs(residual) * cellSize(p) / (2.0 * gn) - diffusivity);
  }

  void LeftAtrium2DViscousLogImplicit::setupSpecies()
  {
    const Real Dth = m_cfg.thrombosis.diffusivityThrombin;
    const Real Dfg = m_cfg.thrombosis.diffusivityFibrinogen;
    const Real Dfn = m_cfg.thrombosis.diffusivityFibrin;
    const Real keff = m_cfg.thrombosis.reactionRate;
    const Real scale = m_cfg.thrombosis.stabilizationScale;
    const Real dt = m_cfg.dt;
    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;
    const Real thIn = m_cfg.inletThrombin;
    const Real fnIn = m_cfg.inletFibrin;

    // Cell Peclet is of order 1e6 here: these are essentially pure advection
    // problems and their tau is set by the advective and transient scales, not
    // by the flow's tau.
    RealFunction tauThFn = [this, Dth, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dth, 0.0, dt);
    };
    RealFunction tauFgFn = [this, Dfg, keff, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale *
        supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dfg,
          keff * std::abs(m_thCur.getValue(p)), dt);
    };
    RealFunction tauFnFn = [this, Dfn, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dfn, 0.0, dt);
    };

    // Codina crosswind: SUPG leaves the crosswind direction free, which is what
    // produces negative concentration bands across shear layers.
    RealFunction invSpeedSqFn = [this](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      const Real u2 = Math::dot(uc, uc);
      return (u2 > 1.0e-20) ? 1.0 / u2 : 0.0;
    };
    RealFunction kdcThFn = [this, Dth](const Point& p) {
      return crosswind(p, m_thCur.getValue(p), m_thPrev.getValue(p),
        Grad(m_thCur).getValue(p), Dth, 0.0);
    };
    RealFunction kdcFgFn = [this, Dfg, keff](const Point& p) {
      return crosswind(p, m_fgCur.getValue(p), m_fgPrev.getValue(p),
        Grad(m_fgCur).getValue(p), Dfg, keff * m_thCur.getValue(p) * m_fgCur.getValue(p));
    };
    RealFunction kdcFnFn = [this, Dfn, keff](const Point& p) {
      return crosswind(p, m_fnCur.getValue(p), m_fnPrev.getValue(p),
        Grad(m_fnCur).getValue(p), Dfn,
        -keff * m_thCur.getValue(p) * m_fgCur.getValue(p));
    };

    const auto pTh = Dot(m_uOld, Grad(m_vth));
    const auto pFg = Dot(m_uOld, Grad(m_vfg));
    const auto pFn = Dot(m_uOld, Grad(m_vfn));

    m_species = (1.0 / dt) * Integral(m_th, m_vth) -
      (1.0 / dt) * Integral(m_thCur, m_vth) + Dth * Integral(Grad(m_th), Grad(m_vth)) +
      Integral(Dot(m_uOld, Grad(m_th)), m_vth)

      + (1.0 / dt) * Integral(m_fg, m_vfg) - (1.0 / dt) * Integral(m_fgCur, m_vfg) +
      Dfg * Integral(Grad(m_fg), Grad(m_vfg)) + Integral(Dot(m_uOld, Grad(m_fg)), m_vfg) +
      keff * Integral(m_thCur * m_fg, m_vfg)

      + (1.0 / dt) * Integral(m_fn, m_vfn) - (1.0 / dt) * Integral(m_fnCur, m_vfn) +
      Dfn * Integral(Grad(m_fn), Grad(m_vfn)) + Integral(Dot(m_uOld, Grad(m_fn)), m_vfn) -
      keff * Integral(m_thCur * m_fg, m_vfn)

      // Endothelial thrombin flux: a surface flux, armed by the cycle indices.
      // The activation field is identically zero until the first cycle has
      // closed, so this term contributes nothing before the OSI exists.
      - m_cfg.thrombosis.thrombinWallFlux *
        BoundaryIntegral(m_activation * m_vth).over(m_wallSet)

      // SUPG. Only fibrinogen has a sink proportional to its own unknown.
      + (1.0 / dt) * Integral(tauThFn * m_th, pTh) -
      (1.0 / dt) * Integral(tauThFn * m_thCur, pTh) +
      Integral(tauThFn * Dot(m_uOld, Grad(m_th)), pTh)

      + (1.0 / dt) * Integral(tauFgFn * m_fg, pFg) -
      (1.0 / dt) * Integral(tauFgFn * m_fgCur, pFg) +
      Integral(tauFgFn * Dot(m_uOld, Grad(m_fg)), pFg) +
      keff * Integral(tauFgFn * m_thCur * m_fg, pFg)

      + (1.0 / dt) * Integral(tauFnFn * m_fn, pFn) -
      (1.0 / dt) * Integral(tauFnFn * m_fnCur, pFn) +
      Integral(tauFnFn * Dot(m_uOld, Grad(m_fn)), pFn) -
      keff * Integral(tauFnFn * m_thCur * m_fg, pFn)

      // Codina crosswind, one pair per species.
      + Integral(kdcThFn * Grad(m_th), Grad(m_vth)) -
      Integral(kdcThFn * invSpeedSqFn * Dot(m_uOld, Grad(m_th)), Dot(m_uOld, Grad(m_vth)))

      + Integral(kdcFgFn * Grad(m_fg), Grad(m_vfg)) -
      Integral(kdcFgFn * invSpeedSqFn * Dot(m_uOld, Grad(m_fg)), Dot(m_uOld, Grad(m_vfg)))

      + Integral(kdcFnFn * Grad(m_fn), Grad(m_vfn)) -
      Integral(kdcFnFn * invSpeedSqFn * Dot(m_uOld, Grad(m_fn)), Dot(m_uOld, Grad(m_vfn)))

      // Pure advection needs a Dirichlet condition on the inflow boundary. The
      // veins carry plasma fibrinogen, the circulating thrombin level, and no
      // fibrin. Imposing it on the whole PV patch rather than only where
      // u.n < 0 is what makes this a *boundary* condition and not a switch:
      // the mitral outlet is left free, so the reverse-flow terms there are
      // the ones that have to hold the transport together during backflow.
      + DirichletBC(m_th, RealFunction(Real(thIn))).on(m_inletSet) +
      DirichletBC(m_fn, RealFunction(Real(fnIn))).on(m_inletSet) +
      DirichletBC(m_fg, RealFunction(Real(fg0))).on(m_inletSet);
  }

  void LeftAtrium2DViscousLogImplicit::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const Real etaS = m_cfg.oldroydB.etaS;
    const auto& sigma = m_sigma;
    const auto nx = Component(normal, 0);
    const auto ny = Component(normal, 1);

    // t = (2 eta_s eps(u) + sigma) n: the recovered rows of 2 eps(u) plus the
    // elastic stress, Voigt (xx, xy, yy), which is nodal already.
    const auto traction = VectorFunction(
      etaS * Dot(m_symRec0, normal) + Component(sigma, 0) * nx + Component(sigma, 1) * ny,
      etaS * Dot(m_symRec1, normal) + Component(sigma, 1) * nx + Component(sigma, 2) * ny);
    const auto wallStress = traction - Dot(traction, normal) * normal;

    // An L2 projection restricted to the wall, regularised in the interior so
    // the mass matrix stays invertible off it. It is used here, and not nodal
    // interpolation, for one reason: a wall node belongs to two facets with
    // different normals, so tau_w has two nodal values and the projection is
    // what averages them by facet measure. Everything downstream of this --
    // |tau_w|, TAWSS, OSI, the activation weight -- is nodal.
    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  void LeftAtrium2DViscousLogImplicit::updateStabilization()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;
    const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
    const Real gradDiv = m_cfg.useVMS ? m_cfg.gradDivScale : 0.0;

    // Cellwise parameters at u^n, Eqs. (40)-(42); tau_K is the dynamic alpha1.
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

    // Pi_h of every split term, from (u^n, p^n, sigma^n) (Sec. 7.1).
    m_piConv.project(Mult(Jacobian(m_uOld), m_uOld));
    m_sub.project(VectorFunction(dim,
      [this, rho, dt, dim](const Point& p) -> Math::SpatialVector<Real> {
        // u' = tau_K rho (u'^n/dt - ((grad u^n) u^n - Pi))
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
    // sigma^n = s (exp(psi^n) - I), with the operators of setupFlow(), all at
    // (psi^n, u^n).
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

  bool LeftAtrium2DViscousLogImplicit::solveFlow()
  {
    const bool trace = isRoot() && m_step < 3;
    const auto phase = [trace](const char* what) {
      if (trace)
      {
        ThreeDInfo() << what << " ..." << Alert::Raise;
        std::cout.flush();
      }
    };

    const auto vmsStart = CoronaryClock::now();

    phase("VMS: updating the stabilisation (alphas and projections)");
    updateStabilization();

    m_timing.vms = secondsSince(vmsStart);

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Assembling the flow system ("
                   << (m_vh.getSize() + m_sh.getSize() + m_tauh.getSize())
                   << " unknowns) ..." << Alert::Raise;
      std::cout.flush();
    }

    // Newton about (u^k, psi^k), from (u^n, psi^n); every iteration
    // re-assembles everything and solves once.
    ::KSPConvergedReason reason = KSP_CONVERGED_ITS;
    PetscErrorCode ierr = PETSC_SUCCESS;
    ::Vec increment = PETSC_NULLPTR;
    ierr = VecDuplicate(m_psiIt.getData(), &increment);
    assert(ierr == PETSC_SUCCESS);
    ::Vec uIncrement = PETSC_NULLPTR;
    ierr = VecDuplicate(m_uIt.getData(), &uIncrement);
    assert(ierr == PETSC_SUCCESS);
    m_timing.assembly = 0.0;
    m_timing.solve = 0.0;
    m_conformationIts = 0;
    m_uIt.setData(m_uOld.getData());
    for (int k = 0; k < std::max(1, m_cfg.conformationIterations); ++k)
    {
      const auto assemblyStart = CoronaryClock::now();
      m_flow.assemble();
      m_timing.assembly += secondsSince(assemblyStart);

      if (!m_flowFieldSplitsSet)
      {
        if (isRoot())
          ThreeDInfo() << "Creating field splits ..." << Alert::Raise;
        m_flow.setFieldSplits();
        m_flowFieldSplitsSet = true;
      }

      const auto solveStart = CoronaryClock::now();
      m_flow.solve(m_flowKSP);
      m_timing.solve += secondsSince(solveStart);

      ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
      assert(ierr == PETSC_SUCCESS);
      if (reason <= 0)
        break;

      // psi^{k+1} = psi^k + omega dpsi, omega < 1 only to cap max|dpsi| at
      // newtonMaxStep; the increment is reported relative to max(1, |psi^k|).
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

      if (isRoot() && m_step < 3)
        KSPInfo() << "Newton " << m_conformationIts << ": max|dpsi| = " << stepSize
                  << (omega < 1.0 ? " (damped)" : "") << Alert::Raise;
      if (!std::isfinite(m_psiIncrement) ||
          (omega == 1.0 && m_psiIncrement < m_cfg.conformationTolerance))
        break;
    }
    ierr = VecDestroy(&increment);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDestroy(&uIncrement);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    // Only a real stall is worth a line.
    if (isRoot() && m_conformationIts >= m_cfg.conformationIterations &&
        m_psiIncrement >= 100.0 * m_cfg.conformationTolerance)
      Alert::Warning() << "[log] Newton stopped after " << m_conformationIts
                       << " iterations with max|dpsi| = " << m_psiIncrement << Alert::Raise;

    // The accepted state is the (damped) iterate pair, so u and psi stay
    // consistent with each other.
    m_uOld.setData(m_uIt.getData());
    m_u.getSolution().setData(m_uIt.getData());
    m_psiOld.setData(m_psiIt.getData());
    m_psi.getSolution().setData(m_psiIt.getData());
    ++m_psiRevision;
    updateConformation();
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.get().getData());

    // A direct solve always "succeeds": check the physics, not the solver.
    m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
    m_stress = std::max(std::abs(m_sigma.max()), std::abs(m_sigma.min()));

    return reason > 0 && std::isfinite(m_speed) && m_speed <= m_cfg.maxVelocity;
  }

  void LeftAtrium2DViscousLogImplicit::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Recovering the wall shear stress ..." << Alert::Raise;
      std::cout.flush();
    }

    // The two rows of 2 eps(u) = grad u + grad u^T, recovered onto the nodes:
    //   row 0 = (2 du_x/dx,            du_x/dy + du_y/dx)
    //   row 1 = (du_x/dy + du_y/dx,    2 du_y/dy)
    // grad u_h is elementwise constant on P1, so a recovery is unavoidable and
    // an L2 projection is the right one here -- it is a linear functional of
    // the solution, unlike the indices built from it further down.
    const auto jac = Jacobian(uSol);
    const auto offDiagonal = Component(jac, 0, 1) + Component(jac, 1, 0);

    projectVector(VectorFunction(2.0 * Component(jac, 0, 0), offDiagonal), m_symRec0);
    projectVector(VectorFunction(offDiagonal, 2.0 * Component(jac, 1, 1)), m_symRec1);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

    // Keep only the wall. tau_w is an L2 projection, not a nodal quantity, and
    // the operator is (M_wall + reg M_vol). ON the wall the boundary mass
    // dominates by reg*h ~ 5e-7, so the recovered traction there is exact --
    // that part is sound. OFF the wall the only equation a node has is
    // reg M x = 0, and for a consistent P1 mass matrix that forces
    // x_i = -(1/M_ii) sum_j M_ij x_j with every M_ij > 0: the tail ALTERNATES
    // IN SIGN from node to node, decaying by about a quarter per layer. A
    // sign-alternating field is exactly what renders as dots rather than as a
    // field, and it is what put a 0/0 inside OSI. It is not a solver failure:
    // after Jacobi the operator has a condition number of about 3.
    //
    // Zeroing it leaves the wall untouched, makes the two cycle accumulators
    // identically zero off the wall, and costs one nodal pass.
    const size_t dim = m_mesh.getSpaceDimension();
    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    const auto& wssSol = m_wssTrial.getSolution();
    m_wss = Math::SpatialVector<Real>{{0.0, 0.0}};
    m_wss.project(Region::Boundary,
      VectorFunction(dim,
        [&wssSol, dim](const Point& p) -> Math::SpatialVector<Real> {
          const auto w = wssSol.getValue(p);
          Math::SpatialVector<Real> out(dim);
          for (Index c = 0; c < static_cast<Index>(dim); ++c)
            out(c) = w(c);
          return out;
        }),
      onWall);

    // Cycle accumulators. Two are needed and they are not interchangeable: the
    // vector integral measures how much net direction survives, the scalar one
    // how much shear was applied regardless of direction. Both are formed by
    // VecAXPY, i.e. degree of freedom by degree of freedom.
    axpy(m_cfg.dt, m_wss.getData(), m_netShear.getData());

    // Nodal interpolation, NOT an L2 projection. |tau_w| is finite on the wall
    // and zero in the interior; the L2 projection of that onto a continuous P1
    // field undershoots, and a single negative node makes int|tau_w| dt -- and
    // therefore TAWSS, which is non-negative by definition -- come out
    // negative. Interpolating the magnitude at the nodes cannot.
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_cfg.dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void LeftAtrium2DViscousLogImplicit::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // The three indices live ON THE WALL and nowhere else.
    //
    // They used to be projected over the whole domain, which is not merely
    // untidy: off the wall tau_w decays to zero, so TAWSS -> 0, the logistic
    // SATURATES at its maximum 1/(1+exp(-tau_a/w)) = 0.935, and OSI becomes
    // |int tau dt| / int |tau| dt with both accumulators vanishing -- a 0/0
    // that drifts to 1/2. The product is 0.935 * 1 = 0.933, and that is
    // exactly the number the run reported as maxActivation, with maxOSI =
    // 0.499 beside it. The field maximum was sitting in the middle of the
    // cavity, where there is no endothelium at all, and the picture was of
    // the interior rather than of the wall.
    //
    // The wall flux itself was never affected -- BoundaryIntegral(...) over
    // the wall only ever reads wall nodes, whose TAWSS and OSI are genuine --
    // so this changes the diagnostics and the XDMF fields, not the physics of
    // any run already completed. What it does change is the reported maxima,
    // which were interior artefacts and are now wall values.
    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    // Zeroed first, so a node that is not on the wall carries 0 rather than
    // whatever the previous cycle left there.
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);

    // TAWSS = (1/T) int |tau_w| dt. Still an exact scaling of the accumulator,
    // node by node: evaluating a P1 field at its own node returns the nodal
    // value, so nothing here can change its sign.
    m_tawss.project(Region::Boundary,
      RealFunction([this, elapsed](const Point& p) -> Real {
        return m_absShear.getValue(p) / elapsed;
      }),
      onWall);

    // OSI = (1/2)[1 - |int tau_w dt| / int |tau_w| dt]. Both accumulators are
    // built from the same nodal values, so the triangle inequality holds node
    // by node and OSI lands in [0, 1/2] on its own; the clamp is insurance.
    // Evaluated at the nodes: a ratio of two fields taken through quadrature
    // is not the ratio of the two nodal fields, and it is not bounded by 1/2
    // either.
    m_osi.project(Region::Boundary, RealFunction([this](const Point& p) -> Real {
      const Real abs = m_absShear.getValue(p);
      if (abs <= 0.0)
        return 0.0;
      const auto net = m_netShear.getValue(p);
      const Real mag = std::sqrt(Math::dot(net, net));
      return std::clamp<Real>(0.5 * (1.0 - mag / abs), 0.0, 0.5);
    }),
      onWall);

    // Smooth, bounded activation: a logistic in the measured shear threshold
    // times the oscillatory index mapped onto [0,1]. Clamped to [0,1] so the
    // endothelial thrombin flux can never turn into a sink -- a negative
    // activation is what drives thrombin, and with it fibrin, negative.
    // Reads the two fields written just above, so it must come last.
    const Real tauA = m_cfg.thrombosis.activationShearStress;
    const Real width = std::max<Real>(m_cfg.thrombosis.activationShearWidth, 1e-12);
    m_activation.project(Region::Boundary,
      RealFunction([this, tauA, width](const Point& p) -> Real {
        const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
        const Real osi = std::clamp<Real>(2.0 * m_osi.getValue(p), 0.0, 1.0);
        return std::clamp<Real>(low * osi, 0.0, 1.0);
      }),
      onWall);

    // The ghost entries must be zeroed too, or the next accumulation reads a
    // stale halo.
    const auto zeroWithGhosts = [](::Vec& v) {
      PetscErrorCode e = VecZeroEntries(v);
      assert(e == PETSC_SUCCESS);
      e = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(e == PETSC_SUCCESS);
      e = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(e == PETSC_SUCCESS);
      (void)e;
    };

    zeroWithGhosts(m_netShear.getData());
    zeroWithGhosts(m_absShear.getData());

    m_indicesReady = true;
  }

  void LeftAtrium2DViscousLogImplicit::solveSpecies()
  {
    const auto speciesStart = CoronaryClock::now();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Assembling and solving the coagulation kinetics ..."
                   << Alert::Raise;
      std::cout.flush();
    }

    m_species.assemble();
    m_species.solve(m_speciesKSP);

    // History rotates before the current level is overwritten: the crosswind
    // residual needs both c^n and c^{n-1}.
    m_thPrev.setData(m_thCur.getData());
    m_fgPrev.setData(m_fgCur.getData());
    m_fnPrev.setData(m_fnCur.getData());
    m_thCur.setData(m_th.getSolution().getData());
    m_fgCur.setData(m_fg.getSolution().getData());
    m_fnCur.setData(m_fn.getSolution().getData());

    m_timing.species = secondsSince(speciesStart);
  }

  void LeftAtrium2DViscousLogImplicit::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    // In 2D these are flow rates per unit depth, m^2/s. n is outward, so qIn
    // is negative while the veins fill the atrium.
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_qIn = m_flux(m_one);

    for (size_t i = 0; i < m_cfg.labels.inlets.size(); ++i)
    {
      m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.inlets[i]);
      m_flux.assemble();
      m_qInPatch[i] = m_flux(m_one);
    }

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);

    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_outletPressure = m_flux(m_one) / m_outletMeasure;
  }

  void LeftAtrium2DViscousLogImplicit::writeCSVHeader()
  {
    m_csv << "t,cycle,pPV,pMV,dp,qIn,qOut,";
    for (const auto tag : m_cfg.labels.inlets)
      m_csv << "qIn" << tag << ',';
    m_csv << "pOutletMean,maxU,uScale,maxSigma,"
          << "maxTAWSS,maxOSI,maxActivation,maxThrombin,minFibrinogen,maxFibrin\n";
  }

  void LeftAtrium2DViscousLogImplicit::writeCSVRow(int cycle)
  {
    // Every one of these reduces over the communicator, so they are taken on
    // all ranks before the root-only write.
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();
    const Real activation = m_activation.max();
    const Real th = m_thCur.max();
    const Real fg = m_fgCur.min();
    const Real fn = m_fnCur.max();

    if (!isRoot())
      return;

    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << m_pOut << ','
          << (m_pIn - m_pOut) << ',' << m_qIn << ',' << m_qOut << ',';
    for (const auto q : m_qInPatch)
      m_csv << q << ',';
    m_csv << m_outletPressure << ',' << m_speed << ',' << m_velocityScale << ',' << m_stress
          << ',' << tawss
          << ',' << osi << ',' << activation << ',' << th << ',' << fg << ',' << fn
          << '\n';
    m_csv.flush();
  }

  int LeftAtrium2DViscousLogImplicit::run()
  {
    if (!m_initialized)
      initialize();

    const int stepsPerCycle = static_cast<int>(m_cfg.period / m_cfg.dt + 0.5);
    const int totalCycles = m_cfg.flowCycles + m_cfg.speciesCycles;
    const int totalSteps = totalCycles * stepsPerCycle;

    if (isRoot())
      Alert::Info() << "[run] " << totalCycles << " cycles of " << stepsPerCycle
                    << " steps (" << totalSteps << " total); species from cycle "
                    << (m_cfg.flowCycles + 1)
                    << " on, and never before the first OSI; XDMF every "
                    << m_cfg.outputEvery << " steps and at every cycle boundary"
                    << Alert::Raise;

    Real cycleElapsed = 0.0;
    const auto runStart = CoronaryClock::now();

    for (int step = 0; step < totalSteps; ++step)
    {
      m_step = step;
      m_timing = Timing{};
      const auto stepStart = CoronaryClock::now();

      m_t += m_cfg.dt;
      const int cycle = step / stepsPerCycle;
      const bool endOfCycle = (step % stepsPerCycle == stepsPerCycle - 1);

      // Both tractions come from the tabulated waveforms. p_pv - p_mv is
      // identically zero while the mitral valve is shut, so the pair carries
      // the valve and no diode is imposed on top of it.
      m_pIn = m_inletWave(m_t) + m_cfg.pressureOffset;
      m_pOut = m_outletWave(m_t) + m_cfg.pressureOffset;

      if (!solveFlow())
      {
        Alert::Exception() << "[flow] diverged at step " << (step + 1)
                           << ": max|u| = " << m_speed
                           << " m/s = " << (m_speed / m_velocityScale)
                           << "x sqrt(2 max|dp|/rho), i.e. a dynamic head of "
                           << (0.5 * m_cfg.rho * m_speed * m_speed)
                           << " Pa against a driving head of at most "
                           << (0.5 * m_cfg.rho * m_velocityScale * m_velocityScale)
                           << " Pa" << Alert::Raise;
        return 1;
      }

      computeWallShear();
      cycleElapsed += m_cfg.dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      m_timing.fluxes = secondsSince(fluxStart);

      // The indices are closed BEFORE the species are advanced, so that within
      // this step the activation field the wall flux reads is the one that
      // belongs to the cycle just finished.
      if (endOfCycle)
      {
        closeCycle(cycleElapsed);
        cycleElapsed = 0.0;

        // VecMax is collective: every rank must reach it, so the reductions
        // happen outside the root-only print.
        const Real maxTawss = m_tawss.max();
        const Real maxOsi = m_osi.max();
        const Real maxActivation = m_activation.max();

        if (isRoot())
          Alert::Info() << "[cycle " << (cycle + 1) << "/" << totalCycles
                        << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                        << "  maxActivation=" << maxActivation << Alert::Raise;
      }

      // Kinetics: only once the flow is periodic AND the first OSI exists.
      if (m_cfg.solveKinetics && m_indicesReady && cycle >= m_cfg.flowCycles)
        solveSpecies();

      writeCSVRow(cycle);

      // The flow is written from the very first step, so the warm-up cycles are
      // available before the species are switched on.
      const bool output =
        endOfCycle || (m_cfg.outputEvery > 0 && step % m_cfg.outputEvery == 0);
      if (output)
      {
        const auto outputStart = CoronaryClock::now();
        m_xdmf.write(m_t).flush();
        m_timing.output = secondsSince(outputStart);
      }

      m_timing.total = secondsSince(stepStart);

      if (isRoot() && (m_step < 3 || step % 20 == 0 || endOfCycle))
      {
        const Real elapsed = secondsSince(runStart);
        const Real perStep = elapsed / static_cast<Real>(step + 1);

        Alert::Info() << "---- Step " << (step + 1) << "/" << totalSteps << "  cycle "
                      << (cycle + 1) << "/" << totalCycles << "  t = " << m_t << " s"
                      << (cycle < m_cfg.flowCycles ? "  (warm-up)" : "") << Alert::Raise;

        Alert::Info() << "[2D] max|u|=" << m_speed << " m/s ("
                      << (m_speed / m_velocityScale) << "x scale)  max|sigma|=" << m_stress
                      << " Pa  newton=" << m_conformationIts << " (|dpsi|=" << m_psiIncrement
                      << ")  p_pv=" << m_pIn
                      << " Pa  p_mv=" << m_pOut << " Pa  dp=" << (m_pIn - m_pOut)
                      << " Pa  qIn=" << m_qIn << "  qOut=" << m_qOut << " m^2/s"
                      << Alert::Raise;

        Alert::Info() << "[PV] per ostium (m^2/s): " << m_cfg.labels.inlets[0] << "="
                      << m_qInPatch[0] << "  " << m_cfg.labels.inlets[1] << "="
                      << m_qInPatch[1] << "  " << m_cfg.labels.inlets[2] << "="
                      << m_qInPatch[2] << "  " << m_cfg.labels.inlets[3] << "="
                      << m_qInPatch[3] << Alert::Raise;

        Alert::Info() << "[timing] vms=" << m_timing.vms << "  asm=" << m_timing.assembly
                      << "  ksp=" << m_timing.solve << "  wss=" << m_timing.shear
                      << "  flux=" << m_timing.fluxes << "  species=" << m_timing.species
                      << "  out=" << m_timing.output << "  total=" << m_timing.total
                      << " s  |  ETA " << (perStep * (totalSteps - step - 1) / 60.0)
                      << " min" << Alert::Raise;
        std::cout.flush();
      }
    }

    m_xdmf.close();
    if (isRoot())
      m_csv.close();

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

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;

    {
      Rodin::Examples::Heart::LeftAtrium2DViscousLogImplicit::Config cfg;

      char buffer[512];
      PetscBool got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mesh", buffer, sizeof(buffer), &got);
      if (got)
        cfg.meshPath = buffer;

      got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_pv", buffer, sizeof(buffer), &got);
      if (got)
        cfg.inletPressurePath = buffer;

      got = PETSC_FALSE;
      PetscOptionsGetString(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mv", buffer, sizeof(buffer), &got);
      if (got)
        cfg.outletPressurePath = buffer;

      PetscReal real = 0.0;
      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mesh_scale", &real, &got);
      if (got)
        cfg.meshScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_dt", &real, &got);
      if (got)
        cfg.dt = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_period", &real, &got);
      if (got)
        cfg.period = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_th_in", &real, &got);
      if (got)
        cfg.inletThrombin = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_inlet_impedance", &real, &got);
      if (got)
        cfg.inletImpedance = real;

      PetscInt integer = 0;
      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_flow_cycles", &integer, &got);
      if (got)
        cfg.flowCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_species_cycles", &integer, &got);
      if (got)
        cfg.speciesCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_output_every", &integer, &got);
      if (got)
        cfg.outputEvery = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_vms_scale", &real, &got);
      if (got)
        cfg.vmsScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_graddiv_scale", &real, &got);
      if (got)
        cfg.gradDivScale = real;

      const auto getReal = [&](const char* name, Rodin::Real& out) {
        PetscBool set = PETSC_FALSE;
        PetscReal value = 0.0;
        PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &set);
        if (set)
          out = value;
      };
      getReal("-la2d_pressure_scale", cfg.pressureScale);
      getReal("-la2d_stress_div_scale", cfg.stressDivScale);
      getReal("-la2d_stress_scale", cfg.stressScale);
      getReal("-la2d_eta_s", cfg.oldroydB.etaS);
      getReal("-la2d_eta_p", cfg.oldroydB.etaP);
      getReal("-la2d_lambda", cfg.oldroydB.lambda);
      getReal("-la2d_lambda0_factor", cfg.oldroydB.lambda0Factor);
      getReal("-la2d_lambda0_min", cfg.oldroydB.lambda0Min);
      getReal("-la2d_ptt_epsilon", cfg.oldroydB.pttEpsilon);
      getReal("-la2d_conformation_tol", cfg.conformationTolerance);
      getReal("-la2d_newton_max_step", cfg.newtonMaxStep);
      {
        PetscInt its = cfg.conformationIterations;
        PetscBool set = PETSC_FALSE;
        PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_conformation_its", &its, &set);
        if (set)
          cfg.conformationIterations = static_cast<int>(its);
      }

      PetscBool flag = PETSC_FALSE;
      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_vms", &flag, &got);
      if (got)
        cfg.useVMS = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_kinetics", &flag, &got);
      if (got)
        cfg.solveKinetics = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(
        PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_inlet_conformation", &flag, &got);
      if (got)
        cfg.inletConformation = (flag == PETSC_TRUE);

      Rodin::Examples::Heart::LeftAtrium2DViscousLogImplicit simulation(context, cfg);
      status = simulation.initialize().run();
    }

    PetscFinalize();
    return status;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }
}
