/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ArterialLesionAxi_viscous_logarithmic_implicit.cpp
 * @brief Driver for pulsatile sPTT/Oldroyd-B flow through an idealised
 *        stenosis or aneurysm of a pipe, axisymmetric, fully-implicit
 *        log-conformation; see the header.
 *
 * Run (from the build directory):
 *   mpirun -n 5 ./examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit \
 *     -al_mesh ../resources/examples/viscoelastic_fluids/S75_axi_medium.mesh \
 *     -al_re 300 -al_wo 4 -al_amplitude 0.5 -al_wi 1
 *
 * Options (all prefixed -al_): those of ArterialLesion2D, plus axis_penalty,
 * fluid (sptt | cy | gsptt), eta_ref, cy_eta_zero, cy_eta_inf, cy_lambda,
 * cy_n, cy_a, inlet_developed_conformation, trace_newton (see the header).
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

#include "ArterialLesionAxi_viscous_logarithmic_implicit.h"
#include "CoronaryArtery/CoronaryArteryAlerts.h"
#include "CoronaryArtery/CoronaryArteryTiming.h"

namespace Rodin::Examples::ViscoelasticFluids
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  using Heart::CoronaryClock;
  using Heart::KSPInfo;
  using Heart::ThreeDInfo;
  using Heart::secondsSince;

  // The in-plane operators below are those of ArterialLesion2D (and of
  // Heart/LeftAtrium2D_viscous_logarithmic_implicit); the hoop block psi_t is
  // the only addition.
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

    // ---- in-plane tensors in Voigt order (xx, xr, rr) -------------------------

    /// @brief E0 = e_x e_x, E1 = e_x e_r + e_r e_x, E2 = e_r e_r.
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

    /// @brief e_r e_r: selects u_r in the axis penalty.
    MatrixFunction<Math::Matrix<Real>> radialProjector()
    {
      Math::Matrix<Real> P = Math::Matrix<Real>::Zero(2, 2);
      P(1, 1) = 1;
      return MatrixFunction(P);
    }

    template <class S>
    auto tensor(const S& s)
    {
      return voigtBasis(0) * Component(s, 0) + voigtBasis(1) * Component(s, 1) +
        voigtBasis(2) * Component(s, 2);
    }

    template <class M>
    auto voigt(const M& m)
    {
      return VectorFunction{ Component(m, 0, 0), Component(m, 0, 1), Component(m, 1, 1) };
    }

    /// @brief (u . grad) sigma, componentwise: no swirl, no Christoffel terms.
    template <class S, class U>
    auto advection(const S& s, const U& u)
    {
      return tensor(Mult(Jacobian(s), u));
    }

    /// @brief In-plane (div sigma)_i = d sigma_ij / d x_j (stabilisation only).
    template <class S>
    auto divergence(const S& s)
    {
      const auto dx = Mult(Jacobian(s), unit(0));
      const auto dy = Mult(Jacobian(s), unit(1));
      return unit(0) * (Component(dx, 0) + Component(dy, 1)) +
        unit(1) * (Component(dx, 1) + Component(dy, 2));
    }

    // ---- log-conformation ---------------------------------------------------

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

    struct LogModel
    {
        Real lambda = 0.0;
        Real lambda0 = 0.0;
        Real epsilon = 0.0;   ///< sPTT
        Real c = 0.0;         ///< 2 (1 - lambda0/lambda), on eps(u)
    };

    /// @brief f = 1 + epsilon (lambda/lambda0)(tr exp(psi_P) + e^{psi_t} - 3).
    Real pttFactorOf(Real traceP, Real psiT, const LogModel& m)
    {
      return 1.0 + m.epsilon * (m.lambda / m.lambda0) * (traceP + std::exp(psiT) - 3.0);
    }

    /**
     * @brief Everything the forms need at one quadrature point, from
     *        (psi_P^k, psi_t^k, L_P^k = grad u^k).
     *
     * @details In-plane block, as in the planar driver: exp, Dexp[B_j], N_P =
     *          -Dexp^{-1}[L T + T L^T - c eps] + (f/lambda)(I - exp(-psi_P)),
     *          JN_j = dN_P/dpsi_j, Q_ab = dN_P/dL_ab; plus JNT = dN_P/dpsi_t
     *          (through f). Hoop block: N_t = a(psi_t) L_t + g(psi_P, psi_t),
     *          L_t = u_r/r, with a, da = a', g and its derivatives gT, gP_j.
     *          All derivatives in psi by central differences, as in 2D.
     */
    struct LogPoint
    {
        Eigen::Matrix2d exp;
        std::array<Eigen::Matrix2d, 3> D;
        Real f = 1.0;
        Eigen::Matrix2d N;
        std::array<Eigen::Matrix2d, 3> JN;
        Eigen::Matrix2d JNT;
        std::array<std::array<Eigen::Matrix2d, 2>, 2> Q;
        Real psiT = 0.0;
        Real expT = 1.0;
        Real a = 0.0;
        Real da = 0.0;
        Real g = 0.0;
        Real gT = 0.0;
        std::array<Real, 3> gP{{0.0, 0.0, 0.0}};
        /// @brief V_m = sum_j D2exp_{psi^k}[d_j psi^k, B_m] e_j, the derivative
        ///        in psi of div exp(psi) = sum_j Dexp_psi[d_j psi] e_j beyond
        ///        the one carried by d_j psi, and Vk = sum_m V_m psi^k_m.
        std::array<Eigen::Vector2d, 3> V;
        Eigen::Vector2d Vk;
    };

    Eigen::Matrix2d nonlinearMap(const Eigen::Matrix2d& psi, Real psiT,
      const Eigen::Matrix2d& L, const LogModel& m, Real* fOut = nullptr)
    {
      const LogEig e = logEig(psi);
      const Eigen::Matrix2d T = e.exp();
      const Eigen::Matrix2d H = L * T + T * L.transpose() - 0.5 * m.c * (L + L.transpose());
      const Real f = pttFactorOf(T.trace(), psiT, m);
      if (fOut)
        *fOut = f;
      return -e.dexpInverse(H) +
        (f / m.lambda) * (Eigen::Matrix2d::Identity() - e.expMinus());
    }

    /// @brief g = (f/lambda)(1 - e^{-psi_t}), the relaxation of the hoop block.
    Real hoopRelaxation(const Eigen::Matrix2d& psi, Real psiT, const LogModel& m)
    {
      const Real f = pttFactorOf(logEig(psi).exp().trace(), psiT, m);
      return (f / m.lambda) * (1.0 - std::exp(-psiT));
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

    /// @brief A scalar field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseScalar final
      : public RealFunctionBase<PointwiseScalar<F>>
    {
      public:
        using Parent = RealFunctionBase<PointwiseScalar<F>>;

        explicit PointwiseScalar(F f, size_t order = 2) : m_f(std::move(f)), m_order(order) {}
        PointwiseScalar(const PointwiseScalar& other) : Parent(other), m_f(other.m_f), m_order(other.m_order) {}
        PointwiseScalar(PointwiseScalar&& other)
          : Parent(std::move(other)), m_f(std::move(other.m_f)), m_order(other.m_order) {}

        Real getValue(const Point& p) const { return m_f(p); }
        /// Not a polynomial; quadratic is enough against P1 x P1. Order 0 is
        /// used for a function known to be constant, so the quadrature of the
        /// parent formulation is kept.
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept { return m_order; }
        PointwiseScalar* copy() const noexcept override { return new PointwiseScalar(*this); }

      private:
        F m_f;
        size_t m_order;
    };

    /// @brief A 2-vector field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseVector final
      : public VectorFunctionBase<Real, PointwiseVector<F>>
    {
      public:
        using Parent = VectorFunctionBase<Real, PointwiseVector<F>>;

        explicit PointwiseVector(F f) : m_f(std::move(f)) {}
        PointwiseVector(const PointwiseVector& other) : Parent(other), m_f(other.m_f) {}
        PointwiseVector(PointwiseVector&& other) : Parent(std::move(other)), m_f(std::move(other.m_f)) {}

        Math::Vector<Real> getValue(const Point& p) const
        {
          const Eigen::Vector2d v = m_f(p);
          Math::Vector<Real> out(2);
          out << v(0), v(1);
          return out;
        }
        size_t getDimension() const { return 2; }
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept { return 2; }
        PointwiseVector* copy() const noexcept override { return new PointwiseVector(*this); }

      private:
        F m_f;
    };

    template <class GF, class TF, class UF, class Pick>
    auto vectorFromPoint(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, std::uint16_t need, Pick pick)
    {
      return PointwiseVector([&psi, &psiT, &u, m, &revision, need, pick](const Point& p) -> Eigen::Vector2d {
        return pick(logPointAt(psi, psiT, u, m, revision, p, need)); });
    }

    /// @brief 1/r, at the (interior) quadrature points.
    auto inverseRadius()
    {
      return PointwiseScalar([](const Point& p) -> Real {
        return 1.0 / std::max<Real>(p.getPhysicalCoordinates()(1), 1.0e-300); });
    }

    /// @brief The parts of a LogPoint, each computed only when a form asks
    ///        for it: the forms are assembled one integrator at a time over
    ///        the whole mesh, so each pass recomputes the point and should pay
    ///        only for what it reads. psi_t, e^{psi_t}, a and da are always
    ///        set (no eigenproblem).
    enum LogPart : std::uint16_t
    {
      LogHoop = 0,    ///< psiT, expT, a, da
      LogExp = 1,     ///< exp(psi_P)
      LogD = 2,       ///< Dexp[B_j]
      LogN = 4,       ///< N_P, and f
      LogJN = 8,      ///< dN_P/dpsi_j
      LogJNT = 16,    ///< dN_P/dpsi_t
      LogQ = 32,      ///< Q_ab
      LogG = 64,      ///< g
      LogGT = 128,    ///< dg/dpsi_t
      LogGP = 256,    ///< dg/dpsi_j
      LogD2 = 512     ///< V_m, Vk (second derivative of exp along d_j psi^k)
    };

    /// @brief logPoint(psi_P(x), psi_t(x), grad u(x)), cached per quadrature
    ///        point as in the planar driver; the revision is bumped whenever
    ///        any of (psi_P, psi_t, u) changes. Only the parts in `need` are
    ///        guaranteed.
    ///
    /// @details In-plane block, as in the planar driver: exp, Dexp[B_j], N_P =
    ///          -Dexp^{-1}[L T + T L^T - c eps] + (f/lambda)(I - exp(-psi_P)),
    ///          JN_j = dN_P/dpsi_j, Q_ab = dN_P/dL_ab; plus JNT = dN_P/dpsi_t
    ///          (through f). Hoop block: N_t = a(psi_t) L_t + g(psi_P, psi_t),
    ///          L_t = u_r/r, with a, da = a', g and its derivatives gT, gP_j.
    ///          All derivatives in psi by central differences, as in 2D.
    template <class GF, class TF, class UF>
    const LogPoint& logPointAt(const GF& psi, const TF& psiT, const UF& u,
      const LogModel& m, std::uint64_t revision, const Point& p, std::uint16_t need)
    {
      struct Entry
      {
        const void* field = nullptr;
        std::uint64_t revision = 0;
        Index cell = 0;
        Real r0 = 0, r1 = 0;
        std::uint16_t have = 0;
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
        // Hoop block: Dexp^{-1} on a scalar divides by e^{psi_t}, so
        // N_t = -(2 L_t e^{psi_t} - c L_t) e^{-psi_t} + g = a L_t + g.
        LogPoint& out = hit->value;
        out.psiT = psiT.getValue(p);
        out.expT = std::exp(out.psiT);
        out.a = -2.0 + m.c / out.expT;
        out.da = -m.c / out.expT;
      }

      Entry& e = *hit;
      const std::uint16_t missing = need & ~e.have;
      if (!missing)
        return e.value;

      LogPoint& out = e.value;
      const Real pT = out.psiT;
      const auto& B = voigtBasisEigen();
      const Real h = 1.0e-6 *
        std::max<Real>({ 1.0, e.psi.cwiseAbs().maxCoeff(), std::abs(pT) });
      if ((missing & (LogN | LogJN | LogJNT)) && !e.hasL)
      {
        const auto J = Jacobian(u).getValue(p);
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
        for (size_t j = 0; j < 3; ++j)
          out.D[j] = e.eig.dexp(B[j]);
      if (missing & LogN)
        out.N = nonlinearMap(e.psi, pT, e.L, m, &out.f);
      if (missing & LogJN)
        for (size_t j = 0; j < 3; ++j)
          out.JN[j] = (nonlinearMap(e.psi + h * B[j], pT, e.L, m) -
            nonlinearMap(e.psi - h * B[j], pT, e.L, m)) / (2.0 * h);
      if (missing & LogJNT)
        out.JNT = (nonlinearMap(e.psi, pT + h, e.L, m) - nonlinearMap(e.psi, pT - h, e.L, m)) / (2.0 * h);
      if (missing & LogQ)
        for (size_t a = 0; a < 2; ++a)
          for (size_t b = 0; b < 2; ++b)
          {
            Eigen::Matrix2d E = Eigen::Matrix2d::Zero();
            E(a, b) = 1.0;
            out.Q[a][b] = -e.eig.dexpInverse(E * out.exp + out.exp * E.transpose() -
              0.5 * m.c * (E + E.transpose()));
          }
      if (missing & LogG)
        out.g = hoopRelaxation(e.psi, pT, m);
      if (missing & LogGT)
        out.gT = (hoopRelaxation(e.psi, pT + h, m) - hoopRelaxation(e.psi, pT - h, m)) / (2.0 * h);
      if (missing & LogGP)
        for (size_t j = 0; j < 3; ++j)
          out.gP[j] = (hoopRelaxation(e.psi + h * B[j], pT, m) -
            hoopRelaxation(e.psi - h * B[j], pT, m)) / (2.0 * h);
      if (missing & LogD2)
      {
        // d_j psi^k as 2x2 tensors, from the Voigt gradient (3 x 2).
        const auto G = Jacobian(psi).getValue(p);
        std::array<Eigen::Matrix2d, 2> g;
        for (size_t j = 0; j < 2; ++j)
        {
          g[j] << G(0, j), G(1, j), G(1, j), G(2, j);
        }
        out.Vk.setZero();
        const Real pk[3] = { e.psi(0, 0), e.psi(0, 1), e.psi(1, 1) };
        for (size_t mm = 0; mm < 3; ++mm)
        {
          const LogEig ep = logEig(e.psi + h * B[mm]);
          const LogEig em = logEig(e.psi - h * B[mm]);
          out.V[mm].setZero();
          for (size_t j = 0; j < 2; ++j)
          {
            const Eigen::Matrix2d d2 = (ep.dexp(g[j]) - em.dexp(g[j])) / (2.0 * h);
            out.V[mm] += d2.col(j);
          }
          out.Vk += out.V[mm] * pk[mm];
        }
      }
      e.have |= need;
      return out;
    }

    template <class GF, class TF, class UF, class Pick>
    auto fromPoint(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, std::uint16_t need, Pick pick)
    {
      return PointwiseTensor([&psi, &psiT, &u, m, &revision, need, pick](const Point& p) -> Eigen::Matrix2d {
        return pick(logPointAt(psi, psiT, u, m, revision, p, need)); });
    }

    template <class GF, class TF, class UF, class Pick>
    auto scalarFromPoint(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, std::uint16_t need, Pick pick)
    {
      return PointwiseScalar([&psi, &psiT, &u, m, &revision, need, pick](const Point& p) -> Real {
        return pick(logPointAt(psi, psiT, u, m, revision, p, need)); });
    }

    /// @brief Dexp_{psi_P(x)}[w] = sum_j w_j Dexp[B_j] for a Voigt w.
    template <class GF, class TF, class UF, class W>
    auto dexp(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, const W& w)
    {
      const auto D = [&](size_t j) {
        return fromPoint(psi, psiT, u, m, revision, LogD, [j](const LogPoint& q) { return q.D[j]; });
      };
      return D(0) * Component(w, 0) + D(1) * Component(w, 1) + D(2) * Component(w, 2);
    }

    /// @brief JN_{psi(x)}[w] = sum_j w_j dN_P/dpsi_j.
    template <class GF, class TF, class UF, class W>
    auto jacobianN(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, const W& w)
    {
      const auto J = [&](size_t j) {
        return fromPoint(psi, psiT, u, m, revision, LogJN, [j](const LogPoint& q) { return q.JN[j]; });
      };
      return J(0) * Component(w, 0) + J(1) * Component(w, 1) + J(2) * Component(w, 2);
    }

    /// @brief sum_ab L_ab(v) Q_ab: the exact in-plane velocity dependence.
    template <class GF, class TF, class UF, class V>
    auto velocityMap(const GF& psi, const TF& psiT, const UF& u, const LogModel& m,
      const std::uint64_t& revision, const V& v)
    {
      const auto Q = [&](size_t a, size_t b) {
        return fromPoint(psi, psiT, u, m, revision, LogQ, [a, b](const LogPoint& q) { return q.Q[a][b]; });
      };
      const auto col0 = Mult(Jacobian(v), unit(0));
      const auto col1 = Mult(Jacobian(v), unit(1));
      return Q(0, 0) * Component(col0, 0) + Q(1, 0) * Component(col0, 1) +
        Q(0, 1) * Component(col1, 0) + Q(1, 1) * Component(col1, 1);
    }

    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  // ==========================================================================
  // PulsatilePipeInflow
  // ==========================================================================
  PulsatilePipeInflow::Complex PulsatilePipeInflow::besselIScaled(int nu, Complex z)
  {
    assert(nu == 0 || nu == 1);
    if (std::abs(z) <= 17.0)
    {
      // sum_m (z/2)^{2m+nu}/(m! (m+nu)!), times e^{-z}.
      const Complex q = 0.25 * z * z;
      Complex term = (nu == 0) ? Complex(1.0) : 0.5 * z;
      Complex sum = term;
      for (int m = 1; m < 200; ++m)
      {
        term *= q / (static_cast<Real>(m) * static_cast<Real>(m + nu));
        sum += term;
        if (std::abs(term) < 1e-17 * std::abs(sum))
          break;
      }
      return sum * std::exp(-z);
    }
    // Hankel expansion, -pi/2 < arg z < 3pi/2, both exponentials kept since
    // arg k -> pi/2 as the fluid becomes elastic:
    //   I_nu(z) ~ e^z/sqrt(2 pi z) sum_k (-1)^k a_k/z^k
    //           + e^{-z + i pi (nu + 1/2)}/sqrt(2 pi z) sum_k a_k/z^k.
    const Real mu = 4.0 * nu * nu;
    Complex a(1.0), s1(1.0), s2(1.0);
    for (int k = 1; k < 40; ++k)
    {
      const Complex next = a * (mu - static_cast<Real>((2 * k - 1) * (2 * k - 1))) /
        (static_cast<Real>(8 * k) * z);
      if (std::abs(next) > std::abs(a))
        break;   // asymptotic series: stop at the smallest term
      a = next;
      s1 += ((k % 2) ? -1.0 : 1.0) * a;
      s2 += a;
      if (std::abs(a) < 1e-17)
        break;
    }
    const Complex phase = std::exp(Complex(0.0, M_PI * (nu + 0.5)));
    return (s1 + std::exp(-2.0 * z) * phase * s2) / std::sqrt(2.0 * M_PI * z);
  }

  PulsatilePipeInflow& PulsatilePipeInflow::setSinusoidal(Real amplitude)
  {
    // A sin(w t) = 2 Re[(-i A/2) e^{i w t}].
    m_c.assign({ Complex(1.0, 0.0), Complex(0.0, -0.5 * amplitude) });
    m_k.clear();
    return *this;
  }

  PulsatilePipeInflow& PulsatilePipeInflow::load(const std::string& path, size_t harmonics)
  {
    std::ifstream file(path);
    if (!file)
      throw std::runtime_error("Failed to open the flow-rate file " + path);

    std::vector<Real> t, q;
    std::string line;
    while (std::getline(file, line))
    {
      const auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || line[first] == '#')
        continue;
      std::istringstream row(line.substr(first));
      Real ti = 0.0, qi = 0.0;
      if (!(row >> ti >> qi))
        throw std::runtime_error("Malformed line in " + path + ": " + line);
      if (!t.empty() && !(ti > t.back()))
        throw std::runtime_error("Non-increasing time column in " + path);
      t.push_back(ti);
      q.push_back(qi);
    }
    if (t.size() < 3)
      throw std::runtime_error("Fewer than three samples in " + path);

    // Trapezoidal Fourier coefficients over one period, time normalised by
    // the file's own period: the file fixes the shape, the run fixes T.
    const Real period = t.back() - t.front();
    Real mean = 0.0;
    for (size_t i = 0; i + 1 < t.size(); ++i)
      mean += 0.5 * (q[i] + q[i + 1]) * (t[i + 1] - t[i]);
    mean /= period;
    if (!(mean > 0.0))
      throw std::runtime_error("The flow rate in " + path + " has no positive mean.");

    m_c.assign(harmonics + 1, Complex(0.0, 0.0));
    for (size_t n = 0; n <= harmonics; ++n)
    {
      Complex c(0.0, 0.0);
      for (size_t i = 0; i + 1 < t.size(); ++i)
      {
        const auto e = [&](size_t j) {
          const Real theta = -2.0 * M_PI * static_cast<Real>(n) * (t[j] - t.front()) / period;
          return (q[j] / mean) * Complex(std::cos(theta), std::sin(theta));
        };
        c += 0.5 * (e(i) + e(i + 1)) * (t[i + 1] - t[i]);
      }
      m_c[n] = c / period;
    }
    m_c[0] = Complex(1.0, 0.0);
    m_k.clear();
    return *this;
  }

  PulsatilePipeInflow& PulsatilePipeInflow::setFluid(
    Real rho, Real etaS, Real etaP, Real lambda, Real omega, Real radius)
  {
    m_omega = omega;
    m_radius = radius;
    m_k.assign(m_c.size(), Complex(0.0, 0.0));
    for (size_t n = 1; n < m_c.size(); ++n)
    {
      // Maxwell complex viscosity of the polymer plus the solvent.
      const Real w = static_cast<Real>(n) * omega;
      const Complex eta = etaS + etaP / Complex(1.0, w * lambda);
      m_k[n] = std::sqrt(Complex(0.0, w * rho) / eta);
    }
    return *this;
  }

  PulsatilePipeInflow::Complex PulsatilePipeInflow::mode(size_t n, Real r) const
  {
    assert(n < m_k.size());
    const Complex k = m_k[n];
    const Complex kR = k * m_radius;
    // Quasi-steady limit, where the ratio below is 0/0.
    if (std::abs(kR) < 1.0e-3)
      return Complex(2.0 * (1.0 - r * r / (m_radius * m_radius)), 0.0);
    // I_0(kr)/I_0(kR) = e^{k(r - R)} Ihat_0(kr)/Ihat_0(kR), Ihat = e^{-z} I:
    // nothing overflows at large Wo (Re k > 0, r <= R).
    const Complex i0R = besselIScaled(0, kR);
    const Complex ratio = std::exp(k * (r - m_radius)) * besselIScaled(0, k * r) / i0R;
    const Complex i1Ratio = besselIScaled(1, kR) / i0R;
    return (1.0 - ratio) / (1.0 - 2.0 * i1Ratio / kR);
  }

  PulsatilePipeInflow::Real PulsatilePipeInflow::flowRate(Real t) const
  {
    Real out = m_c[0].real();
    for (size_t n = 1; n < m_c.size(); ++n)
    {
      const Real theta = static_cast<Real>(n) * m_omega * t;
      out += 2.0 * (m_c[n] * Complex(std::cos(theta), std::sin(theta))).real();
    }
    return out;
  }

  PulsatilePipeInflow::Real PulsatilePipeInflow::velocity(Real r, Real t) const
  {
    assert(m_k.size() == m_c.size());
    r = std::min(std::abs(r), m_radius);
    Real out = m_c[0].real() * 2.0 * (1.0 - r * r / (m_radius * m_radius));
    for (size_t n = 1; n < m_c.size(); ++n)
    {
      const Real theta = static_cast<Real>(n) * m_omega * t;
      out += 2.0 * (m_c[n] * mode(n, r) * Complex(std::cos(theta), std::sin(theta))).real();
    }
    return out;
  }

  PulsatilePipeInflow::Real PulsatilePipeInflow::getPeak() const
  {
    if (!(m_omega > 0.0))
      return flowRate(0.0);
    const int samples = 2000;
    Real peak = flowRate(0.0);
    for (int i = 1; i < samples; ++i)
      peak = std::max(peak, flowRate(2.0 * M_PI * static_cast<Real>(i) / (samples * m_omega)));
    return peak;
  }

  // ==========================================================================
  // ArterialLesionAxiViscousLogImplicit
  // ==========================================================================
  ArterialLesionAxiViscousLogImplicit::AttributeSet
  ArterialLesionAxiViscousLogImplicit::makeWallSet(const Config& cfg)
  {
    return AttributeSet(cfg.labels.wall.begin(), cfg.labels.wall.end());
  }

  ArterialLesionAxiViscousLogImplicit::ArterialLesionAxiViscousLogImplicit(
    const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_wallSet(makeWallSet(m_cfg)),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_tauh(std::integral_constant<size_t, 1>{}, m_mesh, size_t(3)),
      m_ch(m_mesh),
      m_radius(m_sh),
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
      m_psiT(m_sh),
      m_chiT(m_sh),
      m_psiTOld(m_sh),
      m_psiTIt(m_sh),
      m_sigma(m_tauh),
      m_sigmaT(m_sh),
      m_wTrial(m_vh),
      m_wTest(m_vh),
      m_tauK(m_ch),
      m_alpha1(m_ch),
      m_alpha2(m_ch),
      m_alpha3(m_ch),
      m_alphaPsi(m_ch),
      m_piConv(m_vh, "al_proj_"),
      m_sub(m_vh, "al_proj_"),
      m_subOld(m_vh),
      m_piDiv(m_sh, "al_proj_"),
      m_piGradP(m_vh, "al_proj_"),
      m_piDivSigma(m_vh, "al_proj_"),
      m_piEps(m_tauh, "al_proj_"),
      m_piAdvPsi(m_tauh, "al_proj_"),
      m_piAdvPsiT(m_sh, "al_proj_"),
      m_wss(m_vh),
      m_symRec0(m_vh),
      m_symRec1(m_vh),
      m_netShear(m_vh),
      m_absShear(m_sh),
      m_shearMagnitude(m_sh),
      m_tawss(m_sh),
      m_osi(m_sh),
      m_qFlux(m_sh),
      m_one(m_sh),
      m_flux(m_qFlux),
      m_flow(m_u, m_p, m_psi, m_psiT, m_v, m_q, m_chi, m_chiT),
      m_flowKSP(m_flow),
      m_flowGN(m_u, m_p, m_v, m_q),
      m_flowGNKSP(m_flowGN),
      m_vectorProjection(m_wTrial, m_wTest),
      m_vectorProjectionKSP(m_vectorProjection),
      m_wssTrial(m_vh),
      m_wssTest(m_vh),
      m_wssProjection(m_wssTrial, m_wssTest),
      m_wssKSP(m_wssProjection)
  {
    deriveParameters();

    const auto cellCount = m_mesh.getCellCount();
    const auto vertexCount = m_mesh.getVertexCount();
    if (isRoot())
      Alert::Info() << "[mesh] " << m_cfg.meshPath << "  cells=" << cellCount
                    << " vertices=" << vertexCount << " velocity DOFs=" << m_vh.getSize()
                    << " pressure DOFs=" << m_sh.getSize() << Alert::Raise;
  }

  ArterialLesionAxiViscousLogImplicit::~ArterialLesionAxiViscousLogImplicit() = default;

  bool ArterialLesionAxiViscousLogImplicit::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  void ArterialLesionAxiViscousLogImplicit::deriveParameters()
  {
    auto& ob = m_cfg.oldroydB;
    const Real D = m_cfg.diameter;
    const Real rho = m_cfg.rho;

    if (m_cfg.fluid != "sptt" && m_cfg.fluid != "cy" && m_cfg.fluid != "gsptt")
      throw std::runtime_error("fluid must be sptt, cy or gsptt.");
    if ((m_cfg.fluid == "gsptt" || m_cfg.inletDevelopedConformation) &&
        (ob.lambda0Factor != 1.0 || ob.lambda0Min > 0.0))
      throw std::runtime_error("gsptt and inlet_developed_conformation require lambda0 = lambda "
                               "(lambda0_factor 1, lambda0_min 0).");
    if (m_cfg.fluid == "cy")
      ob.etaP = 0.0;   // no polymer: the (u, p) problem is solved
    const Real eta0 = referenceViscosity();

    if (!(m_cfg.reynolds > 0.0))
      throw std::runtime_error("Re must be positive.");
    if (!(m_cfg.womersley > 0.0))
      throw std::runtime_error("Wo must be positive; steady inflow is amplitude 0.");
    if (!(m_cfg.stepsPerCycle > 0))
      throw std::runtime_error("stepsPerCycle must be positive.");

    m_meanVelocity = m_cfg.reynolds * eta0 / (rho * D);
    const Real omega = 4.0 * m_cfg.womersley * m_cfg.womersley * eta0 / (rho * D * D);
    m_period = 2.0 * M_PI / omega;
    m_dt = m_period / static_cast<Real>(m_cfg.stepsPerCycle);

    if (m_cfg.deborah > 0.0)
      ob.lambda = m_cfg.deborah * m_period;
    else if (m_cfg.weissenberg > 0.0)
      ob.lambda = m_cfg.weissenberg * D / m_meanVelocity;
    if (!(ob.lambda > 0.0))
      throw std::runtime_error("lambda must be positive (set lambda, Wi or De).");

    if (m_cfg.flowWaveformPath.empty())
      m_inflow.setSinusoidal(m_cfg.amplitude);
    else
      m_inflow.load(m_cfg.flowWaveformPath, static_cast<size_t>(std::max(1, m_cfg.harmonics)));
    // Developed pulsatile profile of the inlet: the Maxwell/Oldroyd-B modes
    // of the solvent at its high-shear viscosity and the polymer. For "cy" and
    // "gsptt" the exact developed profile is not available in closed form;
    // the inlet lies 10 D upstream of the lesion, where it adjusts.
    if (m_cfg.fluid == "sptt")
      m_inflow.setFluid(rho, ob.etaS, ob.etaP, ob.lambda, omega, 0.5 * D);
    else
      m_inflow.setFluid(rho, m_cfg.carreauYasuda.etaInf, ob.etaP, ob.lambda, omega, 0.5 * D);

    if (m_cfg.fluid != "sptt")
    {
      // The solvent of "gsptt" must stay positive on the whole shear-rate range.
      for (Real gd = 1.0e-3; gd < 1.0e5; gd *= 1.5)
        if (!(solventViscosity(gd) > 0.0))
          throw std::runtime_error("gsptt: eta_CY(gd) - eta_p/f(gd) <= 0 at gd = " +
            std::to_string(gd) + " 1/s; lower eta_p or epsilon.");
    }

    if (isRoot())
    {
      const Real wi = ob.lambda * m_meanVelocity / D;
      const Real de = ob.lambda / m_period;
      Alert::Info() << "[groups] Re=" << m_cfg.reynolds << "  Wo=" << m_cfg.womersley
                    << "  Wi=lambda U/D=" << wi << "  De=lambda/T=" << de
                    << "  El=Wi/Re=" << (wi / m_cfg.reynolds)
                    << "  UT/D=" << (m_meanVelocity * m_period / D)
                    << "  beta=" << (ob.etaS / (ob.etaS + ob.etaP)) << "  epsilon=" << ob.pttEpsilon
                    << "  fluid=" << m_cfg.fluid << "  eta_ref=" << eta0 << Alert::Raise;
      if (m_cfg.fluid != "sptt")
      {
        const auto& cy = m_cfg.carreauYasuda;
        Alert::Info() << "[fluid] Carreau-Yasuda eta_0=" << cy.etaZero << " eta_inf=" << cy.etaInf
                      << " lambda=" << cy.lambda << " n=" << cy.n << " a=" << cy.a
                      << "  solvent eta_s(gd) at 1, 10, 100, 1000 1/s = " << solventViscosity(1.0)
                      << ", " << solventViscosity(10.0) << ", " << solventViscosity(100.0) << ", "
                      << solventViscosity(1000.0) << " Pa s" << Alert::Raise;
      }
      Alert::Info() << "[inflow] D=" << D << " m  Ubar=" << m_meanVelocity
                    << " m/s  T=" << m_period << " s  lambda=" << ob.lambda
                    << " s  lambda0=" << lambda0() << " s  dt=" << m_dt << " s  q_peak="
                    << m_inflow.getPeak() << "  modes="
                    << (m_cfg.flowWaveformPath.empty() ? std::string("sinusoid")
                                                       : m_cfg.flowWaveformPath)
                    << Alert::Raise;
      if (m_cfg.flowWaveformPath.empty() && std::abs(m_cfg.amplitude) > 1.0)
        Alert::Warning() << "[inflow] |A| > 1: the inlet reverses, and psi = 0 is then "
                            "imposed on an outflow boundary." << Alert::Raise;
    }
  }

  ArterialLesionAxiViscousLogImplicit::MeshType ArterialLesionAxiViscousLogImplicit::makeMesh(
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
          "ArterialLesionAxi expects a triangular MEDIT mesh written as "
          "\"Dimension 2\"; generate it with "
          "examples/viscoelastic_fluids/make_lesion_mesh.py --axi.");

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
    mesh.scale(cfg.diameter);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.reconcile(1);

    return mesh;
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  bool ArterialLesionAxiViscousLogImplicit::isGeneralizedNewtonian() const noexcept
  {
    return m_cfg.fluid == "cy";
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::referenceViscosity() const
  {
    if (m_cfg.etaRef > 0.0)
      return m_cfg.etaRef;
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS + m_cfg.oldroydB.etaP;
    return m_cfg.carreauYasuda.etaInf;
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::carreauYasuda(Real gammaDot) const
  {
    const auto& cy = m_cfg.carreauYasuda;
    return cy.etaInf + (cy.etaZero - cy.etaInf) *
      std::pow(1.0 + std::pow(cy.lambda * gammaDot, cy.a), (cy.n - 1.0) / cy.a);
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::steadyPTTFactor(Real gammaDot) const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real c = 2.0 * ob.pttEpsilon * ob.lambda * ob.lambda * gammaDot * gammaDot;
    if (c <= 0.0)
      return 1.0;
    // f^3 - f^2 - c = 0 has one root f >= 1; Newton from above converges
    // monotonically.
    Real f = 1.0 + std::cbrt(c);
    for (int it = 0; it < 50; ++it)
    {
      const Real g = f * f * f - f * f - c;
      const Real dg = 3.0 * f * f - 2.0 * f;
      const Real step = g / dg;
      f -= step;
      if (std::abs(step) < 1.0e-13 * f)
        break;
    }
    return std::max<Real>(f, 1.0);
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::solventViscosity(Real gammaDot) const
  {
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS;
    const Real etaCY = carreauYasuda(gammaDot);
    if (m_cfg.fluid == "cy")
      return etaCY;
    return etaCY - m_cfg.oldroydB.etaP / steadyPTTFactor(gammaDot);
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::shearRateAt(const Point& p) const
  {
    const auto J = Jacobian(m_uOld).getValue(p);
    const Real exx = J(0, 0), err = J(1, 1), exr = 0.5 * (J(0, 1) + J(1, 0));
    const Real r = p.getPhysicalCoordinates()(1);
    const Real ett = r > 1.0e-300 ? m_uOld.getValue(p)(1) / r : err;
    return std::sqrt(2.0 * (exx * exx + err * err + ett * ett + 2.0 * exr * exr));
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::solventViscosityAt(const Point& p) const
  {
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS;
    return solventViscosity(shearRateAt(p));
  }

  Math::SpatialVector<ArterialLesionAxiViscousLogImplicit::Real>
  ArterialLesionAxiViscousLogImplicit::steadyShearLogConformation(Real dudr) const
  {
    const Real f = steadyPTTFactor(std::abs(dudr));
    const Real sxr = m_cfg.oldroydB.lambda * dudr / f;
    Eigen::Matrix2d cm;
    cm << 1.0 + 2.0 * sxr * sxr, sxr, sxr, 1.0;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(cm);
    const Eigen::Matrix2d lg = eig.eigenvectors() *
      eig.eigenvalues().array().log().matrix().asDiagonal() * eig.eigenvectors().transpose();
    return pack(lg);
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = solventViscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::pttFactor(Real traceExp) const
  {
    const auto& ob = m_cfg.oldroydB;
    if (ob.pttEpsilon == 0.0)
      return 1.0;
    return 1.0 + ob.pttEpsilon * (ob.lambda / lambda0()) * (traceExp - 3.0);
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::alpha3At(const Point& p) const
  {
    const auto& ob = m_cfg.oldroydB;
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real speed = std::sqrt(Math::dot(uc, uc));
    const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
    const Real f = pttFactor(logEig(unpack(m_psiOld.getValue(p))).exp().trace() +
      std::exp(m_psiTOld.getValue(p)));
    return 1.0 / (4.0 * f / (2.0 * ob.etaP) +
      0.25 * (ob.lambda * speed / (2.0 * ob.etaP * h) + ob.lambda * gradNorm / ob.etaP));
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::lambda0() const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = std::max(ob.lambda0Factor * ob.lambda, ob.lambda0Min);
    if (!(l0 > 0.0))
      throw std::runtime_error("lambda0 = max(k lambda, lambda0_min) must be positive.");
    return l0;
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::ramp(Real t) const
  {
    const Real tr = m_cfg.rampCycles * m_period;
    if (!(tr > 0.0) || t >= tr)
      return 1.0;
    return 0.5 * (1.0 - std::cos(M_PI * t / tr));
  }

  ArterialLesionAxiViscousLogImplicit::Real ArterialLesionAxiViscousLogImplicit::inflowFactor(Real t) const
  {
    return ramp(t) * m_inflow.flowRate(t);
  }

  void ArterialLesionAxiViscousLogImplicit::updateConformation()
  {
    const Real s = m_cfg.oldroydB.etaP / lambda0();
    const auto field = [](auto f) { return VectorFunction(size_t(3), f); };

    m_sigma.project(field([&](const Point& p) -> Math::SpatialVector<Real> {
      return pack(s * (logEig(unpack(m_psiIt.getValue(p))).exp() - Eigen::Matrix2d::Identity()));
    }));
    m_sigmaT.project(RealFunction([&](const Point& p) -> Real {
      return s * std::expm1(m_psiTIt.getValue(p));
    }));
  }

  void ArterialLesionAxiViscousLogImplicit::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  ArterialLesionAxiViscousLogImplicit& ArterialLesionAxiViscousLogImplicit::initialize()
  {
    setupSpaces();
    if (isGeneralizedNewtonian())
      setupFlowGN();
    else
      setupFlow();
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

  void ArterialLesionAxiViscousLogImplicit::setupSpaces()
  {
    const auto zeroVector = Math::SpatialVector<Real>{{0.0, 0.0}};
    const auto zeroTensor = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_radius.project(RealFunction([](const Point& p) -> Real {
      return std::max<Real>(p.getPhysicalCoordinates()(1), 0.0); }));

    m_uOld = zeroVector;
    m_uIt = zeroVector;
    m_subOld = zeroVector;
    m_sub.get() = zeroVector;
    m_psiOld = zeroTensor;
    m_psiIt = zeroTensor;
    m_psiTOld = Real(0);
    m_psiTIt = Real(0);
    ++m_psiRevision;
    m_wss = zeroVector;
    m_netShear = zeroVector;
    m_symRec0 = zeroVector;
    m_symRec1 = zeroVector;

    m_absShear = Real(0);
    m_shearMagnitude = Real(0);
    m_tawss = Real(0);
    m_osi = Real(0);
    m_one = Real(1);

    configureMassSolver(m_vectorProjectionKSP, "al_vproj_");
    configureMassSolver(m_wssKSP, "al_wss_");

    updateConformation();  // psi = 0: tau = I, sigma = 0

    m_u.setName("velocity");
    m_p.setName("pressure");
    m_psi.setName("logConformation");
    m_psiT.setName("logConformationTheta");
    m_sigma.setName("elasticStress");
    m_sigmaT.setName("elasticStressTheta");
    m_tawss.setName("TAWSS");
    m_osi.setName("OSI");
    m_wss.setName("shearStress");

    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("logConformation", m_psi.getSolution());
    m_xdmf.add("logConformationTheta", m_psiT.getSolution());
    m_xdmf.add("elasticStress", m_sigma);
    m_xdmf.add("elasticStressTheta", m_sigmaT);
    m_xdmf.add("TAWSS", m_tawss);
    m_xdmf.add("OSI", m_osi);
    m_xdmf.add("shearStress", m_wss);

    // int r ds over the inlet is R^2/2 = D^2/8 for the meridian half-plane;
    // the planar channel gives 0 (y in [-R, R]), any other scale is wrong.
    m_inletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.inlet});
    m_outletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.outlet});
    const Real expected = m_cfg.diameter * m_cfg.diameter / 8.0;

    if (isRoot())
      Alert::Info() << "[boundary] int r ds: inlet = " << m_inletMeasure
                    << "  outlet = " << m_outletMeasure << " m^2  (D^2/8 = " << expected
                    << ")" << Alert::Raise;
    if (std::abs(m_inletMeasure / expected - 1.0) > 1.0e-2)
      throw std::runtime_error(
        "The inlet does not satisfy int r ds = D^2/8: the mesh is not the "
        "axisymmetric half-plane in units of D (make_lesion_mesh.py --axi).");
  }

  ArterialLesionAxiViscousLogImplicit::Real
  ArterialLesionAxiViscousLogImplicit::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_radius, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-300);
  }

  void ArterialLesionAxiViscousLogImplicit::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const auto& ob = m_cfg.oldroydB;
    const Real eta0 = referenceViscosity();
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;                  // sigma = s (exp(psi) - I)
    // eta_s(gd(u^n)); constant (order 0) for "sptt", so that assembly is the
    // parent one.
    const auto etaS = PointwiseScalar([this](const Point& p) -> Real {
      return solventViscosityAt(p); }, m_cfg.fluid == "sptt" ? size_t(0) : size_t(2));
    const Real chiScale = 2.0;
    const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };
    const auto& r = m_radius;
    const auto invR = inverseRadius();

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto ur = Component(m_u, 1);
    const auto vr = Component(m_v, 1);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    // r div u^n = r div_P u^n + u^n_r, in the skew-symmetric (Temam) term.
    const auto temam = (r * Div(m_uOld) + Dot(m_uOld, unit(1))) * Dot(m_u, m_v);

    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);
    const Real axisGamma = m_cfg.axisPenalty * eta0 / m_cfg.diameter;

    RealFunction pOutFn = [this](const Point&) { return m_cfg.outletPressure; };
    VectorFunction inletVelocity(dim, [this](const Point& p) -> Math::SpatialVector<Real> {
      const Real y = p.getPhysicalCoordinates()(1);
      return Math::SpatialVector<Real>{{
        m_meanVelocity * ramp(m_t) * m_inflow.velocity(y, m_t), 0.0 }};
    });

    const auto streamV = Mult(Jacobian(m_v), m_uOld);
    const auto X = tensor(m_chi);
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(2, 2));

    // Everything pointwise is at (psi_P^k, psi_t^k, u^k).
    const auto at = [this, model](std::uint16_t need, auto pick) {
      return fromPoint(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, need, pick); };
    const auto atS = [this, model](std::uint16_t need, auto pick) {
      return scalarFromPoint(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, need, pick); };
    const auto Dexp = [this, model](const auto& w) {
      return dexp(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, w); };
    const auto JN = [this, model](const auto& w) {
      return jacobianN(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, w); };
    const auto Gu = [this, model](const auto& v) {
      return velocityMap(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, v); };
    const auto divExp = [&](const auto& psi) {
      return Mult(Dexp(Mult(Jacobian(psi), unit(0))), unit(0)) +
        Mult(Dexp(Mult(Jacobian(psi), unit(1))), unit(1));
    };
    // Newton on div exp(psi): divExp carries Dexp_{psi^k}[d_j psi]; the
    // dependence of Dexp on psi adds sum_m V_m (psi_m - psi^k_m). Without it
    // the Newton iterations converge only linearly once psi is large.
    const auto atV = [this, model](auto pick) {
      return vectorFromPoint(m_psiIt, m_psiTIt, m_uIt, model, m_psiRevision, LogD2, pick); };
    const auto divExpCorrection = [&](const auto& psi) {
      return atV([](const LogPoint& q) { return q.V[0]; }) * Component(psi, 0) +
        atV([](const LogPoint& q) { return q.V[1]; }) * Component(psi, 1) +
        atV([](const LogPoint& q) { return q.V[2]; }) * Component(psi, 2);
    };
    const auto divExpCorrectionKnown = atV([](const LogPoint& q) { return q.Vk; });

    // In-plane momentum: exp(psi_P) ~ T + C at psi^k, as in 2D.
    const auto T = Dexp(m_psi);
    const auto Tn = at(LogExp, [](const LogPoint& q) { return q.exp; });
    const auto C = Tn - Dexp(m_psiIt);
    // Hoop stress: e^{psi_t} ~ e^{psi_t^k}(1 + psi_t - psi_t^k).
    const auto expT = atS(LogHoop, [](const LogPoint& q) { return q.expT; });
    const auto sigmaTKnown = atS(LogHoop, [s](const LogPoint& q) {
      return s * (q.expT * (1.0 - q.psiT) - 1.0); });

    // In-plane psi equation, Newton about (psi^k, u^k), plus dN_P/dpsi_t.
    const auto N = at(LogN, [](const LogPoint& q) { return q.N; });
    const auto JNT = at(LogJNT, [](const LogPoint& q) { return q.JNT; });
    const auto known = N - JN(m_psiIt) - JNT * m_psiTIt - Gu(m_uIt)
      - (1.0 / dt) * tensor(m_psiOld) - advection(m_psiIt, m_uIt);

    // Hoop psi equation: r [d_t psi_t + u . grad psi_t + g] + a u_r = 0,
    // linearised about (psi^k, u^k); a L_t is bilinear in (psi_t, u).
    const auto aT = atS(LogHoop, [](const LogPoint& q) { return q.a; });
    const auto daT = atS(LogHoop, [](const LogPoint& q) { return q.da; });
    const auto gV = atS(LogG, [](const LogPoint& q) { return q.g; });
    const auto gT = atS(LogGT, [](const LogPoint& q) { return q.gT; });
    const auto gP = [&](size_t j) {
      return atS(LogGP, [j](const LogPoint& q) { return q.gP[j]; }); };
    const auto gPsi = [&](const auto& w) {
      return gP(0) * Component(w, 0) + gP(1) * Component(w, 1) + gP(2) * Component(w, 2); };
    const auto knownT =
      r * (gV - gT * m_psiTIt - gPsi(m_psiIt) - (1.0 / dt) * m_psiTOld
        - Dot(m_uIt, Grad(m_psiTIt)))
      - daT * Dot(m_uIt, unit(1)) * m_psiTIt;

    // Full divergence for the grad-div term: div u = div_P u + u_r/r.
    const auto divU = Div(m_u) + invR * ur;
    const auto divV = Div(m_v) + invR * vr;
    const auto advT = Dot(m_uOld, Grad(m_chiT));   // u^n . grad chi_t

    const auto body =
      // ---- Momentum (weighted by r) ----------------------------------------
        (rho / dt) * Integral(r * m_u, m_v) - (rho / dt) * Integral(r * m_uOld, m_v)
      + rho * Integral(r * convU, m_v) + 0.5 * rho * Integral(temam)
      + 2.0 * Integral(r * (etaS * symU), symV) + 2.0 * Integral((etaS * invR) * ur, vr)
      + s * Integral(r * T, symV) + s * Integral(r * C, symV) - s * Integral(r * I, symV)
      + s * Integral(expT * m_psiT, vr) + Integral(sigmaTKnown, vr)
      - Integral(r * m_p, Div(m_v)) - Integral(m_p, vr)

      // ---- Continuity: r div u = r div_P u + u_r ----------------------------
      + Integral(r * Div(m_u), m_q) + Integral(ur, m_q)
      + m_cfg.pressurePenalty * Integral(r * m_p, m_q)

      // ---- In-plane psi equation ----------------------------------------------
      + Integral(r * ((1.0 / dt) * tensor(m_psi) + JN(m_psi) + advection(m_psi, m_uIt)), X)
      + Integral(r * (JNT * m_psiT), X)
      + Integral(r * (advection(m_psiIt, m_u) + Gu(m_u)), X)
      + Integral(r * known, X)

      // ---- Hoop psi equation --------------------------------------------------
      + Integral(r * ((1.0 / dt) * m_psiT + gT * m_psiT + Dot(m_uIt, Grad(m_psiT))), m_chiT)
      + Integral(daT * Dot(m_uIt, unit(1)) * m_psiT, m_chiT)
      + Integral(r * gPsi(m_psi), m_chiT)
      + Integral(r * Dot(m_u, Grad(m_psiTIt)), m_chiT)
      + Integral(aT * ur, m_chiT)
      + Integral(knownT, m_chiT)

      // ---- S1, momentum ---------------------------------------------------
      + rho * rho * Integral(r * (m_tauK * convU), streamV)
      - rho * rho * Integral(r * (m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get())), streamV)
      + m_cfg.pressureScale * Integral(r * (m_alpha1 * Grad(m_p)), Grad(m_q))
      - m_cfg.pressureScale * Integral(r * (m_alpha1 * m_piGradP.get()), Grad(m_q))
      + chiScale * m_cfg.stressDivScale * s * Integral(r * (m_alpha1 * divExp(m_psi)), divergence(m_chi))
      + chiScale * m_cfg.stressDivScale * s * Integral(r * (m_alpha1 * divExpCorrection(m_psi)), divergence(m_chi))
      - chiScale * m_cfg.stressDivScale * s * Integral(r * (m_alpha1 * divExpCorrectionKnown), divergence(m_chi))
      - chiScale * m_cfg.stressDivScale * Integral(r * (m_alpha1 * m_piDivSigma.get()), divergence(m_chi))

      // ---- S2, continuity (full axisymmetric divergence) -----------------------
      + Integral(r * (m_alpha2 * divU), divV)
      - Integral(r * (m_alpha2 * m_piDiv.get()), divV)

      // ---- S3: u-psi compatibility, and the psi advection -------------------
      + Integral(r * (m_alpha3 * symU), symV)
      - Integral(r * (m_alpha3 * tensor(m_piEps.get())), symV)
      + Integral(r * (m_alphaPsi * advection(m_psi, m_uOld)), advection(m_chi, m_uOld))
      - Integral(r * (m_alphaPsi * tensor(m_piAdvPsi.get())), advection(m_chi, m_uOld))
      + Integral(r * (m_alphaPsi * Dot(m_uOld, Grad(m_psiT))), advT)
      - Integral(r * (m_alphaPsi * m_piAdvPsiT.get()), advT)

      // ---- Outlet: traction p_out n and directional do-nothing ---------------
      + BoundaryIntegral(r * pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)
      + BoundaryIntegral(r * outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Axis: u_r = 0 by penalty -------------------------------------------
      + BoundaryIntegral(RealFunction(axisGamma) * Dot(Mult(radialProjector(), m_u), m_v)).over(m_cfg.labels.axis)

      // ---- Inlet: developed pulsatile pipe profile; walls: no slip ------------
      + DirichletBC(m_u, inletVelocity).on(m_cfg.labels.inlet)
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // Inlet conformation: relaxed (psi = 0), or the steady simple-shear state
    // of the mean Poiseuille profile, du_x/dr = -16 Ubar r/D^2.
    VectorFunction inletConformation(size_t(3), [this](const Point& p) -> Math::SpatialVector<Real> {
      if (!m_cfg.inletDevelopedConformation)
        return Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
      const Real r = p.getPhysicalCoordinates()(1);
      const Real D = m_cfg.diameter;
      return steadyShearLogConformation(-16.0 * m_meanVelocity * r / (D * D));
    });
    if (m_cfg.inletConformation)
      m_flow = body + DirichletBC(m_psi, inletConformation).on(m_cfg.labels.inlet)
        + DirichletBC(m_psiT, RealFunction(0.0)).on(m_cfg.labels.inlet);
    else
      m_flow = body;
  }

  void ArterialLesionAxiViscousLogImplicit::setupFlowGN()
  {
    // The momentum, continuity, S1 (convection, pressure), S2, outlet, axis
    // and Dirichlet terms of setupFlow, with sigma = 0 and the Carreau-Yasuda
    // viscosity at u^n. No psi, no S3.
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const Real eta0 = referenceViscosity();
    const auto& r = m_radius;
    const auto invR = inverseRadius();
    const auto etaS = PointwiseScalar([this](const Point& p) -> Real {
      return solventViscosityAt(p); });

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto ur = Component(m_u, 1);
    const auto vr = Component(m_v, 1);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = (r * Div(m_uOld) + Dot(m_uOld, unit(1))) * Dot(m_u, m_v);

    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);
    const Real axisGamma = m_cfg.axisPenalty * eta0 / m_cfg.diameter;

    RealFunction pOutFn = [this](const Point&) { return m_cfg.outletPressure; };
    VectorFunction inletVelocity(dim, [this](const Point& p) -> Math::SpatialVector<Real> {
      const Real y = p.getPhysicalCoordinates()(1);
      return Math::SpatialVector<Real>{{
        m_meanVelocity * ramp(m_t) * m_inflow.velocity(y, m_t), 0.0 }};
    });

    const auto streamV = Mult(Jacobian(m_v), m_uOld);
    const auto divU = Div(m_u) + invR * ur;
    const auto divV = Div(m_v) + invR * vr;

    m_flowGN =
        (rho / dt) * Integral(r * m_u, m_v) - (rho / dt) * Integral(r * m_uOld, m_v)
      + rho * Integral(r * convU, m_v) + 0.5 * rho * Integral(temam)
      + 2.0 * Integral(r * (etaS * symU), symV) + 2.0 * Integral((etaS * invR) * ur, vr)
      - Integral(r * m_p, Div(m_v)) - Integral(m_p, vr)

      + Integral(r * Div(m_u), m_q) + Integral(ur, m_q)
      + m_cfg.pressurePenalty * Integral(r * m_p, m_q)

      + rho * rho * Integral(r * (m_tauK * convU), streamV)
      - rho * rho * Integral(r * (m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get())), streamV)
      + m_cfg.pressureScale * Integral(r * (m_alpha1 * Grad(m_p)), Grad(m_q))
      - m_cfg.pressureScale * Integral(r * (m_alpha1 * m_piGradP.get()), Grad(m_q))

      + Integral(r * (m_alpha2 * divU), divV)
      - Integral(r * (m_alpha2 * m_piDiv.get()), divV)

      + BoundaryIntegral(r * pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)
      + BoundaryIntegral(r * outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)
      + BoundaryIntegral(RealFunction(axisGamma) * Dot(Mult(radialProjector(), m_u), m_v)).over(m_cfg.labels.axis)

      + DirichletBC(m_u, inletVelocity).on(m_cfg.labels.inlet)
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);
  }

  void ArterialLesionAxiViscousLogImplicit::setupWallShear()
  {
    // The wall normal lies in the meridian plane (n_theta = 0), so the wall
    // traction is the in-plane one of the planar driver.
    const auto normal = BoundaryNormal(m_mesh);
    const auto etaS = PointwiseScalar([this](const Point& p) -> Real {
      return solventViscosityAt(p); }, m_cfg.fluid == "sptt" ? size_t(0) : size_t(2));
    const auto& sigma = m_sigma;
    const auto nx = Component(normal, 0);
    const auto ny = Component(normal, 1);

    const auto traction = VectorFunction(
      etaS * Dot(m_symRec0, normal) + Component(sigma, 0) * nx + Component(sigma, 1) * ny,
      etaS * Dot(m_symRec1, normal) + Component(sigma, 1) * nx + Component(sigma, 2) * ny);
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  void ArterialLesionAxiViscousLogImplicit::updateStabilization()
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
    if (!isGeneralizedNewtonian())
    {
      m_alpha3.project(RealFunction([this](const Point& p) {
        return m_cfg.stressScale * alpha3At(p); }));
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
    m_piDiv.project(Div(m_uOld) + inverseRadius() * Dot(m_uOld, unit(1)));
    m_piGradP.project(Grad(m_p.getSolution()));
    if (isGeneralizedNewtonian())
      return;
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;
    const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };
    const auto Dexp = [this, model](const auto& w) {
      return dexp(m_psiOld, m_psiTOld, m_uOld, model, m_psiRevision, w); };
    const auto twoEpsN = Jacobian(m_uOld) + Transpose(Jacobian(m_uOld));
    m_piDivSigma.project(s * (Mult(Dexp(Mult(Jacobian(m_psiOld), unit(0))), unit(0)) +
      Mult(Dexp(Mult(Jacobian(m_psiOld), unit(1))), unit(1))));
    m_piEps.project(voigt(0.5 * twoEpsN));
    m_piAdvPsi.project(Mult(Jacobian(m_psiOld), m_uOld));
    m_piAdvPsiT.project(Dot(m_uOld, Grad(m_psiTOld)));
  }

  bool ArterialLesionAxiViscousLogImplicit::solveFlow()
  {
    const bool trace = isRoot() && (m_step < 3 || m_cfg.traceNewton);
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

    ::KSPConvergedReason reason = KSP_CONVERGED_ITS;
    if (isGeneralizedNewtonian())
    {
      // One linear solve per step: convection and viscosity at u^n (Picard),
      // as the momentum equation of the viscoelastic problem.
      if (trace)
      {
        ThreeDInfo() << "Assembling the (u, p) system ("
                     << (m_vh.getSize() + m_sh.getSize()) << " unknowns) ..." << Alert::Raise;
        std::cout.flush();
      }
      const auto assemblyStart = CoronaryClock::now();
      m_flowGN.assemble();
      m_timing.assembly = secondsSince(assemblyStart);
      if (!m_flowGNFieldSplitsSet)
      {
        if (isRoot())
          ThreeDInfo() << "Creating field splits ..." << Alert::Raise;
        m_flowGN.setFieldSplits();
        m_flowGNFieldSplitsSet = true;
      }
      const auto solveStart = CoronaryClock::now();
      m_flowGN.solve(m_flowGNKSP);
      m_timing.solve = secondsSince(solveStart);
      PetscErrorCode ierr = KSPGetConvergedReason(m_flowGNKSP.getHandle(), &reason);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      m_conformationIts = 1;
      m_psiIncrement = 0.0;

      m_uOld.setData(m_u.getSolution().getData());
      m_uIt.setData(m_uOld.getData());
      if (m_cfg.useVMS)
        m_subOld.setData(m_sub.get().getData());
      m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
      m_stress = 0.0;
      const Real guard =
        m_cfg.maxVelocityFactor * m_meanVelocity * std::max<Real>(1.0, m_inflow.getPeak());
      return reason > 0 && std::isfinite(m_speed) && m_speed <= guard;
    }

    if (trace)
    {
      ThreeDInfo() << "Assembling the flow system ("
                   << (m_vh.getSize() + 2 * m_sh.getSize() + m_tauh.getSize())
                   << " unknowns) ..." << Alert::Raise;
      std::cout.flush();
    }

    PetscErrorCode ierr = PETSC_SUCCESS;
    ::Vec increment = PETSC_NULLPTR;
    ierr = VecDuplicate(m_psiIt.getData(), &increment);
    assert(ierr == PETSC_SUCCESS);
    ::Vec incrementT = PETSC_NULLPTR;
    ierr = VecDuplicate(m_psiTIt.getData(), &incrementT);
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

      // The Newton step is measured on the whole psi = (psi_P, psi_t).
      ierr = VecWAXPY(increment, -1.0, m_psiIt.getData(), m_psi.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(incrementT, -1.0, m_psiTIt.getData(), m_psiT.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      Real stepP = 0.0, stepT = 0.0, scaleP = 0.0, scaleT = 0.0;
      ierr = VecNorm(increment, NORM_INFINITY, &stepP);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(incrementT, NORM_INFINITY, &stepT);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_psiIt.getData(), NORM_INFINITY, &scaleP);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_psiTIt.getData(), NORM_INFINITY, &scaleT);
      assert(ierr == PETSC_SUCCESS);
      const Real stepSize = std::max(stepP, stepT);
      const Real psiScale = std::max(scaleP, scaleT);
      const Real omega = (std::isfinite(stepSize) && stepSize > m_cfg.newtonMaxStep)
        ? m_cfg.newtonMaxStep / stepSize : 1.0;
      m_psiIncrement = stepSize / std::max<Real>(1.0, psiScale);
      axpy(omega, increment, m_psiIt.getData());
      axpy(omega, incrementT, m_psiTIt.getData());
      ierr = VecWAXPY(uIncrement, -1.0, m_uIt.getData(), m_u.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      axpy(omega, uIncrement, m_uIt.getData());
      ++m_psiRevision;
      ++m_conformationIts;

      if (trace)
        KSPInfo() << "Newton " << m_conformationIts << ": max|dpsi| = " << stepSize
                  << " (hoop " << stepT << ")" << (omega < 1.0 ? " (damped)" : "")
                  << Alert::Raise;
      if (!std::isfinite(m_psiIncrement) ||
          (omega == 1.0 && m_psiIncrement < m_cfg.conformationTolerance))
        break;
    }
    ierr = VecDestroy(&increment);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDestroy(&incrementT);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDestroy(&uIncrement);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    if (isRoot() && m_conformationIts >= m_cfg.conformationIterations &&
        m_psiIncrement >= 100.0 * m_cfg.conformationTolerance)
      Alert::Warning() << "[log] Newton stopped after " << m_conformationIts
                       << " iterations with max|dpsi| = " << m_psiIncrement << Alert::Raise;

    m_uOld.setData(m_uIt.getData());
    m_u.getSolution().setData(m_uIt.getData());
    m_psiOld.setData(m_psiIt.getData());
    m_psi.getSolution().setData(m_psiIt.getData());
    m_psiTOld.setData(m_psiTIt.getData());
    m_psiT.getSolution().setData(m_psiTIt.getData());
    ++m_psiRevision;
    updateConformation();
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.get().getData());

    m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
    m_stress = std::max({ std::abs(m_sigma.max()), std::abs(m_sigma.min()),
      std::abs(m_sigmaT.max()), std::abs(m_sigmaT.min()) });

    const Real guard =
      m_cfg.maxVelocityFactor * m_meanVelocity * std::max<Real>(1.0, m_inflow.getPeak());
    return reason > 0 && std::isfinite(m_speed) && m_speed <= guard;
  }

  void ArterialLesionAxiViscousLogImplicit::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    const auto jac = Jacobian(uSol);
    const auto offDiagonal = Component(jac, 0, 1) + Component(jac, 1, 0);
    projectVector(VectorFunction(2.0 * Component(jac, 0, 0), offDiagonal), m_symRec0);
    projectVector(VectorFunction(offDiagonal, 2.0 * Component(jac, 1, 1)), m_symRec1);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

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

    axpy(m_dt, m_wss.getData(), m_netShear.getData());
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void ArterialLesionAxiViscousLogImplicit::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    m_tawss = Real(0);
    m_osi = Real(0);

    m_tawss.project(Region::Boundary,
      RealFunction([this, elapsed](const Point& p) -> Real {
        return m_absShear.getValue(p) / elapsed;
      }),
      onWall);

    m_osi.project(Region::Boundary, RealFunction([this](const Point& p) -> Real {
      const Real abs = m_absShear.getValue(p);
      if (abs <= 0.0)
        return 0.0;
      const auto net = m_netShear.getValue(p);
      const Real mag = std::sqrt(Math::dot(net, net));
      return std::clamp<Real>(0.5 * (1.0 - mag / abs), 0.0, 0.5);
    }),
      onWall);

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
  }

  void ArterialLesionAxiViscousLogImplicit::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();
    const auto& r = m_radius;

    // Q = 2 pi int u . n r ds (m^3/s); n is outward, so the inlet is negated.
    m_flux = BoundaryIntegral(r * Dot(uSol, normal), m_qFlux).over(m_cfg.labels.inlet);
    m_flux.assemble();
    m_qIn = -2.0 * M_PI * m_flux(m_one);

    m_flux = BoundaryIntegral(r * Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = 2.0 * M_PI * m_flux(m_one);

    // Area-averaged pressures: int p r ds / int r ds.
    m_flux = BoundaryIntegral(r * m_p.getSolution(), m_qFlux).over(m_cfg.labels.inlet);
    m_flux.assemble();
    m_inletPressure = m_flux(m_one) / m_inletMeasure;

    m_flux = BoundaryIntegral(r * m_p.getSolution(), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_outletMeanPressure = m_flux(m_one) / m_outletMeasure;
  }

  void ArterialLesionAxiViscousLogImplicit::writeCSVHeader()
  {
    m_csv << "t,cycle,qTarget,qIn,qOut,pInletMean,pOutletMean,dp,maxU,maxSigma,"
          << "newtonIts,dpsi,maxTAWSS,maxOSI\n";
  }

  void ArterialLesionAxiViscousLogImplicit::writeCSVRow(int cycle)
  {
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();

    if (!isRoot())
      return;

    const Real area = 0.25 * M_PI * m_cfg.diameter * m_cfg.diameter;
    const Real qTarget = m_meanVelocity * area * m_inflowNow;
    m_csv << m_t << ',' << cycle << ',' << qTarget << ',' << m_qIn << ',' << m_qOut << ','
          << m_inletPressure << ',' << m_outletMeanPressure << ','
          << (m_inletPressure - m_outletMeanPressure) << ',' << m_speed << ',' << m_stress
          << ',' << m_conformationIts << ',' << m_psiIncrement << ',' << tawss << ',' << osi
          << '\n';
    m_csv.flush();
  }

  int ArterialLesionAxiViscousLogImplicit::run()
  {
    if (!m_initialized)
      initialize();

    const int stepsPerCycle = m_cfg.stepsPerCycle;
    const int totalSteps = m_cfg.cycles * stepsPerCycle;

    if (isRoot())
      Alert::Info() << "[run] " << m_cfg.cycles << " cycles of " << stepsPerCycle
                    << " steps (" << totalSteps << " total), inflow ramped over "
                    << m_cfg.rampCycles << " cycle(s); XDMF every " << m_cfg.outputEvery
                    << " steps and at every cycle boundary" << Alert::Raise;

    Real cycleElapsed = 0.0;
    const auto runStart = CoronaryClock::now();

    for (int step = 0; step < totalSteps; ++step)
    {
      m_step = step;
      m_timing = Timing{};
      const auto stepStart = CoronaryClock::now();

      m_t += m_dt;
      const int cycle = step / stepsPerCycle;
      const bool endOfCycle = (step % stepsPerCycle == stepsPerCycle - 1);
      m_inflowNow = inflowFactor(m_t);

      if (!solveFlow())
      {
        Alert::Exception() << "[flow] diverged at step " << (step + 1)
                           << ": max|u| = " << m_speed << " m/s = "
                           << (m_speed / m_meanVelocity) << " Ubar" << Alert::Raise;
        return 1;
      }

      computeWallShear();
      cycleElapsed += m_dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      m_timing.fluxes = secondsSince(fluxStart);

      if (endOfCycle)
      {
        closeCycle(cycleElapsed);
        cycleElapsed = 0.0;

        const Real maxTawss = m_tawss.max();
        const Real maxOsi = m_osi.max();
        if (isRoot())
          Alert::Info() << "[cycle " << (cycle + 1) << "/" << m_cfg.cycles
                        << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                        << Alert::Raise;
      }

      writeCSVRow(cycle);

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
                      << (cycle + 1) << "/" << m_cfg.cycles << "  t = " << m_t << " s"
                      << Alert::Raise;
        Alert::Info() << "[axi] q=" << m_inflowNow << "  qIn=" << m_qIn << "  qOut=" << m_qOut
                      << " m^3/s  max|u|=" << m_speed << " m/s  max|sigma|=" << m_stress
                      << " Pa  dp=" << (m_inletPressure - m_outletMeanPressure)
                      << " Pa  newton=" << m_conformationIts << " (|dpsi|=" << m_psiIncrement
                      << ")" << Alert::Raise;
        Alert::Info() << "[timing] vms=" << m_timing.vms << "  asm=" << m_timing.assembly
                      << "  ksp=" << m_timing.solve << "  wss=" << m_timing.shear
                      << "  flux=" << m_timing.fluxes << "  out=" << m_timing.output
                      << "  total=" << m_timing.total << " s  |  ETA "
                      << (perStep * (totalSteps - step - 1) / 60.0) << " min" << Alert::Raise;
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
      using Simulation = Rodin::Examples::ViscoelasticFluids::ArterialLesionAxiViscousLogImplicit;
      Simulation::Config cfg;

      const auto getString = [](const char* name, std::string& out) {
        char buffer[1024];
        PetscBool set = PETSC_FALSE;
        PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, name, buffer, sizeof(buffer), &set);
        if (set)
          out = buffer;
      };
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

      getString("-al_mesh", cfg.meshPath);
      getString("-al_xdmf", cfg.xdmfBasename);
      getString("-al_csv", cfg.csvPath);
      getString("-al_waveform", cfg.flowWaveformPath);

      getReal("-al_diameter", cfg.diameter);
      getReal("-al_rho", cfg.rho);
      getReal("-al_eta_s", cfg.oldroydB.etaS);
      getReal("-al_eta_p", cfg.oldroydB.etaP);
      getReal("-al_lambda", cfg.oldroydB.lambda);
      getReal("-al_lambda0_factor", cfg.oldroydB.lambda0Factor);
      getReal("-al_lambda0_min", cfg.oldroydB.lambda0Min);
      getReal("-al_ptt_epsilon", cfg.oldroydB.pttEpsilon);

      getReal("-al_re", cfg.reynolds);
      getReal("-al_wo", cfg.womersley);
      getReal("-al_amplitude", cfg.amplitude);
      getReal("-al_wi", cfg.weissenberg);
      getReal("-al_de", cfg.deborah);
      getReal("-al_ramp_cycles", cfg.rampCycles);
      getInt("-al_harmonics", cfg.harmonics);
      getInt("-al_steps_per_cycle", cfg.stepsPerCycle);
      getInt("-al_cycles", cfg.cycles);
      getInt("-al_output_every", cfg.outputEvery);

      getReal("-al_outlet_pressure", cfg.outletPressure);
      getReal("-al_axis_penalty", cfg.axisPenalty);
      getReal("-al_vms_scale", cfg.vmsScale);
      getReal("-al_graddiv_scale", cfg.gradDivScale);
      getReal("-al_pressure_scale", cfg.pressureScale);
      getReal("-al_stress_div_scale", cfg.stressDivScale);
      getReal("-al_stress_scale", cfg.stressScale);
      getInt("-al_conformation_its", cfg.conformationIterations);
      getReal("-al_conformation_tol", cfg.conformationTolerance);
      getReal("-al_newton_max_step", cfg.newtonMaxStep);
      getReal("-al_max_velocity_factor", cfg.maxVelocityFactor);
      getString("-al_fluid", cfg.fluid);
      getReal("-al_eta_ref", cfg.etaRef);
      getReal("-al_cy_eta_zero", cfg.carreauYasuda.etaZero);
      getReal("-al_cy_eta_inf", cfg.carreauYasuda.etaInf);
      getReal("-al_cy_lambda", cfg.carreauYasuda.lambda);
      getReal("-al_cy_n", cfg.carreauYasuda.n);
      getReal("-al_cy_a", cfg.carreauYasuda.a);
      getBool("-al_vms", cfg.useVMS);
      getBool("-al_inlet_conformation", cfg.inletConformation);
      getBool("-al_inlet_developed_conformation", cfg.inletDevelopedConformation);
      getBool("-al_trace_newton", cfg.traceNewton);

      Simulation simulation(context, cfg);
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
