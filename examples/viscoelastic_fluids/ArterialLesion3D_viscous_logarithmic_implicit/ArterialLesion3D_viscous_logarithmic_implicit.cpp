/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ArterialLesion3D_viscous_logarithmic_implicit.cpp
 * @brief Driver for pulsatile sPTT/Oldroyd-B flow through an idealised
 *        stenosis or aneurysm in a 3D pipe, fully-implicit log-conformation;
 *        see the header.
 *
 * Run (from the build directory):
 *   mpirun -n 8 ./examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit \
 *     -al_mesh ../resources/examples/viscoelastic_fluids/S75_pipe_coarse.mesh \
 *     -al_re 300 -al_wo 4 -al_amplitude 0.5 -al_wi 1
 *
 * Options (all prefixed -al_), those of the 2D driver plus dt:
 *   mesh, xdmf, csv, waveform                      strings
 *   diameter, rho, eta_s, eta_p, lambda, ptt_epsilon,
 *   lambda0_factor, lambda0_min                    fluid and pipe (SI)
 *   re, wo, amplitude, wi, de, ramp_cycles         dimensionless inflow
 *   harmonics, steps_per_cycle, cycles, output_every
 *   dt                                             step (s); overrides steps_per_cycle
 *   outlet_pressure, vms_scale, graddiv_scale, pressure_scale,
 *   stress_div_scale, stress_scale, conformation_its, conformation_tol,
 *   newton_max_step, max_velocity_factor
 *   vms, inlet_conformation, trace_newton          booleans
 *   inlet_perturbation, perturbation_cycles        m = 1 seed (see the header)
 *   newton_d2 (bool), psi_order (int)              Newton D2exp term; psi quadrature degree
 *   adaptive_dt (bool), max_step_retries, dt_shrink, dt_grow, dt_grow_after,
 *   newton_easy, newton_hard                       step control (header, Config)
 *   newton_predict (bool), newton_predict_safety, newton_calibrate_every
 *                                                  Newton acceptance by C d_k^2
 *   pc_reuse_solves, pc_reuse_growth               lagged preconditioner (Krylov only)
 *   outlet_conformation_inflow (bool), psi_skew (bool)  psi outlet inflow; skew advection
 *
 * The CSV ends with maxTAWSSLesion, maxOSILesion (tag 5 only) and uyMean,
 * uzMean, the volume means of u_y and u_z over Ubar: both vanish for an
 * axisymmetric flow, so they measure the m = 1 (deflected or tilted jet) mode.
 */
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <unordered_map>
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

#include "ArterialLesion3D_viscous_logarithmic_implicit.h"
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

  // The block below is the anonymous namespace of
  // ArterialLesion2D_viscous_logarithmic_implicit.cpp with 2 -> 3 and psi
  // split into its diagonal (D) and off-diagonal (O) blocks.
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

    // ---- psi in Voigt order, 3D -------------------------------------------
    // Block D holds (xx, yy, zz), block O holds (xy, xz, yz). tensor<B>()
    // turns one block into its part of the symmetric tensor, so every pairing
    // below is psi:chi restricted to a pair of blocks. Since E_i : E_{3+k} = 0,
    // tensor<D>(a) : tensor<O>(b) = 0 pointwise.

    enum class Block { D, O };

    constexpr size_t voigtIndex(Block b, size_t c)
    {
      return (b == Block::D ? 0 : 3) + c;
    }

    /// @brief E_0..2 = e_i e_i; E_3 = e_x e_y + e_y e_x, E_4 = e_x e_z + e_z
    ///        e_x, E_5 = e_y e_z + e_z e_y.
    const std::array<Eigen::Matrix3d, 6>& voigtBasisEigen()
    {
      static const std::array<Eigen::Matrix3d, 6> B = [] {
        std::array<Eigen::Matrix3d, 6> b;
        for (auto& m : b)
          m.setZero();
        b[0](0, 0) = 1;
        b[1](1, 1) = 1;
        b[2](2, 2) = 1;
        b[3](0, 1) = b[3](1, 0) = 1;
        b[4](0, 2) = b[4](2, 0) = 1;
        b[5](1, 2) = b[5](2, 1) = 1;
        return b;
      }();
      return B;
    }

    MatrixFunction<Math::Matrix<Real>> voigtBasis(size_t k)
    {
      const Math::Matrix<Real> E = voigtBasisEigen()[k];
      return MatrixFunction(E);
    }

    /// @brief e_j, valued in Math::SpatialVector (no allocation per evaluation).
    VectorFunction<Real, Real, Real> unit(size_t j)
    {
      return VectorFunction{ j == 0 ? 1.0 : 0.0, j == 1 ? 1.0 : 0.0, j == 2 ? 1.0 : 0.0 };
    }

    /// @brief tensor<B>(a) : tensor<C>(b) = w_B delta_BC a . b, w_D = 1, w_O = 2.
    constexpr Real voigtWeight(Block b)
    {
      return b == Block::D ? 1.0 : 2.0;
    }

    /// @brief The block B of the tensor: sum_c E_{B,c} s_c.
    template <Block B, class S>
    auto tensor(const S& s)
    {
      return voigtBasis(voigtIndex(B, 0)) * Component(s, 0) +
        voigtBasis(voigtIndex(B, 1)) * Component(s, 1) +
        voigtBasis(voigtIndex(B, 2)) * Component(s, 2);
    }

    /// @brief The D and O blocks of a symmetric tensor function, for projections.
    template <class M>
    auto voigtD(const M& m)
    {
      return VectorFunction{ Component(m, 0, 0), Component(m, 1, 1), Component(m, 2, 2) };
    }

    template <class M>
    auto voigtO(const M& m)
    {
      return VectorFunction{ Component(m, 0, 1), Component(m, 0, 2), Component(m, 1, 2) };
    }

    /// @brief Block B of (u . grad) psi.
    template <Block B, class S, class U>
    auto advection(const S& s, const U& u)
    {
      return tensor<B>(Mult(Jacobian(s), u));
    }

    /// @brief Block B's contribution to (div chi)_i = d chi_ij / d x_j.
    template <Block B, class S>
    auto divergence(const S& s)
    {
      const auto dx = Mult(Jacobian(s), unit(0));
      const auto dy = Mult(Jacobian(s), unit(1));
      const auto dz = Mult(Jacobian(s), unit(2));
      if constexpr (B == Block::D)
      {
        return unit(0) * Component(dx, 0) + unit(1) * Component(dy, 1) +
          unit(2) * Component(dz, 2);
      }
      else
      {
        // chi_xy = s_0, chi_xz = s_1, chi_yz = s_2.
        return unit(0) * (Component(dy, 0) + Component(dz, 1)) +
          unit(1) * (Component(dx, 0) + Component(dz, 2)) +
          unit(2) * (Component(dx, 1) + Component(dy, 2));
      }
    }

    // ---- log-conformation ---------------------------------------------------

    /// @brief Eigen-decomposition of a symmetric 3x3 psi, with exp, Dexp and
    ///        Dexp^{-1} by Daleckii-Krein: in the eigenbasis R, R^T Dexp[H] R =
    ///        F o (R^T H R), F_ij = (e^li - e^lj)/(li - lj), F_ii = e^li; F > 0,
    ///        so Dexp^{-1} needs no special case for equal eigenvalues.
    struct LogEig
    {
        Eigen::Matrix3d R;
        Eigen::Vector3d l;
        Eigen::Matrix3d F;

        Eigen::Matrix3d exp() const
        {
          return R * l.array().exp().matrix().asDiagonal() * R.transpose();
        }
        Eigen::Matrix3d expMinus() const
        {
          return R * (-l).array().exp().matrix().asDiagonal() * R.transpose();
        }
        Eigen::Matrix3d dexp(const Eigen::Matrix3d& H) const
        {
          return R * F.cwiseProduct(R.transpose() * H * R) * R.transpose();
        }
        Eigen::Matrix3d dexpInverse(const Eigen::Matrix3d& H) const
        {
          return R * (R.transpose() * H * R).cwiseQuotient(F) * R.transpose();
        }
    };

    LogEig logEig(const Eigen::Matrix3d& psi)
    {
      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(psi);
      LogEig out;
      out.R = eig.eigenvectors();
      out.l = eig.eigenvalues();
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          const Real ei = std::exp(out.l(i)), d = out.l(j) - out.l(i);
          out.F(i, j) = (i == j) ? ei
            : (std::abs(d) > 1e-12 ? ei * std::expm1(d) / d : ei * (1.0 + 0.5 * d));
        }
      return out;
    }

    template <class VD, class VO>
    Eigen::Matrix3d unpack(const VD& d, const VO& o)
    {
      Eigen::Matrix3d m;
      m << d(0), o(0), o(1),
           o(0), d(1), o(2),
           o(1), o(2), d(2);
      return m;
    }

    Math::SpatialVector<Real> packD(const Eigen::Matrix3d& m)
    {
      return Math::SpatialVector<Real>{{ m(0, 0), m(1, 1), m(2, 2) }};
    }

    Math::SpatialVector<Real> packO(const Eigen::Matrix3d& m)
    {
      return Math::SpatialVector<Real>{{ m(0, 1), m(0, 2), m(1, 2) }};
    }

    /// @brief The constants of the psi equation.
    struct LogModel
    {
        Real lambda = 0.0;
        Real lambda0 = 0.0;
        Real epsilon = 0.0;   ///< sPTT
        Real c = 0.0;         ///< 2 (1 - lambda0/lambda), on eps(u)
    };

    /// @brief What the forms read at one quadrature point, from (psi^k, L^k =
    ///        grad u^k): exp(psi) and Dexp[E_j] (momentum, S1), and the Newton
    ///        data of the psi equation paired against the test blocks.
    ///
    /// @details With N = -Dexp^{-1}[L T + T L^T - c eps] + (f/lambda)(I -
    ///          exp(-psi)), its FD Jacobian JN_j = dN/dpsi_j, the exact
    ///          (linear) velocity dependence Q_ab = -Dexp^{-1}[E_ab T + T
    ///          E_ab^T - c sym(E_ab)] and R = N - JN[psi^k] - sum_ab L^k_ab
    ///          Q_ab, the pairings against block C (0 = D, 1 = O) are
    ///            A[C][B](k, c) = JN_{B,c} : E_{C,k},
    ///            P[C][b](k, a) = Q_ab : E_{C,k},
    ///            r[C](k) = R : E_{C,k}.
    struct LogPoint
    {
        Eigen::Matrix3d exp;
        std::array<Eigen::Matrix3d, 6> D;
        std::array<std::array<Eigen::Matrix3d, 2>, 2> A;
        std::array<std::array<Eigen::Matrix3d, 3>, 2> P;
        std::array<Eigen::Vector3d, 2> r;
        /// @brief V_m = sum_j (D2exp_{psi^k}[d_j psi^k, E_m]) e_j, m = 0..5:
        ///        the derivative of div exp(psi) = sum_j Dexp_psi[d_j psi] e_j
        ///        in psi beyond the one carried by d_j psi; Vk = sum_m V_m psi^k_m.
        std::array<Eigen::Vector3d, 6> V;
        Eigen::Vector3d Vk;
    };

    Eigen::Matrix3d nonlinearMap(
      const Eigen::Matrix3d& psi, const Eigen::Matrix3d& L, const LogModel& m, Real* fOut = nullptr)
    {
      const LogEig e = logEig(psi);
      const Eigen::Matrix3d T = e.exp();
      const Eigen::Matrix3d H = L * T + T * L.transpose() - 0.5 * m.c * (L + L.transpose());
      const Real f = 1.0 + m.epsilon * (m.lambda / m.lambda0) * (T.trace() - 3.0);
      if (fOut)
        *fOut = f;
      return -e.dexpInverse(H) +
        (f / m.lambda) * (Eigen::Matrix3d::Identity() - e.expMinus());
    }

    /// @brief The parts of a LogPoint. Dexp (exp and Dexp[E_j]) costs one
    ///        eigenproblem; Newton (A, P, r) thirteen more, twelve of them for
    ///        the FD Jacobian. The step-value store (stabilisation) reads only
    ///        Dexp.
    enum LogPart : std::uint8_t
    {
      LogDexp = 1,
      LogNewton = 2,
      LogD2 = 4       ///< V_m, Vk (second derivative of exp along d_j psi^k)
    };

    /// @brief exp(psi) and Dexp[E_j], from the eigen-decomposition of psi.
    void logPointDexp(LogPoint& out, const LogEig& e)
    {
      out.exp = e.exp();
      const auto& B = voigtBasisEigen();
      for (size_t j = 0; j < 6; ++j)
        out.D[j] = e.dexp(B[j]);
    }

    /// @brief A, P and r; reads out.exp, so logPointDexp() comes first.
    void logPointNewton(LogPoint& out, const LogEig& e, const Eigen::Matrix3d& psi,
      const Eigen::Matrix3d& L, const LogModel& m)
    {
      const auto& B = voigtBasisEigen();
      // Central differences in the six Voigt directions (twelve 3x3
      // eigenproblems); the map is smooth in psi.
      const Real h = 1.0e-6 * std::max<Real>(1.0, psi.cwiseAbs().maxCoeff());
      std::array<Eigen::Matrix3d, 6> JN;
      for (size_t j = 0; j < 6; ++j)
        JN[j] = (nonlinearMap(psi + h * B[j], L, m) - nonlinearMap(psi - h * B[j], L, m)) / (2.0 * h);
      std::array<std::array<Eigen::Matrix3d, 3>, 3> Q;
      for (size_t a = 0; a < 3; ++a)
        for (size_t b = 0; b < 3; ++b)
        {
          Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
          E(a, b) = 1.0;
          Q[a][b] = -e.dexpInverse(E * out.exp + out.exp * E.transpose() -
            0.5 * m.c * (E + E.transpose()));
        }
      // psi = sum_j psi_j E_j, psi_j = psi : E_j / |E_j|^2.
      Eigen::Matrix3d R = nonlinearMap(psi, L, m);
      for (size_t j = 0; j < 6; ++j)
        R -= (psi.cwiseProduct(B[j]).sum() / B[j].squaredNorm()) * JN[j];
      for (size_t a = 0; a < 3; ++a)
        for (size_t b = 0; b < 3; ++b)
          R -= L(a, b) * Q[a][b];
      for (size_t C = 0; C < 2; ++C)
        for (size_t k = 0; k < 3; ++k)
        {
          const Eigen::Matrix3d& Ek = B[3 * C + k];
          out.r[C](k) = R.cwiseProduct(Ek).sum();
          for (size_t Bk = 0; Bk < 2; ++Bk)
            for (size_t c = 0; c < 3; ++c)
              out.A[C][Bk](k, c) = JN[3 * Bk + c].cwiseProduct(Ek).sum();
          for (size_t a = 0; a < 3; ++a)
            for (size_t b = 0; b < 3; ++b)
              out.P[C][b](k, a) = Q[a][b].cwiseProduct(Ek).sum();
        }
    }

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

    /// @brief A 3x3 tensor field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseTensor final
      : public MatrixFunctionBase<Real, PointwiseTensor<F>>
    {
      public:
        using Parent = MatrixFunctionBase<Real, PointwiseTensor<F>>;

        explicit PointwiseTensor(F f) : m_f(std::move(f)) {}
        PointwiseTensor(const PointwiseTensor& other) : Parent(other), m_f(other.m_f) {}
        PointwiseTensor(PointwiseTensor&& other) : Parent(std::move(other)), m_f(std::move(other.m_f)) {}

        Math::SpatialMatrix<Real> getValue(const Point& p) const { return Math::SpatialMatrix<Real>(m_f(p)); }
        size_t getRows() const { return 3; }
        size_t getColumns() const { return 3; }
        /// Not a polynomial; quadratic is enough against P1 x P1.
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept { return 2; }
        PointwiseTensor* copy() const noexcept override { return new PointwiseTensor(*this); }

      private:
        F m_f;
    };

    /// @brief The fields a LogPoint is computed from: (psi_D, psi_O, u) and
    ///        the revision that keys the cache.
    template <class GF>
    struct LogFields
    {
        const GF* psiD;
        const GF* psiO;
        const GF* u;
        LogModel model;
        const std::uint64_t* revision;
    };

    /// @brief logPoint(psi(x), grad u(x)) at every quadrature point of the
    ///        current Newton iterate, computed once and kept until the
    ///        iterate changes.
    ///
    /// @details Rodin assembles integrator by integrator over the whole mesh
    ///          so a small recency cache recomputes the twelve eigenproblems
    ///          of every point once per integrator that reads them. One store per
    ///          psi_D field (the iterate psi^k and the step value psi^n),
    ///          cleared when the revision changes; keyed by the polytope
    ///          (dimension, index) and the reference coordinates. Every
    ///          integral reading it uses one rule (psiOrder in setupFlow), so
    ///          a cell holds one set of points.
    /// @brief A quadrature point: polytope (dimension, index), reference
    ///        coordinates.
    struct PointKey
    {
        size_t dimension;
        Index index;
        Real r0, r1, r2;
        bool operator==(const PointKey&) const = default;

        static PointKey of(const Point& p)
        {
          const auto& polytope = p.getPolytope();
          const auto& rc = p.getReferenceCoordinates();
          return PointKey{ polytope.getDimension(), polytope.getIndex(),
            rc(0), rc.size() > 1 ? rc(1) : 0.0, rc.size() > 2 ? rc(2) : 0.0 };
        }
    };

    struct PointKeyHash
    {
        size_t operator()(const PointKey& k) const noexcept
        {
          size_t h = std::hash<size_t>()(k.dimension);
          for (const size_t v : { std::hash<Index>()(k.index), std::hash<Real>()(k.r0),
                 std::hash<Real>()(k.r1), std::hash<Real>()(k.r2) })
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
          return h;
        }
    };

    class LogStore
    {
      public:
        /// @brief A stored point and the parts already computed for it.
        struct Entry
        {
            LogPoint value;
            std::uint8_t have = 0;
        };

        template <class GF>
        const Entry& at(const LogFields<GF>& s, const Point& p, const PointKey& key,
          std::uint8_t need)
        {
          const std::uint64_t revision = *s.revision;

          std::lock_guard<std::mutex> lock(m_mutex);
          if (revision != m_revision)
          {
            m_points.clear();
            m_revision = revision;
          }
          Entry& e = m_points[key];
          if (!(need & ~e.have))
            return e;

          // The Newton part reads exp, so it brings the Dexp part along.
          const Eigen::Matrix3d psi = unpack(s.psiD->getValue(p), s.psiO->getValue(p));
          const LogEig eig = logEig(psi);
          if (!(e.have & LogDexp))
          {
            logPointDexp(e.value, eig);
            e.have |= LogDexp;
          }
          if (need & LogNewton & ~e.have)
          {
            const auto J = Jacobian(*s.u).getValue(p);   // Math::SpatialMatrix, not Eigen
            Eigen::Matrix3d L;
            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                L(a, b) = J(a, b);
            logPointNewton(e.value, eig, psi, L, s.model);
            e.have |= LogNewton;
          }
          if (need & LogD2 & ~e.have)
          {
            // d_j psi^k = sum_c GD(c, j) E_c + GO(c, j) E_{3+c}.
            const auto GD = Jacobian(*s.psiD).getValue(p);
            const auto GO = Jacobian(*s.psiO).getValue(p);
            const auto& B = voigtBasisEigen();
            std::array<Eigen::Matrix3d, 3> g;
            for (int j = 0; j < 3; ++j)
            {
              g[j].setZero();
              for (int k = 0; k < 3; ++k)
                g[j] += GD(k, j) * B[k] + GO(k, j) * B[3 + k];
            }
            const Real h = 1.0e-6 * std::max<Real>(1.0, psi.cwiseAbs().maxCoeff());
            const Real pk[6] = { psi(0, 0), psi(1, 1), psi(2, 2), psi(0, 1), psi(0, 2), psi(1, 2) };
            e.value.Vk.setZero();
            for (size_t m = 0; m < 6; ++m)
            {
              const LogEig ep = logEig(psi + h * B[m]);
              const LogEig em = logEig(psi - h * B[m]);
              e.value.V[m].setZero();
              for (int j = 0; j < 3; ++j)
                e.value.V[m] += ((ep.dexp(g[j]) - em.dexp(g[j])) / (2.0 * h)).col(j);
              e.value.Vk += pk[m] * e.value.V[m];
            }
            e.have |= LogD2;
          }
          return e;
        }

        /// @brief Frees every stored point.
        void release()
        {
          std::lock_guard<std::mutex> lock(m_mutex);
          std::unordered_map<PointKey, Entry, PointKeyHash>().swap(m_points);
          m_revision = std::numeric_limits<std::uint64_t>::max();
        }

      private:
        std::mutex m_mutex;
        std::uint64_t m_revision = std::numeric_limits<std::uint64_t>::max();
        std::unordered_map<PointKey, Entry, PointKeyHash> m_points;
    };

    /// @brief One store per psi_D field, and a counter bumped by every
    ///        release so that no thread reuses a pointer into a freed store.
    struct LogStores
    {
        std::map<const void*, LogStore> stores;
        std::mutex mutex;
        std::atomic<std::uint64_t> epoch{0};
    };

    LogStores& logStores()
    {
      static LogStores all;
      return all;
    }

    /// @brief Frees the points stored for `field` (a psi_D field) once the
    ///        forms that read them have been assembled.
    void releaseLogPoints(const void* field)
    {
      auto& all = logStores();
      LogStore* store = nullptr;
      {
        std::lock_guard<std::mutex> lock(all.mutex);
        const auto it = all.stores.find(field);
        if (it == all.stores.end())
          return;
        store = &it->second;
      }
      all.epoch.fetch_add(1);
      store->release();
    }

    template <class GF>
    const LogPoint& logPointAt(const LogFields<GF>& s, const Point& p, std::uint8_t need)
    {
      // Consecutive evaluations are almost always at the same point: reuse
      // the last entry without locking. References into an unordered_map
      // survive insertions, and a store is cleared only when the revision
      // changes, between assemblies, or by a release, which bumps the epoch.
      struct Last
      {
          const void* field = nullptr;
          std::uint64_t revision = 0;
          std::uint64_t epoch = 0;
          PointKey key{};
          const LogStore::Entry* entry = nullptr;
      };
      thread_local Last last;
      auto& all = logStores();
      const PointKey key = PointKey::of(p);
      const std::uint64_t revision = *s.revision;
      const std::uint64_t epoch = all.epoch.load();
      if (last.entry && last.field == s.psiD && last.revision == revision &&
          last.epoch == epoch && last.key == key && (last.entry->have & need) == need)
        return last.entry->value;

      LogStore* store = nullptr;
      {
        std::lock_guard<std::mutex> lock(all.mutex);
        store = &all.stores[s.psiD];
      }
      const LogStore::Entry& entry = store->at(s, p, key, need);
      last = Last{ s.psiD, revision, epoch, key, &entry };
      return entry.value;
    }

    /// @brief A pointwise matrix drawn from the cached LogPoint by `pick`,
    ///        which reads only the parts in `need`.
    template <class GF, class Pick>
    auto fromPoint(const LogFields<GF>& s, std::uint8_t need, Pick pick)
    {
      return PointwiseTensor([s, need, pick](const Point& p) -> Eigen::Matrix3d {
        return pick(logPointAt(s, p, need)); });
    }

    /// @brief Dexp_{psi(x)}[w] restricted to block B: sum_c w_c Dexp[E_{B,c}],
    ///        for w a block of psi or of a derivative of it; trial or field.
    template <Block B, class GF, class W>
    auto dexp(const LogFields<GF>& s, const W& w)
    {
      const auto D = [&](size_t c) {
        const size_t j = voigtIndex(B, c);
        return fromPoint(s, LogDexp, [j](const LogPoint& q) { return q.D[j]; });
      };
      return D(0) * Component(w, 0) + D(1) * Component(w, 1) + D(2) * Component(w, 2);
    }

    constexpr size_t blockIndex(Block b)
    {
      return b == Block::D ? 0 : 1;
    }

    /// @brief A pointwise 3-vector drawn from the cached LogPoint by `pick`,
    ///        which reads only the parts in `need`.
    template <class GF, class Pick>
    auto fromPointVector(const LogFields<GF>& s, std::uint8_t need, Pick pick)
    {
      return VectorFunction(size_t(3), [s, need, pick](const Point& p) -> Math::SpatialVector<Real> {
        const Eigen::Vector3d v = pick(logPointAt(s, p, need));
        return Math::SpatialVector<Real>{{ v(0), v(1), v(2) }};
      });
    }

    // The psi equation is tested against each block chi_C as a plain vector:
    // for a tensor Y, Y : tensor<C>(chi_C) = pair<C>(Y) . chi_C, pair<C>(Y)_k =
    // Y : E_{C,k}. The pairing is folded into the pointwise data, so the test
    // side carries no coefficient.

    /// @brief pair<C>(JN[w] + mass w) for w = psi_B: (A[C][B] + mass I) w.
    template <Block C, Block B, class GF, class W>
    auto pairedJacobianN(const LogFields<GF>& s, Real mass, const W& w)
    {
      constexpr size_t c = blockIndex(C), b = blockIndex(B);
      return Mult(fromPoint(s, LogNewton, [mass](const LogPoint& q) -> Eigen::Matrix3d {
        return q.A[c][b] + mass * Eigen::Matrix3d::Identity(); }), w);
    }

    /// @brief pair<C>(sum_ab L_ab(v) Q_ab) = sum_b P[C][b] (J(v) e_b), for a
    ///        vector v (trial or field): the exact velocity dependence of the
    ///        psi equation at psi^k, L_ab = d v_a / d x_b.
    template <Block C, class GF, class V>
    auto pairedVelocityMap(const LogFields<GF>& s, const V& v)
    {
      constexpr size_t c = blockIndex(C);
      const auto P = [&](size_t b) {
        return fromPoint(s, LogNewton, [b](const LogPoint& q) -> Eigen::Matrix3d { return q.P[c][b]; });
      };
      return Mult(P(0), Mult(Jacobian(v), unit(0))) + Mult(P(1), Mult(Jacobian(v), unit(1))) +
        Mult(P(2), Mult(Jacobian(v), unit(2)));
    }

    /// @brief pair<C>(N - JN[psi^k] - sum_ab L^k_ab Q_ab), grade 0.
    template <Block C, class GF>
    auto pairedResidual(const LogFields<GF>& s)
    {
      constexpr size_t c = blockIndex(C);
      return fromPointVector(s, LogNewton, [](const LogPoint& q) -> Eigen::Vector3d { return q.r[c]; });
    }

    /// @brief Block B's contribution to div exp(psi) = sum_j Dexp[d_j psi] e_j
    ///        (chain rule), Dexp at the point of `s`.
    template <Block B, class GF, class S>
    auto divExp(const LogFields<GF>& s, const S& psi)
    {
      return Mult(dexp<B>(s, Mult(Jacobian(psi), unit(0))), unit(0)) +
        Mult(dexp<B>(s, Mult(Jacobian(psi), unit(1))), unit(1)) +
        Mult(dexp<B>(s, Mult(Jacobian(psi), unit(2))), unit(2));
    }

    /// @brief divExp<B> for a trial function psi_B: sum_c Dexp[E_{B,c}]
    ///        grad psi_c. A basis function of psi_B has one nonzero component
    ///        c, where this is Dexp[E_{B,c}] grad phi, the value divExp<B>
    ///        gives (the other terms are exact zeros), with three pointwise
    ///        matrices per point instead of nine.
    template <Block B, class GF, class S>
    auto divExpTrial(const LogFields<GF>& s, const S& psi)
    {
      const auto D = [&](size_t c) {
        const size_t j = voigtIndex(B, c);
        return fromPoint(s, LogDexp, [j](const LogPoint& q) { return q.D[j]; });
      };
      const auto gradOf = [&](size_t c) { return Mult(Transpose(Jacobian(psi)), unit(c)); };
      return Mult(D(0), gradOf(0)) + Mult(D(1), gradOf(1)) + Mult(D(2), gradOf(2));
    }

    /// @brief Newton on div exp(psi): divExpTrial carries Dexp_{psi^k}[d_j
    ///        psi]; the dependence of Dexp on psi adds sum_m V_m psi_m - Vk
    ///        (block B's part of the sum here, Vk in divExpKnownD2). Without
    ///        it Newton converges only linearly once psi is large.
    template <Block B, class GF, class S>
    auto divExpTrialD2(const LogFields<GF>& s, const S& psi)
    {
      const auto V = [&](size_t c) {
        const size_t m = voigtIndex(B, c);
        return fromPointVector(s, LogD2, [m](const LogPoint& q) -> Eigen::Vector3d { return q.V[m]; });
      };
      return V(0) * Component(psi, 0) + V(1) * Component(psi, 1) + V(2) * Component(psi, 2);
    }

    template <class GF>
    auto divExpKnownD2(const LogFields<GF>& s)
    {
      return fromPointVector(s, LogD2, [](const LogPoint& q) -> Eigen::Vector3d { return q.Vk; });
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

  PulsatilePipeInflow& PulsatilePipeInflow::setMeanProfile(std::vector<Real> r, std::vector<Real> phi)
  {
    assert(r.size() == phi.size() && r.size() >= 2);
    m_meanR = std::move(r);
    m_meanPhi = std::move(phi);
    return *this;
  }

  PulsatilePipeInflow::Real PulsatilePipeInflow::velocity(Real r, Real t) const
  {
    assert(m_k.size() == m_c.size());
    r = std::min(std::abs(r), m_radius);
    Real phi0 = 2.0 * (1.0 - r * r / (m_radius * m_radius));
    if (!m_meanR.empty())
    {
      const auto it = std::upper_bound(m_meanR.begin(), m_meanR.end(), r);
      const size_t j = std::min<size_t>(
        std::max<std::ptrdiff_t>(it - m_meanR.begin(), 1), m_meanR.size() - 1);
      const Real w = (r - m_meanR[j - 1]) / (m_meanR[j] - m_meanR[j - 1]);
      phi0 = (1.0 - w) * m_meanPhi[j - 1] + w * m_meanPhi[j];
    }
    Real out = m_c[0].real() * phi0;
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
  // ArterialLesion3DViscousLogImplicit
  // ==========================================================================
  ArterialLesion3DViscousLogImplicit::AttributeSet
  ArterialLesion3DViscousLogImplicit::makeWallSet(const Config& cfg)
  {
    return AttributeSet(cfg.labels.wall.begin(), cfg.labels.wall.end());
  }

  ArterialLesion3DViscousLogImplicit::ArterialLesion3DViscousLogImplicit(
    const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_wallSet(makeWallSet(m_cfg)),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_tauh(std::integral_constant<size_t, 1>{}, m_mesh, 3),
      m_ch(m_mesh),
      m_u(m_vh),
      m_p(m_sh),
      m_v(m_vh),
      m_q(m_sh),
      m_uOld(m_vh),
      m_uIt(m_vh),
      m_uNm1(m_vh),
      m_uHist(m_vh),
      m_uAdv(m_vh),
      m_pN(m_sh),
      m_pNm1(m_sh),
      m_pIt(m_sh),
      m_psiD(m_tauh),
      m_psiO(m_tauh),
      m_chiD(m_tauh),
      m_chiO(m_tauh),
      m_psiOldD(m_tauh),
      m_psiOldO(m_tauh),
      m_psiNm1D(m_tauh),
      m_psiNm1O(m_tauh),
      m_psiHistD(m_tauh),
      m_psiHistO(m_tauh),
      m_psiItD(m_tauh),
      m_psiItO(m_tauh),
      m_sigmaD(m_tauh),
      m_sigmaO(m_tauh),
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
      m_piEpsD(m_tauh, "al_proj_"),
      m_piEpsO(m_tauh, "al_proj_"),
      m_piAdvPsiD(m_tauh, "al_proj_"),
      m_piAdvPsiO(m_tauh, "al_proj_"),
      m_wss(m_vh),
      m_symRec0(m_vh),
      m_symRec1(m_vh),
      m_symRec2(m_vh),
      m_netShear(m_vh),
      m_absShear(m_sh),
      m_shearMagnitude(m_sh),
      m_tawss(m_sh),
      m_osi(m_sh),
      m_tawssLesion(m_sh),
      m_osiLesion(m_sh),
      m_qFlux(m_sh),
      m_one(m_sh),
      m_flux(m_qFlux),
      m_flow(m_u, m_p, m_psiD, m_psiO, m_v, m_q, m_chiD, m_chiO),
      m_flowKSP(m_flow),
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
    const auto velocityDOFs = m_vh.getSize();
    const auto pressureDOFs = m_sh.getSize();
    const auto stressDOFs = 2 * m_tauh.getSize();
    if (isRoot())
      Alert::Info() << "[mesh] " << m_cfg.meshPath << "  cells=" << cellCount
                    << " vertices=" << vertexCount << " velocity DOFs=" << velocityDOFs
                    << " pressure DOFs=" << pressureDOFs << " psi DOFs=" << stressDOFs
                    << Alert::Raise;
  }

  ArterialLesion3DViscousLogImplicit::~ArterialLesion3DViscousLogImplicit() = default;

  bool ArterialLesion3DViscousLogImplicit::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  void ArterialLesion3DViscousLogImplicit::deriveParameters()
  {
    auto& ob = m_cfg.oldroydB;
    const Real D = m_cfg.diameter;
    const Real rho = m_cfg.rho;

    if (m_cfg.fluid != "sptt" && m_cfg.fluid != "cy" && m_cfg.fluid != "gsptt")
      throw std::runtime_error("fluid must be sptt, cy or gsptt.");
    if (m_cfg.fluid == "cy")
    {
      // Carreau-Yasuda as "gsptt" with a vanishing polymer (file header).
      ob.etaP = 1.0e-6 * m_cfg.carreauYasuda.etaInf;
      m_cfg.deborah = 1.0e-3;
      m_cfg.weissenberg = 0.0;
      ob.pttEpsilon = 0.0;
      ob.lambda0Factor = 1.0;
      ob.lambda0Min = 0.0;
    }
    if ((m_cfg.fluid == "gsptt" || m_cfg.inletDevelopedConformation) &&
        (ob.lambda0Factor != 1.0 || ob.lambda0Min > 0.0))
      throw std::runtime_error("gsptt and inlet_developed_conformation require lambda0 = lambda "
                               "(lambda0_factor 1, lambda0_min 0).");
    const Real eta0 = referenceViscosity();

    if (!(m_cfg.reynolds > 0.0))
      throw std::runtime_error("Re must be positive.");
    if (!(m_cfg.womersley > 0.0))
      throw std::runtime_error("Wo must be positive; steady inflow is amplitude 0.");
    // Re = rho Ubar D/eta_0 and Wo = (D/2) sqrt(omega rho/eta_0).
    m_meanVelocity = m_cfg.reynolds * eta0 / (rho * D);
    const Real omega = 4.0 * m_cfg.womersley * m_cfg.womersley * eta0 / (rho * D * D);
    m_period = 2.0 * M_PI / omega;
    if (m_cfg.timeStep > 0.0)
      m_cfg.stepsPerCycle = std::max(1, static_cast<int>(std::lround(m_period / m_cfg.timeStep)));
    if (!(m_cfg.stepsPerCycle > 0))
      throw std::runtime_error("stepsPerCycle must be positive.");
    m_dt = m_period / static_cast<Real>(m_cfg.stepsPerCycle);
    m_dtNominal = m_dt;
    if (m_cfg.timeOrder != 1 && m_cfg.timeOrder != 2)
      throw std::runtime_error("timeOrder must be 1 (BDF1) or 2 (BDF2).");

    // De = lambda/T wins over Wi = lambda Ubar/D, which wins over lambda.
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
    if (m_cfg.fluid == "sptt")
      m_inflow.setFluid(rho, ob.etaS, ob.etaP, ob.lambda, omega, 0.5 * D);
    else
    {
      // Mean mode: the developed steady Carreau-Yasuda profile (also that of
      // "gsptt"); oscillatory modes: Womersley at eta_inf, as in the
      // axisymmetric driver.
      m_inflow.setFluid(rho, m_cfg.carreauYasuda.etaInf, ob.etaP, ob.lambda, omega, 0.5 * D);
      computeDevelopedProfile();
      m_inflow.setMeanProfile(m_inletR, m_inletPhi);
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
                      << "  eta_p=" << ob.etaP << "  solvent eta_s(gd) at 1, 10, 100, 1000 1/s = "
                      << solventViscosity(1.0) << ", " << solventViscosity(10.0) << ", "
                      << solventViscosity(100.0) << ", " << solventViscosity(1000.0) << " Pa s"
                      << Alert::Raise;
      }
      Alert::Info() << "[inflow] D=" << D << " m  Ubar=" << m_meanVelocity
                    << " m/s  T=" << m_period << " s  lambda=" << ob.lambda
                    << " s  lambda0=" << lambda0() << " s  dt=" << m_dt << " s  BDF" << m_cfg.timeOrder << "  q_peak="
                    << m_inflow.getPeak() << "  modes="
                    << (m_cfg.flowWaveformPath.empty() ? std::string("sinusoid")
                                                       : m_cfg.flowWaveformPath)
                    << Alert::Raise;
      if (m_cfg.flowWaveformPath.empty() && std::abs(m_cfg.amplitude) > 1.0)
        Alert::Warning() << "[inflow] |A| > 1: the inlet reverses, and psi = 0 is then "
                            "imposed on an outflow boundary." << Alert::Raise;
    }
  }

  ArterialLesion3DViscousLogImplicit::MeshType ArterialLesion3DViscousLogImplicit::makeMesh(
    const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3)
        throw std::runtime_error(
          "ArterialLesion3D expects a tetrahedral MEDIT mesh written as "
          "\"Dimension 3\"; generate it with "
          "examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/"
          "make_lesion_mesh_3d.py.");

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);
      mesh.getConnectivity().compute(D - 1, 1);
      mesh.getConnectivity().compute(1, 0);

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
    // The mesh is written in units of D.
    mesh.scale(cfg.diameter);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.getConnectivity().compute(D - 1, 1);
    mesh.getConnectivity().compute(1, 0);
    // Faces, then edges, as in Heart/CoronaryArtery_FSI_PETSc_MPI.
    mesh.reconcile(2);
    mesh.reconcile(1);

    return mesh;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = solventViscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::referenceViscosity() const
  {
    if (m_cfg.etaRef > 0.0)
      return m_cfg.etaRef;
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS + m_cfg.oldroydB.etaP;
    return m_cfg.carreauYasuda.etaInf;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::carreauYasuda(Real gammaDot) const
  {
    const auto& cy = m_cfg.carreauYasuda;
    return cy.etaInf + (cy.etaZero - cy.etaInf) *
      std::pow(1.0 + std::pow(cy.lambda * gammaDot, cy.a), (cy.n - 1.0) / cy.a);
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::steadyPTTFactor(Real gammaDot) const
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

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::solventViscosity(Real gammaDot) const
  {
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS;
    const Real etaCY = carreauYasuda(gammaDot);
    if (m_cfg.fluid == "cy")
      return etaCY;
    return etaCY - m_cfg.oldroydB.etaP / steadyPTTFactor(gammaDot);
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::shearRateAt(const Point& p) const
  {
    // At u^n (BDF1) or u* (BDF2).
    const auto J = Jacobian(m_uAdv).getValue(p);
    Real s2 = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        const Real e = 0.5 * (J(i, j) + J(j, i));
        s2 += e * e;
      }
    return std::sqrt(2.0 * s2);
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::solventViscosityAt(const Point& p) const
  {
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS;
    return solventViscosity(shearRateAt(p));
  }

  void ArterialLesion3DViscousLogImplicit::computeDevelopedProfile()
  {
    const Real R = 0.5 * m_cfg.diameter;
    const Real etaInf = m_cfg.carreauYasuda.etaInf;
    const size_t n = 801;
    std::vector<Real> r(n), gd(n), u(n);
    for (size_t i = 0; i < n; ++i)
      r[i] = R * static_cast<Real>(i) / static_cast<Real>(n - 1);

    // eta_CY(gd) gd is increasing in gd, and eta_CY >= eta_inf bounds gd.
    const auto shearRate = [this, etaInf](Real tau) {
      if (tau <= 0.0)
        return Real(0);
      Real lo = 0.0, hi = tau / etaInf;
      for (int it = 0; it < 200; ++it)
      {
        const Real mid = 0.5 * (lo + hi);
        if (carreauYasuda(mid) * mid < tau) lo = mid; else hi = mid;
        if (hi - lo <= 1.0e-14 * hi) break;
      }
      return 0.5 * (lo + hi);
    };
    const auto meanVelocity = [&](Real G) {
      for (size_t i = 0; i < n; ++i)
        gd[i] = shearRate(0.5 * G * r[i]);
      u[n - 1] = 0.0;
      for (size_t i = n - 1; i-- > 0;)
        u[i] = u[i + 1] + 0.5 * (gd[i] + gd[i + 1]) * (r[i + 1] - r[i]);
      Real q = 0.0;   // int u r dr
      for (size_t i = 0; i + 1 < n; ++i)
        q += 0.5 * (u[i] * r[i] + u[i + 1] * r[i + 1]) * (r[i + 1] - r[i]);
      return 2.0 * q / (R * R);
    };

    // Bisection on log G: the Newtonian gradients at eta_0 and eta_inf bracket G.
    Real lo = std::log(8.0 * etaInf * m_meanVelocity / (R * R));
    Real hi = std::log(8.0 * m_cfg.carreauYasuda.etaZero * m_meanVelocity / (R * R)) + 1.0;
    for (int it = 0; it < 200; ++it)
    {
      const Real mid = 0.5 * (lo + hi);
      if (meanVelocity(std::exp(mid)) < m_meanVelocity) lo = mid; else hi = mid;
      if (hi - lo < 1.0e-13) break;
    }
    const Real G = std::exp(0.5 * (lo + hi));
    const Real U = meanVelocity(G);

    m_inletR = r;
    m_inletPhi.resize(n);
    m_inletShear = gd;
    for (size_t i = 0; i < n; ++i)
      m_inletPhi[i] = u[i] / U;   // unit mean, exactly on this grid

    if (isRoot())
      Alert::Info() << "[inflow] developed Carreau-Yasuda profile: dp/dx=" << G
                    << " Pa/m  u_max/Ubar=" << m_inletPhi[0] << "  wall shear rate="
                    << gd[n - 1] << " 1/s  tau_w=" << 0.5 * G * R << " Pa" << Alert::Raise;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::interpolateInlet(const std::vector<Real>& table, Real r) const
  {
    r = std::clamp<Real>(std::abs(r), m_inletR.front(), m_inletR.back());
    const auto it = std::upper_bound(m_inletR.begin(), m_inletR.end(), r);
    const size_t j = std::min<size_t>(
      std::max<std::ptrdiff_t>(it - m_inletR.begin(), 1), m_inletR.size() - 1);
    const Real w = (r - m_inletR[j - 1]) / (m_inletR[j] - m_inletR[j - 1]);
    return (1.0 - w) * table[j - 1] + w * table[j];
  }

  Math::SpatialVector<ArterialLesion3DViscousLogImplicit::Real> ArterialLesion3DViscousLogImplicit::inletLogConformation(
    const Point& p, bool offDiagonal) const
  {
    Math::SpatialVector<Real> out{{0.0, 0.0, 0.0}};
    if (!m_cfg.inletDevelopedConformation)
      return out;
    const auto& x = p.getPhysicalCoordinates();
    const Real r = std::sqrt(x(1) * x(1) + x(2) * x(2));
    if (r <= 1.0e-14 * m_cfg.diameter)
      return out;
    const Real D = m_cfg.diameter;
    // du_x/dr < 0 of the developed mean profile.
    const Real dudr = m_inletShear.empty() ? -16.0 * m_meanVelocity * r / (D * D)
                                           : -interpolateInlet(m_inletShear, r);
    // Steady simple shear of the sPTT mode (lambda0 = lambda), in (x, r):
    // c = [[1 + 2 a^2, a], [a, 1]], a = lambda du_x/dr / f; c_tt = 1.
    const Real a = m_cfg.oldroydB.lambda * dudr / steadyPTTFactor(std::abs(dudr));
    Eigen::Matrix2d cm;
    cm << 1.0 + 2.0 * a * a, a, a, 1.0;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(cm);
    const Eigen::Matrix2d lg = eig.eigenvectors() *
      eig.eigenvalues().array().log().matrix().asDiagonal() * eig.eigenvectors().transpose();
    const Real ny = x(1) / r, nz = x(2) / r;
    // log c = Lxx e_x e_x + Lxr (e_x e_r + e_r e_x) + Lrr e_r e_r (L_tt = 0).
    if (offDiagonal)
      out = Math::SpatialVector<Real>{{ lg(0, 1) * ny, lg(0, 1) * nz, lg(1, 1) * ny * nz }};
    else
      out = Math::SpatialVector<Real>{{ lg(0, 0), lg(1, 1) * ny * ny, lg(1, 1) * nz * nz }};
    return out;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::pttFactor(Real traceExp) const
  {
    const auto& ob = m_cfg.oldroydB;
    if (ob.pttEpsilon == 0.0)
      return 1.0;
    return 1.0 + ob.pttEpsilon * (ob.lambda / lambda0()) * (traceExp - 3.0);
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::alpha3At(const Point& p) const
  {
    const auto& ob = m_cfg.oldroydB;
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real speed = std::sqrt(Math::dot(uc, uc));
    const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
    const Real f = pttFactor(
      logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp().trace());
    return 1.0 / (4.0 * f / (2.0 * ob.etaP) +
      0.25 * (ob.lambda * speed / (2.0 * ob.etaP * h) + ob.lambda * gradNorm / ob.etaP));
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::lambda0() const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = std::max(ob.lambda0Factor * ob.lambda, ob.lambda0Min);
    if (!(l0 > 0.0))
      throw std::runtime_error("lambda0 = max(k lambda, lambda0_min) must be positive.");
    return l0;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::ramp(Real t) const
  {
    const Real tr = m_cfg.rampCycles * m_period;
    if (!(tr > 0.0) || t >= tr)
      return 1.0;
    return 0.5 * (1.0 - std::cos(M_PI * t / tr));
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::inflowFactor(Real t) const
  {
    return ramp(t) * m_inflow.flowRate(t);
  }

  void ArterialLesion3DViscousLogImplicit::updateConformation()
  {
    const Real s = m_cfg.oldroydB.etaP / lambda0();
    const auto field = [](auto f) { return VectorFunction(size_t(3), f); };
    const auto sigmaAt = [this, s](const Point& p) -> Eigen::Matrix3d {
      return s * (logEig(unpack(m_psiItD.getValue(p), m_psiItO.getValue(p))).exp() -
        Eigen::Matrix3d::Identity());
    };

    m_sigmaD.project(field([&](const Point& p) -> Math::SpatialVector<Real> {
      return packD(sigmaAt(p)); }));
    m_sigmaO.project(field([&](const Point& p) -> Math::SpatialVector<Real> {
      return packO(sigmaAt(p)); }));
  }

  void ArterialLesion3DViscousLogImplicit::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  void ArterialLesion3DViscousLogImplicit::combine(
    Real a, const ::Vec& x, Real b, const ::Vec& y, ::Vec& z)
  {
    // z = b y, then z += a x; z may alias neither x nor y.
    PetscErrorCode ierr = VecCopy(y, z);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecScale(z, b);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecAXPY(z, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(z, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(z, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::bdfRatio() const
  {
    // The first step takes x^{-1} = x^0, for which w does not enter the
    // history; w = 1 keeps its mass step at 2 dt/3.
    if (m_cfg.timeOrder != 2 || m_history < 1 || !(m_dtPrev > 0.0))
      return 1.0;
    return m_dt / m_dtPrev;
  }

  ArterialLesion3DViscousLogImplicit::Real ArterialLesion3DViscousLogImplicit::massStep() const
  {
    if (m_cfg.timeOrder != 2)
      return m_dt;
    // (1 + 2 w)/(1 + w) x^{n+1} - (1 + w) x^n + w^2/(1 + w) x^{n-1}, over dt,
    // is (x^{n+1} - x_hist)/dtM with dtM = dt (1 + w)/(1 + 2 w). At w = 1 the
    // expression is evaluated as 2 dt/3, as before the variable step.
    const Real w = bdfRatio();
    return (1.0 + w) * m_dt / (1.0 + 2.0 * w);
  }

  void ArterialLesion3DViscousLogImplicit::prepareStep()
  {
    if (m_cfg.timeOrder == 2)
    {
      // x_hist = ((1 + w)^2 x^n - w^2 x^{n-1})/(1 + 2 w), u* = (1 + w) u^n -
      // w u^{n-1}; at w = 1, (4 x^n - x^{n-1})/3 and 2 u^n - u^{n-1}.
      const Real w = bdfRatio();
      const Real a = (1.0 + w) * (1.0 + w) / (1.0 + 2.0 * w);
      const Real b = -w * w / (1.0 + 2.0 * w);
      combine(a, m_uOld.getData(), b, m_uNm1.getData(), m_uHist.getData());
      combine(1.0 + w, m_uOld.getData(), -w, m_uNm1.getData(), m_uAdv.getData());
      combine(a, m_psiOldD.getData(), b, m_psiNm1D.getData(), m_psiHistD.getData());
      combine(a, m_psiOldO.getData(), b, m_psiNm1O.getData(), m_psiHistO.getData());
    }
    else
    {
      m_uHist.setData(m_uOld.getData());
      m_uAdv.setData(m_uOld.getData());
      m_psiHistD.setData(m_psiOldD.getData());
      m_psiHistO.setData(m_psiOldO.getData());
    }
  }

  void ArterialLesion3DViscousLogImplicit::advanceHistory()
  {
    if (m_cfg.timeOrder == 2)
    {
      // x^{n-1} <- x^n, x^n <- x^{n+1}; the history of the next step is
      // prepareStep's, once its dt is known.
      m_uNm1.setData(m_uOld.getData());
      m_uOld.setData(m_uIt.getData());
      m_pNm1.setData(m_pN.getData());
      m_pN.setData(m_p.getSolution().getData());
      m_psiNm1D.setData(m_psiOldD.getData());
      m_psiNm1O.setData(m_psiOldO.getData());
      m_psiOldD.setData(m_psiItD.getData());
      m_psiOldO.setData(m_psiItO.getData());
    }
    else
    {
      m_uOld.setData(m_uIt.getData());
      m_pN.setData(m_p.getSolution().getData());
      m_psiOldD.setData(m_psiItD.getData());
      m_psiOldO.setData(m_psiItO.getData());
    }
  }

  void ArterialLesion3DViscousLogImplicit::rejectStep()
  {
    // x^n, x^{n-1} and the history are untouched by a failed attempt; the
    // trial solutions (read as p^n by the BDF1 projections, and as the Krylov
    // initial guess) and the iterates go back to x^n.
    m_u.getSolution().setData(m_uOld.getData());
    m_p.getSolution().setData(m_pN.getData());
    m_psiD.getSolution().setData(m_psiOldD.getData());
    m_psiO.getSolution().setData(m_psiOldO.getData());
    m_uIt.setData(m_uOld.getData());
    m_psiItD.setData(m_psiOldD.getData());
    m_psiItO.setData(m_psiOldO.getData());
    ++m_psiRevision;
    // The matrix of the retry differs (dt): rebuild the preconditioner.
    m_pcStale = true;
  }

  ArterialLesion3DViscousLogImplicit& ArterialLesion3DViscousLogImplicit::initialize()
  {
    setupSpaces();
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

  void ArterialLesion3DViscousLogImplicit::setupSpaces()
  {
    const auto zero = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_uOld = zero;
    m_uIt = zero;
    m_uNm1 = zero;
    m_uHist = zero;
    m_uAdv = zero;
    m_pN = Real(0);
    m_pNm1 = Real(0);
    m_pIt = Real(0);
    m_subOld = zero;
    m_sub.get() = zero;
    m_psiOldD = zero;
    m_psiOldO = zero;
    m_psiNm1D = zero;
    m_psiNm1O = zero;
    m_psiHistD = zero;
    m_psiHistO = zero;
    m_psiItD = zero;
    m_psiItO = zero;
    ++m_psiRevision;
    m_wss = zero;
    m_netShear = zero;
    m_symRec0 = zero;
    m_symRec1 = zero;
    m_symRec2 = zero;

    m_absShear = Real(0);
    m_shearMagnitude = Real(0);
    m_tawss = Real(0);
    m_osi = Real(0);
    m_tawssLesion = Real(0);
    m_osiLesion = Real(0);
    m_one = Real(1);

    configureMassSolver(m_vectorProjectionKSP, "al_vproj_");
    configureMassSolver(m_wssKSP, "al_wss_");

    updateConformation();  // psi = 0: tau = I, sigma = 0

    m_u.setName("velocity");
    m_p.setName("pressure");
    m_psiD.setName("logConformationDiag");
    m_psiO.setName("logConformationOff");
    m_sigmaD.setName("elasticStressDiag");
    m_sigmaO.setName("elasticStressOff");
    m_tawss.setName("TAWSS");
    m_osi.setName("OSI");
    m_wss.setName("shearStress");

    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    // (xx, yy, zz) and (xy, xz, yz).
    m_xdmf.add("logConformationDiag", m_psiD.getSolution());
    m_xdmf.add("logConformationOff", m_psiO.getSolution());
    m_xdmf.add("elasticStressDiag", m_sigmaD);
    m_xdmf.add("elasticStressOff", m_sigmaO);
    m_xdmf.add("TAWSS", m_tawss);
    m_xdmf.add("OSI", m_osi);
    m_xdmf.add("shearStress", m_wss);

    // In 3D these are areas. The inlet must measure pi D^2/4, up to the
    // polygonal approximation of the circle.
    m_inletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.inlet});
    m_outletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.outlet});
    m_flux = Integral(m_one, m_qFlux);
    m_flux.assemble();
    m_volume = std::max<Real>(m_flux(m_one), 1e-30);

    if (isRoot())
    {
      const Real area = 0.25 * M_PI * m_cfg.diameter * m_cfg.diameter;
      Alert::Info() << "[boundary] inlet area = " << m_inletMeasure
                    << " m^2  outlet area = " << m_outletMeasure << " m^2  (pi D^2/4 = "
                    << area << " m^2)" << Alert::Raise;
      if (std::abs(m_inletMeasure / area - 1.0) > 2.0e-2)
        Alert::Warning() << "[boundary] the inlet is not the disk of diameter D: the "
                            "mesh is not the pipe in units of D, and the imposed flow "
                            "rate is not Ubar pi D^2/4 q(t)." << Alert::Raise;
    }
  }

  ArterialLesion3DViscousLogImplicit::Real
  ArterialLesion3DViscousLogImplicit::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-18);
  }

  void ArterialLesion3DViscousLogImplicit::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    // BDF mass terms (x - x_hist)/dtM: dtM = dt (BDF1) or 2 dt/3 (BDF2).
    const Real dtM = massStep();
    m_formDt = dt;
    m_formMassStep = dtM;
    const auto& ob = m_cfg.oldroydB;
    // eta_s(gd(u^n or u*)); constant (order 0) for "sptt", so that assembly is
    // the parent one.
    const auto etaS = PointwiseScalar([this](const Point& p) -> Real {
      return solventViscosityAt(p); }, m_cfg.fluid == "sptt" ? size_t(0) : size_t(2));
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;                  // sigma = s (exp(psi) - I)
    /// The psi equation is 2x the sigma equation of the exp form (2D header).
    const Real chiScale = 2.0;
    const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };
    constexpr Block D = Block::D;
    constexpr Block O = Block::O;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    // Galerkin convection with u^n (BDF1) or u* (BDF2).
    const auto convU = Mult(Jacobian(m_u), m_uAdv);
    const auto temam = Dot(Div(m_uAdv) * m_u, m_v);

    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);

    // Time-dependent data enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pOutFn = [this](const Point&) { return m_cfg.outletPressure; };
    VectorFunction inletVelocity(dim, [this](const Point& p) -> Math::SpatialVector<Real> {
      const auto& x = p.getPhysicalCoordinates();
      const Real r = std::sqrt(x(1) * x(1) + x(2) * x(2));
      const Real R = 0.5 * m_cfg.diameter;
      const Real seed = (m_cfg.inletPerturbation != 0.0 &&
          m_t < m_cfg.perturbationCycles * m_period) ? 1.0 + m_cfg.inletPerturbation * x(1) / R : 1.0;
      return Math::SpatialVector<Real>{{
        m_meanVelocity * ramp(m_t) * m_inflow.velocity(r, m_t) * seed, 0.0, 0.0 }};
    });

    const auto streamV = Mult(Jacobian(m_v), m_uOld);  // (grad v) u^n

    // (A u^n) . (B u^n) = (A u^n (u^n)^T) : B. Bilinear streamline terms are
    // written this way so that u^n sits on the trial side.
    const auto uu = PointwiseTensor([this](const Point& p) -> Eigen::Matrix3d {
      const auto u = m_uOld.getValue(p);
      const Eigen::Vector3d w(u(0), u(1), u(2));
      return w * w.transpose();
    });
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(3, 3));

    // Everything pointwise is at (psi^k, u^k), the Newton linearisation point.
    const LogFields<VectorGridFunctionType> at{ &m_psiItD, &m_psiItO, &m_uIt, model, &m_psiRevision };

    // Momentum: exp(psi) ~ exp(psi^k) + Dexp[psi - psi^k], all at psi^k.
    const auto Tn = fromPoint(at, LogDexp, [](const LogPoint& q) { return q.exp; });
    const auto C = Tn - dexp<D>(at, m_psiItD) - dexp<O>(at, m_psiItO);

    // psi equation, Newton about (psi^k, u^k); see the 2D driver:
    //   psi/dt + JN[psi] + u^k.grad psi + u.grad psi^k + G[u]
    //     = psi^n/dt - N + JN[psi^k] + G[u^k] + u^k.grad psi^k,
    // with (dt, psi^n) -> (dtM, psi history) for BDF2.
    // G[v] = sum_ab L_ab(v) Q_ab, tested against chi_D and chi_O. Mass and
    // advection are block-diagonal, tensor<B>(a) : tensor<B>(b) = w_B a.b.
    const Real wD = voigtWeight(D), wO = voigtWeight(O);
    const auto known = [&]<Block B>(const auto& psiOld, const auto& psiIt) {
      return pairedResidual<B>(at) - (voigtWeight(B) / dtM) * psiOld
        - voigtWeight(B) * Mult(Jacobian(psiIt), m_uIt);
    };

    // One quadrature rule for every integral that reads the pointwise log
    // data. Newton needs JN[psi^k], G[u^k] and Dexp[psi^k] inside the known
    // data to cancel exactly against their trial counterparts at the fixed
    // point, which requires the same rule on both sides; the pointwise data
    // has no polynomial degree to infer an order from (a callable reports
    // none, which gives the known terms a one-point rule). One rule also lets
    // LogStore compute each point once. 4 integrates the order-2 pointwise
    // data against P1 x P1.
    const size_t psiOrder = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));
    const auto withOrder = [psiOrder](auto integral) {
      integral.setOrder(psiOrder);
      return integral;
    };

    const auto psiDD = withOrder(Integral(pairedJacobianN<D, D>(at, wD / dtM, m_psiD)
      + wD * Mult(Jacobian(m_psiD), m_uIt), m_chiD));
    const auto psiOO = withOrder(Integral(pairedJacobianN<O, O>(at, wO / dtM, m_psiO)
      + wO * Mult(Jacobian(m_psiO), m_uIt), m_chiO));
    const auto psiOD = withOrder(Integral(pairedJacobianN<D, O>(at, 0.0, m_psiO), m_chiD));
    const auto psiDO = withOrder(Integral(pairedJacobianN<O, D>(at, 0.0, m_psiD), m_chiO));
    const auto uD = withOrder(Integral(
      wD * Mult(Jacobian(m_psiItD), m_u) + pairedVelocityMap<D>(at, m_u), m_chiD));
    const auto uO = withOrder(Integral(
      wO * Mult(Jacobian(m_psiItO), m_u) + pairedVelocityMap<O>(at, m_u), m_chiO));
    const auto knownD = withOrder(Integral(known.template operator()<D>(m_psiHistD, m_psiItD), m_chiD));
    const auto knownO = withOrder(Integral(known.template operator()<O>(m_psiHistO, m_psiItO), m_chiO));

    const auto body =
      // ---- Momentum, sigma = s (exp(psi) - I) --------------------------------
        (rho / dtM) * Integral(m_u, m_v) - (rho / dtM) * Integral(m_uHist, m_v)
      + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)
      + 2.0 * Integral(etaS * symU, symV)
      + s * withOrder(Integral(dexp<D>(at, m_psiD), symV))
      + s * withOrder(Integral(dexp<O>(at, m_psiO), symV))
      + s * withOrder(Integral(C, symV)) - s * Integral(I, symV)
      - Integral(m_p, Div(m_v))

      // ---- Continuity -----------------------------------------------------
      + Integral(Div(m_u), m_q) + m_cfg.pressurePenalty * Integral(m_p, m_q)

      // ---- Constitutive, psi form, BDF1/BDF2, Newton ------------------------
      + psiDD + psiOO + psiOD + psiDO   // psi_B : chi_C
      + uD + uO                         // u : chi_C
      + knownD + knownO                 // known

      // ---- S1, momentum ---------------------------------------------------
      + rho * rho * Integral(m_tauK * Mult(Jacobian(m_u), uu), Jacobian(m_v))
      - rho * rho * Integral(m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get()), streamV)
      + m_cfg.pressureScale * Integral(m_alpha1 * Grad(m_p), Grad(m_q))
      - m_cfg.pressureScale * Integral(m_alpha1 * m_piGradP.get(), Grad(m_q))
      + chiScale * m_cfg.stressDivScale * s * withOrder(Integral(m_alpha1 * divExpTrial<D>(at, m_psiD), divergence<D>(m_chiD)))
      + chiScale * m_cfg.stressDivScale * s * withOrder(Integral(m_alpha1 * divExpTrial<O>(at, m_psiO), divergence<D>(m_chiD)))
      + chiScale * m_cfg.stressDivScale * s * withOrder(Integral(m_alpha1 * divExpTrial<D>(at, m_psiD), divergence<O>(m_chiO)))
      + chiScale * m_cfg.stressDivScale * s * withOrder(Integral(m_alpha1 * divExpTrial<O>(at, m_psiO), divergence<O>(m_chiO)))
      - chiScale * m_cfg.stressDivScale * Integral(m_alpha1 * m_piDivSigma.get(), divergence<D>(m_chiD))
      - chiScale * m_cfg.stressDivScale * Integral(m_alpha1 * m_piDivSigma.get(), divergence<O>(m_chiO))

      // ---- S2, continuity -------------------------------------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v))
      - Integral(m_alpha2 * m_piDiv.get(), Div(m_v))

      // ---- S3: u-psi compatibility, and the psi advection (block-diagonal) ---
      + Integral(m_alpha3 * symU, symV)
      - Integral(m_alpha3 * (tensor<D>(m_piEpsD.get()) + tensor<O>(m_piEpsO.get())), symV)
      + wD * Integral(m_alphaPsi * Mult(Jacobian(m_psiD), uu), Jacobian(m_chiD))
      + wO * Integral(m_alphaPsi * Mult(Jacobian(m_psiO), uu), Jacobian(m_chiO))
      - wD * Integral(m_alphaPsi * m_piAdvPsiD.get(), Mult(Jacobian(m_chiD), m_uOld))
      - wO * Integral(m_alphaPsi * m_piAdvPsiO.get(), Mult(Jacobian(m_chiO), m_uOld))

      // ---- Outlet: traction p_out n and directional do-nothing ---------------
      + BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)
      + BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Inlet: developed pulsatile profile; walls: no slip -----------------
      + DirichletBC(m_u, inletVelocity).on(m_cfg.labels.inlet)
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // The psi equation is hyperbolic in psi: relaxed fluid enters at the inlet.
    // Relaxed fluid (psi = 0), or the steady simple-shear state of the
    // developed mean profile (inletDevelopedConformation).
    VectorFunction inletPsiD(size_t(3), [this](const Point& p) -> Math::SpatialVector<Real> {
      return inletLogConformation(p, false); });
    VectorFunction inletPsiO(size_t(3), [this](const Point& p) -> Math::SpatialVector<Real> {
      return inletLogConformation(p, true); });
    const auto finish = [&](const auto& form) {
      if (m_cfg.inletConformation)
        m_flow = form
          + DirichletBC(m_psiD, inletPsiD).on(m_cfg.labels.inlet)
          + DirichletBC(m_psiO, inletPsiO).on(m_cfg.labels.inlet);
      else
        m_flow = form;
    };
    // Optional psi terms. Each changes the type of the form, so they are
    // appended through generic lambdas, in the order outlet, skew, finish.
    //  - Outlet (Config::outletConformationInflow): where u*.n < 0 the
    //    upwind term |u*.n|_- (psi - psi^n) : chi imposes weakly that the
    //    entering fluid carries the conformation of step n there; it vanishes
    //    identically while the outlet is an outflow.
    //  - Skew (Config::psiSkewAdvection): + (1/2)(div u*) psi : chi.
    const auto inflowOut = Max(-Dot(m_uAdv, normal), 0.0);
    const auto outlet = m_cfg.labels.outlet;
    const auto withSkew = [&](const auto& form) {
      if (m_cfg.psiSkewAdvection)
        finish(form
          + 0.5 * wD * Integral(Dot(Div(m_uAdv) * m_psiD, m_chiD))
          + 0.5 * wO * Integral(Dot(Div(m_uAdv) * m_psiO, m_chiO)));
      else
        finish(form);
    };
    const auto withOutlet = [&](const auto& form) {
      if (m_cfg.outletConformationInflow)
        withSkew(form
          + wD * BoundaryIntegral(inflowOut * Dot(m_psiD, m_chiD)).over(outlet)
          - wD * BoundaryIntegral(inflowOut * Dot(m_psiOldD, m_chiD)).over(outlet)
          + wO * BoundaryIntegral(inflowOut * Dot(m_psiO, m_chiO)).over(outlet)
          - wO * BoundaryIntegral(inflowOut * Dot(m_psiOldO, m_chiO)).over(outlet));
      else
        withSkew(form);
    };
    if (m_cfg.newtonD2)
    {
      // S1, the D2exp part of the Newton linearisation of div exp(psi).
      const Real k1 = chiScale * m_cfg.stressDivScale * s;
      withOutlet(body
        + k1 * withOrder(Integral(m_alpha1 * divExpTrialD2<D>(at, m_psiD), divergence<D>(m_chiD)))
        + k1 * withOrder(Integral(m_alpha1 * divExpTrialD2<O>(at, m_psiO), divergence<D>(m_chiD)))
        + k1 * withOrder(Integral(m_alpha1 * divExpTrialD2<D>(at, m_psiD), divergence<O>(m_chiO)))
        + k1 * withOrder(Integral(m_alpha1 * divExpTrialD2<O>(at, m_psiO), divergence<O>(m_chiO)))
        - k1 * withOrder(Integral(m_alpha1 * divExpKnownD2(at), divergence<D>(m_chiD)))
        - k1 * withOrder(Integral(m_alpha1 * divExpKnownD2(at), divergence<O>(m_chiO))));
    }
    else
      withOutlet(body);
  }

  void ArterialLesion3DViscousLogImplicit::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    // eta_s at the recovered rate of strain: gd^2 = 2 eps:eps = |2 eps|^2/2,
    // the rows of 2 eps being m_symRec0..2.
    RealFunction etaS([this](const Point& p) -> Real {
      if (m_cfg.fluid == "sptt")
        return m_cfg.oldroydB.etaS;
      const auto a = m_symRec0.getValue(p), b = m_symRec1.getValue(p), c = m_symRec2.getValue(p);
      const Real n2 = Math::dot(a, a) + Math::dot(b, b) + Math::dot(c, c);
      return solventViscosity(std::sqrt(0.5 * n2));
    });
    const auto nx = Component(normal, 0);
    const auto ny = Component(normal, 1);
    const auto nz = Component(normal, 2);
    // sigma_xx, yy, zz = D_0, D_1, D_2; sigma_xy, xz, yz = O_0, O_1, O_2.
    const auto& sigmaD = m_sigmaD;
    const auto& sigmaO = m_sigmaO;
    const auto sD = [&sigmaD](size_t c) { return Component(sigmaD, c); };
    const auto sO = [&sigmaO](size_t c) { return Component(sigmaO, c); };

    // t = (2 eta_s eps(u) + sigma) n.
    const auto traction = VectorFunction(
      etaS * Dot(m_symRec0, normal) + sD(0) * nx + sO(0) * ny + sO(1) * nz,
      etaS * Dot(m_symRec1, normal) + sO(0) * nx + sD(1) * ny + sO(2) * nz,
      etaS * Dot(m_symRec2, normal) + sO(1) * nx + sO(2) * ny + sD(2) * nz);
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  void ArterialLesion3DViscousLogImplicit::updateStabilization()
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
    // BDF2 projects in solveFlow, from the extrapolation (Config::bdf2Extrapolate).
    if (m_cfg.timeOrder == 1)
      projectConsistent(false, 0u);
  }

  void ArterialLesion3DViscousLogImplicit::projectConsistent(bool iterate, unsigned mask)
  {
    // Source of each projection: step n, or the iterate (the extrapolation
    // (u*, p*, psi*) at the first Newton iterate of a BDF2 step).
    const auto pick = [iterate, mask](unsigned bit) { return iterate && (mask & bit); };
    const VectorGridFunctionType& uDiv = pick(1u) ? m_uIt : m_uOld;
    const ScalarGridFunctionType& p = pick(2u) ? m_pIt : m_p.getSolution();
    const bool sigmaIt = pick(4u);
    const VectorGridFunctionType& psiSD = sigmaIt ? m_psiItD : m_psiOldD;
    const VectorGridFunctionType& psiSO = sigmaIt ? m_psiItO : m_psiOldO;
    const VectorGridFunctionType& uS = sigmaIt ? m_uIt : m_uOld;
    const VectorGridFunctionType& uE = pick(8u) ? m_uIt : m_uOld;
    const VectorGridFunctionType& psiAD = pick(16u) ? m_psiItD : m_psiOldD;
    const VectorGridFunctionType& psiAO = pick(16u) ? m_psiItO : m_psiOldO;
    m_piDiv.project(Div(uDiv));
    m_piGradP.project(Grad(p));
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;
    const LogModel model{ ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda) };
    // At the iterate the pointwise data are those of the assembly (same key).
    const LogFields<VectorGridFunctionType> atX{ &psiSD, &psiSO, &uS, model, &m_psiRevision };
    m_piDivSigma.project(s * (divExp<Block::D>(atX, psiSD) + divExp<Block::O>(atX, psiSO)));
    // Nothing else reads psi^n pointwise this step: free its points.
    if (!sigmaIt)
      releaseLogPoints(&m_psiOldD);
    const auto epsN = 0.5 * (Jacobian(uE) + Transpose(Jacobian(uE)));
    m_piEpsD.project(voigtD(epsN));
    m_piEpsO.project(voigtO(epsN));
    // The advecting velocity of the S3 psi term is u^n on both sides.
    m_piAdvPsiD.project(Mult(Jacobian(psiAD), m_uOld));
    m_piAdvPsiO.project(Mult(Jacobian(psiAO), m_uOld));
  }

  bool ArterialLesion3DViscousLogImplicit::solveFlow()
  {
    const bool trace = isRoot() && (m_step < 3 || m_cfg.traceNewton);
    const auto phase = [trace](const char* what) {
      if (trace)
      {
        ThreeDInfo() << what << " ..." << Alert::Raise;
        std::cout.flush();
      }
    };

    m_failure.clear();
    m_newtonPredicted = false;

    const auto vmsStart = CoronaryClock::now();
    phase("VMS: updating the stabilisation (alphas and projections)");
    updateStabilization();
    m_timing.vms += secondsSince(vmsStart);

    if (trace)
    {
      ThreeDInfo() << "Assembling the flow system ("
                   << (m_vh.getSize() + m_sh.getSize() + 2 * m_tauh.getSize())
                   << " unknowns) ..." << Alert::Raise;
      std::cout.flush();
    }

    // Newton about (u^k, psi^k), from (u^n, psi^n), as in the 2D driver.
    ::KSPConvergedReason reason = KSP_CONVERGED_ITS;
    PetscErrorCode ierr = PETSC_SUCCESS;
    ::Vec incrementD = PETSC_NULLPTR;
    ierr = VecDuplicate(m_psiItD.getData(), &incrementD);
    assert(ierr == PETSC_SUCCESS);
    ::Vec incrementO = PETSC_NULLPTR;
    ierr = VecDuplicate(m_psiItO.getData(), &incrementO);
    assert(ierr == PETSC_SUCCESS);
    ::Vec uIncrement = PETSC_NULLPTR;
    ierr = VecDuplicate(m_uIt.getData(), &uIncrement);
    assert(ierr == PETSC_SUCCESS);
    // The Krylov initial guess of a solve with a reused preconditioner, kept
    // in case that solve has to be repeated with a fresh one.
    ::Vec guess = PETSC_NULLPTR;
    m_conformationIts = 0;
    m_newtonUndamped = 0;

    // Newton starts from x^n (BDF1) or from the extrapolation (1 + w) x^n - w
    // x^{n-1} (BDF2). The extrapolation of psi is skipped on the second step:
    // the start from rest with a non-zero inlet conformation makes psi^1 -
    // psi^0 a jump, and 2 psi^1 a poor guess (Newton diverged for the "gsptt"
    // fluid). The iterates are set explicitly, since a rejected attempt may
    // have left them anywhere.
    const Real w = bdfRatio();
    m_uIt.setData(m_uAdv.getData());
    if (m_cfg.timeOrder == 2)
      combine(1.0 + w, m_pN.getData(), -w, m_pNm1.getData(), m_pIt.getData());
    if (m_cfg.timeOrder == 2 && m_history >= 2)
    {
      combine(1.0 + w, m_psiOldD.getData(), -w, m_psiNm1D.getData(), m_psiItD.getData());
      combine(1.0 + w, m_psiOldO.getData(), -w, m_psiNm1O.getData(), m_psiItO.getData());
    }
    else
    {
      m_psiItD.setData(m_psiOldD.getData());
      m_psiItO.setData(m_psiOldO.getData());
    }
    ++m_psiRevision;

    // The quadratic prediction d_{k+1} ~ C d_k^2 uses the largest recent
    // estimate of C; it is not used on calibration steps, so that C keeps
    // being measured.
    Real contraction = -1.0;
    for (const Real c : m_contraction)
      contraction = std::max(contraction, c);
    const bool calibrate =
      m_cfg.newtonCalibrateEvery > 0 && m_history % m_cfg.newtonCalibrateEvery == 0;
    const bool predict = m_cfg.newtonPredict && !calibrate && contraction > 0.0;
    const auto recordContraction = [this, &contraction](Real previous, Real increment) {
      // A zero increment is a Krylov solve that stopped on its initial
      // residual: below the resolution of the linear solver, it says nothing
      // about C.
      if (!(previous > 0.0) || !(increment > 0.0))
        return;
      // The first step (start from rest, inlet jump) is not representative.
      if (m_history == 0)
        return;
      constexpr size_t ContractionSamples = 16;  ///< recent estimates of C kept
      const Real c = increment / (previous * previous);
      if (m_contraction.size() < ContractionSamples)
        m_contraction.push_back(c);
      else
        m_contraction[m_contractionNext] = c;
      m_contractionNext = (m_contractionNext + 1) % ContractionSamples;
      // A measurement of this step counts at once: away from the quadratic
      // regime (a large dt, a hard step) it is the one that matters.
      contraction = std::max(contraction, c);
    };

    bool converged = false;
    Real previous = -1.0;  // the last undamped increment, < 0 if none
    const int maxIterations = std::max(1, m_cfg.conformationIterations);
    for (int k = 0; k < maxIterations; ++k)
    {
      if (m_cfg.timeOrder == 2 && k == 0)
      {
        // The projections of the step, from the extrapolation (u*, p*, psi*).
        const auto projStart = CoronaryClock::now();
        projectConsistent(true, static_cast<unsigned>(m_cfg.bdf2Extrapolate));
        m_timing.vms += secondsSince(projStart);
      }

      // The Krylov initial guess is gathered from the trial solutions by
      // assemble(): the iterate (u^k, p^k, psi^k), not x^n. The pressure of
      // the BDF1 projections is read above, before this.
      m_u.getSolution().setData(m_uIt.getData());
      m_psiD.getSolution().setData(m_psiItD.getData());
      m_psiO.getSolution().setData(m_psiItO.getData());
      if (k == 0 && m_cfg.timeOrder == 2)
        m_p.getSolution().setData(m_pIt.getData());

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

      // Lagged preconditioner (Config::pcReuseSolves): reused while the Krylov
      // solver needs at most pcReuseGrowth times the iterations it needed
      // with a fresh one; a solve that fails with it is repeated, from the
      // same guess, with a fresh one.
      const auto solveStart = CoronaryClock::now();
      ::KSP ksp = m_flowKSP.getHandle();
      ::Vec x = m_flow.getLinearSystem().getSolution();
      // A lagged preconditioner needs a Krylov method: with preonly it would be
      // the solution operator of an earlier matrix.
      ::KSPType kspType = nullptr;
      ierr = KSPGetType(ksp, &kspType);
      assert(ierr == PETSC_SUCCESS);
      const bool krylov = kspType && std::strcmp(kspType, KSPPREONLY) != 0;
      bool reuse = krylov && m_cfg.pcReuseSolves > 1 && m_pcBuilt && !m_pcStale &&
        m_pcSolves < m_cfg.pcReuseSolves;
      if (reuse)
      {
        if (!guess)
        {
          ierr = VecDuplicate(x, &guess);
          assert(ierr == PETSC_SUCCESS);
        }
        ierr = VecCopy(x, guess);
        assert(ierr == PETSC_SUCCESS);
      }
      ierr = KSPSetReusePreconditioner(ksp, reuse ? PETSC_TRUE : PETSC_FALSE);
      assert(ierr == PETSC_SUCCESS);
      m_flow.solve(m_flowKSP);
      PetscInt its = 0;
      ierr = KSPGetConvergedReason(ksp, &reason);
      assert(ierr == PETSC_SUCCESS);
      ierr = KSPGetIterationNumber(ksp, &its);
      assert(ierr == PETSC_SUCCESS);
      if (reuse && reason <= 0)
      {
        m_timing.kspIterations += static_cast<int>(its);
        ierr = VecCopy(guess, x);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPSetReusePreconditioner(ksp, PETSC_FALSE);
        assert(ierr == PETSC_SUCCESS);
        m_flow.solve(m_flowKSP);
        ierr = KSPGetConvergedReason(ksp, &reason);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPGetIterationNumber(ksp, &its);
        assert(ierr == PETSC_SUCCESS);
        reuse = false;
      }
      if (!reuse)
      {
        m_pcBuilt = true;
        m_pcStale = false;
        m_pcSolves = 0;
        // The reference count: a first Newton solve, which starts from the
        // predictor; later ones start next to the solution and need few.
        m_pcReferenceIts = k == 0 ? static_cast<int>(its) : -1;
        ++m_timing.pcBuilds;
      }
      else if (m_pcReferenceIts < 0 && k == 0)
        m_pcReferenceIts = static_cast<int>(its);
      else if (m_pcReferenceIts >= 0 &&
               static_cast<Real>(its) > m_cfg.pcReuseGrowth * static_cast<Real>(m_pcReferenceIts))
        m_pcStale = true;  // this solve is kept; the next one rebuilds
      ++m_pcSolves;
      m_timing.kspIterations += static_cast<int>(its);
      m_timing.solve += secondsSince(solveStart);

      if (reason <= 0)
      {
        m_failure = "linear solver failed (KSPConvergedReason " +
          std::to_string(static_cast<int>(reason)) + ", " + std::to_string(its) + " iterations)";
        break;
      }

      // max |psi^{k+1} - psi^k| over both blocks, relative to max(1, |psi^k|).
      ierr = VecWAXPY(incrementD, -1.0, m_psiItD.getData(), m_psiD.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(incrementO, -1.0, m_psiItO.getData(), m_psiO.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      Real stepD = 0.0, stepO = 0.0, scaleD = 0.0, scaleO = 0.0;
      ierr = VecNorm(incrementD, NORM_INFINITY, &stepD);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(incrementO, NORM_INFINITY, &stepO);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_psiItD.getData(), NORM_INFINITY, &scaleD);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(m_psiItO.getData(), NORM_INFINITY, &scaleO);
      assert(ierr == PETSC_SUCCESS);
      // std::max would drop a NaN in its second argument.
      const Real stepSize = std::isfinite(stepD) && std::isfinite(stepO)
        ? std::max(stepD, stepO) : std::numeric_limits<Real>::quiet_NaN();
      const Real psiScale = std::max(scaleD, scaleO);
      const Real omega = (std::isfinite(stepSize) && stepSize > m_cfg.newtonMaxStep)
        ? m_cfg.newtonMaxStep / stepSize : 1.0;
      m_psiIncrement = stepSize / std::max<Real>(1.0, psiScale);
      axpy(omega, incrementD, m_psiItD.getData());
      axpy(omega, incrementO, m_psiItO.getData());
      ierr = VecWAXPY(uIncrement, -1.0, m_uIt.getData(), m_u.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      axpy(omega, uIncrement, m_uIt.getData());
      ++m_psiRevision;
      ++m_conformationIts;

      const bool undamped = omega == 1.0;
      if (undamped)
        ++m_newtonUndamped;
      if (undamped && m_psiIncrement >= m_cfg.conformationTolerance)
        recordContraction(previous, m_psiIncrement);
      const bool predicted = undamped && predict && contraction > 0.0 &&
        m_psiIncrement >= m_cfg.conformationTolerance &&
        contraction * m_psiIncrement * m_psiIncrement <
          m_cfg.newtonPredictSafety * m_cfg.conformationTolerance;
      if (trace)
        KSPInfo() << "Newton " << m_conformationIts << ": max|dpsi| = " << stepSize
                  << " (relative " << m_psiIncrement << ")" << (omega < 1.0 ? " (damped)" : "")
                  << "  ksp its " << its << (reuse ? " (lagged PC)" : "")
                  << (predicted ? "  accepted: C d^2 = " + std::to_string(contraction *
                        m_psiIncrement * m_psiIncrement) : std::string())
                  << Alert::Raise;
      if (!std::isfinite(m_psiIncrement))
      {
        m_failure = "non-finite Newton increment";
        break;
      }
      if (undamped && m_psiIncrement < m_cfg.conformationTolerance)
      {
        recordContraction(previous, m_psiIncrement);
        converged = true;
        break;
      }
      // A growing undamped increment: Newton is outside its basin, and the
      // remaining iterations only cost time (rejected with step control).
      if (m_cfg.adaptiveTimeStep && undamped && previous > 0.0 && m_psiIncrement > previous)
      {
        m_failure = "Newton increment grew (" + std::to_string(previous) + " -> " +
          std::to_string(m_psiIncrement) + ")";
        break;
      }
      if (predicted)
      {
        m_newtonPredicted = true;
        converged = true;
        break;
      }
      previous = undamped ? m_psiIncrement : -1.0;
    }
    ierr = VecDestroy(&incrementD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDestroy(&incrementO);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecDestroy(&uIncrement);
    assert(ierr == PETSC_SUCCESS);
    if (guess)
    {
      ierr = VecDestroy(&guess);
      assert(ierr == PETSC_SUCCESS);
    }
    (void)ierr;

    if (!converged && m_failure.empty())
    {
      if (m_cfg.adaptiveTimeStep)
        m_failure = "Newton did not converge in " + std::to_string(m_conformationIts) +
          " iterations (|dpsi| = " + std::to_string(m_psiIncrement) + ")";
      else if (isRoot() && m_psiIncrement >= 100.0 * m_cfg.conformationTolerance)
        // Without step control the step is accepted, as before.
        Alert::Warning() << "[log] Newton stopped after " << m_conformationIts
                         << " iterations with max|dpsi| = " << m_psiIncrement << Alert::Raise;
    }
    if (!m_failure.empty())
      return false;

    // The divergence guard, on u^{n+1} before anything is advanced.
    m_speed = std::max(std::abs(m_uIt.max()), std::abs(m_uIt.min()));
    const Real guard =
      m_cfg.maxVelocityFactor * m_meanVelocity * std::max<Real>(1.0, m_inflow.getPeak());
    if (!std::isfinite(m_speed) || m_speed > guard)
    {
      m_failure = "divergence guard: max|u| = " + std::to_string(m_speed) + " m/s = " +
        std::to_string(m_speed / m_meanVelocity) + " Ubar";
      return false;
    }

    advanceHistory();
    m_u.getSolution().setData(m_uIt.getData());
    m_psiD.getSolution().setData(m_psiItD.getData());
    m_psiO.getSolution().setData(m_psiItO.getData());
    ++m_psiRevision;
    updateConformation();
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.get().getData());

    m_stress = std::max({ std::abs(m_sigmaD.max()), std::abs(m_sigmaD.min()),
      std::abs(m_sigmaO.max()), std::abs(m_sigmaO.min()) });
    return true;
  }

  void ArterialLesion3DViscousLogImplicit::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    // Rows of 2 eps(u), recovered onto the nodes (see the 2D driver).
    const auto jac = Jacobian(uSol);
    const auto e01 = Component(jac, 0, 1) + Component(jac, 1, 0);
    const auto e02 = Component(jac, 0, 2) + Component(jac, 2, 0);
    const auto e12 = Component(jac, 1, 2) + Component(jac, 2, 1);
    projectVector(VectorFunction(2.0 * Component(jac, 0, 0), e01, e02), m_symRec0);
    projectVector(VectorFunction(e01, 2.0 * Component(jac, 1, 1), e12), m_symRec1);
    projectVector(VectorFunction(e02, e12, 2.0 * Component(jac, 2, 2)), m_symRec2);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

    // Keep only the wall: off it the regularised projection carries a
    // sign-alternating tail (2D driver, computeWallShear).
    const size_t dim = m_mesh.getSpaceDimension();
    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    const auto& wssSol = m_wssTrial.getSolution();
    m_wss = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
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

  void ArterialLesion3DViscousLogImplicit::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // TAWSS and OSI on the wall only (2D driver, closeCycle).
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

    // The same indices restricted to the lesion wall: the global maxima can
    // sit at the inlet, where the imposed profile meets the discrete flow.
    const auto onLesion = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && *a == m_cfg.labels.lesion;
    };
    m_tawssLesion = Real(0);
    m_osiLesion = Real(0);
    m_tawssLesion.project(Region::Boundary,
      RealFunction([this](const Point& p) -> Real { return m_tawss.getValue(p); }), onLesion);
    m_osiLesion.project(Region::Boundary,
      RealFunction([this](const Point& p) -> Real { return m_osi.getValue(p); }), onLesion);
  }

  void ArterialLesion3DViscousLogImplicit::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    // Volumetric flow rates (m^3/s). n is outward, so the inlet flux is
    // negated to read as the inflow.
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.inlet);
    m_flux.assemble();
    m_qIn = -m_flux(m_one);

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);

    // Mean pressures: their difference is the pressure drop across the lesion
    // and the two straight segments.
    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.labels.inlet);
    m_flux.assemble();
    m_inletPressure = m_flux(m_one) / m_inletMeasure;

    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_outletMeanPressure = m_flux(m_one) / m_outletMeasure;

    // Volume means of the transverse velocity (m = 1 content), over Ubar.
    const Real norm = m_volume * std::max<Real>(m_meanVelocity, 1e-30);
    m_flux = Integral(Component(uSol, 1), m_qFlux);
    m_flux.assemble();
    m_uyMean = m_flux(m_one) / norm;
    m_flux = Integral(Component(uSol, 2), m_qFlux);
    m_flux.assemble();
    m_uzMean = m_flux(m_one) / norm;
  }

  void ArterialLesion3DViscousLogImplicit::writeCSVHeader()
  {
    m_csv << "t,cycle,qTarget,qIn,qOut,pInletMean,pOutletMean,dp,maxU,maxSigma,"
          << "newtonIts,dpsi,maxTAWSS,maxOSI,maxTAWSSLesion,maxOSILesion,uyMean,uzMean,dt,retries\n";
  }

  void ArterialLesion3DViscousLogImplicit::writeCSVRow(int cycle)
  {
    // Reductions over the communicator, taken on all ranks first.
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();
    const Real tawssLesion = m_tawssLesion.max();
    const Real osiLesion = m_osiLesion.max();

    if (!isRoot())
      return;

    const Real qTarget =
      m_meanVelocity * 0.25 * M_PI * m_cfg.diameter * m_cfg.diameter * m_inflowNow;
    m_csv << m_t << ',' << cycle << ',' << qTarget << ',' << m_qIn << ',' << m_qOut << ','
          << m_inletPressure << ',' << m_outletMeanPressure << ','
          << (m_inletPressure - m_outletMeanPressure) << ',' << m_speed << ',' << m_stress
          << ',' << m_conformationIts << ',' << m_psiIncrement << ',' << tawss << ',' << osi
          << ',' << tawssLesion << ',' << osiLesion << ',' << m_uyMean << ',' << m_uzMean
          << ',' << m_dt << ',' << m_timing.retries << '\n';
    m_csv.flush();
  }

  int ArterialLesion3DViscousLogImplicit::run()
  {
    if (!m_initialized)
      initialize();

    const int stepsPerCycle = m_cfg.stepsPerCycle;
    const int nominalSteps = m_cfg.cycles * stepsPerCycle;
    // A step whose end lies within TimeRoundoff dt0 of a cycle boundary or an
    // output time reaches it: far below any step the control takes (dt0
    // dtShrink^maxStepRetries), far above the roundoff of t summed over a run.
    constexpr Real TimeRoundoff = 1.0e-6;
    const Real dtMin = m_dtNominal * std::pow(m_cfg.dtShrink, m_cfg.maxStepRetries);
    // XDMF output, as before at the nominal step: after the first step and
    // every outputEvery nominal steps from there; cycle boundaries always.
    const Real outputInterval =
      m_cfg.outputEvery > 0 ? static_cast<Real>(m_cfg.outputEvery) * m_dtNominal : 0.0;
    Real nextOutput = m_dtNominal;
    const Real endTime = static_cast<Real>(m_cfg.cycles) * m_period;

    if (isRoot())
    {
      Alert::Info() << "[run] " << m_cfg.cycles << " cycles of " << stepsPerCycle
                    << " steps (" << nominalSteps << " total), inflow ramped over "
                    << m_cfg.rampCycles << " cycle(s); XDMF every " << m_cfg.outputEvery
                    << " steps and at every cycle boundary" << Alert::Raise;
      if (m_cfg.adaptiveTimeStep)
        Alert::Info() << "[run] step control: retry with dt x " << m_cfg.dtShrink << " (at most "
                      << m_cfg.maxStepRetries << " times, dt >= " << dtMin << " s), growth x "
                      << m_cfg.dtGrow << " after " << m_cfg.dtGrowAfter << " steps of <= "
                      << m_cfg.newtonEasy << " Newton iterations, up to dt0 = " << m_dtNominal
                      << " s; Newton prediction " << (m_cfg.newtonPredict ? "on" : "off")
                      << "; lagged preconditioner over " << m_cfg.pcReuseSolves << " solves"
                      << Alert::Raise;
    }

    Real cycleElapsed = 0.0;
    const auto runStart = CoronaryClock::now();
    Real dtTarget = m_dtNominal;
    int easySteps = 0;
    int cycle = 0;
    int step = 0;

    while (cycle < m_cfg.cycles)
    {
      m_step = step;
      m_timing = Timing{};
      const auto stepStart = CoronaryClock::now();
      const Real tStart = m_t;
      const Real cycleEnd = static_cast<Real>(cycle + 1) * m_period;

      // The target step, clipped to end on the cycle boundary. Within two
      // targets of it the remainder is split in two equal steps, so that no
      // sliver of a step is left; a remainder that is a whole number of
      // targets up to roundoff keeps the target.
      Real dt = dtTarget;
      {
        const Real remaining = cycleEnd - tStart;
        const Real ratio = remaining / dtTarget;
        if (ratio <= 1.0 + TimeRoundoff)
          dt = ratio >= 1.0 - TimeRoundoff ? dtTarget : remaining;
        else if (ratio < 2.0 - TimeRoundoff)
          dt = 0.5 * remaining;
      }

      for (int retries = 0;; ++retries)
      {
        m_dt = dt;
        // setupFlow bakes dt and the BDF mass step into the forms.
        if (m_dt != m_formDt || massStep() != m_formMassStep)
        {
          setupFlow();
          m_pcStale = true;
        }
        m_t = tStart + m_dt;
        m_inflowNow = inflowFactor(m_t);
        prepareStep();
        if (solveFlow())
          break;

        rejectStep();
        const Real next = m_cfg.dtShrink * m_dt;
        if (!m_cfg.adaptiveTimeStep || retries >= m_cfg.maxStepRetries ||
            next < (1.0 - TimeRoundoff) * dtMin)
        {
          Alert::Exception() << "[flow] diverged at step " << (step + 1) << ", t = " << m_t
                             << " s, dt = " << m_dt << " s: " << m_failure << Alert::Raise;
          return 1;
        }
        if (isRoot())
          Alert::Warning() << "[step] " << (step + 1) << " rejected at t = " << m_t
                           << " s, dt = " << m_dt << " s: " << m_failure
                           << "; retrying with dt = " << next << " s" << Alert::Raise;
        ++m_timing.retries;
        dt = next;
      }
      m_dtPrev = m_dt;
      ++m_history;
      ++step;

      computeWallShear();
      cycleElapsed += m_dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      m_timing.fluxes = secondsSince(fluxStart);

      const bool endOfCycle = cycleEnd - m_t <= TimeRoundoff * m_dtNominal;
      if (endOfCycle)
      {
        closeCycle(cycleElapsed);
        cycleElapsed = 0.0;

        const Real maxTawss = m_tawss.max();
        const Real maxOsi = m_osi.max();
        const Real lesionTawss = m_tawssLesion.max();
        const Real lesionOsi = m_osiLesion.max();
        if (isRoot())
          Alert::Info() << "[cycle " << (cycle + 1) << "/" << m_cfg.cycles
                        << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                        << "  lesion: maxTAWSS=" << lesionTawss << " Pa  maxOSI=" << lesionOsi
                        << Alert::Raise;
      }

      writeCSVRow(cycle);

      bool output = endOfCycle;
      if (outputInterval > 0.0 && m_t >= nextOutput - TimeRoundoff * m_dtNominal)
      {
        output = true;
        while (nextOutput <= m_t + TimeRoundoff * m_dtNominal)
          nextOutput += outputInterval;
      }
      if (output)
      {
        const auto outputStart = CoronaryClock::now();
        m_xdmf.write(m_t).flush();
        m_timing.output = secondsSince(outputStart);
      }

      m_timing.total = secondsSince(stepStart);

      if (isRoot() && (m_step < 3 || m_step % 20 == 0 || endOfCycle || m_timing.retries > 0))
      {
        const Real elapsed = secondsSince(runStart);
        const Real remainingTime = std::max<Real>(endTime - m_t, 0.0);

        Alert::Info() << "---- Step " << step << "/" << nominalSteps << "  cycle "
                      << (cycle + 1) << "/" << m_cfg.cycles << "  t = " << m_t << " s"
                      << "  dt = " << m_dt << " s" << Alert::Raise;
        Alert::Info() << "[3D] q=" << m_inflowNow << "  qIn=" << m_qIn << "  qOut=" << m_qOut
                      << " m^3/s  max|u|=" << m_speed << " m/s  max|sigma|=" << m_stress
                      << " Pa  dp=" << (m_inletPressure - m_outletMeanPressure)
                      << " Pa  newton=" << m_conformationIts << (m_newtonPredicted ? "p" : "")
                      << " (|dpsi|=" << m_psiIncrement
                      << ")  <uy,uz>/Ubar=(" << m_uyMean << ", " << m_uzMean << ")" << Alert::Raise;
        Alert::Info() << "[timing] vms=" << m_timing.vms << "  asm=" << m_timing.assembly
                      << "  ksp=" << m_timing.solve << " (" << m_timing.kspIterations << " its, "
                      << m_timing.pcBuilds << " PC builds)  wss=" << m_timing.shear
                      << "  flux=" << m_timing.fluxes << "  out=" << m_timing.output
                      << "  retries=" << m_timing.retries
                      << "  total=" << m_timing.total << " s  |  ETA "
                      << (elapsed / std::max<Real>(m_t, 1e-30) * remainingTime / 60.0) << " min"
                      << Alert::Raise;
        std::cout.flush();
      }

      if (endOfCycle)
        ++cycle;

      // Step size of the next step (Config::adaptiveTimeStep).
      if (m_cfg.adaptiveTimeStep)
      {
        if (m_timing.retries > 0)
        {
          // Go on at the step that worked.
          dtTarget = m_dt;
          easySteps = 0;
        }
        else if (m_newtonUndamped >= m_cfg.newtonHard)
        {
          dtTarget = std::max(m_cfg.dtShrink * dtTarget, dtMin);
          easySteps = 0;
        }
        else if (m_newtonUndamped <= m_cfg.newtonEasy && dtTarget < m_dtNominal)
        {
          if (++easySteps >= m_cfg.dtGrowAfter)
          {
            dtTarget = std::min(m_cfg.dtGrow * dtTarget, m_dtNominal);
            easySteps = 0;
          }
        }
        else
          easySteps = 0;
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
      using Simulation = Rodin::Examples::ViscoelasticFluids::ArterialLesion3DViscousLogImplicit;
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
      getString("-al_fluid", cfg.fluid);

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
      getReal("-al_dt", cfg.timeStep);
      getInt("-al_cycles", cfg.cycles);
      getInt("-al_time_order", cfg.timeOrder);
      getReal("-al_eta_ref", cfg.etaRef);
      getReal("-al_cy_eta_zero", cfg.carreauYasuda.etaZero);
      getReal("-al_cy_eta_inf", cfg.carreauYasuda.etaInf);
      getReal("-al_cy_lambda", cfg.carreauYasuda.lambda);
      getReal("-al_cy_n", cfg.carreauYasuda.n);
      getReal("-al_cy_a", cfg.carreauYasuda.a);
      getBool("-al_inlet_developed_conformation", cfg.inletDevelopedConformation);
      getInt("-al_bdf2_extrapolate", cfg.bdf2Extrapolate);
      getInt("-al_output_every", cfg.outputEvery);

      getReal("-al_outlet_pressure", cfg.outletPressure);
      getReal("-al_vms_scale", cfg.vmsScale);
      getReal("-al_graddiv_scale", cfg.gradDivScale);
      getReal("-al_pressure_scale", cfg.pressureScale);
      getReal("-al_stress_div_scale", cfg.stressDivScale);
      getReal("-al_stress_scale", cfg.stressScale);
      getInt("-al_conformation_its", cfg.conformationIterations);
      getReal("-al_conformation_tol", cfg.conformationTolerance);
      getReal("-al_newton_max_step", cfg.newtonMaxStep);
      getReal("-al_max_velocity_factor", cfg.maxVelocityFactor);
      getBool("-al_vms", cfg.useVMS);
      getBool("-al_inlet_conformation", cfg.inletConformation);
      getBool("-al_trace_newton", cfg.traceNewton);
      getBool("-al_newton_d2", cfg.newtonD2);
      getInt("-al_psi_order", cfg.psiQuadratureOrder);
      getReal("-al_inlet_perturbation", cfg.inletPerturbation);
      getReal("-al_perturbation_cycles", cfg.perturbationCycles);
      getBool("-al_adaptive_dt", cfg.adaptiveTimeStep);
      getInt("-al_max_step_retries", cfg.maxStepRetries);
      getReal("-al_dt_shrink", cfg.dtShrink);
      getReal("-al_dt_grow", cfg.dtGrow);
      getInt("-al_dt_grow_after", cfg.dtGrowAfter);
      getInt("-al_newton_easy", cfg.newtonEasy);
      getInt("-al_newton_hard", cfg.newtonHard);
      getBool("-al_newton_predict", cfg.newtonPredict);
      getReal("-al_newton_predict_safety", cfg.newtonPredictSafety);
      getInt("-al_newton_calibrate_every", cfg.newtonCalibrateEvery);
      getInt("-al_pc_reuse_solves", cfg.pcReuseSolves);
      getReal("-al_pc_reuse_growth", cfg.pcReuseGrowth);
      getBool("-al_outlet_conformation_inflow", cfg.outletConformationInflow);
      getBool("-al_psi_skew", cfg.psiSkewAdvection);

      // The XDMF and CSV writers do not create directories.
      if (world.rank() == 0)
      {
        for (const std::string& path : { cfg.xdmfBasename, cfg.csvPath })
        {
          const auto dir = std::filesystem::path(path).parent_path();
          if (!dir.empty())
            std::filesystem::create_directories(dir);
        }
      }
      world.barrier();

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
