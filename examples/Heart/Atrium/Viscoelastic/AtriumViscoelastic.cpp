/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file AtriumViscoelastic.cpp
 * @brief Driver for Newtonian, generalised Newtonian, viscoelastic and
 *        thixotropic viscoelastic blood in the fluid region of the Atrium
 *        mesh, fully-implicit log-conformation; see the header and README.md.
 *
 * Run (from the build directory):
 *   mpirun -n 60 ./examples/Heart/Atrium/Viscoelastic/AtriumViscoelastic \
 *     -av_mesh ../resources/examples/Heart/Atrium/Coarse.mesh \
 *     -av_fluid owens -av_cycles 6 -av_steps_per_cycle 850 -av_time_order 2
 *
 * Options (all prefixed -av_):
 *   mesh, xdmf, csv                                strings
 *   mesh_scale, reference_length, reference_velocity, rho, period
 *   p_in, p_in_amplitude, ramp_cycles              venous pressure (Pa) and its ramp
 *   inlet_impedance, inlet_tangential_damping, inlet_backflow, outlet_backflow
 *   lv_model (bool), lv_feedback (bool), mitral_r_open, mitral_r_closed,
 *   p_out, r_out                                   mitral plane (header)
 *   fluid (newtonian|cy|owens|gsptt|sptt), eta_ref,
 *   cy_eta_zero, cy_eta_inf, cy_lambda, cy_n, cy_a
 *   owens_equilibrium (bool), owens_lambda_h, owens_eta_plasma, owens_a0,
 *   owens_a_crit, owens_gd_crit, owens_gd_max, owens_b_floor (bool),
 *   owens_epsilon, owens_n0, owens_n_inlet, owens_reset_cycle, inlet_shear_rate
 *   eta_s, eta_p, lambda, de, ptt_epsilon, lambda0_factor, lambda0_min   sPTT fluids
 *   steps_per_cycle, dt, cycles, time_order, output_every
 *   appendage_center x,y,z  appendage_radius       ball of the appendage (m)
 *   metrics (bool), probes "x,y,z;x,y,z", probe_radius
 *   and the numerical options of the arterial driver: formulation, c1..c4,
 *   dynamic_subscales, vms_scale, graddiv_scale, pressure_scale,
 *   stress_div_scale, stress_scale, conformation_its, conformation_tol,
 *   psi_order, newton_d2, adaptive_dt, max_step_retries, dt_shrink, dt_grow,
 *   dt_grow_after, newton_easy, newton_hard, newton_predict,
 *   newton_predict_safety, newton_calibrate_every, pc_reuse_solves,
 *   pc_reuse_growth, wall_kinematics, inlet_conformation,
 *   outlet_conformation_inflow, psi_skew, two_field, residence_time,
 *   jacobian_check, max_velocity_factor, dump_stress, dump_max, trace_newton.
 *
 * Output: the XDMF fields of the arterial driver, <csv> with one row per
 * accepted step (openings, fluxes, Newton, structure, age, and the volume
 * diagnostics of Config::metrics) and <csv stem>_probes.csv.
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
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "AtriumViscoelastic.h"
#include "../../CoronaryArtery/CoronaryArteryAlerts.h"
#include "../../CoronaryArtery/CoronaryArteryTiming.h"

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
  // viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit.cpp,
  // unchanged: psi split into its diagonal (D) and off-diagonal (O) blocks.
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

    enum class Block
    {
      D,
      O
    };

    constexpr size_t voigtIndex(Block b, size_t c)
    {
      return (b == Block::D ? 0 : 3) + c;
    }

    /**
     * @brief E_0..2 = e_i e_i; E_3 = e_x e_y + e_y e_x, E_4 = e_x e_z + e_z
     *        e_x, E_5 = e_y e_z + e_z e_y.
     */
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
      return VectorFunction{j == 0 ? 1.0 : 0.0, j == 1 ? 1.0 : 0.0, j == 2 ? 1.0 : 0.0};
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
      return VectorFunction{Component(m, 0, 0), Component(m, 1, 1), Component(m, 2, 2)};
    }

    template <class M>
    auto voigtO(const M& m)
    {
      return VectorFunction{Component(m, 0, 1), Component(m, 0, 2), Component(m, 1, 2)};
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

    /**
     * @brief Eigen-decomposition of a symmetric 3x3 psi, with exp, Dexp and
     *        Dexp^{-1} by Daleckii-Krein: in the eigenbasis R, R^T Dexp[H] R =
     *        F o (R^T H R), F_ij = (e^li - e^lj)/(li - lj), F_ii = e^li; F > 0,
     *        so Dexp^{-1} needs no special case for equal eigenvalues.
     */
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
      {
        for (int j = 0; j < 3; ++j)
        {
          const Real ei = std::exp(out.l(i)), d = out.l(j) - out.l(i);
          out.F(i, j) = (i == j)
            ? ei
            : (std::abs(d) > 1e-12 ? ei * std::expm1(d) / d : ei * (1.0 + 0.5 * d));
        }
      }
      return out;
    }

    template <class VD, class VO>
    Eigen::Matrix3d unpack(const VD& d, const VO& o)
    {
      Eigen::Matrix3d m;
      m << d(0), o(0), o(1), o(0), d(1), o(2), o(1), o(2), d(2);
      return m;
    }

    Math::SpatialVector<Real> packD(const Eigen::Matrix3d& m)
    {
      return Math::SpatialVector<Real>{{m(0, 0), m(1, 1), m(2, 2)}};
    }

    Math::SpatialVector<Real> packO(const Eigen::Matrix3d& m)
    {
      return Math::SpatialVector<Real>{{m(0, 1), m(0, 2), m(1, 2)}};
    }

    /// @brief The constants of the psi equation.
    struct LogModel
    {
      Real lambda = 0.0;
      Real lambda0 = 0.0;
      Real epsilon = 0.0;   ///< sPTT
      Real c = 0.0;         ///< 2 (1 - lambda0/lambda), on eps(u)

      /**
       * @brief Exp form (Config::formulation "exp"): the BDF history of the
       *        conformation, E_hist = bdfA exp(psi^n) + bdfB exp(psi^{n-1}),
       *        and the mass step, d/dt E ~ (E - E_hist)/dtM.
       */
      Real bdfA = 1.0;
      Real bdfB = 0.0;
      Real dtM = 1.0;
    };

    /**
     * @brief The model at a point: lambda (= lambda0) from the Owens field
     *        by cell when there is one, the constants otherwise.
     */
    LogModel pointModel(const LogModel& m,
      const AtriumViscoelastic::CellGridFunctionType* lambda, const Point& p)
    {
      if (!lambda)
        return m;
      LogModel out = m;
      out.lambda = lambda->getValue(p);
      out.lambda0 = out.lambda;
      return out;
    }

    /**
     * @brief Config::wallKinematics: the wall-node indicator and the nodal
     *        wall normal (by component), set by the driver; null when off.
     */
    struct WallKinematics
    {
      const AtriumViscoelastic::ScalarGridFunctionType* flag = nullptr;
      std::array<const AtriumViscoelastic::ScalarGridFunctionType*, 3> normal{};
    };

    WallKinematics& wallKinematics()
    {
      static WallKinematics w;
      return w;
    }

    /**
     * @brief The wall weight w in [0, 1] at p and the unit wall normal n; 0
     *        off the first element layer or when the correction is off.
     */
    Real wallWeight(const Point& p, Eigen::Vector3d& n)
    {
      const auto& wk = wallKinematics();
      if (!wk.flag)
        return 0.0;
      const Real w = wk.flag->getValue(p);
      if (!(w > 1.0e-12))
        return 0.0;
      n = Eigen::Vector3d(
        wk.normal[0]->getValue(p), wk.normal[1]->getValue(p), wk.normal[2]->getValue(p));
      const Real nn = n.norm();
      if (!(nn > 0.0))
        return 0.0;
      n /= nn;
      return std::min<Real>(w, 1.0);
    }

    /**
     * @brief L' = (1 - w) L + w (I - n n^T) L n n^T: the pure shear a n^T,
     *        a = (I - n n^T) L n, that a no-slip wall allows, blended in.
     */
    Eigen::Matrix3d wallGradient(
      const Eigen::Matrix3d& L, Real w, const Eigen::Vector3d& n)
    {
      if (w <= 0.0)
        return L;
      const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
      return (1.0 - w) * L + w * (P * L * n) * n.transpose();
    }

    /**
     * @brief What the forms read at one quadrature point, from (psi^k, L^k =
     *        grad u^k): exp(psi) and Dexp[E_j] (momentum, S1), and the Newton
     *        data of the psi equation paired against the test blocks.
     *
     * @details With N = -Dexp^{-1}[L T + T L^T - c eps] + (f/lambda)(I -
     *          exp(-psi)), its FD Jacobian JN_j = dN/dpsi_j, the exact
     *          (linear) velocity dependence Q_ab = -Dexp^{-1}[E_ab T + T
     *          E_ab^T - c sym(E_ab)] and R = N - JN[psi^k] - sum_ab L^k_ab
     *          Q_ab, the pairings against block C (0 = D, 1 = O) are
     *            A[C][B](k, c) = JN_{B,c} : E_{C,k},
     *            P[C][b](k, a) = Q_ab : E_{C,k},
     *            r[C](k) = R : E_{C,k}.
     */
    struct LogPoint
    {
      Eigen::Matrix3d exp;
      std::array<Eigen::Matrix3d, 6> D;
      std::array<std::array<Eigen::Matrix3d, 2>, 2> A;
      std::array<std::array<Eigen::Matrix3d, 3>, 2> P;
      std::array<Eigen::Vector3d, 2> r;

      /**
       * @brief V_m = sum_j (D2exp_{psi^k}[d_j psi^k, E_m]) e_j, m = 0..5:
       *        the derivative of div exp(psi) = sum_j Dexp_psi[d_j psi] e_j
       *        in psi beyond the one carried by d_j psi; Vk = sum_m V_m psi^k_m.
       */
      std::array<Eigen::Vector3d, 6> V;
      Eigen::Vector3d Vk;

      /**
       * @brief Exp form (LogExpForm), at (psi^k, u^k), with L = grad u^k,
       *        E = exp(psi^k), k_a = lambda/(2 lambda0) and the conformation
       *        residual
       *          R = f (E - I)/(2 lambda0) - eps(u)
       *              + k_a [(E - E_hist)/dtM + (u.grad) E - (L E + E L^T) + 2 eps(u)],
       *        the sigma equation over 2 eta_p (Moreno et al. 2019, eq. 2.10):
       *        G_j = d_j E = Dexp[d_j psi]; W_m = sum_j u_j D2exp[d_j psi, E_m];
       *        M_m = dR/dpsi_m but for the k_a (u.grad) Dexp[.] part,
       *          M_m = (f D_m + (df/dpsi_m)(E - I))/(2 lambda0)
       *                + k_a (D_m/dtM + W_m - (L D_m + D_m L^T));
       *        K = R - sum_m psi_m M_m - k_a sum_j u_j G_j - R_u[u], the part
       *        of R that is not linear in (psi, u) at (psi^k, u^k), R_u the
       *        exact (linear) velocity dependence; R itself.
       */
      std::array<Eigen::Matrix3d, 3> G;
      std::array<Eigen::Matrix3d, 6> W;
      std::array<Eigen::Matrix3d, 6> M;
      Eigen::Matrix3d K;
      Eigen::Matrix3d R;

      /**
       * @brief The BDF part of R, k_a (E - E_hist)/dtM, and its part that
       *        is not linear in psi, Kt = k_a (E - E_hist - sum_m psi_m D_m)/dtM
       *        (Config::subscaleTimeResidual).
       */
      Eigen::Matrix3d rate;
      Eigen::Matrix3d Kt;
    };

    Eigen::Matrix3d nonlinearMap(const Eigen::Matrix3d& psi, const Eigen::Matrix3d& L,
      const LogModel& m, Real* fOut = nullptr)
    {
      const LogEig e = logEig(psi);
      const Eigen::Matrix3d T = e.exp();
      const Eigen::Matrix3d H =
        L * T + T * L.transpose() - 0.5 * m.c * (L + L.transpose());
      const Real f = 1.0 + m.epsilon * (m.lambda / m.lambda0) * (T.trace() - 3.0);
      if (fOut)
        *fOut = f;
      return -e.dexpInverse(H) +
        (f / m.lambda) * (Eigen::Matrix3d::Identity() - e.expMinus());
    }

    /**
     * @brief The parts of a LogPoint. Dexp (exp and Dexp[E_j]) costs one
     *        eigenproblem; Newton (A, P, r) thirteen more, twelve of them for
     *        the FD Jacobian. The step-value store (stabilisation) reads only
     *        Dexp.
     */
    enum LogPart : std::uint8_t
    {
      LogDexp = 1,
      LogNewton = 2,
      LogD2 = 4, ///< V_m, Vk (second derivative of exp along d_j psi^k)
      LogExpForm = 8 ///< G, W, M, K, R of the exp form
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
      const Eigen::Matrix3d& Lraw, const LogModel& m, Real w = 0.0,
      const Eigen::Vector3d& n = Eigen::Vector3d::Zero())
    {
      // Config::wallKinematics: N is evaluated at L' = wallGradient(L), linear
      // in L, so its velocity derivative Q'_ab = dN/dL_ab follows by the chain
      // rule (below) and the residual keeps the same form in terms of L.
      const Eigen::Matrix3d L = wallGradient(Lraw, w, n);
      const auto& B = voigtBasisEigen();
      // Central differences in the six Voigt directions (twelve 3x3
      // eigenproblems); the map is smooth in psi.
      const Real h = 1.0e-6 * std::max<Real>(1.0, psi.cwiseAbs().maxCoeff());
      std::array<Eigen::Matrix3d, 6> JN;
      for (size_t j = 0; j < 6; ++j)
      {
        JN[j] =
          (nonlinearMap(psi + h * B[j], L, m) - nonlinearMap(psi - h * B[j], L, m)) /
          (2.0 * h);
      }
      std::array<std::array<Eigen::Matrix3d, 3>, 3> Q;
      for (size_t a = 0; a < 3; ++a)
      {
        for (size_t b = 0; b < 3; ++b)
        {
          Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
          E(a, b) = 1.0;
          Q[a][b] = -e.dexpInverse(
            E * out.exp + out.exp * E.transpose() - 0.5 * m.c * (E + E.transpose()));
        }
      }
      if (w > 0.0)
      {
        // dL'_cd/dL_ab = (1 - w) d_ca d_db + w P_ca n_b n_d, P = I - n n^T:
        // Q'_ab = (1 - w) Q_ab + w n_b sum_c P_ca (sum_d Q_cd n_d).
        const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
        std::array<Eigen::Matrix3d, 3> Qn;
        for (size_t c = 0; c < 3; ++c)
        {
          Qn[c].setZero();
          for (size_t dd = 0; dd < 3; ++dd)
            Qn[c] += n(dd) * Q[c][dd];
        }
        std::array<std::array<Eigen::Matrix3d, 3>, 3> Qw;
        for (size_t a = 0; a < 3; ++a)
        {
          Eigen::Matrix3d sa = Eigen::Matrix3d::Zero();
          for (size_t c = 0; c < 3; ++c)
            sa += P(c, a) * Qn[c];
          for (size_t b = 0; b < 3; ++b)
            Qw[a][b] = (1.0 - w) * Q[a][b] + w * n(b) * sa;
        }
        Q = Qw;
      }
      // psi = sum_j psi_j E_j, psi_j = psi : E_j / |E_j|^2.
      Eigen::Matrix3d R = nonlinearMap(psi, L, m);
      for (size_t j = 0; j < 6; ++j)
        R -= (psi.cwiseProduct(B[j]).sum() / B[j].squaredNorm()) * JN[j];
      for (size_t a = 0; a < 3; ++a)
      {
        for (size_t b = 0; b < 3; ++b)
          R -= Lraw(a, b) * Q[a][b];
      }
      for (size_t C = 0; C < 2; ++C)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          const Eigen::Matrix3d& Ek = B[3 * C + k];
          out.r[C](k) = R.cwiseProduct(Ek).sum();
          for (size_t Bk = 0; Bk < 2; ++Bk)
          {
            for (size_t c = 0; c < 3; ++c)
              out.A[C][Bk](k, c) = JN[3 * Bk + c].cwiseProduct(Ek).sum();
          }
          for (size_t a = 0; a < 3; ++a)
          {
            for (size_t b = 0; b < 3; ++b)
              out.P[C][b](k, a) = Q[a][b].cwiseProduct(Ek).sum();
          }
        }
      }
    }

    /// @brief A scalar field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseScalar final : public RealFunctionBase<PointwiseScalar<F>>
    {
      public:
        using Parent = RealFunctionBase<PointwiseScalar<F>>;

        explicit PointwiseScalar(F f, size_t order = 2)
          : m_f(std::move(f)),
            m_order(order)
        {}
        PointwiseScalar(const PointwiseScalar& other)
          : Parent(other),
            m_f(other.m_f),
            m_order(other.m_order)
        {}
        PointwiseScalar(PointwiseScalar&& other)
          : Parent(std::move(other)),
            m_f(std::move(other.m_f)),
            m_order(other.m_order)
        {}

        Real getValue(const Point& p) const
        {
          return m_f(p);
        }

        /**
         * Not a polynomial; quadratic is enough against P1 x P1. Order 0 is
         * used for a function known to be constant, so the quadrature of the
         * parent formulation is kept.
         */
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
        {
          return m_order;
        }
        PointwiseScalar* copy() const noexcept override
        {
          return new PointwiseScalar(*this);
        }

      private:
        F m_f;
        size_t m_order;
    };

    /// @brief A 3x3 tensor field given by a callable, evaluated pointwise.
    template <class F>
    class PointwiseTensor final : public MatrixFunctionBase<Real, PointwiseTensor<F>>
    {
      public:
        using Parent = MatrixFunctionBase<Real, PointwiseTensor<F>>;

        explicit PointwiseTensor(F f)
          : m_f(std::move(f))
        {}
        PointwiseTensor(const PointwiseTensor& other)
          : Parent(other),
            m_f(other.m_f)
        {}
        PointwiseTensor(PointwiseTensor&& other)
          : Parent(std::move(other)),
            m_f(std::move(other.m_f))
        {}

        Math::SpatialMatrix<Real> getValue(const Point& p) const
        {
          return Math::SpatialMatrix<Real>(m_f(p));
        }
        size_t getRows() const
        {
          return 3;
        }
        size_t getColumns() const
        {
          return 3;
        }

        /// Not a polynomial; quadratic is enough against P1 x P1.
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
        {
          return 2;
        }
        PointwiseTensor* copy() const noexcept override
        {
          return new PointwiseTensor(*this);
        }

      private:
        F m_f;
    };

    /**
     * @brief A constant 3x3 tensor in fixed-size storage (Math::Matrix is
     *        dynamic: products with it allocate at every evaluation).
     */
    class ConstantTensor final : public MatrixFunctionBase<Real, ConstantTensor>
    {
      public:
        using Parent = MatrixFunctionBase<Real, ConstantTensor>;

        explicit ConstantTensor(const Eigen::Matrix3d& m)
          : m_m(m)
        {}
        ConstantTensor(const ConstantTensor& other)
          : Parent(other),
            m_m(other.m_m)
        {}
        ConstantTensor(ConstantTensor&& other)
          : Parent(std::move(other)),
            m_m(other.m_m)
        {}

        Math::SpatialMatrix<Real> getValue(const Point&) const
        {
          return Math::SpatialMatrix<Real>(m_m);
        }
        size_t getRows() const
        {
          return 3;
        }
        size_t getColumns() const
        {
          return 3;
        }
        Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
        {
          return 0;
        }
        ConstantTensor* copy() const noexcept override
        {
          return new ConstantTensor(*this);
        }

      private:
        Eigen::Matrix3d m_m;
    };

    /// @brief tensor<B>, with fixed-size basis tensors (exp formulation).
    template <Block B, class S>
    auto tensorFixed(const S& s)
    {
      const auto& E = voigtBasisEigen();
      return ConstantTensor(E[voigtIndex(B, 0)]) * Component(s, 0) +
        ConstantTensor(E[voigtIndex(B, 1)]) * Component(s, 1) +
        ConstantTensor(E[voigtIndex(B, 2)]) * Component(s, 2);
    }

    /**
     * @brief The fields a LogPoint is computed from: (psi_D, psi_O, u) and
     *        the revision that keys the cache.
     */
    template <class GF>
    struct LogFields
    {
      const GF* psiD;
      const GF* psiO;
      const GF* u;
      LogModel model;
      const std::uint64_t* revision;

      /// @brief psi^n and psi^{n-1}, for E_hist (LogExpForm only).
      const GF* psiND = nullptr;
      const GF* psiNO = nullptr;
      const GF* psiNm1D = nullptr;
      const GF* psiNm1O = nullptr;

      /// @brief Owens: lambda by cell (nullptr: LogModel::lambda).
      const AtriumViscoelastic::CellGridFunctionType* lambda = nullptr;
    };

    /**
     * @brief logPoint(psi(x), grad u(x)) at every quadrature point of the
     *        current Newton iterate, computed once and kept until the
     *        iterate changes.
     *
     * @details Rodin assembles integrator by integrator over the whole mesh
     *          so a small recency cache recomputes the twelve eigenproblems
     *          of every point once per integrator that reads them. One store per
     *          psi_D field (the iterate psi^k and the step value psi^n),
     *          cleared when the revision changes; keyed by the polytope
     *          (dimension, index) and the reference coordinates. Every
     *          integral reading it uses one rule (psiOrder in setupFlow), so
     *          a cell holds one set of points.
     * @brief A quadrature point: polytope (dimension, index), reference
     *        coordinates.
     */
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
        return PointKey{polytope.getDimension(), polytope.getIndex(), rc(0),
          rc.size() > 1 ? rc(1) : 0.0, rc.size() > 2 ? rc(2) : 0.0};
      }
    };

    struct PointKeyHash
    {
      size_t operator()(const PointKey& k) const noexcept
      {
        size_t h = std::hash<size_t>()(k.dimension);
        for (const size_t v : {std::hash<Index>()(k.index), std::hash<Real>()(k.r0),
               std::hash<Real>()(k.r1), std::hash<Real>()(k.r2)})
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
        const Entry& at(
          const LogFields<GF>& s, const Point& p, const PointKey& key, std::uint8_t need)
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
          const LogModel md = pointModel(s.model, s.lambda, p);
          if (!(e.have & LogDexp))
          {
            logPointDexp(e.value, eig);
            e.have |= LogDexp;
          }
          if (need & LogNewton & ~e.have)
          {
            const auto J = Jacobian(*s.u).getValue(p); // Math::SpatialMatrix, not Eigen
            Eigen::Matrix3d L;
            for (int a = 0; a < 3; ++a)
            {
              for (int b = 0; b < 3; ++b)
                L(a, b) = J(a, b);
            }
            Eigen::Vector3d nw = Eigen::Vector3d::Zero();
            const Real ww = wallWeight(p, nw);
            logPointNewton(e.value, eig, psi, L, md, ww, nw);
            e.have |= LogNewton;
          }
          const bool wantD2 = need & LogD2 & ~e.have;
          const bool wantExp = need & LogExpForm & ~e.have;
          if (wantD2 || wantExp)
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
            // The exp form also needs u^k, L^k and (u^k . grad) psi^k.
            Eigen::Vector3d u = Eigen::Vector3d::Zero();
            Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d gu = Eigen::Matrix3d::Zero();
            if (wantExp)
            {
              const auto uv = s.u->getValue(p);
              const auto J = Jacobian(*s.u).getValue(p);
              for (int a = 0; a < 3; ++a)
              {
                u(a) = uv(a);
                for (int b = 0; b < 3; ++b)
                  L(a, b) = J(a, b);
              }
              for (int j = 0; j < 3; ++j)
                gu += u(j) * g[j];
            }
            // Central differences of Dexp in the six Voigt directions (twelve
            // eigenproblems), shared by V_m (div exp, S1) and W_m (exp form).
            const Real h = 1.0e-6 * std::max<Real>(1.0, psi.cwiseAbs().maxCoeff());
            const Real pk[6] = {
              psi(0, 0), psi(1, 1), psi(2, 2), psi(0, 1), psi(0, 2), psi(1, 2)};
            if (wantD2)
              e.value.Vk.setZero();
            for (size_t m = 0; m < 6; ++m)
            {
              const LogEig ep = logEig(psi + h * B[m]);
              const LogEig em = logEig(psi - h * B[m]);
              if (wantD2)
              {
                e.value.V[m].setZero();
                for (int j = 0; j < 3; ++j)
                  e.value.V[m] += ((ep.dexp(g[j]) - em.dexp(g[j])) / (2.0 * h)).col(j);
                e.value.Vk += pk[m] * e.value.V[m];
              }
              if (wantExp)
                e.value.W[m] = (ep.dexp(gu) - em.dexp(gu)) / (2.0 * h);
            }
            if (wantD2)
              e.have |= LogD2;
            if (wantExp)
            {
              const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
              const Eigen::Matrix3d& E = e.value.exp;
              const Real ka = md.lambda / (2.0 * md.lambda0);
              const Real ptt = md.epsilon * (md.lambda / md.lambda0);
              const Real fk = 1.0 + ptt * (E.trace() - 3.0);
              // E_hist from psi^n and psi^{n-1} at this point.
              const auto expAt = [&p](const GF* d, const GF* o) -> Eigen::Matrix3d {
                return logEig(unpack(d->getValue(p), o->getValue(p))).exp();
              };
              Eigen::Matrix3d Eh = md.bdfA * expAt(s.psiND, s.psiNO);
              if (md.bdfB != 0.0)
                Eh += md.bdfB * expAt(s.psiNm1D, s.psiNm1O);
              for (int j = 0; j < 3; ++j)
                e.value.G[j] = eig.dexp(g[j]);
              Eigen::Matrix3d adv = Eigen::Matrix3d::Zero(); // (u . grad) E
              for (int j = 0; j < 3; ++j)
                adv += u(j) * e.value.G[j];
              const Eigen::Matrix3d symL = 0.5 * (L + L.transpose());
              Eigen::Matrix3d sumPsiM = Eigen::Matrix3d::Zero();
              for (size_t m = 0; m < 6; ++m)
              {
                const Eigen::Matrix3d& Dm = e.value.D[m];
                // df/dpsi_m = ptt tr(Dexp[E_m]) (= ptt tr(E E_m)).
                e.value.M[m] =
                  (fk * Dm + ptt * Dm.trace() * (E - Id)) / (2.0 * md.lambda0) +
                  ka * (Dm / md.dtM + e.value.W[m] - (L * Dm + Dm * L.transpose()));
                sumPsiM += pk[m] * e.value.M[m];
              }
              const Eigen::Matrix3d relax = fk * (E - Id) / (2.0 * md.lambda0);
              const Eigen::Matrix3d rate = ka * (E - Eh) / md.dtM;
              e.value.R = relax - symL + rate +
                ka * (adv - (L * E + E * L.transpose()) + 2.0 * symL);
              // R = sum_m psi_m M_m + k_a (u . grad) Dexp[psi] + R_u[u] + K, with
              // R_u[u] = (k_a 2 - 1) eps(u) + k_a (sum_j u_j G_j - (L E + E L^T)).
              e.value.K = relax + rate - sumPsiM - ka * adv;
              e.value.rate = rate;
              e.value.Kt = rate;
              for (size_t m = 0; m < 6; ++m)
                e.value.Kt -= (ka / md.dtM) * pk[m] * e.value.D[m];
              e.have |= LogExpForm;
            }
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

    /**
     * @brief One store per psi_D field, and a counter bumped by every
     *        release so that no thread reuses a pointer into a freed store.
     */
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

    /**
     * @brief Frees the points stored for `field` (a psi_D field) once the
     *        forms that read them have been assembled.
     */
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
      last = Last{s.psiD, revision, epoch, key, &entry};
      return entry.value;
    }

    /**
     * @brief A pointwise matrix drawn from the cached LogPoint by `pick`,
     *        which reads only the parts in `need`.
     */
    template <class GF, class Pick>
    auto fromPoint(const LogFields<GF>& s, std::uint8_t need, Pick pick)
    {
      return PointwiseTensor([s, need, pick](const Point& p) -> Eigen::Matrix3d {
        return pick(logPointAt(s, p, need));
      });
    }

    /**
     * @brief Dexp_{psi(x)}[w] restricted to block B: sum_c w_c Dexp[E_{B,c}],
     *        for w a block of psi or of a derivative of it; trial or field.
     */
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

    /**
     * @brief A pointwise 3-vector drawn from the cached LogPoint by `pick`,
     *        which reads only the parts in `need`.
     */
    template <class GF, class Pick>
    auto fromPointVector(const LogFields<GF>& s, std::uint8_t need, Pick pick)
    {
      return VectorFunction(
        size_t(3), [s, need, pick](const Point& p) -> Math::SpatialVector<Real> {
          const Eigen::Vector3d v = pick(logPointAt(s, p, need));
          return Math::SpatialVector<Real>{{v(0), v(1), v(2)}};
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
      return Mult(fromPoint(s, LogNewton,
                    [mass](const LogPoint& q) -> Eigen::Matrix3d {
                      return q.A[c][b] + mass * Eigen::Matrix3d::Identity();
                    }),
        w);
    }

    /**
     * @brief pair<C>(sum_ab L_ab(v) Q_ab) = sum_b P[C][b] (J(v) e_b), for a
     *        vector v (trial or field): the exact velocity dependence of the
     *        psi equation at psi^k, L_ab = d v_a / d x_b.
     */
    template <Block C, class GF, class V>
    auto pairedVelocityMap(const LogFields<GF>& s, const V& v)
    {
      constexpr size_t c = blockIndex(C);
      const auto P = [&](size_t b) {
        return fromPoint(
          s, LogNewton, [b](const LogPoint& q) -> Eigen::Matrix3d { return q.P[c][b]; });
      };
      return Mult(P(0), Mult(Jacobian(v), unit(0))) +
        Mult(P(1), Mult(Jacobian(v), unit(1))) + Mult(P(2), Mult(Jacobian(v), unit(2)));
    }

    /// @brief pair<C>(N - JN[psi^k] - sum_ab L^k_ab Q_ab), grade 0.
    template <Block C, class GF>
    auto pairedResidual(const LogFields<GF>& s)
    {
      constexpr size_t c = blockIndex(C);
      return fromPointVector(
        s, LogNewton, [](const LogPoint& q) -> Eigen::Vector3d { return q.r[c]; });
    }

    /**
     * @brief Block B's contribution to div exp(psi) = sum_j Dexp[d_j psi] e_j
     *        (chain rule), Dexp at the point of `s`.
     */
    template <Block B, class GF, class S>
    auto divExp(const LogFields<GF>& s, const S& psi)
    {
      return Mult(dexp<B>(s, Mult(Jacobian(psi), unit(0))), unit(0)) +
        Mult(dexp<B>(s, Mult(Jacobian(psi), unit(1))), unit(1)) +
        Mult(dexp<B>(s, Mult(Jacobian(psi), unit(2))), unit(2));
    }

    /**
     * @brief divExp<B> for a trial function psi_B: sum_c Dexp[E_{B,c}]
     *        grad psi_c. A basis function of psi_B has one nonzero component
     *        c, where this is Dexp[E_{B,c}] grad phi, the value divExp<B>
     *        gives (the other terms are exact zeros), with three pointwise
     *        matrices per point instead of nine.
     */
    template <Block B, class GF, class S>
    auto divExpTrial(const LogFields<GF>& s, const S& psi)
    {
      const auto D = [&](size_t c) {
        const size_t j = voigtIndex(B, c);
        return fromPoint(s, LogDexp, [j](const LogPoint& q) { return q.D[j]; });
      };
      const auto gradOf = [&](
                            size_t c) { return Mult(Transpose(Jacobian(psi)), unit(c)); };
      return Mult(D(0), gradOf(0)) + Mult(D(1), gradOf(1)) + Mult(D(2), gradOf(2));
    }

    /**
     * @brief Newton on div exp(psi): divExpTrial carries Dexp_{psi^k}[d_j
     *        psi]; the dependence of Dexp on psi adds sum_m V_m psi_m - Vk
     *        (block B's part of the sum here, Vk in divExpKnownD2). Without
     *        it Newton converges only linearly once psi is large.
     */
    template <Block B, class GF, class S>
    auto divExpTrialD2(const LogFields<GF>& s, const S& psi)
    {
      const auto V = [&](size_t c) {
        const size_t m = voigtIndex(B, c);
        return fromPointVector(
          s, LogD2, [m](const LogPoint& q) -> Eigen::Vector3d { return q.V[m]; });
      };
      return V(0) * Component(psi, 0) + V(1) * Component(psi, 1) +
        V(2) * Component(psi, 2);
    }

    template <class GF>
    auto divExpKnownD2(const LogFields<GF>& s)
    {
      return fromPointVector(
        s, LogD2, [](const LogPoint& q) -> Eigen::Vector3d { return q.Vk; });
    }

    // ---- Exp form (Config::formulation "exp") ----------------------------------
    //
    // The constitutive equation is the sigma equation over 2 eta_p written in
    // E = exp(psi) (Moreno, Codina, Baiges, Castillo 2019, eq. 2.10), tested
    // with a stress-like chi, and linearised by Newton about (psi^k, u^k):
    //   R(psi, u) ~ sum_m M_m psi_m + k_a Dexp[(u^k . grad) psi] + R_u[u] + K.

    /// @brief Trial part of R in the block B of psi.
    template <Block B, class GF, class S, class U>
    auto expFormTrialPsi(const LogFields<GF>& s, Real ka, const S& psi, const U& uk)
    {
      const auto Mc = [&](size_t c) {
        const size_t m = voigtIndex(B, c);
        return fromPoint(s, LogExpForm, [m](const LogPoint& q) { return q.M[m]; });
      };
      const auto Dc = [&](size_t c) {
        const size_t m = voigtIndex(B, c);
        return fromPoint(s, LogDexp, [m](const LogPoint& q) { return q.D[m]; });
      };
      const auto adv = Mult(Jacobian(psi), uk);
      return Mc(0) * Component(psi, 0) + Mc(1) * Component(psi, 1) +
        Mc(2) * Component(psi, 2) + ka * Dc(0) * Component(adv, 0) +
        ka * Dc(1) * Component(adv, 1) + ka * Dc(2) * Component(adv, 2);
    }

    /**
     * @brief Trial part of R in u: R_u[u] = (2 k_a - 1) eps(u) + k_a (sum_j
     *        u_j d_j E - (grad u E + E grad u^T)), E and d_j E at psi^k.
     */
    template <class GF, class U>
    auto expFormTrialU(const LogFields<GF>& s, Real ka, const U& u)
    {
      const auto Gj = [&](size_t j) {
        return fromPoint(s, LogExpForm, [j](const LogPoint& q) { return q.G[j]; });
      };
      const auto E = fromPoint(s, LogDexp, [](const LogPoint& q) { return q.exp; });
      return (ka - 0.5) * (Jacobian(u) + Transpose(Jacobian(u))) +
        ka * Gj(0) * Component(u, 0) + ka * Gj(1) * Component(u, 1) +
        ka * Gj(2) * Component(u, 2) + (-ka) * Mult(Jacobian(u), E) +
        (-ka) * Mult(E, Transpose(Jacobian(u)));
    }

    /**
     * @brief The adjoint test operator of the constitutive subscale, block C
     *        of chi (Moreno, Codina, Baiges 2020, eq. 27):
     *          A(chi) = f chi/(2 eta_p) - (lambda/(2 eta_p)) ((u^n . grad) chi
     *                   + chi L^n + (L^n)^T chi),  L^n = grad u^n;
     *        the eps(v) part is paired separately.
     */
    template <Block C, class T, class F, class U, class L>
    auto adjointTest(
      const T& chi, const F& fN, Real kA, Real kB, const U& un, const L& ln)
    {
      const auto X = tensorFixed<C>(chi);
      return kA * fN * X + (-kB) * tensorFixed<C>(Mult(Jacobian(chi), un)) +
        (-kB) * Mult(X, ln) + (-kB) * Mult(Transpose(ln), X);
    }

    /**
     * @brief Dynamic subgrid scales at the quadrature points, step n: the
     *        velocity parts u1 (convection), u2 (pressure gradient), u3 (div
     *        sigma) of the term-by-term momentum subscale, and the stress
     *        subscale sig. Written for step n+1 once the step is accepted, then
     *        committed; keyed like LogStore, so every integral that reads it
     *        uses the psi quadrature rule. Missing points read zero (start).
     */
    class SubscaleStore
    {
      public:
        struct Entry
        {
          Eigen::Vector3d u1 = Eigen::Vector3d::Zero();
          Eigen::Vector3d u2 = Eigen::Vector3d::Zero();
          Eigen::Vector3d u3 = Eigen::Vector3d::Zero();
          Eigen::Matrix3d sig = Eigen::Matrix3d::Zero();
        };

        const Entry& at(const Point& p) const
        {
          static const Entry zero{};
          std::lock_guard<std::mutex> lock(m_mutex);
          const auto it = m_current.find(PointKey::of(p));
          return it == m_current.end() ? zero : it->second;
        }

        void stage(const Point& p, const Entry& e)
        {
          std::lock_guard<std::mutex> lock(m_mutex);
          m_next[PointKey::of(p)] = e;
        }

        /// @brief The staged values become step n.
        void commit()
        {
          std::lock_guard<std::mutex> lock(m_mutex);
          m_current.swap(m_next);
          m_next.clear();
        }

        void discard()
        {
          std::lock_guard<std::mutex> lock(m_mutex);
          m_next.clear();
        }

      private:
        mutable std::mutex m_mutex;
        std::unordered_map<PointKey, Entry, PointKeyHash> m_current;
        std::unordered_map<PointKey, Entry, PointKeyHash> m_next;
    };

    /// @brief One simulation per process: the store of its subscales.
    SubscaleStore& subscales()
    {
      static SubscaleStore store;
      return store;
    }

    /// @brief A velocity subscale of step n, as a vector function.
    template <class Pick>
    auto subscaleVector(Pick pick)
    {
      return VectorFunction(
        size_t(3), [pick](const Point& p) -> Math::SpatialVector<Real> {
          const Eigen::Vector3d v = pick(subscales().at(p));
          return Math::SpatialVector<Real>{{v(0), v(1), v(2)}};
        });
    }

    /**
     * @brief The projections invert mass matrices: CG + Jacobi, never the
     *        global direct solver used by the coupled flow system.
     */
    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  // ==========================================================================
  // AtriumViscoelastic
  // ==========================================================================
  AtriumViscoelastic::AttributeSet AtriumViscoelastic::makeWallSet(const Config& cfg)
  {
    return AttributeSet(cfg.labels.wall.begin(), cfg.labels.wall.end());
  }

  AtriumViscoelastic::AtriumViscoelastic(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_input(makeVentricleInput(m_cfg, &m_outletMeanPressure)),
      m_model(m_input),
      m_parent(makeParentMesh(context, m_cfg)),
      m_mesh(makeFluidMesh(m_parent, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_wallSet(makeWallSet(m_cfg)),
      m_inletSet(m_cfg.labels.inlets.begin(), m_cfg.labels.inlets.end()),
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
      m_piConv(m_vh, "av_proj_"),
      m_sub(m_vh, "av_proj_"),
      m_subOld(m_vh),
      m_piDiv(m_sh, "av_proj_"),
      m_piGradP(m_vh, "av_proj_"),
      m_piDivSigma(m_vh, "av_proj_"),
      m_piEpsD(m_tauh, "av_proj_"),
      m_piEpsO(m_tauh, "av_proj_"),
      m_piAdvPsiD(m_tauh, "av_proj_"),
      m_piAdvPsiO(m_tauh, "av_proj_"),
      m_piRD(m_tauh, "av_proj_"),
      m_piRO(m_tauh, "av_proj_"),
      m_alpha3Dyn(m_ch),
      m_pttN(m_ch),
      m_n(m_sh),
      m_nTest(m_sh),
      m_nOld(m_sh),
      m_nNm1(m_sh),
      m_nHist(m_sh),
      m_lambdaCell(m_ch),
      m_structure(m_n, m_nTest),
      m_structureKSP(m_structure),
      m_wss(m_vh),
      m_symRec0(m_vh),
      m_symRec1(m_vh),
      m_symRec2(m_vh),
      m_netShear(m_vh),
      m_absShear(m_sh),
      m_shearMagnitude(m_sh),
      m_tawss(m_sh),
      m_osi(m_sh),
      m_tawssAppendage(m_sh),
      m_osiAppendage(m_sh),
      m_wallArea(m_sh),
      m_appendageArea(m_sh),
      m_wallFlag(m_sh),
      m_wallN0(m_sh),
      m_wallN1(m_sh),
      m_wallN2(m_sh),
      m_wssRec(m_vh),
      m_qFlux(m_sh),
      m_one(m_sh),
      m_flux(m_qFlux),
      m_wiCell(m_ch),
      m_aSum(m_sh),
      m_bSum(m_sh),
      m_rateWork(m_sh),
      m_subscaleSweep(m_qFlux),
      m_flow(m_u, m_p, m_psiD, m_psiO, m_v, m_q, m_chiD, m_chiO),
      m_flowKSP(m_flow),
      m_flow2(m_u, m_p, m_v, m_q),
      m_flow2KSP(m_flow2),
      m_age(m_sh),
      m_ageTest(m_sh),
      m_ageOld(m_sh),
      m_ageNm1(m_sh),
      m_ageHist(m_sh),
      m_ageProblem(m_age, m_ageTest),
      m_ageKSP(m_ageProblem),
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
    {
      Alert::Info() << "[mesh] " << m_cfg.meshPath << "  cells=" << cellCount
                    << " vertices=" << vertexCount << " velocity DOFs=" << velocityDOFs
                    << " pressure DOFs=" << pressureDOFs << " psi DOFs=" << stressDOFs
                    << Alert::Raise;
    }

    // Partition quality: owned cells per rank. The internal partitioner is
    // much worse than Scotch (RODIN_USE_SCOTCH): the slowest rank sets the
    // pace of assembly and solve.
    {
      const auto& comm = m_mesh.getContext().getCommunicator();
      const auto& shard = m_mesh.getShard();
      const size_t D = shard.getDimension();
      size_t owned = 0, ghost = 0;
      for (Index i = 0; i < shard.getCellCount(); ++i)
        (shard.isOwned(D, i) ? owned : ghost) += 1;
      size_t minOwned = 0, maxOwned = 0, sumOwned = 0, maxGhost = 0;
      boost::mpi::all_reduce(comm, owned, minOwned, boost::mpi::minimum<size_t>());
      boost::mpi::all_reduce(comm, owned, maxOwned, boost::mpi::maximum<size_t>());
      boost::mpi::all_reduce(comm, owned, sumOwned, std::plus<size_t>());
      boost::mpi::all_reduce(comm, ghost, maxGhost, boost::mpi::maximum<size_t>());
      const Real mean = static_cast<Real>(sumOwned) / static_cast<Real>(comm.size());
      const Real imbalance = mean > 0.0 ? static_cast<Real>(maxOwned) / mean : 1.0;
#ifdef RODIN_USE_SCOTCH
      const char* partitioner = "Scotch";
#else
      const char* partitioner =
        "BalancedCompact (built without Scotch, RODIN_USE_SCOTCH=OFF)";
#endif
      if (isRoot())
      {
        Alert::Info() << "[partition] " << partitioner << ", " << comm.size()
                      << " ranks: owned cells min/mean/max = " << minOwned << "/"
                      << static_cast<size_t>(mean) << "/" << maxOwned
                      << " (max/mean = " << imbalance
                      << "), ghost cells max = " << maxGhost << Alert::Raise;
        if (comm.size() > 1 && imbalance > 1.2)
        {
          Alert::Warning() << "[partition] the largest shard has " << imbalance
                           << " times the mean: rebuild Rodin with -DRODIN_USE_SCOTCH=ON"
                           << Alert::Raise;
        }
      }
    }
  }

  AtriumViscoelastic::~AtriumViscoelastic() = default;

  bool AtriumViscoelastic::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  void AtriumViscoelastic::deriveParameters()
  {
    auto& ob = m_cfg.oldroydB;
    const Real D = m_cfg.referenceLength;
    const Real rho = m_cfg.rho;

    if (m_cfg.formulation != "exp" && m_cfg.formulation != "psi" &&
      m_cfg.formulation != "psi_rb")
      throw std::runtime_error("formulation must be psi_rb, exp or psi.");
    if (m_cfg.formulation != "psi" && m_cfg.psiSkewAdvection)
      throw std::runtime_error("psi_skew is a term of the psi formulation only.");
    if (m_cfg.fluid != "sptt" && m_cfg.fluid != "cy" && m_cfg.fluid != "gsptt" &&
      m_cfg.fluid != "owens" && m_cfg.fluid != "newtonian")
      throw std::runtime_error("fluid must be sptt, cy, gsptt, newtonian or owens.");
    if (m_cfg.fluid == "owens")
    {
      // Plasma solvent; the polymer at infinite shear, eta_p,inf = eta_CY(inf)
      // - eta_plasma, relaxes in lambda_H: G = eta_p/lambda = eta_p,inf/lambda_H
      // is the (constant) modulus, so sigma = G (exp(psi) - I) as for gsptt
      // and only lambda varies, by cell (Config::Owens).
      const auto& ow = m_cfg.owens;
      if (!(ow.lambdaH > 0.0))
        throw std::runtime_error("owens: lambda_H must be positive.");
      ob.etaS = ow.etaPlasma;
      ob.etaP = m_cfg.carreauYasuda.etaInf - ow.etaPlasma;
      if (!(ob.etaP > 0.0))
        throw std::runtime_error("owens: eta_CY(inf) - eta_plasma must be positive.");
      ob.lambda = ow.lambdaH;
      ob.pttEpsilon = ow.epsilon;
      ob.lambda0Factor = 1.0;
      ob.lambda0Min = 0.0;
      if (isRoot() && m_cfg.deborah > 0.0)
      {
        Alert::Warning() << "[fluid] owens: De is ignored; lambda_H = " << ow.lambdaH
                         << " s is the single-cell time." << Alert::Raise;
      }
      m_cfg.deborah = 0.0;
      if (!(ow.gdCrit > 0.0 && ow.gdMax > ow.gdCrit && ow.a0 >= 0.0))
        throw std::runtime_error("owens: need 0 < gd_crit < gd_max and a0 >= 0.");
      // The peak of b(gd) (bFloor), on a log grid up to gd_max.
      m_owensBPeak = 0.0;
      m_owensGdPeak = std::numeric_limits<Real>::infinity();
      for (Real gd = 1.0e-2; gd <= ow.gdMax; gd *= 1.02)
      {
        const Real b = owensBRaw(gd);
        if (b > m_owensBPeak)
        {
          m_owensBPeak = b;
          m_owensGdPeak = gd;
        }
      }
    }
    if (m_cfg.fluid == "cy")
    {
      // Carreau-Yasuda as "gsptt" with a vanishing polymer (file header).
      // The polymer must stay passive: with De = 1e-3 (lambda ~ 1 ms) and
      // epsilon = 0 it was an Oldroyd-B with lambda eps_dot > 1/2 in the
      // throat at peak systole, psi grew without bound (run 11: psi ~ 27 at
      // t = 1.22 s) and s exp(psi) reached the momentum. De = 1e-5 keeps
      // lambda gd < 0.1 up to gd = 1e4 1/s, and epsilon = 0.5 bounds psi.
      ob.etaP = 1.0e-6 * m_cfg.carreauYasuda.etaInf;
      m_cfg.deborah = 1.0e-5;
      ob.pttEpsilon = 0.5;
      ob.lambda0Factor = 1.0;
      ob.lambda0Min = 0.0;
    }
    if (m_cfg.fluid == "newtonian")
    {
      // Constant viscosity eta_ref (eta_inf by default). The polymer is the
      // passive one of "cy" (unused by the two-field system).
      ob.etaS = m_cfg.etaRef > 0.0 ? m_cfg.etaRef : m_cfg.carreauYasuda.etaInf;
      ob.etaP = 1.0e-6 * ob.etaS;
      m_cfg.deborah = 1.0e-5;
      ob.pttEpsilon = 0.5;
      ob.lambda0Factor = 1.0;
      ob.lambda0Min = 0.0;
    }
    if (twoField() && m_cfg.formulation != "psi_rb")
    {
      throw std::runtime_error(
        "the two-field system (cy, newtonian) is written for formulation psi_rb.");
    }
    if (m_cfg.fluid == "gsptt" && (ob.lambda0Factor != 1.0 || ob.lambda0Min > 0.0))
    {
      throw std::runtime_error(
        "gsptt requires lambda0 = lambda (lambda0_factor 1, lambda0_min 0).");
    }
    if (m_cfg.owens.equilibrium && !isOwens())
      throw std::runtime_error("owens_equilibrium is an option of fluid owens.");
    const Real eta0 = referenceViscosity();

    if (!(m_cfg.period > 0.0))
      throw std::runtime_error("The cardiac period must be positive.");
    if (!(m_cfg.meshScale > 0.0) || !(m_cfg.referenceLength > 0.0) ||
      !(m_cfg.referenceVelocity > 0.0))
    {
      throw std::runtime_error(
        "mesh_scale, reference_length and reference_velocity must be positive.");
    }
    if (m_cfg.useVentricularModel && m_cfg.ventricularFeedback &&
      !(m_cfg.lv.Kat > 0.0 && m_cfg.lv.Kp > 0.0))
    {
      throw std::runtime_error(
        "The coupled mitral valve needs positive conductances Kat and Kp.");
    }
    if (!m_cfg.useVentricularModel && !(m_cfg.prescribedOutletResistance >= 0.0))
      throw std::runtime_error("The prescribed outlet resistance must not be negative.");
    m_meanVelocity = m_cfg.referenceVelocity;
    m_period = m_cfg.period;
    if (m_cfg.timeStep > 0.0)
    {
      m_cfg.stepsPerCycle =
        std::max(1, static_cast<int>(std::lround(m_period / m_cfg.timeStep)));
    }
    if (!(m_cfg.stepsPerCycle > 0))
      throw std::runtime_error("stepsPerCycle must be positive.");
    m_dt = m_period / static_cast<Real>(m_cfg.stepsPerCycle);
    m_dtNominal = m_dt;
    if (m_cfg.timeOrder != 1 && m_cfg.timeOrder != 2)
      throw std::runtime_error("timeOrder must be 1 (BDF1) or 2 (BDF2).");

    // De = lambda/T wins over lambda.
    if (m_cfg.deborah > 0.0)
      ob.lambda = m_cfg.deborah * m_period;
    if (!(ob.lambda > 0.0))
      throw std::runtime_error("lambda must be positive (set lambda or De).");

    if (m_cfg.fluid != "sptt" && m_cfg.fluid != "newtonian")
    {
      if (!isOwens())
      {
        // The solvent must be positive AND its shear stress eta_s(gd) gd
        // increasing: the solvent carries the instantaneous response, and a
        // falling solvent stress is a negative differential viscosity, an
        // ill-posed (Hadamard-unstable) momentum equation whose blow-up is
        // faster on finer meshes (R with eta_p 0.0104, lambda 0.1 s: eta_s
        // down to 1.6e-5 Pa s and falling stress for 1.75 < gd < 16 1/s).
        Real gdPrev = 0.0, stressPrev = 0.0;
        for (Real gd = 1.0e-3; gd < 1.0e5; gd *= 1.02)
        {
          const Real es = solventViscosity(gd);
          if (!(es > 0.0))
          {
            throw std::runtime_error("gsptt: eta_CY(gd) - eta_p/f(gd) <= 0 at gd = " +
              std::to_string(gd) + " 1/s; lower eta_p or epsilon.");
          }
          if (gdPrev > 0.0 && !(es * gd > stressPrev))
          {
            throw std::runtime_error(
              "gsptt: the solvent shear stress (eta_CY - eta_p/f) gd "
              "decreases between gd = " +
              std::to_string(gdPrev) + " and " + std::to_string(gd) +
              " 1/s (negative differential viscosity, ill-posed); lower eta_p (at lambda "
              "0.1 s, "
              "epsilon 0.2: eta_p <= 7.2e-3 Pa s).");
          }
          gdPrev = gd;
          stressPrev = es * gd;
        }
      }
    }

    if (isRoot())
    {
      // Groups at the reference scales (U, L): orders of magnitude only, the
      // flow sets its own speed. For Owens lambda is lambda_H, the high-shear
      // value; its range is in the table below.
      const Real wi = ob.lambda * m_meanVelocity / D;
      const Real de = ob.lambda / m_period;
      const Real re = rho * m_meanVelocity * D / eta0;
      Alert::Info() << "[groups] at U=" << m_meanVelocity << " m/s, L=" << D
                    << " m: Re=" << re << "  Wo="
                    << (0.5 * D * std::sqrt(2.0 * M_PI * rho / (m_period * eta0)))
                    << "  Wi=lambda U/L=" << wi << "  De=lambda/T=" << de
                    << "  El=Wi/Re=" << (wi / re)
                    << "  beta=" << (ob.etaS / (ob.etaS + ob.etaP))
                    << "  epsilon=" << ob.pttEpsilon << "  fluid=" << m_cfg.fluid
                    << (m_cfg.owens.equilibrium ? " (equilibrium structure)" : "")
                    << "  eta_ref=" << eta0 << Alert::Raise;
      if (m_cfg.fluid != "sptt" && m_cfg.fluid != "newtonian")
      {
        const auto& cy = m_cfg.carreauYasuda;
        Alert::Info() << "[fluid] Carreau-Yasuda eta_0=" << cy.etaZero
                      << " eta_inf=" << cy.etaInf << " lambda=" << cy.lambda
                      << " n=" << cy.n << " a=" << cy.a << "  eta_p=" << ob.etaP
                      << "  solvent eta_s(gd) at 1, 10, 100, 1000 1/s = "
                      << solventViscosity(1.0) << ", " << solventViscosity(10.0) << ", "
                      << solventViscosity(100.0) << ", " << solventViscosity(1000.0)
                      << " Pa s" << Alert::Raise;
        if (isOwens())
        {
          const auto& ow = m_cfg.owens;
          const Real G = ob.etaP / ob.lambda;
          Alert::Info() << "[fluid] Owens: lambda_H=" << ow.lambdaH
                        << " s  eta_plasma=" << ob.etaS << "  eta_p,inf=" << ob.etaP
                        << " Pa s  G=" << G << " Pa  epsilon=" << ob.pttEpsilon
                        << "  aN0: a0=" << ow.a0 << " ac=" << owensA(ow.gdCrit)
                        << " 1/s at gd_c=" << ow.gdCrit << ", 0 above gd_max=" << ow.gdMax
                        << " 1/s  b floor=" << (ow.bFloor ? m_owensBPeak : 0.0)
                        << " 1/s above " << m_owensGdPeak << " 1/s  n(t=0)="
                        << (ow.nInitial > 0.0       ? ow.nInitial
                               : ow.nInitial == 0.0 ? inletAggregateSize()
                                                    : owensNst(0.0))
                        << "  n_in=" << inletAggregateSize()
                        << (ow.equilibrium ? "  equilibrium structure, n = n_st(gd)" : "")
                        << Alert::Raise;
          std::ostringstream table;
          table << "[fluid] Owens steady shear, gd: n_st  mu_st (s)  G mu_st/(eta_CY - "
                   "eta_s)  "
                   "relaxation time of n 2/(b (2 n_st - 1)) (s) |";
          for (const Real gd : {0.0, 1.0, 10.0, 100.0, 300.0, 1000.0, 5000.0})
          {
            const Real nst = owensNst(gd);
            const Real lam = owensLambda(nst, gd);
            const Real b = owensB(gd);
            table << "  " << gd << ": " << nst << "  " << lam << "  "
                  << (G * lam / std::max<Real>(carreauYasuda(gd) - ob.etaS, 1.0e-300))
                  << "  "
                  << (b > 0.0 ? 2.0 / (b * (2.0 * nst - 1.0))
                              : std::numeric_limits<Real>::infinity())
                  << " |";
          }
          Alert::Info() << table.str() << Alert::Raise;
          // The steady Weissenberg number lambda_st gd, and that of blood that
          // keeps the structure of rest when sheared (n = n_st(0)): the range
          // spanned by the thixotropic fluid.
          std::ostringstream wiTable;
          wiTable << "[fluid] Owens Wi = lambda gd, gd: at n_st(gd)  at n_st(0) |";
          for (const Real gd : {1.0, 10.0, 50.0, 100.0, 300.0})
          {
            wiTable << "  " << gd << ": " << (owensLambda(owensNst(gd), gd) * gd) << "  "
                    << (owensLambda(owensNst(0.0), gd) * gd) << " |";
          }
          Alert::Info() << wiTable.str() << Alert::Raise;
        }
      }
      Alert::Info() << "[openings] p_pv = " << m_cfg.inletPressureMean << " + "
                    << m_cfg.inletPressureAmplitude << " sin(2 pi t/T) Pa, ramped over "
                    << m_cfg.rampCycles << " cycle(s);  mitral plane: "
                    << (m_cfg.useVentricularModel
                           ? (m_cfg.ventricularFeedback
                                 ? "lumped ventricle, coupled valve"
                                 : "lumped ventricle, one-way")
                           : "prescribed pressure and resistance")
                    << ";  T=" << m_period << " s  lambda=" << ob.lambda
                    << " s  lambda0=" << lambda0() << " s  dt=" << m_dt << " s  BDF"
                    << m_cfg.timeOrder << Alert::Raise;
    }
  }

  AtriumViscoelastic::MeshType AtriumViscoelastic::makeParentMesh(
    const Context::MPI& context, const Config& cfg)
  {
    // The combined mesh, loaded on the root rank, partitioned once and
    // distributed, as in Atrium::makeMesh.
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3 ||
        mesh.getCellCount() == 0)
      {
        throw std::runtime_error(
          "AtriumViscoelastic expects a nonempty 3D MEDIT volume mesh: " + cfg.meshPath +
          " (a Git LFS pointer?).");
      }

      for (size_t d = 1; d <= 3; ++d)
      {
        mesh.getConnectivity().compute(d, 0);
        mesh.getConnectivity().compute(d, d);
        mesh.getConnectivity().compute(d, d - 1);
        mesh.getConnectivity().compute(d - 1, d);
      }

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
    for (size_t d = 1; d <= 3; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(d, d);
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
    }
    // Faces, then edges, as in Atrium::makeMesh.
    mesh.reconcile(2);
    mesh.reconcile(1);

    return mesh;
  }

  AtriumViscoelastic::FluidMeshType AtriumViscoelastic::makeFluidMesh(
    const MeshType& parent, const Config& cfg)
  {
    // The owned cells of the fluid material, as in Atrium/Fluid: the submesh
    // keeps the partition and the boundary attributes of the parent. A mesh
    // with no cell of that material is a mesh of the fluid alone (Coarse.mesh,
    // the volume of the PTAG atrium) and is used whole.
    const auto& comm = parent.getContext().getCommunicator();
    const Attribute material = cfg.labels.fluidMaterial;
    size_t owned = 0, ownedFluid = 0;
    for (auto cell = parent.getCell(); cell; ++cell)
    {
      if (!parent.getShard().isOwned(3, cell->getIndex()))
        continue;
      ++owned;
      if (cell->getAttribute() == material)
        ++ownedFluid;
    }
    const size_t total = boost::mpi::all_reduce(comm, owned, std::plus<size_t>());
    const size_t fluid = boost::mpi::all_reduce(comm, ownedFluid, std::plus<size_t>());
    if (total == 0)
      throw std::runtime_error("No cells in " + cfg.meshPath + ".");
    const bool whole = fluid == 0;
    if (comm.rank() == RootRank)
    {
      if (whole)
      {
        Alert::Info() << "[mesh] no cell of material " << material
                      << ": the whole mesh is the fluid (" << total << " cells)"
                      << Alert::Raise;
      }
      else
      {
        Alert::Info() << "[mesh] fluid material " << material << ": " << fluid << " of "
                      << total << " cells" << Alert::Raise;
      }
    }
    FluidMeshType::Builder builder;
    builder.initialize(parent);
    for (auto cell = parent.getCell(); cell; ++cell)
    {
      if (parent.getShard().isOwned(3, cell->getIndex()) &&
        (whole || cell->getAttribute() == material))
        builder.include(3, cell->getIndex());
    }
    auto mesh = builder.finalize();
    size_t invalid = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getGeometry() != Geometry::Polytope::Type::Tetrahedron &&
        mesh.getShard().isOwned(3, cell->getIndex()))
        ++invalid;
    }
    invalid = boost::mpi::all_reduce(
      parent.getContext().getCommunicator(), invalid, std::plus<size_t>());
    if (invalid)
      throw std::runtime_error("The fluid region must be tetrahedral (P1 spaces).");
    mesh.scale(cfg.meshScale);
    for (size_t d = 1; d <= 3; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
    }
    return mesh;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  AtriumViscoelastic::Real AtriumViscoelastic::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = solventViscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::referenceViscosity() const
  {
    if (m_cfg.etaRef > 0.0)
      return m_cfg.etaRef;
    if (m_cfg.fluid == "sptt")
      return m_cfg.oldroydB.etaS + m_cfg.oldroydB.etaP;
    if (m_cfg.fluid == "newtonian")
      return m_cfg.oldroydB.etaS;
    return m_cfg.carreauYasuda.etaInf;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::carreauYasuda(Real gammaDot) const
  {
    const auto& cy = m_cfg.carreauYasuda;
    return cy.etaInf +
      (cy.etaZero - cy.etaInf) *
      std::pow(1.0 + std::pow(cy.lambda * gammaDot, cy.a), (cy.n - 1.0) / cy.a);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::steadyPTTFactor(Real gammaDot) const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real lam = steadyLambda(gammaDot);
    const Real c = 2.0 * ob.pttEpsilon * lam * lam * gammaDot * gammaDot;
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

  AtriumViscoelastic::Real AtriumViscoelastic::solventViscosity(Real gammaDot) const
  {
    if (constantSolvent())
      return m_cfg.oldroydB.etaS;
    const Real etaCY = carreauYasuda(gammaDot);
    if (m_cfg.fluid == "cy")
      return etaCY;
    return etaCY - m_cfg.oldroydB.etaP / steadyPTTFactor(gammaDot);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::shearRateAt(const Point& p) const
  {
    // At u^n (BDF1) or u* (BDF2).
    const auto J = Jacobian(m_uAdv).getValue(p);
    Real s2 = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const Real e = 0.5 * (J(i, j) + J(j, i));
        s2 += e * e;
      }
    }
    return std::sqrt(2.0 * s2);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::solventViscosityAt(const Point& p) const
  {
    if (constantSolvent())
      return m_cfg.oldroydB.etaS;
    return solventViscosity(shearRateAt(p));
  }

  AtriumViscoelastic::Real AtriumViscoelastic::pttFactor(Real traceExp) const
  {
    const auto& ob = m_cfg.oldroydB;
    if (ob.pttEpsilon == 0.0)
      return 1.0;
    return 1.0 + ob.pttEpsilon * (ob.lambda / lambda0()) * (traceExp - 3.0);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::alpha3At(const Point& p) const
  {
    const auto& ob = m_cfg.oldroydB;
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real speed = std::sqrt(Math::dot(uc, uc));
    const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
    const Real f = (isOwens() ? ob.lambda / m_lambdaCell.getValue(p) : 1.0) *
      pttFactor(
        logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp().trace());
    return 1.0 /
      (4.0 * f / (2.0 * ob.etaP) +
        0.25 *
          (ob.lambda * speed / (2.0 * ob.etaP * h) + ob.lambda * gradNorm / ob.etaP));
  }

  AtriumViscoelastic::Real AtriumViscoelastic::steadyLambda(Real gammaDot) const
  {
    const auto& ob = m_cfg.oldroydB;
    if (!isOwens())
      return ob.lambda;
    // G lambda_st = eta_CY - eta_plasma: the steady shear viscosity of a
    // Hookean mode is eta_p = G lambda whatever the rate.
    const Real G = ob.etaP / ob.lambda;
    return std::max<Real>(carreauYasuda(gammaDot) - ob.etaS, 1.0e-12 * ob.etaP) / G;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::owensA(Real gammaDot) const
  {
    // C1 cubics with zero slope at 0, gd_c and gd_max (Config::Owens).
    const auto& ow = m_cfg.owens;
    const Real gd = std::abs(gammaDot);
    const Real ac = ow.aCrit >= 0.0 ? ow.aCrit : 2.0 * ow.a0;
    const auto h = [](Real s) { return s * s * (3.0 - 2.0 * s); };
    if (gd <= ow.gdCrit)
      return ow.a0 + (ac - ow.a0) * h(gd / ow.gdCrit);
    if (gd <= ow.gdMax)
      return ac * (1.0 - h((gd - ow.gdCrit) / (ow.gdMax - ow.gdCrit)));
    return 0.0;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::owensBRaw(Real gammaDot) const
  {
    // Eq. (28), b = aN0/(n_st (n_st - 1)).
    const Real a = owensA(gammaDot);
    const Real nst = owensNst(gammaDot);
    return (a > 0.0 && nst > 1.0) ? a / (nst * (nst - 1.0)) : 0.0;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::owensB(Real gammaDot) const
  {
    const Real b = owensBRaw(gammaDot);
    if (m_cfg.owens.bFloor && std::abs(gammaDot) >= m_owensGdPeak)
      return std::max(b, m_owensBPeak);
    return b;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::owensNst(Real gammaDot) const
  {
    // Eq. (27) with eta(gd) - eta_s for the polymer: n_st >= 1.
    const auto& ob = m_cfg.oldroydB;
    const Real r = std::max<Real>((carreauYasuda(gammaDot) - ob.etaS) / ob.etaP, 1.0);
    return r * (1.0 + 1.5 * owensA(gammaDot) * m_cfg.owens.lambdaH);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::owensLambda(Real n, Real gammaDot) const
  {
    // Eq. (23): mu = n lambda_H/(1 + g_n n lambda_H), g_n n = b n (n - 1)/2 + aN0.
    const auto& ow = m_cfg.owens;
    n = std::max<Real>(n, 1.0);
    const Real gn = 0.5 * owensB(gammaDot) * n * (n - 1.0) + owensA(gammaDot);
    return n * ow.lambdaH / (1.0 + gn * ow.lambdaH);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::shearRateOf(
    const VectorGridFunctionType& u, const Point& p) const
  {
    const auto J = Jacobian(u).getValue(p);
    Real s2 = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const Real e = 0.5 * (J(i, j) + J(j, i));
        s2 += e * e;
      }
    }
    return std::sqrt(2.0 * s2);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::inletAggregateSize() const
  {
    // Owens::nInlet, or the steady value at the venous shear rate.
    if (m_cfg.owens.nInlet > 0.0)
      return m_cfg.owens.nInlet;
    return owensNst(m_cfg.inletShearRate);
  }

  void AtriumViscoelastic::updateStructureFields()
  {
    // lambda(n^n, gd(u^n)) by cell; with Owens::equilibrium n = n_st(gd), so
    // lambda = (eta_CY(gd) - eta_plasma)/G.
    if (m_cfg.owens.equilibrium)
    {
      m_lambdaCell.project(RealFunction(
        [this](const Point& p) { return steadyLambda(shearRateOf(m_uOld, p)); }));
    }
    else
    {
      m_lambdaCell.project(RealFunction([this](const Point& p) {
        return owensLambda(m_nOld.getValue(p), shearRateOf(m_uOld, p));
      }));
    }
    m_nMin = m_nOld.min();
    m_nMax = m_nOld.max();
    m_lambdaMax = m_lambdaCell.max();
  }

  void AtriumViscoelastic::setupStructure()
  {
    // D n/Dt = -(1/2) b (n - n_st)(n + n_st - 1): BDF with the mass step of the
    // flow (n_hist as x_hist in prepareStep), advected by u^{n+1}, the reaction
    // linearised about n^n, c = (1/2) b(gd) (n^n + n_st - 1) > 0:
    //   (n - n_hist)/dtM + u.grad n + c n = c n_st,
    // SUPG with tau = [(2/dtM)^2 + (2 |u|/h)^2 + c^2]^-1/2. Pure advection: the
    // entering blood carries n_in (veins) or n^n (mitral plane), imposed by
    // the upwind term |u.n|_- (n - n_in) on the openings.
    // The coefficients read dt and the fields when assembled, so the form is
    // built once. It is solved after the flow step is accepted, with u^{n+1}
    // in m_uOld and n^n in m_nOld.
    const auto reaction = [this](const Point& p) -> Real {
      const Real gd = shearRateOf(m_uOld, p);
      return 0.5 * owensB(gd) * (m_nOld.getValue(p) + owensNst(gd) - 1.0);
    };
    const auto tau = [this, reaction](const Point& p) -> Real {
      const auto uc = m_uOld.getValue(p);
      const Real a = 2.0 / massStep();
      const Real b = 2.0 * std::sqrt(Math::dot(uc, uc)) / cellSize(p);
      const Real c = reaction(p);
      return 1.0 / std::sqrt(a * a + b * b + c * c);
    };
    const auto known = [this, reaction](const Point& p) -> Real {
      return m_nHist.getValue(p) / massStep() +
        reaction(p) * owensNst(shearRateOf(m_uOld, p));
    };
    RealFunction massFn = [this](const Point&) { return 1.0 / massStep(); };
    RealFunction reactionFn = [reaction](const Point& p) { return reaction(p); };
    RealFunction tauFn = [tau](const Point& p) { return tau(p); };
    RealFunction tauMassFn = [this, tau](const Point& p) { return tau(p) / massStep(); };
    RealFunction tauReactionFn = [tau, reaction](
                                   const Point& p) { return tau(p) * reaction(p); };
    RealFunction knownFn = [known](const Point& p) { return known(p); };
    RealFunction tauKnownFn = [tau, known](const Point& p) { return tau(p) * known(p); };
    // Entering blood carries n_in on the veins and n^n on the mitral plane.
    const auto normal = BoundaryNormal(m_mesh);
    const auto inflow = Max(-Dot(m_uOld, normal), 0.0);
    RealFunction inletN = [this](const Point&) { return inletAggregateSize(); };

    const auto adv = Dot(m_uOld, Grad(m_n));
    const auto advTest = Dot(m_uOld, Grad(m_nTest));
    m_structure = Integral(massFn * m_n, m_nTest) + Integral(adv, m_nTest) +
      Integral(reactionFn * m_n, m_nTest) - Integral(knownFn, m_nTest) +
      Integral(tauMassFn * m_n, advTest) + Integral(tauFn * adv, advTest) +
      Integral(tauReactionFn * m_n, advTest) - Integral(tauKnownFn, advTest) +
      BoundaryIntegral(inflow * m_n, m_nTest).over(m_inletSet) -
      BoundaryIntegral(inflow * inletN, m_nTest).over(m_inletSet) +
      BoundaryIntegral(inflow * m_n, m_nTest).over(m_cfg.labels.outlet) -
      BoundaryIntegral(inflow * m_nOld, m_nTest).over(m_cfg.labels.outlet);
  }

  void AtriumViscoelastic::solveStructure()
  {
    if (m_cfg.owens.equilibrium)
    {
      // n = n_st(gd(u^{n+1})), for the output; lambda follows in the update.
      m_nOld.project(RealFunction(
        [this](const Point& p) { return owensNst(shearRateOf(m_uOld, p)); }));
      m_nNm1.setData(m_nOld.getData());
      updateStructureFields();
      return;
    }
    // n_hist as x_hist in prepareStep, with this step's dt ratio.
    if (m_cfg.timeOrder == 2)
    {
      const Real w = bdfRatio();
      const Real a = (1.0 + w) * (1.0 + w) / (1.0 + 2.0 * w);
      const Real b = -w * w / (1.0 + 2.0 * w);
      combine(a, m_nOld.getData(), b, m_nNm1.getData(), m_nHist.getData());
    }
    else
      m_nHist.setData(m_nOld.getData());

    m_structure.assemble();
    m_structure.solve(m_structureKSP);

    // SUPG is not monotone: n stays within [1, n_st(0)], the single cell and
    // the rest state.
    {
      ::Vec v = m_n.getSolution().getData();
      PetscInt size = 0;
      PetscScalar* data = nullptr;
      PetscErrorCode ierr = VecGetLocalSize(v, &size);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecGetArray(v, &data);
      assert(ierr == PETSC_SUCCESS);
      const Real nMax = owensNst(0.0);
      for (PetscInt i = 0; i < size; ++i)
        data[i] = std::isfinite(data[i]) ? std::clamp<Real>(data[i], 1.0, nMax) : 1.0;
      ierr = VecRestoreArray(v, &data);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
    m_nNm1.setData(m_nOld.getData());
    m_nOld.setData(m_n.getSolution().getData());
    updateStructureFields();
  }

  AtriumViscoelastic::Real AtriumViscoelastic::lambda0() const
  {
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = std::max(ob.lambda0Factor * ob.lambda, ob.lambda0Min);
    if (!(l0 > 0.0))
      throw std::runtime_error("lambda0 = max(k lambda, lambda0_min) must be positive.");
    return l0;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::ramp(Real t) const
  {
    const Real tr = m_cfg.rampCycles * m_period;
    if (!(tr > 0.0) || t >= tr)
      return 1.0;
    return 0.5 * (1.0 - std::cos(M_PI * t / tr));
  }

  void AtriumViscoelastic::updateConformation()
  {
    if (twoField())
    {
      const auto zero = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
      m_sigmaD = zero;
      m_sigmaO = zero;
      return;
    }
    const Real s = m_cfg.oldroydB.etaP / lambda0();
    const auto field = [](auto f) { return VectorFunction(size_t(3), f); };
    const auto sigmaAt = [this, s](const Point& p) -> Eigen::Matrix3d {
      return s *
        (logEig(unpack(m_psiItD.getValue(p), m_psiItO.getValue(p))).exp() -
          Eigen::Matrix3d::Identity());
    };

    m_sigmaD.project(field(
      [&](const Point& p) -> Math::SpatialVector<Real> { return packD(sigmaAt(p)); }));
    m_sigmaO.project(field(
      [&](const Point& p) -> Math::SpatialVector<Real> { return packO(sigmaAt(p)); }));
  }

  void AtriumViscoelastic::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  void AtriumViscoelastic::combine(
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

  AtriumViscoelastic::Real AtriumViscoelastic::bdfRatio() const
  {
    // The first step takes x^{-1} = x^0, for which w does not enter the
    // history; w = 1 keeps its mass step at 2 dt/3.
    if (m_cfg.timeOrder != 2 || m_history < 1 || !(m_dtPrev > 0.0))
      return 1.0;
    return m_dt / m_dtPrev;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::massStep() const
  {
    if (m_cfg.timeOrder != 2)
      return m_dt;
    // (1 + 2 w)/(1 + w) x^{n+1} - (1 + w) x^n + w^2/(1 + w) x^{n-1}, over dt,
    // is (x^{n+1} - x_hist)/dtM with dtM = dt (1 + w)/(1 + 2 w). At w = 1 the
    // expression is evaluated as 2 dt/3, as before the variable step.
    const Real w = bdfRatio();
    return (1.0 + w) * m_dt / (1.0 + 2.0 * w);
  }

  void AtriumViscoelastic::prepareStep()
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

  void AtriumViscoelastic::advanceHistory()
  {
    ++m_stepRevision; // the step data (psi^n, u^n) change
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

  void AtriumViscoelastic::rejectStep()
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

  AtriumViscoelastic& AtriumViscoelastic::initialize()
  {
    setupSpaces();
    setupFlow();
    if (isOwens())
      setupStructure();
    if (m_cfg.residenceTime)
      setupResidence();
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

  void AtriumViscoelastic::setupSpaces()
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
    m_piRD.get() = zero;
    m_piRO.get() = zero;
    m_pttN = Real(1);
    m_alpha3Dyn = Real(0);
    if (isOwens())
    {
      const auto& ow = m_cfg.owens;
      m_nOld = Real(ow.nInitial > 0.0 ? ow.nInitial
          : ow.nInitial == 0.0        ? inletAggregateSize()
                                      : owensNst(0.0));
      m_nNm1.setData(m_nOld.getData());
      m_nHist.setData(m_nOld.getData());
      m_n.getSolution().setData(m_nOld.getData());
      updateStructureFields();
    }
    else
    {
      m_nOld = Real(1);
      m_nNm1 = Real(1);
      m_nHist = Real(1);
      m_lambdaCell = Real(m_cfg.oldroydB.lambda);
      m_lambdaMax = m_cfg.oldroydB.lambda;
    }
    m_aSum = Real(0);
    m_bSum = Real(0);
    m_rateWork = Real(0);
    m_wiCell = Real(0);
    m_ageOld = Real(0);
    m_ageNm1 = Real(0);
    m_ageHist = Real(0);
    m_age.getSolution() = Real(0);
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
    m_tawssAppendage = Real(0);
    m_osiAppendage = Real(0);
    m_one = Real(1);

    configureMassSolver(m_vectorProjectionKSP, "av_vproj_");
    configureMassSolver(m_wssKSP, "av_wss_");
    // The aggregate-size transport is not symmetric: GMRES + block Jacobi (ILU).
    setPrefixedDefault("av_struct_", "ksp_type", "gmres");
    setPrefixedDefault("av_struct_", "pc_type", "bjacobi");
    setPrefixedDefault("av_struct_", "ksp_rtol", "1e-10");
    m_structureKSP.setPrefix("av_struct_");

    updateConformation(); // psi = 0: tau = I, sigma = 0

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
    if (isOwens())
    {
      m_nOld.setName("aggregateSize");
      m_xdmf.add("aggregateSize", m_nOld);
    }
    if (m_cfg.residenceTime)
    {
      m_ageOld.setName("residenceTime");
      m_xdmf.add("residenceTime", m_ageOld);
      // The age transport is not symmetric: GMRES + block Jacobi (ILU).
      setPrefixedDefault("av_age_", "ksp_type", "gmres");
      setPrefixedDefault("av_age_", "pc_type", "bjacobi");
      setPrefixedDefault("av_age_", "ksp_rtol", "1e-10");
      m_ageKSP.setPrefix("av_age_");
    }

    // Areas of the openings, volume of the fluid and of the appendage.
    m_inletMeasure = boundaryMeasure(m_inletSet);
    m_outletMeasure = boundaryMeasure(AttributeSet{m_cfg.labels.outlet});
    m_flux = Integral(m_one, m_qFlux);
    m_flux.assemble();
    m_volume = std::max<Real>(m_flux(m_one), 1e-30);
    m_appendageVolume =
      integrate([this](const Point& p) -> Real { return inAppendage(p) ? 1.0 : 0.0; });
    if (!(m_appendageVolume > 1.0e-9 * m_volume))
    {
      if (m_cfg.appendage.radius > 0.0 && isRoot())
      {
        Alert::Warning()
          << "[appendage] no fluid inside the ball: its centre is in metres, "
             "the coordinates of the XDMF output."
          << Alert::Raise;
      }
      m_appendageVolume = 0.0;
    }

    // The ventricle starts from the atrial pressure of the fluid at rest
    // (Atrium/Fluid): zero under a pressure ramp, the venous mean otherwise.
    m_outletMeanPressure = m_cfg.rampCycles > 0.0 ? Real(0) : m_cfg.inletPressureMean;
    {
      VentricleType::State s0;
      s0.t = 0.0;
      s0.y = m_cfg.lv.initialY;
      s0.v = m_cfg.lv.initialV;
      s0.pv = m_input.pAt(0.0) + m_cfg.lv.initialPvOffset;
      s0.par = m_cfg.lv.initialPar;
      s0.pd = m_cfg.lv.initialPd;
      s0.w = m_input.m0(s0.ec);
      m_model.setMaxIterations(m_cfg.lv.maxIterations)
        .setAbsoluteTolerance(m_cfg.lv.absoluteTolerance)
        .setRelativeTolerance(m_cfg.lv.relativeTolerance)
        .setStepTolerance(m_cfg.lv.stepTolerance)
        .setDampingFactor(m_cfg.lv.dampingFactor);
      m_model.initialize(s0);
      m_lvTime0 = 0.0;
      m_lvTime1 = 0.0;
      m_lvPressure0 = m_model.getState().pv;
      m_lvPressure1 = m_lvPressure0;
    }
    updateBoundaryData(0.0);

    // Nodal wall/appendage indicators (see m_wallArea). Interpolating a P1 field
    // at its own node leaves O(1e-21) round-off where it is zero, hence the
    // thresholds relative to the maximum.
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(m_wallSet);
    m_flux.assemble();
    m_wallArea.setData(m_flux.getVector());
    m_wallAreaThreshold = 1.0e-6 * m_wallArea.max();
    m_flux = BoundaryIntegral(
      RealFunction([this](const Point& p) -> Real { return inAppendage(p) ? 1.0 : 0.0; }),
      m_qFlux)
               .over(m_wallSet);
    m_flux.assemble();
    m_appendageArea.setData(m_flux.getVector());
    m_appendageAreaThreshold = 1.0e-6 * m_appendageArea.max();

    // Config::wallKinematics: wall-node indicator (0/1, projected over the
    // cells) and the area-weighted nodal wall normal int_Gw n phi_i dA.
    m_wallFlag.project(RealFunction([this](const Point& p) -> Real {
      return m_wallArea.getValue(p) > m_wallAreaThreshold ? 1.0 : 0.0;
    }));
    {
      const auto normal = BoundaryNormal(m_mesh);
      ScalarGridFunctionType* comp[3] = {&m_wallN0, &m_wallN1, &m_wallN2};
      for (size_t c = 0; c < 3; ++c)
      {
        m_flux = BoundaryIntegral(Component(normal, c), m_qFlux).over(m_wallSet);
        m_flux.assemble();
        comp[c]->setData(m_flux.getVector());
      }
    }
    if (m_cfg.wallKinematics)
    {
      auto& wk = wallKinematics();
      wk.flag = &m_wallFlag;
      wk.normal = {&m_wallN0, &m_wallN1, &m_wallN2};
    }
    else
      wallKinematics() = WallKinematics{};
    if (isRoot())
    {
      Alert::Info() << "[wall] no-slip kinematics in the psi equation: "
                    << (m_cfg.wallKinematics ? "on" : "off") << Alert::Raise;
    }

    // Queried on every rank, printed on one.
    const auto fluidCells = m_mesh.getCellCount();
    if (isRoot())
    {
      Alert::Info() << "[boundary] venous area = " << m_inletMeasure
                    << " m^2  mitral area = " << m_outletMeasure
                    << " m^2  fluid volume = " << (1.0e6 * m_volume)
                    << " mL  appendage volume = " << (1.0e6 * m_appendageVolume)
                    << " mL  mean cell size (V/cells)^(1/3) = "
                    << (1.0e3 *
                         std::cbrt(
                           m_volume / static_cast<Real>(std::max<size_t>(fluidCells, 1))))
                    << " mm" << Alert::Raise;
    }
    if (!(m_inletMeasure > 1.0e-12) || !(m_outletMeasure > 1.0e-12))
    {
      throw std::runtime_error(
        "No venous or mitral faces: check the labels of the mesh "
        "(veins 1, 2, 3, 5; mitral plane 6; wall 4) and -av_mesh_scale.");
    }
    {
      // A boundary face with none of these labels would be a traction-free
      // opening: the area of such faces is reported.
      const Real wallMeasure = boundaryMeasure(m_wallSet);
      m_flux = BoundaryIntegral(m_one, m_qFlux);
      m_flux.assemble();
      const Real boundary = m_flux(m_one);
      const Real other = boundary - m_inletMeasure - m_outletMeasure - wallMeasure;
      if (isRoot())
      {
        Alert::Info() << "[boundary] wall area = " << wallMeasure
                      << " m^2  total boundary = " << boundary << " m^2" << Alert::Raise;
        if (std::abs(other) > 1.0e-6 * boundary)
        {
          Alert::Warning() << "[boundary] " << other
                           << " m^2 of the boundary carries none of the labels 1, 2, 3, "
                              "5, 6, 4: those faces are traction-free openings"
                           << Alert::Raise;
        }
      }
    }
    setupProbes();
  }

  AtriumViscoelastic::Real AtriumViscoelastic::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-18);
  }

  AtriumViscoelastic::Real AtriumViscoelastic::streamlineLength(
    const Point& p, const Math::SpatialVector<Real>& u)
  {
    const Real speed = std::sqrt(Math::dot(u, u));
    if (!(speed > 0.0))
      return cellSize(p);
    // The extent of the cell along u/|u|.
    Real lo = std::numeric_limits<Real>::max();
    Real hi = std::numeric_limits<Real>::lowest();
    for (auto it = p.getPolytope().getVertex(); it; ++it)
    {
      const auto x = it->getCoordinates();
      Real d = 0.0;
      for (Index c = 0; c < static_cast<Index>(u.size()); ++c)
        d += x(c) * u(c);
      d /= speed;
      lo = std::min(lo, d);
      hi = std::max(hi, d);
    }
    return hi > lo ? hi - lo : cellSize(p);
  }

  void AtriumViscoelastic::setupFlowExp()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    // BDF mass terms (x - x_hist)/dtM, for u in the momentum and for E in the
    // constitutive equation.
    const Real dtM = massStep();
    m_formDt = dt;
    m_formMassStep = dtM;
    const auto& ob = m_cfg.oldroydB;
    const auto etaS =
      PointwiseScalar([this](const Point& p) -> Real { return solventViscosityAt(p); },
        constantSolvent() ? size_t(0) : size_t(2));
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0; // sigma = s (E - I)
    const Real ka = ob.lambda / (2.0 * l0); // lambda/(2 lambda0)
    const Real kA = 1.0 / (2.0 * ob.etaP); // adjoint: reaction
    const Real kB = ob.lambda / (2.0 * ob.etaP); // adjoint: transport
    // Memory of the dynamic subscales (zero for quasi-static ones):
    // rho/dt in the velocity subscales, lambda/(2 eta_p dt) in the stress one.
    const Real memU = m_cfg.dynamicSubscales ? rho / dt : 0.0;
    const Real memS = m_cfg.dynamicSubscales ? kB / dt : 0.0;
    const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
    // E_hist = a E^n + b E^{n-1}, as x_hist in prepareStep.
    Real bdfA = 1.0, bdfB = 0.0;
    if (m_cfg.timeOrder == 2)
    {
      const Real w = bdfRatio();
      bdfA = (1.0 + w) * (1.0 + w) / (1.0 + 2.0 * w);
      bdfB = -w * w / (1.0 + 2.0 * w);
    }
    const LogModel model{
      ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda), bdfA, bdfB, dtM};
    constexpr Block D = Block::D;
    constexpr Block O = Block::O;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto convU = Mult(Jacobian(m_u), m_uAdv);
    const auto temam = Dot(Div(m_uAdv) * m_u, m_v);
    // Both kinds of opening may reverse: with the outward normal, the entering
    // flux carries the convective boundary term of either sign (Atrium/Fluid).
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;
    // Time-dependent data enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };

    // The stabilisation advects with u^n, the velocity of its parameters.
    const auto streamV = Mult(Jacobian(m_v), m_uOld);
    const auto uu = PointwiseTensor([this](const Point& p) -> Eigen::Matrix3d {
      const auto u = m_uOld.getValue(p);
      const Eigen::Vector3d w(u(0), u(1), u(2));
      return w * w.transpose();
    });
    const auto Ln = Jacobian(m_uOld);
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(3, 3));

    // Pointwise data at the Newton point (psi^k, u^k); E_hist from psi^n, psi^{n-1}.
    LogFields<VectorGridFunctionType> at{&m_psiItD, &m_psiItO, &m_uIt, model,
      &m_psiRevision, &m_psiOldD, &m_psiOldO, &m_psiNm1D, &m_psiNm1O};
    at.lambda = lambdaField();

    // Momentum: exp(psi) ~ exp(psi^k) + Dexp[psi - psi^k].
    const auto Tn = fromPoint(at, LogDexp, [](const LogPoint& q) { return q.exp; });
    const auto C = Tn - dexp<D>(at, m_psiItD) - dexp<O>(at, m_psiItO);

    // One quadrature rule for every integral that reads pointwise data (the
    // log data and the subscales): see setupFlow.
    const size_t psiOrder = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));
    const auto withOrder = [psiOrder](auto integral) {
      integral.setOrder(psiOrder);
      return integral;
    };

    // ---- Constitutive equation, exp form: R(psi, u) : chi = 0 --------------
    const auto RpsiD = expFormTrialPsi<D>(at, ka, m_psiD, m_uIt);
    const auto RpsiO = expFormTrialPsi<O>(at, ka, m_psiO, m_uIt);
    const auto Ru = expFormTrialU(at, ka, m_u);
    const auto Kexp = fromPoint(at, LogExpForm, [](const LogPoint& q) { return q.K; });
    const auto chiD = tensorFixed<D>(m_chiD);
    const auto chiO = tensorFixed<O>(m_chiO);

    // ---- Residual-based OSS of the constitutive equation --------------------
    // <s~, A(chi) + eps(v)>, s~ = alpha3' (memS s~^n - P[R]), P[R] = R - Pi[R^n].
    const auto AD = adjointTest<D>(m_chiD, m_pttN, kA, kB, m_uOld, Ln);
    const auto AO = adjointTest<O>(m_chiO, m_pttN, kA, kB, m_uOld, Ln);
    const auto stressKnown =
      PointwiseTensor([this, at, memS](const Point& p) -> Eigen::Matrix3d {
        const Eigen::Matrix3d K = logPointAt(at, p, LogExpForm).K;
        const Eigen::Matrix3d piR =
          unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p));
        const Real a3 = m_alpha3Dyn.getValue(p);
        return -a3 * (K - piR) + a3 * memS * subscales().at(p).sig;
      });
    const Real ma3 = -1.0;
    // Galerkin and OSS test functions of the constitutive equation, merged so
    // that each (trial, test) pair is assembled once: chi_B - alpha3' A(chi_B).
    const auto testD = chiD + ma3 * m_alpha3Dyn * AD;
    const auto testO = chiO + ma3 * m_alpha3Dyn * AO;
    // d_t E_h leaves the residual of the stress subscale (Config::
    // subscaleTimeResidual off): + alpha3' <k_a (E - E_hist)/dtM, A(chi) +
    // eps(v)>, E linearised as in the momentum, takes it back.
    const Real cT = m_cfg.subscaleTimeResidual ? 0.0 : 1.0;
    const auto rateTrial =
      (cT * ka / dtM) * m_alpha3Dyn * (dexp<D>(at, m_psiD) + dexp<O>(at, m_psiO));
    const auto rateKnown = cT * m_alpha3Dyn *
      fromPoint(at, LogExpForm, [](const LogPoint& q) { return q.Kt; });

    // ---- Subscales of step n --------------------------------------------------
    const auto u1n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u1; });
    const auto u2n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u2; });
    const auto u3n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u3; });
    const Real sd = m_cfg.stressDivScale;
    const Real ps = m_cfg.pressureScale;

    const auto body =
      // ---- Momentum, sigma = s (exp(psi) - I) --------------------------------
      (rho / dtM) * Integral(m_u, m_v) - (rho / dtM) * Integral(m_uHist, m_v) +
      rho * Integral(Dot(convU, m_v)) +
      0.5 * rho * Integral(temam)
      // (merged with the eps(v) part of the constitutive OSS, see below)
      + withOrder(Integral(2.0 * etaS * symU + ma3 * m_alpha3Dyn * Ru, symV)) +
      withOrder(Integral(s * dexp<D>(at, m_psiD) + ma3 * m_alpha3Dyn * RpsiD, symV)) +
      withOrder(Integral(s * dexp<O>(at, m_psiO) + ma3 * m_alpha3Dyn * RpsiO, symV)) +
      s * withOrder(Integral(C, symV)) - s * Integral(I, symV) -
      Integral(m_p, Div(m_v))

      // ---- Continuity -----------------------------------------------------
      + Integral(Div(m_u), m_q) +
      m_cfg.pressurePenalty * Integral(m_p, m_q)

      // ---- Constitutive, exp form, Newton, with its OSS -----------------------
      // <R, chi> + <s~, A(chi) + eps(v)>, s~ = alpha3' (memS s~^n - R + Pi[R^n]);
      // the eps(v) part of the bilinear terms is in the momentum integrals.
      + withOrder(Integral(RpsiD, testD)) + withOrder(Integral(RpsiO, testD)) +
      withOrder(Integral(RpsiD, testO)) + withOrder(Integral(RpsiO, testO)) +
      withOrder(Integral(Ru, testD)) + withOrder(Integral(Ru, testO)) +
      withOrder(Integral(Kexp, chiD)) + withOrder(Integral(Kexp, chiO)) +
      withOrder(Integral(stressKnown, AD)) + withOrder(Integral(stressKnown, AO)) +
      withOrder(Integral(stressKnown, symV)) + withOrder(Integral(rateTrial, AD)) +
      withOrder(Integral(rateTrial, AO)) + withOrder(Integral(rateTrial, symV)) +
      withOrder(Integral(rateKnown, AD)) + withOrder(Integral(rateKnown, AO)) +
      withOrder(Integral(rateKnown, symV))

      // ---- Momentum, split OSS, dynamic subscales ------------------------------
      // u1: <u1, -rho (grad v) u^n>, u1 = alpha1' (rho u1^n/dt - rho P[(grad u) u^n])
      + vms * rho * rho * Integral(m_tauK * Mult(Jacobian(m_u), uu), Jacobian(m_v)) -
      vms * rho * rho *
        withOrder(Integral(m_tauK * (m_piConv.get() + (memU / rho) * u1n), streamV))
      // u2: <u2, -grad q>, u2 = alpha1' (rho u2^n/dt - P[grad p])
      + ps * Integral(m_tauK * Grad(m_p), Grad(m_q)) -
      ps * withOrder(Integral(m_tauK * (m_piGradP.get() + memU * u2n), Grad(m_q)))
      // u3: <u3, div chi>, u3 = alpha1' (rho u3^n/dt + P[s div E])
      - sd *
        withOrder(
          Integral(m_tauK * (m_piDivSigma.get() - memU * u3n), divergence<D>(m_chiD))) -
      sd *
        withOrder(
          Integral(m_tauK * (m_piDivSigma.get() - memU * u3n), divergence<O>(m_chiO)))

      // ---- Continuity, OSS of div u (algebraic subscale) ----------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v)) -
      Integral(m_alpha2 * m_piDiv.get(), Div(m_v))

      // ---- Openings: tractions -p_pv n and -(p_out + Z u.n) n ---------------
      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet) +
      m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.labels.outlet) +
      m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet) +
      BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Wall: no slip --------------------------------------------------------
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // The psi equation is hyperbolic in psi: relaxed fluid (psi = 0) enters
    // through the veins, imposed by the upwind term |u*.n|_- psi : chi, which
    // vanishes while a vein is an outflow (Config::inletConformation).
    // E - I at the Newton point, minus its linear part: the known side of
    // k_a |u*.n|_- (E - I) : chi on the veins (E = I in the entering fluid).
    const auto inletKnown =
      PointwiseTensor([this, at](const Point& p) -> Eigen::Matrix3d {
        const LogPoint& q = logPointAt(at, p, LogDexp);
        const auto d = m_psiItD.getValue(p);
        const auto o = m_psiItO.getValue(p);
        const Real pk[6] = {d(0), d(1), d(2), o(0), o(1), o(2)};
        Eigen::Matrix3d out = q.exp - Eigen::Matrix3d::Identity();
        for (size_t m = 0; m < 6; ++m)
          out -= pk[m] * q.D[m];
        return out;
      });
    const auto inflowOut = Max(-Dot(m_uAdv, normal), 0.0);
    const auto finish = [&](const auto& form) {
      if (m_cfg.inletConformation)
      {
        m_flow = form +
          ka *
            BoundaryIntegral(inflowOut * Dot(dexp<D>(at, m_psiD), chiD))
              .over(m_inletSet) +
          ka *
            BoundaryIntegral(inflowOut * Dot(dexp<O>(at, m_psiO), chiD))
              .over(m_inletSet) +
          ka *
            BoundaryIntegral(inflowOut * Dot(dexp<D>(at, m_psiD), chiO))
              .over(m_inletSet) +
          ka *
            BoundaryIntegral(inflowOut * Dot(dexp<O>(at, m_psiO), chiO))
              .over(m_inletSet) +
          ka * BoundaryIntegral(inflowOut * Dot(chiD, inletKnown)).over(m_inletSet) +
          ka * BoundaryIntegral(inflowOut * Dot(chiO, inletKnown)).over(m_inletSet);
      }
      else
        m_flow = form;
    };
    // Outlet (Config::outletConformationInflow): where u*.n < 0, the upwind
    // term k_a |u*.n|_- (E - E^n) : chi, E linearised as above.
    const auto outlet = m_cfg.labels.outlet;
    const auto outletKnown =
      PointwiseTensor([this, at](const Point& p) -> Eigen::Matrix3d {
        const LogPoint& q = logPointAt(at, p, LogDexp);
        const auto d = m_psiItD.getValue(p);
        const auto o = m_psiItO.getValue(p);
        const Real pk[6] = {d(0), d(1), d(2), o(0), o(1), o(2)};
        Eigen::Matrix3d out =
          q.exp - logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp();
        for (size_t m = 0; m < 6; ++m)
          out -= pk[m] * q.D[m];
        return out;
      });
    const auto withOutlet = [&](const auto& form) {
      if (m_cfg.outletConformationInflow)
      {
        finish(form +
          ka * BoundaryIntegral(inflowOut * Dot(dexp<D>(at, m_psiD), chiD)).over(outlet) +
          ka * BoundaryIntegral(inflowOut * Dot(dexp<O>(at, m_psiO), chiD)).over(outlet) +
          ka * BoundaryIntegral(inflowOut * Dot(dexp<D>(at, m_psiD), chiO)).over(outlet) +
          ka * BoundaryIntegral(inflowOut * Dot(dexp<O>(at, m_psiO), chiO)).over(outlet) +
          ka * BoundaryIntegral(inflowOut * Dot(chiD, outletKnown)).over(outlet) +
          ka * BoundaryIntegral(inflowOut * Dot(chiO, outletKnown)).over(outlet));
      }
      else
        finish(form);
    };
    // u3 trial: s div(Dexp[psi]) (+ the D2exp part of the Newton linearisation
    // of div exp(psi)), each (trial, test) pair assembled once.
    const Real k1 = sd * s;
    if (m_cfg.newtonD2)
    {
      const auto dD = divExpTrial<D>(at, m_psiD) + divExpTrialD2<D>(at, m_psiD);
      const auto dO = divExpTrial<O>(at, m_psiO) + divExpTrialD2<O>(at, m_psiO);
      withOutlet(body + k1 * withOrder(Integral(m_tauK * dD, divergence<D>(m_chiD))) +
        k1 * withOrder(Integral(m_tauK * dO, divergence<D>(m_chiD))) +
        k1 * withOrder(Integral(m_tauK * dD, divergence<O>(m_chiO))) +
        k1 * withOrder(Integral(m_tauK * dO, divergence<O>(m_chiO))) -
        k1 * withOrder(Integral(m_tauK * divExpKnownD2(at), divergence<D>(m_chiD))) -
        k1 * withOrder(Integral(m_tauK * divExpKnownD2(at), divergence<O>(m_chiO))));
    }
    else
    {
      withOutlet(body +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<D>(at, m_psiD), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<O>(at, m_psiO), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<D>(at, m_psiD), divergence<O>(m_chiO))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<O>(at, m_psiO), divergence<O>(m_chiO))));
    }
  }

  void AtriumViscoelastic::setupFlowPsiRB()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const Real dtM = massStep();
    m_formDt = dt;
    m_formMassStep = dtM;
    const auto& ob = m_cfg.oldroydB;
    const auto etaS =
      PointwiseScalar([this](const Point& p) -> Real { return solventViscosityAt(p); },
        constantSolvent() ? size_t(0) : size_t(2));
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0; // sigma = s (exp(psi) - I)
    const Real kB = ob.lambda / (2.0 * ob.etaP); // tau_psi = kB alpha3
    // Memory of the dynamic subscales (zero for quasi-static ones): rho/dt in
    // the velocity subscales, 1/dt in the log-conformation one.
    const Real memU = m_cfg.dynamicSubscales ? rho / dt : 0.0;
    const Real memP = m_cfg.dynamicSubscales ? 1.0 / dt : 0.0;
    const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
    const LogModel model{ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda)};
    constexpr Block D = Block::D;
    constexpr Block O = Block::O;
    const Real wD = voigtWeight(D), wO = voigtWeight(O);

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto convU = Mult(Jacobian(m_u), m_uAdv);
    const auto temam = Dot(Div(m_uAdv) * m_u, m_v);
    // Both kinds of opening may reverse: with the outward normal, the entering
    // flux carries the convective boundary term of either sign (Atrium/Fluid).
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;
    // Time-dependent data enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };

    // The stabilisation advects with u^n, the velocity of its parameters.
    const auto streamV = Mult(Jacobian(m_v), m_uOld);
    const auto uu = PointwiseTensor([this](const Point& p) -> Eigen::Matrix3d {
      const auto u = m_uOld.getValue(p);
      const Eigen::Vector3d w(u(0), u(1), u(2));
      return w * w.transpose();
    });
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(3, 3));

    // Pointwise data at the Newton point (psi^k, u^k).
    LogFields<VectorGridFunctionType> at{
      &m_psiItD, &m_psiItO, &m_uIt, model, &m_psiRevision};
    at.lambda = lambdaField();

    // Momentum: exp(psi) ~ exp(psi^k) + Dexp[psi - psi^k].
    const auto Tn = fromPoint(at, LogDexp, [](const LogPoint& q) { return q.exp; });
    const auto C = Tn - dexp<D>(at, m_psiItD) - dexp<O>(at, m_psiItO);

    const size_t psiOrder = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));
    const auto withOrder = [psiOrder](auto integral) {
      integral.setOrder(psiOrder);
      return integral;
    };

    // ---- psi equation, Newton about (psi^k, u^k), as in setupFlow -----------
    // R_psi = (psi - psi_hist)/dtM + (u . grad) psi + N(psi, grad u), its block
    // C paired against chi_C: rCB the trial part in psi_B, rCu in u, kC known.
    const auto rDD =
      pairedJacobianN<D, D>(at, wD / dtM, m_psiD) + wD * Mult(Jacobian(m_psiD), m_uIt);
    const auto rOO =
      pairedJacobianN<O, O>(at, wO / dtM, m_psiO) + wO * Mult(Jacobian(m_psiO), m_uIt);
    const auto rDO = pairedJacobianN<D, O>(at, 0.0, m_psiO);
    const auto rOD = pairedJacobianN<O, D>(at, 0.0, m_psiD);
    const auto rDu = wD * Mult(Jacobian(m_psiItD), m_u) + pairedVelocityMap<D>(at, m_u);
    const auto rOu = wO * Mult(Jacobian(m_psiItO), m_u) + pairedVelocityMap<O>(at, m_u);
    const auto kD = pairedResidual<D>(at) - (wD / dtM) * m_psiHistD -
      wD * Mult(Jacobian(m_psiItD), m_uIt);
    const auto kO = pairedResidual<O>(at) - (wO / dtM) * m_psiHistO -
      wO * Mult(Jacobian(m_psiItO), m_uIt);

    // ---- Residual-based OSS of the psi equation ------------------------------
    // psi~ = tau' (memP psi~^n - R_psi + Pi[R_psi^n]), tau = kB alpha3, enters
    //   <psi~, -(u^n . grad) chi + JN^*[chi]>  (the psi equation: the adjoint
    //                                           of its linearised operator)
    //   <s Dexp[psi~], eps(v)>                 (sigma~ = s Dexp[psi~] in the
    //                                           momentum),
    // the second giving the u-psi compatibility. Both carry the full residual,
    // so they vanish where it is resolved. The block-diagonal parts are a
    // Petrov-Galerkin test chi_C + tau' (u^n . grad) chi_C - tau' A[C][C]^T
    // chi_C / w_C on the Galerkin terms; the D-O coupling of JN^* is tested
    // separately (crossD, crossO). Config::adjointReaction (off by default)
    // switches JN^* on.
    const auto tauP = kB * m_alpha3Dyn;
    const auto tauR = ((m_cfg.adjointReaction ? 1.0 : 0.0) * kB) * m_alpha3Dyn;
    // A[C][B]^T / w_B: JN^* restricted to the blocks (C, B), for the residual
    // of block B (paired with weight w_B) tested with chi_C. Like the rest of
    // the adjoint (u^n, alpha^n), JN is taken at step n, (psi^n, grad u^n):
    // the term stays linear in the unknowns and Newton exact.
    LogFields<VectorGridFunctionType> atN{
      &m_psiOldD, &m_psiOldO, &m_uOld, model, &m_stepRevision};
    atN.lambda = lambdaField();
    const auto adjointJN = [&atN](size_t c, size_t b, Real w) {
      return fromPoint(atN, LogNewton, [c, b, w](const LogPoint& q) -> Eigen::Matrix3d {
        return q.A[c][b].transpose() / w;
      });
    };
    const auto stabD = tauP * Mult(Jacobian(m_chiD), m_uOld) +
      (-1.0) * tauR * Mult(adjointJN(0, 0, wD), m_chiD);
    const auto stabO = tauP * Mult(Jacobian(m_chiO), m_uOld) +
      (-1.0) * tauR * Mult(adjointJN(1, 1, wO), m_chiO);
    const auto testD = m_chiD + stabD;
    const auto testO = m_chiO + stabO;
    // d_t psi_h leaves the residual of psi~ (Config::subscaleTimeResidual off):
    // cT (psi - psi_hist) is taken back from the stabilisation terms below
    // (the Petrov-Galerkin test, and sigma~ in the momentum).
    const Real cT = m_cfg.subscaleTimeResidual ? 0.0 : 1.0 / dtM;
    const Real hT = m_cfg.subscaleTimeResidual ? 1.0 : 0.0; // psi_hist in R_known
    const auto crossD =
      (-1.0) * tauR * Mult(adjointJN(0, 1, wO), m_chiD); // residual O, chi_D
    const auto crossO =
      (-1.0) * tauR * Mult(adjointJN(1, 0, wD), m_chiO); // residual D, chi_O
    // tau' pair<C>(JN[memP psi~^n + Pi[R^n]]) = tau' sum_B A[C][B] y_B.
    const auto knownJN = [this, atN, memP](size_t c) {
      return VectorFunction(
        size_t(3), [this, atN, memP, c](const Point& p) -> Math::SpatialVector<Real> {
          const LogPoint& q = logPointAt(atN, p, LogNewton);
          const Eigen::Matrix3d y =
            unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p)) +
            memP * subscales().at(p).sig;
          const Eigen::Vector3d yD(y(0, 0), y(1, 1), y(2, 2)),
            yO(y(0, 1), y(0, 2), y(1, 2));
          const Eigen::Vector3d v = q.A[c][0] * yD + q.A[c][1] * yO;
          return Math::SpatialVector<Real>{{v(0), v(1), v(2)}};
        });
    };
    // The known part of -psi~ paired against block C: w_C (Pi[R^n] + memP psi~^n)_C.
    const auto knownSub = [this, memP](bool off) {
      return VectorFunction(
        size_t(3), [this, memP, off](const Point& p) -> Math::SpatialVector<Real> {
          const Eigen::Matrix3d m =
            unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p)) +
            memP * subscales().at(p).sig;
          return off
            ? Math::SpatialVector<Real>{{2.0 * m(0, 1), 2.0 * m(0, 2), 2.0 * m(1, 2)}}
            : packD(m);
        });
    };
    // s tau' Dexp[memP psi~^n + Pi[R^n] - R_known] against eps(v).
    const auto stressKnown = PointwiseTensor(
      [this, at, kB, s, memP, dtM, hT](const Point& p) -> Eigen::Matrix3d {
        const LogPoint& q = logPointAt(at, p, LogNewton);
        // R_known = R_rest - hT psi_hist/dtM - (u^k . grad) psi^k, as a tensor.
        Eigen::Matrix3d rk;
        rk << q.r[0](0), 0.5 * q.r[1](0), 0.5 * q.r[1](1), 0.5 * q.r[1](0), q.r[0](1),
          0.5 * q.r[1](2), 0.5 * q.r[1](1), 0.5 * q.r[1](2), q.r[0](2);
        rk -= (hT / dtM) * unpack(m_psiHistD.getValue(p), m_psiHistO.getValue(p));
        const auto GD = Jacobian(m_psiItD).getValue(p);
        const auto GO = Jacobian(m_psiItO).getValue(p);
        const auto uk = m_uIt.getValue(p);
        Eigen::Vector3d advD = Eigen::Vector3d::Zero(), advO = Eigen::Vector3d::Zero();
        for (int c = 0; c < 3; ++c)
        {
          for (int j = 0; j < 3; ++j)
          {
            advD(c) += GD(c, j) * uk(j);
            advO(c) += GO(c, j) * uk(j);
          }
        }
        rk -= unpack(advD, advO);
        const Eigen::Matrix3d X = memP * subscales().at(p).sig +
          unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p)) - rk;
        // Dexp[X] = sum_m X_m D_m, X = sum_m X_m E_m.
        const Real xm[6] = {X(0, 0), X(1, 1), X(2, 2), X(0, 1), X(0, 2), X(1, 2)};
        Eigen::Matrix3d out = Eigen::Matrix3d::Zero();
        for (size_t m = 0; m < 6; ++m)
          out += xm[m] * q.D[m];
        return s * kB * m_alpha3Dyn.getValue(p) * out;
      });
    const Real mst = -s; // -s tau' Dexp[R_trial] against eps(v)

    // ---- Subscales of step n --------------------------------------------------
    const auto u1n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u1; });
    const auto u2n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u2; });
    const auto u3n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u3; });
    const Real sd = m_cfg.stressDivScale;
    const Real ps = m_cfg.pressureScale;

    const auto body =
      // ---- Momentum, sigma = s (exp(psi) - I), with sigma~ = s Dexp[psi~] ------
      (rho / dtM) * Integral(m_u, m_v) - (rho / dtM) * Integral(m_uHist, m_v) +
      rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam) +
      withOrder(Integral(2.0 * etaS * symU +
          mst * tauP * (dexp<D>(at, (1.0 / wD) * rDu) + dexp<O>(at, (1.0 / wO) * rOu)),
        symV)) +
      withOrder(Integral(s * dexp<D>(at, m_psiD) +
          mst * tauP * (dexp<D>(at, (1.0 / wD) * rDD) + dexp<O>(at, (1.0 / wO) * rOD)),
        symV)) +
      withOrder(Integral(s * dexp<O>(at, m_psiO) +
          mst * tauP * (dexp<O>(at, (1.0 / wO) * rOO) + dexp<D>(at, (1.0 / wD) * rDO)),
        symV)) +
      withOrder(Integral(stressKnown, symV))
      // d_t psi_h out of sigma~: + s tau' cT Dexp[psi] (mst tauP rDD carries -s tau'/dtM)
      + (s * cT) *
        withOrder(Integral(tauP * (dexp<D>(at, m_psiD) + dexp<O>(at, m_psiO)), symV)) +
      s * withOrder(Integral(C, symV)) - s * Integral(I, symV) -
      Integral(m_p, Div(m_v))

      // ---- Continuity -----------------------------------------------------
      + Integral(Div(m_u), m_q) +
      m_cfg.pressurePenalty * Integral(m_p, m_q)

      // ---- psi equation, Galerkin and the streamline part of its OSS ---------
      + withOrder(Integral(rDD, testD)) + withOrder(Integral(rDO, testD)) +
      withOrder(Integral(rDu, testD)) + withOrder(Integral(rOO, testO)) +
      withOrder(Integral(rOD, testO)) + withOrder(Integral(rOu, testO)) +
      withOrder(Integral(kD, testD)) + withOrder(Integral(kO, testO)) -
      withOrder(Integral(tauP * knownSub(false), Mult(Jacobian(m_chiD), m_uOld))) -
      withOrder(Integral(tauP * knownSub(true), Mult(Jacobian(m_chiO), m_uOld)))
      // d_t psi_h out of the Petrov-Galerkin test: - cT w_C <psi_C - psi_hist,C, stab_C>
      - (cT * wD) * withOrder(Integral(m_psiD, stabD)) +
      (cT * wD) * withOrder(Integral(m_psiHistD, stabD)) -
      (cT * wO) * withOrder(Integral(m_psiO, stabO)) +
      (cT * wO) * withOrder(Integral(m_psiHistO, stabO))
      // the D-O coupling of the adjoint reaction, and its known part
      + withOrder(Integral(rOO, crossD)) + withOrder(Integral(rOD, crossD)) +
      withOrder(Integral(rOu, crossD)) + withOrder(Integral(rDD, crossO)) +
      withOrder(Integral(rDO, crossO)) + withOrder(Integral(rDu, crossO)) +
      withOrder(Integral(kO, crossD)) + withOrder(Integral(kD, crossO)) -
      (cT * wO) * withOrder(Integral(m_psiO, crossD)) +
      (cT * wO) * withOrder(Integral(m_psiHistO, crossD)) -
      (cT * wD) * withOrder(Integral(m_psiD, crossO)) +
      (cT * wD) * withOrder(Integral(m_psiHistD, crossO)) +
      withOrder(Integral(tauR * knownJN(0), m_chiD)) +
      withOrder(Integral(tauR * knownJN(1), m_chiO))

      // ---- Momentum, split OSS, dynamic subscales (as setupFlowExp) ------------
      + vms * rho * rho * Integral(m_tauK * Mult(Jacobian(m_u), uu), Jacobian(m_v)) -
      vms * rho * rho *
        withOrder(Integral(m_tauK * (m_piConv.get() + (memU / rho) * u1n), streamV)) +
      ps * Integral(m_tauK * Grad(m_p), Grad(m_q)) -
      ps * withOrder(Integral(m_tauK * (m_piGradP.get() + memU * u2n), Grad(m_q))) -
      sd *
        withOrder(
          Integral(m_tauK * (m_piDivSigma.get() - memU * u3n), divergence<D>(m_chiD))) -
      sd *
        withOrder(
          Integral(m_tauK * (m_piDivSigma.get() - memU * u3n), divergence<O>(m_chiO)))

      // ---- Continuity, OSS of div u (algebraic subscale) ----------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v)) -
      Integral(m_alpha2 * m_piDiv.get(), Div(m_v))

      // ---- Openings: tractions -p_pv n and -(p_out + Z u.n) n ---------------
      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet) +
      m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.labels.outlet) +
      m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet) +
      BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Wall: no slip --------------------------------------------------------
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // The psi equation is hyperbolic in psi: relaxed fluid (psi = 0) enters
    // through the veins, imposed by the upwind term |u*.n|_- psi : chi, which
    // vanishes while a vein is an outflow (Config::inletConformation).
    const auto inflowOut = Max(-Dot(m_uAdv, normal), 0.0);
    const auto finish = [&](const auto& form) {
      if (m_cfg.inletConformation)
      {
        m_flow = form +
          wD * BoundaryIntegral(inflowOut * Dot(m_psiD, m_chiD)).over(m_inletSet) +
          wO * BoundaryIntegral(inflowOut * Dot(m_psiO, m_chiO)).over(m_inletSet);
      }
      else
        m_flow = form;
    };
    // Outlet (Config::outletConformationInflow): the upwind term of setupFlow.
    const auto outlet = m_cfg.labels.outlet;
    const auto withOutlet = [&](const auto& form) {
      if (m_cfg.outletConformationInflow)
      {
        finish(form +
          wD * BoundaryIntegral(inflowOut * Dot(m_psiD, m_chiD)).over(outlet) -
          wD * BoundaryIntegral(inflowOut * Dot(m_psiOldD, m_chiD)).over(outlet) +
          wO * BoundaryIntegral(inflowOut * Dot(m_psiO, m_chiO)).over(outlet) -
          wO * BoundaryIntegral(inflowOut * Dot(m_psiOldO, m_chiO)).over(outlet));
      }
      else
        finish(form);
    };
    // u3 trial: s div(Dexp[psi]) (+ the D2exp part of its Newton linearisation).
    const Real k1 = sd * s;
    if (m_cfg.newtonD2)
    {
      const auto dD = divExpTrial<D>(at, m_psiD) + divExpTrialD2<D>(at, m_psiD);
      const auto dO = divExpTrial<O>(at, m_psiO) + divExpTrialD2<O>(at, m_psiO);
      withOutlet(body + k1 * withOrder(Integral(m_tauK * dD, divergence<D>(m_chiD))) +
        k1 * withOrder(Integral(m_tauK * dO, divergence<D>(m_chiD))) +
        k1 * withOrder(Integral(m_tauK * dD, divergence<O>(m_chiO))) +
        k1 * withOrder(Integral(m_tauK * dO, divergence<O>(m_chiO))) -
        k1 * withOrder(Integral(m_tauK * divExpKnownD2(at), divergence<D>(m_chiD))) -
        k1 * withOrder(Integral(m_tauK * divExpKnownD2(at), divergence<O>(m_chiO))));
    }
    else
    {
      withOutlet(body +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<D>(at, m_psiD), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<O>(at, m_psiO), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<D>(at, m_psiD), divergence<O>(m_chiO))) +
        k1 *
          withOrder(
            Integral(m_tauK * divExpTrial<O>(at, m_psiO), divergence<O>(m_chiO))));
    }
  }

  void AtriumViscoelastic::updateStabilizationExp()
  {
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const auto& ob = m_cfg.oldroydB;
    const bool dynamic = m_cfg.dynamicSubscales;
    const Real gradDiv = m_cfg.useVMS ? m_cfg.gradDivScale : 0.0;

    // f(psi^n) by cell: the reaction of alpha3 and of the adjoint, and the
    // polymer viscosity eta_p/f of eta_0.
    // Owens: eta_p = G lambda(p), G = eta_p,inf/lambda_H, so f/(2 eta_p) = (f
    // lambda_H/lambda(p))/(2 eta_p,inf) and eta_p/f likewise: the factor rides on f.
    if (twoField())
      m_pttN = Real(1); // psi = 0
    else
    {
      m_pttN.project(RealFunction([this](const Point& p) {
        const Real f = pttFactor(
          logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp().trace());
        return isOwens() ? f * m_cfg.oldroydB.lambda / m_lambdaCell.getValue(p) : f;
      }));
    }
    // alpha1 = [c1 eta_0/h1^2 + c2 rho |u|/h2]^-1, eta_0 = eta_s + eta_p/f.
    m_alpha1.project(RealFunction([this, rho](const Point& p) {
      const auto u = m_uOld.getValue(p);
      const Real h1 = cellSize(p);
      const Real h2 = streamlineLength(p, u);
      const Real eta0 = solventViscosityAt(p) + m_cfg.oldroydB.etaP / m_pttN.getValue(p);
      return 1.0 /
        (m_cfg.c1 * eta0 / (h1 * h1) + m_cfg.c2 * rho * std::sqrt(Math::dot(u, u)) / h2);
    }));
    // alpha1' = (rho/dt + 1/alpha1)^-1 (dynamic) or alpha1.
    m_tauK.project(RealFunction([this, rho, dt, dynamic](const Point& p) {
      const Real a1 = m_alpha1.getValue(p);
      return dynamic ? 1.0 / (rho / dt + 1.0 / a1) : a1;
    }));
    // alpha2 = h1^2/(c1 alpha1).
    m_alpha2.project(RealFunction([this, gradDiv](const Point& p) {
      const Real h1 = cellSize(p);
      return gradDiv * h1 * h1 / (m_cfg.c1 * m_alpha1.getValue(p));
    }));
    // alpha3 = [c3 f/(2 eta_p) + c4 (lambda |u|/(2 eta_p h2) + lambda |grad u|/eta_p)]^-1
    // and alpha3' = (lambda/(2 eta_p dt) + 1/alpha3)^-1 (dynamic) or alpha3.
    m_alpha3.project(RealFunction([this, ob](const Point& p) {
      const auto u = m_uOld.getValue(p);
      const Real h2 = streamlineLength(p, u);
      const Real speed = std::sqrt(Math::dot(u, u));
      const Real gradNorm = Jacobian(m_uOld).getValue(p).norm();
      const Real f = m_pttN.getValue(p);
      return m_cfg.stressScale /
        (m_cfg.c3 * f / (2.0 * ob.etaP) +
          m_cfg.c4 *
            (ob.lambda * speed / (2.0 * ob.etaP * h2) + ob.lambda * gradNorm / ob.etaP));
    }));
    m_alpha3Dyn.project(RealFunction([this, ob, dt, dynamic](const Point& p) {
      const Real a3 = m_alpha3.getValue(p);
      return dynamic ? 1.0 / (ob.lambda / (2.0 * ob.etaP * dt) + 1.0 / a3) : a3;
    }));

    // Pi[(grad u^n) u^n]; Pi[div u^n], Pi[grad p^n], Pi[s div E^n] in
    // projectConsistent. Pi[R^n] is from the end of the last step.
    m_piConv.project(Mult(Jacobian(m_uOld), m_uOld));
    projectConsistent(false, 0u);
  }

  void AtriumViscoelastic::updateSubscalesExp()
  {
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;
    const Real kB = ob.lambda / (2.0 * ob.etaP);
    const Real memU = rho / dt;
    const Real memS = kB / dt;
    Real bdfA = 1.0, bdfB = 0.0;
    if (m_cfg.timeOrder == 2)
    {
      const Real w = bdfRatio();
      bdfA = (1.0 + w) * (1.0 + w) / (1.0 + 2.0 * w);
      bdfB = -w * w / (1.0 + 2.0 * w);
    }
    // The data of the accepted step: (psi^{n+1}, u^{n+1}) in the iterates, its
    // BDF history and mass step, so R is the residual the step solved.
    const LogModel model{
      ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda), bdfA, bdfB, massStep()};
    LogFields<VectorGridFunctionType> at{&m_psiItD, &m_psiItO, &m_uIt, model,
      &m_psiRevision, &m_psiOldD, &m_psiOldO, &m_psiNm1D, &m_psiNm1O};
    at.lambda = lambdaField();
    const size_t psiOrder = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));
    const Real cT =
      m_cfg.subscaleTimeResidual ? 0.0 : 1.0; // d_t E_h out of R (setupFlowExp)

    if (m_cfg.dynamicSubscales)
    {
      // s~^{n+1} = alpha' (mem s~^n + rhs) at every quadrature point of the
      // psi rule, with the same data as the assembled terms: u^n advects, the
      // projections are those of the step.
      const auto& pSol = m_p.getSolution();
      auto sweep = Integral(
        PointwiseScalar(
          [this, at, s, rho, memU, memS, cT, &pSol](const Point& p) -> Real {
            const LogPoint& q = logPointAt(at, p, LogExpForm);
            const SubscaleStore::Entry& old = subscales().at(p);
            const Real a1 = m_tauK.getValue(p);
            const Real a3 = m_alpha3Dyn.getValue(p);
            const auto toEigen = [](const auto& v) {
              return Eigen::Vector3d(v(0), v(1), v(2));
            };
            const auto J = Jacobian(m_uIt).getValue(p);
            const Eigen::Vector3d un = toEigen(m_uOld.getValue(p));
            Eigen::Vector3d conv;
            for (int a = 0; a < 3; ++a)
              conv(a) = J(a, 0) * un(0) + J(a, 1) * un(1) + J(a, 2) * un(2);
            Eigen::Vector3d divE;
            for (int i = 0; i < 3; ++i)
              divE(i) = q.G[0](i, 0) + q.G[1](i, 1) + q.G[2](i, 2);
            SubscaleStore::Entry e;
            e.u1 =
              a1 * (memU * old.u1 - rho * (conv - toEigen(m_piConv.get().getValue(p))));
            e.u2 = a1 *
              (memU * old.u2 -
                (toEigen(Grad(pSol).getValue(p)) - toEigen(m_piGradP.get().getValue(p))));
            e.u3 =
              a1 * (memU * old.u3 + (s * divE - toEigen(m_piDivSigma.get().getValue(p))));
            const Eigen::Matrix3d piR =
              unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p));
            e.sig = a3 * (memS * old.sig - (q.R - cT * q.rate - piR));
            subscales().stage(p, e);
            return 0.0;
          },
          2),
        m_qFlux);
      sweep.setOrder(psiOrder);
      m_subscaleSweep = sweep;
      m_subscaleSweep.assemble();
      subscales().commit();
    }

    // Pi[R^{n+1}], read as Pi[R^n] by the next step.
    const auto R = fromPoint(at, LogExpForm,
      [cT](const LogPoint& q) -> Eigen::Matrix3d { return q.R - cT * q.rate; });
    m_piRD.project(voigtD(R));
    m_piRO.project(voigtO(R));
  }

  void AtriumViscoelastic::updateSubscalesPsiRB()
  {
    if (twoField())
    {
      // u1 and u2 only (u3 and psi~ vanish with sigma); Pi[R_psi] stays 0.
      if (!m_cfg.dynamicSubscales)
        return;
      const Real rho = m_cfg.rho;
      const Real memU = rho / m_dt;
      const auto& pSol = m_p.getSolution();
      auto sweep = Integral(
        PointwiseScalar(
          [this, rho, memU, &pSol](const Point& p) -> Real {
            const SubscaleStore::Entry& old = subscales().at(p);
            const Real a1 = m_tauK.getValue(p);
            const auto toEigen = [](const auto& v) {
              return Eigen::Vector3d(v(0), v(1), v(2));
            };
            const auto J = Jacobian(m_uIt).getValue(p);
            const Eigen::Vector3d un = toEigen(m_uOld.getValue(p));
            Eigen::Vector3d conv;
            for (int a = 0; a < 3; ++a)
              conv(a) = J(a, 0) * un(0) + J(a, 1) * un(1) + J(a, 2) * un(2);
            SubscaleStore::Entry e;
            e.u1 =
              a1 * (memU * old.u1 - rho * (conv - toEigen(m_piConv.get().getValue(p))));
            e.u2 = a1 *
              (memU * old.u2 -
                (toEigen(Grad(pSol).getValue(p)) - toEigen(m_piGradP.get().getValue(p))));
            e.u3.setZero();
            e.sig.setZero();
            subscales().stage(p, e);
            return 0.0;
          },
          2),
        m_qFlux);
      sweep.setOrder(static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder)));
      m_subscaleSweep = sweep;
      m_subscaleSweep.assemble();
      subscales().commit();
      return;
    }
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const Real dtM = massStep();
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;
    const Real kB = ob.lambda / (2.0 * ob.etaP);
    const Real memU = rho / dt;
    const Real memP = 1.0 / dt;
    const LogModel model{ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda)};
    const size_t psiOrder = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));

    // The data of the accepted step, (psi^{n+1}, u^{n+1}) in the iterates:
    // R_psi = (psi - psi_hist)/dtM + (u . grad) psi + N(psi, grad u), and
    // div exp(psi) = sum_j Dexp[d_j psi] e_j.
    struct StepPoint
    {
      Eigen::Matrix3d R;
      Eigen::Vector3d divE;
    };
    const Real hT = m_cfg.subscaleTimeResidual ? 1.0 : 0.0;
    const auto stepPoint = [this, model, dtM, hT](const Point& p) -> StepPoint {
      const Eigen::Matrix3d psi = unpack(m_psiItD.getValue(p), m_psiItO.getValue(p));
      const auto J = Jacobian(m_uIt).getValue(p);
      const auto uv = m_uIt.getValue(p);
      Eigen::Matrix3d L;
      for (int a = 0; a < 3; ++a)
      {
        for (int b = 0; b < 3; ++b)
          L(a, b) = J(a, b);
      }
      const auto GD = Jacobian(m_psiItD).getValue(p);
      const auto GO = Jacobian(m_psiItO).getValue(p);
      const auto& B = voigtBasisEigen();
      const LogEig eig = logEig(psi);
      StepPoint out;
      // d_t psi_h only with Config::subscaleTimeResidual (P^perp[d_t psi_h] = 0).
      Eigen::Vector3d nw = Eigen::Vector3d::Zero();
      const Real ww = wallWeight(p, nw);
      out.R =
        (hT / dtM) * (psi - unpack(m_psiHistD.getValue(p), m_psiHistO.getValue(p))) +
        nonlinearMap(psi, wallGradient(L, ww, nw), pointModel(model, lambdaField(), p));
      out.divE.setZero();
      for (int j = 0; j < 3; ++j)
      {
        Eigen::Matrix3d g = Eigen::Matrix3d::Zero();
        for (int k = 0; k < 3; ++k)
          g += GD(k, j) * B[k] + GO(k, j) * B[3 + k];
        out.R += uv(j) * g;
        out.divE += eig.dexp(g).col(j);
      }
      return out;
    };

    if (m_cfg.dynamicSubscales)
    {
      const auto& pSol = m_p.getSolution();
      auto sweep = Integral(
        PointwiseScalar(
          [this, stepPoint, s, kB, rho, memU, memP, &pSol](const Point& p) -> Real {
            const StepPoint q = stepPoint(p);
            const SubscaleStore::Entry& old = subscales().at(p);
            const Real a1 = m_tauK.getValue(p);
            const Real tp = kB * m_alpha3Dyn.getValue(p);
            const auto toEigen = [](const auto& v) {
              return Eigen::Vector3d(v(0), v(1), v(2));
            };
            const auto J = Jacobian(m_uIt).getValue(p);
            const Eigen::Vector3d un = toEigen(m_uOld.getValue(p));
            Eigen::Vector3d conv;
            for (int a = 0; a < 3; ++a)
              conv(a) = J(a, 0) * un(0) + J(a, 1) * un(1) + J(a, 2) * un(2);
            SubscaleStore::Entry e;
            e.u1 =
              a1 * (memU * old.u1 - rho * (conv - toEigen(m_piConv.get().getValue(p))));
            e.u2 = a1 *
              (memU * old.u2 -
                (toEigen(Grad(pSol).getValue(p)) - toEigen(m_piGradP.get().getValue(p))));
            e.u3 = a1 *
              (memU * old.u3 + (s * q.divE - toEigen(m_piDivSigma.get().getValue(p))));
            const Eigen::Matrix3d piR =
              unpack(m_piRD.get().getValue(p), m_piRO.get().getValue(p));
            e.sig = tp * (memP * old.sig - (q.R - piR)); // psi~ (dimensionless)
            subscales().stage(p, e);
            return 0.0;
          },
          2),
        m_qFlux);
      sweep.setOrder(psiOrder);
      m_subscaleSweep = sweep;
      m_subscaleSweep.assemble();
      subscales().commit();
    }

    // Pi[R_psi^{n+1}], read as Pi[R_psi^n] by the next step.
    const auto R = PointwiseTensor(
      [stepPoint](const Point& p) -> Eigen::Matrix3d { return stepPoint(p).R; });
    m_piRD.project(voigtD(R));
    m_piRO.project(voigtO(R));
  }

  void AtriumViscoelastic::setupFlowTwoField()
  {
    // setupFlowPsiRB with sigma = 0: momentum with the generalised Newtonian
    // viscosity eta(gd(u*)) (eta_ref for "newtonian"), continuity, the split
    // OSS with dynamic subscales (u1 convection, u2 pressure gradient; u3 =
    // Pi^perp[div sigma] vanishes), the algebraic pressure subscale, and the
    // same boundary conditions. Linear in (u, p) for the frozen u* and eta:
    // one Newton iteration is exact.
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const Real dtM = massStep();
    m_formDt = dt;
    m_formMassStep = dtM;
    const auto etaS =
      PointwiseScalar([this](const Point& p) -> Real { return solventViscosityAt(p); },
        constantSolvent() ? size_t(0) : size_t(2));
    const Real memU = m_cfg.dynamicSubscales ? rho / dt : 0.0;
    const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
    const Real ps = m_cfg.pressureScale;
    const size_t order = static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder));
    const auto withOrder = [order](auto integral) {
      integral.setOrder(order);
      return integral;
    };

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto convU = Mult(Jacobian(m_u), m_uAdv);
    const auto temam = Dot(Div(m_uAdv) * m_u, m_v);
    // Both kinds of opening may reverse: with the outward normal, the entering
    // flux carries the convective boundary term of either sign (Atrium/Fluid).
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;
    // Time-dependent data enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };
    const auto streamV = Mult(Jacobian(m_v), m_uOld);
    const auto uu = PointwiseTensor([this](const Point& p) -> Eigen::Matrix3d {
      const auto u = m_uOld.getValue(p);
      const Eigen::Vector3d w(u(0), u(1), u(2));
      return w * w.transpose();
    });
    const auto u1n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u1; });
    const auto u2n = subscaleVector([](const SubscaleStore::Entry& e) { return e.u2; });

    m_flow2 =
      // ---- Momentum -------------------------------------------------------
      (rho / dtM) * Integral(m_u, m_v) - (rho / dtM) * Integral(m_uHist, m_v) +
      rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam) +
      withOrder(Integral(2.0 * etaS * symU, symV)) -
      Integral(m_p, Div(m_v))
      // ---- Continuity -----------------------------------------------------
      + Integral(Div(m_u), m_q) +
      m_cfg.pressurePenalty * Integral(m_p, m_q)
      // ---- Split OSS, dynamic subscales ----------------------------------
      + vms * rho * rho * Integral(m_tauK * Mult(Jacobian(m_u), uu), Jacobian(m_v)) -
      vms * rho * rho *
        withOrder(Integral(m_tauK * (m_piConv.get() + (memU / rho) * u1n), streamV)) +
      ps * Integral(m_tauK * Grad(m_p), Grad(m_q)) -
      ps * withOrder(Integral(m_tauK * (m_piGradP.get() + memU * u2n), Grad(m_q)))
      // ---- Continuity, OSS of div u (algebraic subscale) ------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v)) -
      Integral(m_alpha2 * m_piDiv.get(), Div(m_v))
      // ---- Openings: tractions -p_pv n and -(p_out + Z u.n) n ---------------
      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet) +
      m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.labels.outlet) +
      m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet) +
      BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Wall: no slip --------------------------------------------------------
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);
  }

  void AtriumViscoelastic::setupResidence()
  {
    // Age a of the blood: d_t a + u.grad a = 1, BDF with the mass step of the
    // flow (a_hist as x_hist), advected by u^{n+1}; the blood that enters
    // through the veins has a = 0, and the backflow of the mitral plane brings
    // in the age a^n of the boundary (upwind terms).
    // SUPG with tau = [(2/dtM)^2 + (2 |u|/h)^2]^-1/2. Solved after the flow
    // step is accepted, with u^{n+1} in m_uOld.
    const auto normal = BoundaryNormal(m_mesh);
    const auto tau = [this](const Point& p) -> Real {
      const auto uc = m_uOld.getValue(p);
      const Real a = 2.0 / massStep();
      const Real b = 2.0 * std::sqrt(Math::dot(uc, uc)) / cellSize(p);
      return 1.0 / std::sqrt(a * a + b * b);
    };
    RealFunction massFn = [this](const Point&) { return 1.0 / massStep(); };
    RealFunction tauFn = [tau](const Point& p) { return tau(p); };
    RealFunction tauMassFn = [this, tau](const Point& p) { return tau(p) / massStep(); };
    RealFunction knownFn = [this](const Point& p) {
      return 1.0 + m_ageHist.getValue(p) / massStep();
    };
    RealFunction tauKnownFn = [this, tau](const Point& p) {
      return tau(p) * (1.0 + m_ageHist.getValue(p) / massStep());
    };
    const auto adv = Dot(m_uOld, Grad(m_age));
    const auto advTest = Dot(m_uOld, Grad(m_ageTest));
    const auto inflowOut = Max(-Dot(m_uOld, normal), 0.0);
    m_ageProblem = Integral(massFn * m_age, m_ageTest) + Integral(adv, m_ageTest) -
      Integral(knownFn, m_ageTest) + Integral(tauMassFn * m_age, advTest) +
      Integral(tauFn * adv, advTest) - Integral(tauKnownFn, advTest) +
      BoundaryIntegral(inflowOut * m_age, m_ageTest).over(m_cfg.labels.outlet) -
      BoundaryIntegral(inflowOut * m_ageOld, m_ageTest).over(m_cfg.labels.outlet) +
      BoundaryIntegral(inflowOut * m_age, m_ageTest).over(m_inletSet);
  }

  void AtriumViscoelastic::solveResidence()
  {
    if (m_cfg.timeOrder == 2)
    {
      const Real w = bdfRatio();
      const Real a = (1.0 + w) * (1.0 + w) / (1.0 + 2.0 * w);
      const Real b = -w * w / (1.0 + 2.0 * w);
      combine(a, m_ageOld.getData(), b, m_ageNm1.getData(), m_ageHist.getData());
    }
    else
      m_ageHist.setData(m_ageOld.getData());
    m_ageProblem.assemble();
    m_ageProblem.solve(m_ageKSP);
    // SUPG is not monotone: 0 <= a <= t.
    {
      ::Vec v = m_age.getSolution().getData();
      PetscInt size = 0;
      PetscScalar* data = nullptr;
      PetscErrorCode ierr = VecGetLocalSize(v, &size);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecGetArray(v, &data);
      assert(ierr == PETSC_SUCCESS);
      for (PetscInt i = 0; i < size; ++i)
        data[i] = std::isfinite(data[i]) ? std::clamp<Real>(data[i], 0.0, m_t) : 0.0;
      ierr = VecRestoreArray(v, &data);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }
    m_ageNm1.setData(m_ageOld.getData());
    m_ageOld.setData(m_age.getSolution().getData());
    // Volume mean, and mean over the appendage.
    m_flux = Integral(m_ageOld, m_qFlux);
    m_flux.assemble();
    m_ageMean = m_flux(m_one) / m_volume;
    if (m_appendageVolume > 0.0)
    {
      m_ageAppendage = integrate([this](const Point& p) -> Real {
        return inAppendage(p) ? m_ageOld.getValue(p) : 0.0;
      }) /
        m_appendageVolume;
    }
  }

  void AtriumViscoelastic::setupFlow()
  {
    if (twoField())
    {
      setupFlowTwoField();
      return;
    }
    if (m_cfg.formulation == "exp")
    {
      setupFlowExp();
      return;
    }
    if (m_cfg.formulation == "psi_rb")
    {
      setupFlowPsiRB();
      return;
    }
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
    const auto etaS =
      PointwiseScalar([this](const Point& p) -> Real { return solventViscosityAt(p); },
        constantSolvent() ? size_t(0) : size_t(2));
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0; // sigma = s (exp(psi) - I)

    /// The psi equation is 2x the sigma equation of the exp form (2D header).
    const Real chiScale = 2.0;
    const LogModel model{ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda)};
    constexpr Block D = Block::D;
    constexpr Block O = Block::O;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    // Galerkin convection with u^n (BDF1) or u* (BDF2).
    const auto convU = Mult(Jacobian(m_u), m_uAdv);
    const auto temam = Dot(Div(m_uAdv) * m_u, m_v);

    // Both kinds of opening may reverse: with the outward normal, the entering
    // flux carries the convective boundary term of either sign (Atrium/Fluid).
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uAdv, normal), 0.0);
    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;
    // Time-dependent data enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };

    const auto streamV = Mult(Jacobian(m_v), m_uOld); // (grad v) u^n

    // (A u^n) . (B u^n) = (A u^n (u^n)^T) : B. Bilinear streamline terms are
    // written this way so that u^n sits on the trial side.
    const auto uu = PointwiseTensor([this](const Point& p) -> Eigen::Matrix3d {
      const auto u = m_uOld.getValue(p);
      const Eigen::Vector3d w(u(0), u(1), u(2));
      return w * w.transpose();
    });
    const auto I = MatrixFunction<Math::Matrix<Real>>(Math::Matrix<Real>::Identity(3, 3));

    // Everything pointwise is at (psi^k, u^k), the Newton linearisation point.
    LogFields<VectorGridFunctionType> at{
      &m_psiItD, &m_psiItO, &m_uIt, model, &m_psiRevision};
    at.lambda = lambdaField();

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
      return pairedResidual<B>(at) - (voigtWeight(B) / dtM) * psiOld -
        voigtWeight(B) * Mult(Jacobian(psiIt), m_uIt);
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

    const auto psiDD = withOrder(Integral(
      pairedJacobianN<D, D>(at, wD / dtM, m_psiD) + wD * Mult(Jacobian(m_psiD), m_uIt),
      m_chiD));
    const auto psiOO = withOrder(Integral(
      pairedJacobianN<O, O>(at, wO / dtM, m_psiO) + wO * Mult(Jacobian(m_psiO), m_uIt),
      m_chiO));
    const auto psiOD =
      withOrder(Integral(pairedJacobianN<D, O>(at, 0.0, m_psiO), m_chiD));
    const auto psiDO =
      withOrder(Integral(pairedJacobianN<O, D>(at, 0.0, m_psiD), m_chiO));
    const auto uD = withOrder(Integral(
      wD * Mult(Jacobian(m_psiItD), m_u) + pairedVelocityMap<D>(at, m_u), m_chiD));
    const auto uO = withOrder(Integral(
      wO * Mult(Jacobian(m_psiItO), m_u) + pairedVelocityMap<O>(at, m_u), m_chiO));
    const auto knownD =
      withOrder(Integral(known.template operator()<D>(m_psiHistD, m_psiItD), m_chiD));
    const auto knownO =
      withOrder(Integral(known.template operator()<O>(m_psiHistO, m_psiItO), m_chiO));

    const auto body =
      // ---- Momentum, sigma = s (exp(psi) - I) --------------------------------
      (rho / dtM) * Integral(m_u, m_v) - (rho / dtM) * Integral(m_uHist, m_v) +
      rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam) +
      2.0 * Integral(etaS * symU, symV) +
      s * withOrder(Integral(dexp<D>(at, m_psiD), symV)) +
      s * withOrder(Integral(dexp<O>(at, m_psiO), symV)) +
      s * withOrder(Integral(C, symV)) - s * Integral(I, symV) -
      Integral(m_p, Div(m_v))

      // ---- Continuity -----------------------------------------------------
      + Integral(Div(m_u), m_q) +
      m_cfg.pressurePenalty * Integral(m_p, m_q)

      // ---- Constitutive, psi form, BDF1/BDF2, Newton ------------------------
      + psiDD + psiOO + psiOD + psiDO // psi_B : chi_C
      + uD + uO // u : chi_C
      + knownD +
      knownO // known

      // ---- S1, momentum ---------------------------------------------------
      + rho * rho * Integral(m_tauK * Mult(Jacobian(m_u), uu), Jacobian(m_v)) -
      rho * rho *
        Integral(m_tauK * (m_piConv.get() + (1.0 / dt) * m_sub.get()), streamV) +
      m_cfg.pressureScale * Integral(m_alpha1 * Grad(m_p), Grad(m_q)) -
      m_cfg.pressureScale * Integral(m_alpha1 * m_piGradP.get(), Grad(m_q)) +
      chiScale * m_cfg.stressDivScale * s *
        withOrder(
          Integral(m_alpha1 * divExpTrial<D>(at, m_psiD), divergence<D>(m_chiD))) +
      chiScale * m_cfg.stressDivScale * s *
        withOrder(
          Integral(m_alpha1 * divExpTrial<O>(at, m_psiO), divergence<D>(m_chiD))) +
      chiScale * m_cfg.stressDivScale * s *
        withOrder(
          Integral(m_alpha1 * divExpTrial<D>(at, m_psiD), divergence<O>(m_chiO))) +
      chiScale * m_cfg.stressDivScale * s *
        withOrder(
          Integral(m_alpha1 * divExpTrial<O>(at, m_psiO), divergence<O>(m_chiO))) -
      chiScale * m_cfg.stressDivScale *
        Integral(m_alpha1 * m_piDivSigma.get(), divergence<D>(m_chiD)) -
      chiScale * m_cfg.stressDivScale *
        Integral(m_alpha1 * m_piDivSigma.get(), divergence<O>(m_chiO))

      // ---- S2, continuity -------------------------------------------------
      + Integral(m_alpha2 * Div(m_u), Div(m_v)) -
      Integral(m_alpha2 * m_piDiv.get(), Div(m_v))

      // ---- S3: u-psi compatibility, and the psi advection (block-diagonal) ---
      + Integral(m_alpha3 * symU, symV) -
      Integral(m_alpha3 * (tensor<D>(m_piEpsD.get()) + tensor<O>(m_piEpsO.get())), symV) +
      wD * Integral(m_alphaPsi * Mult(Jacobian(m_psiD), uu), Jacobian(m_chiD)) +
      wO * Integral(m_alphaPsi * Mult(Jacobian(m_psiO), uu), Jacobian(m_chiO)) -
      wD * Integral(m_alphaPsi * m_piAdvPsiD.get(), Mult(Jacobian(m_chiD), m_uOld)) -
      wO * Integral(m_alphaPsi * m_piAdvPsiO.get(), Mult(Jacobian(m_chiO), m_uOld))

      // ---- Openings: tractions -p_pv n and -(p_out + Z u.n) n ---------------
      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet) +
      m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.labels.outlet) +
      m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet) +
      BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

      // ---- Wall: no slip --------------------------------------------------------
      + DirichletBC(m_u, Zero(dim)).on(m_wallSet);

    // The psi equation is hyperbolic in psi: relaxed fluid (psi = 0) enters
    // through the veins, imposed by the upwind term |u*.n|_- psi : chi, which
    // vanishes while a vein is an outflow (Config::inletConformation).
    const auto inflowOut = Max(-Dot(m_uAdv, normal), 0.0);
    const auto finish = [&](const auto& form) {
      if (m_cfg.inletConformation)
      {
        m_flow = form +
          wD * BoundaryIntegral(inflowOut * Dot(m_psiD, m_chiD)).over(m_inletSet) +
          wO * BoundaryIntegral(inflowOut * Dot(m_psiO, m_chiO)).over(m_inletSet);
      }
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
    const auto outlet = m_cfg.labels.outlet;
    const auto withSkew = [&](const auto& form) {
      if (m_cfg.psiSkewAdvection)
      {
        finish(form + 0.5 * wD * Integral(Dot(Div(m_uAdv) * m_psiD, m_chiD)) +
          0.5 * wO * Integral(Dot(Div(m_uAdv) * m_psiO, m_chiO)));
      }
      else
        finish(form);
    };
    const auto withOutlet = [&](const auto& form) {
      if (m_cfg.outletConformationInflow)
      {
        withSkew(form +
          wD * BoundaryIntegral(inflowOut * Dot(m_psiD, m_chiD)).over(outlet) -
          wD * BoundaryIntegral(inflowOut * Dot(m_psiOldD, m_chiD)).over(outlet) +
          wO * BoundaryIntegral(inflowOut * Dot(m_psiO, m_chiO)).over(outlet) -
          wO * BoundaryIntegral(inflowOut * Dot(m_psiOldO, m_chiO)).over(outlet));
      }
      else
        withSkew(form);
    };
    if (m_cfg.newtonD2)
    {
      // S1, the D2exp part of the Newton linearisation of div exp(psi).
      const Real k1 = chiScale * m_cfg.stressDivScale * s;
      withOutlet(body +
        k1 *
          withOrder(
            Integral(m_alpha1 * divExpTrialD2<D>(at, m_psiD), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_alpha1 * divExpTrialD2<O>(at, m_psiO), divergence<D>(m_chiD))) +
        k1 *
          withOrder(
            Integral(m_alpha1 * divExpTrialD2<D>(at, m_psiD), divergence<O>(m_chiO))) +
        k1 *
          withOrder(
            Integral(m_alpha1 * divExpTrialD2<O>(at, m_psiO), divergence<O>(m_chiO))) -
        k1 * withOrder(Integral(m_alpha1 * divExpKnownD2(at), divergence<D>(m_chiD))) -
        k1 * withOrder(Integral(m_alpha1 * divExpKnownD2(at), divergence<O>(m_chiO))));
    }
    else
      withOutlet(body);
  }

  void AtriumViscoelastic::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    // eta_s at the recovered rate of strain: gd^2 = 2 eps:eps = |2 eps|^2/2,
    // the rows of 2 eps being m_symRec0..2.
    RealFunction etaS([this](const Point& p) -> Real {
      if (constantSolvent())
        return m_cfg.oldroydB.etaS;
      const auto a = m_symRec0.getValue(p), b = m_symRec1.getValue(p),
                 c = m_symRec2.getValue(p);
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
    const auto traction =
      VectorFunction(etaS * Dot(m_symRec0, normal) + sD(0) * nx + sO(0) * ny + sO(1) * nz,
        etaS * Dot(m_symRec1, normal) + sO(0) * nx + sD(1) * ny + sO(2) * nz,
        etaS * Dot(m_symRec2, normal) + sO(1) * nx + sO(2) * ny + sD(2) * nz);
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  void AtriumViscoelastic::updateStabilization()
  {
    if (m_cfg.formulation != "psi")
    {
      updateStabilizationExp();
      return;
    }
    const size_t dim = m_mesh.getSpaceDimension();
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const Real vms = m_cfg.useVMS ? m_cfg.vmsScale : 0.0;
    const Real gradDiv = m_cfg.useVMS ? m_cfg.gradDivScale : 0.0;

    m_tauK.project(RealFunction([this, rho, dt, vms](const Point& p) {
      return vms / (rho / dt + rho / tau1At(p));
    }));
    m_alpha1.project(
      RealFunction([this, rho](const Point& p) { return tau1At(p) / rho; }));
    m_alpha2.project(RealFunction([this, rho, gradDiv](const Point& p) {
      const Real h = cellSize(p);
      return gradDiv * rho * h * h / (4.0 * tau1At(p));
    }));
    m_alpha3.project(
      RealFunction([this](const Point& p) { return m_cfg.stressScale * alpha3At(p); }));
    {
      const auto& ob = m_cfg.oldroydB;
      const Real k = ob.lambda / (2.0 * ob.etaP);
      const Real s = ob.etaP / lambda0();
      m_alphaPsi.project(RealFunction([this, k, s](const Point& p) {
        return 2.0 * k * k * s * m_cfg.stressScale * alpha3At(p);
      }));
    }

    m_piConv.project(Mult(Jacobian(m_uOld), m_uOld));
    m_sub.project(VectorFunction(
      dim, [this, rho, dt, dim](const Point& p) -> Math::SpatialVector<Real> {
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

  void AtriumViscoelastic::projectConsistent(bool iterate, unsigned mask)
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
    if (twoField())
      return; // sigma = 0: Pi[div sigma] stays 0
    const auto& ob = m_cfg.oldroydB;
    const Real l0 = lambda0();
    const Real s = ob.etaP / l0;
    const LogModel model{ob.lambda, l0, ob.pttEpsilon, 2.0 * (1.0 - l0 / ob.lambda)};
    // At the iterate the pointwise data are those of the assembly (same key);
    // step data (psi^n, u^n) are keyed by the step revision, as in setupFlowPsiRB.
    LogFields<VectorGridFunctionType> atX{
      &psiSD, &psiSO, &uS, model, sigmaIt ? &m_psiRevision : &m_stepRevision};
    atX.lambda = lambdaField();
    m_piDivSigma.project(
      s * (divExp<Block::D>(atX, psiSD) + divExp<Block::O>(atX, psiSO)));
    // Nothing else reads psi^n pointwise this step: free its points.
    if (!sigmaIt)
      releaseLogPoints(&m_psiOldD);
    // The residual-based formulations have no term-by-term constitutive projections.
    if (m_cfg.formulation != "psi")
      return;
    const auto epsN = 0.5 * (Jacobian(uE) + Transpose(Jacobian(uE)));
    m_piEpsD.project(voigtD(epsN));
    m_piEpsO.project(voigtO(epsN));
    // The advecting velocity of the S3 psi term is u^n on both sides.
    m_piAdvPsiD.project(Mult(Jacobian(psiAD), m_uOld));
    m_piAdvPsiO.project(Mult(Jacobian(psiAO), m_uOld));
  }

  bool AtriumViscoelastic::solveFlow()
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
      constexpr size_t ContractionSamples = 16; ///< recent estimates of C kept
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
    Real previous = -1.0; // the last undamped increment, < 0 if none
    const int maxIterations = std::max(1, m_cfg.conformationIterations);
    for (int k = 0; k < maxIterations; ++k)
    {
      if (m_cfg.timeOrder == 2 && k == 0 && m_cfg.formulation == "psi")
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

      const bool two = twoField();
      const auto assemblyStart = CoronaryClock::now();
      if (two)
        m_flow2.assemble();
      else
        m_flow.assemble();
      m_timing.assembly += secondsSince(assemblyStart);

      if (two ? !m_flow2FieldSplitsSet : !m_flowFieldSplitsSet)
      {
        if (isRoot())
          ThreeDInfo() << "Creating field splits ..." << Alert::Raise;
        if (two)
        {
          m_flow2.setFieldSplits();
          m_flow2FieldSplitsSet = true;
        }
        else
        {
          m_flow.setFieldSplits();
          m_flowFieldSplitsSet = true;
        }
      }
      const auto solveSystem = [&]() {
        if (two)
          m_flow2.solve(m_flow2KSP);
        else
          m_flow.solve(m_flowKSP);
      };

      // Config::jacobianCheck: copies of the system of iterate 0.
      const bool check =
        !two && m_cfg.jacobianCheck > 0 && k == 0 && m_step + 1 == m_cfg.jacobianCheck;
      ::Mat checkA = PETSC_NULLPTR;
      ::Vec checkB = PETSC_NULLPTR, checkX = PETSC_NULLPTR;
      if (check)
      {
        auto& ls = m_flow.getLinearSystem();
        ierr = MatDuplicate(ls.getOperator(), MAT_COPY_VALUES, &checkA);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(ls.getVector(), &checkB);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(ls.getVector(), checkB);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDuplicate(ls.getSolution(), &checkX);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(ls.getSolution(), checkX);
        assert(ierr == PETSC_SUCCESS);
      }

      // Lagged preconditioner (Config::pcReuseSolves): reused while the Krylov
      // solver needs at most pcReuseGrowth times the iterations it needed
      // with a fresh one; a solve that fails with it is repeated, from the
      // same guess, with a fresh one.
      const auto solveStart = CoronaryClock::now();
      ::KSP ksp = two ? m_flow2KSP.getHandle() : m_flowKSP.getHandle();
      ::Vec x = two ? m_flow2.getLinearSystem().getSolution()
                    : m_flow.getLinearSystem().getSolution();
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
      solveSystem();
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
        solveSystem();
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
        static_cast<Real>(its) >
          m_cfg.pcReuseGrowth * static_cast<Real>(m_pcReferenceIts))
        m_pcStale = true; // this solve is kept; the next one rebuilds
      ++m_pcSolves;
      m_timing.kspIterations += static_cast<int>(its);
      m_timing.solve += secondsSince(solveStart);

      if (check && reason > 0)
      {
        checkJacobian(checkA, checkB, checkX, x);
        MatDestroy(&checkA);
        VecDestroy(&checkB);
        VecDestroy(&checkX);
        m_failure = "Jacobian check done (Config::jacobianCheck): stop";
        m_jacobianChecked = true;
        break;
      }
      if (reason <= 0)
      {
        m_failure = "linear solver failed (KSPConvergedReason " +
          std::to_string(static_cast<int>(reason)) + ", " + std::to_string(its) +
          " iterations)";
        break;
      }

      // max |psi^{k+1} - psi^k| over both blocks, relative to max(1, |psi^k|).
      ierr =
        VecWAXPY(incrementD, -1.0, m_psiItD.getData(), m_psiD.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      ierr =
        VecWAXPY(incrementO, -1.0, m_psiItO.getData(), m_psiO.getSolution().getData());
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
        ? std::max(stepD, stepO)
        : std::numeric_limits<Real>::quiet_NaN();
      const Real psiScale = std::max(scaleD, scaleO);
      const Real omega = (std::isfinite(stepSize) && stepSize > m_cfg.newtonMaxStep)
        ? m_cfg.newtonMaxStep / stepSize
        : 1.0;
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
      {
        KSPInfo() << "Newton " << m_conformationIts << ": max|dpsi| = " << stepSize
                  << " (relative " << m_psiIncrement << ")"
                  << (omega < 1.0 ? " (damped)" : "") << "  ksp its " << its
                  << (reuse ? " (lagged PC)" : "")
                  << (predicted ? "  accepted: C d^2 = " +
                           std::to_string(contraction * m_psiIncrement * m_psiIncrement)
                                : std::string())
                  << Alert::Raise;
      }
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
      if (m_cfg.adaptiveTimeStep && undamped && previous > 0.0 &&
        m_psiIncrement > previous)
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
      {
        m_failure = "Newton did not converge in " + std::to_string(m_conformationIts) +
          " iterations (|dpsi| = " + std::to_string(m_psiIncrement) + ")";
      }
      else if (isRoot() && m_psiIncrement >= 100.0 * m_cfg.conformationTolerance)
      {
        // Without step control the step is accepted, as before.
        Alert::Warning() << "[log] Newton stopped after " << m_conformationIts
                         << " iterations with max|dpsi| = " << m_psiIncrement
                         << Alert::Raise;
      }
    }
    if (!m_failure.empty())
      return false;

    // The divergence guard, on u^{n+1} before anything is advanced.
    m_speed = std::max(std::abs(m_uIt.max()), std::abs(m_uIt.min()));
    const Real guard = m_cfg.maxVelocityFactor * m_meanVelocity;
    if (!std::isfinite(m_speed) || m_speed > guard)
    {
      m_failure = "divergence guard: max|u| = " + std::to_string(m_speed) + " m/s > " +
        std::to_string(guard) + " m/s";
      return false;
    }

    if (m_cfg.formulation == "exp")
      updateSubscalesExp();
    else if (m_cfg.formulation == "psi_rb")
      updateSubscalesPsiRB();
    advanceHistory();
    m_u.getSolution().setData(m_uIt.getData());
    m_psiD.getSolution().setData(m_psiItD.getData());
    m_psiO.getSolution().setData(m_psiItO.getData());
    ++m_psiRevision;
    updateConformation();
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.get().getData());

    m_stress = std::max({std::abs(m_sigmaD.max()), std::abs(m_sigmaD.min()),
      std::abs(m_sigmaO.max()), std::abs(m_sigmaO.min())});
    return true;
  }

  void AtriumViscoelastic::checkJacobian(::Mat A0, ::Vec b0, ::Vec x0, ::Vec x1)
  {
    // The system of a Newton iterate is A(x^k) x = b(x^k), the linearisation
    // of F about x^k: F(x) = A(x) x - b(x), F'(x^k) = A(x^k). With d = x1 - x0
    // (the Newton increment), (F(x0 + e d) - F(x0))/e - A0 d = O(e) when A0
    // is the exact Jacobian; a wrong term leaves an O(1) residue.
    PetscErrorCode ierr = PETSC_SUCCESS;
    const auto copy = [&ierr](::Vec v) {
      ::Vec c = PETSC_NULLPTR;
      ierr = VecDuplicate(v, &c);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecCopy(v, c);
      assert(ierr == PETSC_SUCCESS);
      return c;
    };
    // The field parts: x0 is the iterate, x1 the solution of the solve.
    ::Vec u0 = copy(m_uIt.getData()), u1 = copy(m_u.getSolution().getData());
    ::Vec d0 = copy(m_psiItD.getData()), d1 = copy(m_psiD.getSolution().getData());
    ::Vec o0 = copy(m_psiItO.getData()), o1 = copy(m_psiO.getSolution().getData());
    ::Vec p0 = copy(m_cfg.timeOrder == 2 ? m_pIt.getData() : m_pN.getData());
    ::Vec p1 = copy(m_p.getSolution().getData());

    ::Vec F0 = copy(b0), Fe = copy(b0), Jd = copy(b0), dx = copy(x0);
    ierr = MatMult(A0, x0, F0);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecAXPY(F0, -1.0, b0);
    assert(ierr == PETSC_SUCCESS);
    Real normD = 0.0;
    ierr = VecWAXPY(dx, -1.0, x0, x1);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecNorm(dx, NORM_2, &normD);
    assert(ierr == PETSC_SUCCESS);
    if (isRoot())
    {
      Alert::Info()
        << "[jacobian] step " << (m_step + 1) << ", t = " << m_t << " s: |d| = " << normD
        << "; err(e) = |(F(x0 + e d) - F(x0))/e - A0 d| / |A0 d| should fall like e"
        << Alert::Raise;
    }
    for (const Real e : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7})
    {
      combine(1.0 - e, u0, e, u1, m_uIt.getData());
      combine(1.0 - e, d0, e, d1, m_psiItD.getData());
      combine(1.0 - e, o0, e, o1, m_psiItO.getData());
      m_u.getSolution().setData(m_uIt.getData());
      m_psiD.getSolution().setData(m_psiItD.getData());
      m_psiO.getSolution().setData(m_psiItO.getData());
      combine(1.0 - e, p0, e, p1, m_p.getSolution().getData());
      ++m_psiRevision;
      m_flow.assemble();
      auto& ls = m_flow.getLinearSystem();
      ierr = MatMult(ls.getOperator(), ls.getSolution(), Fe);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecAXPY(Fe, -1.0, ls.getVector());
      assert(ierr == PETSC_SUCCESS);
      // FD = (Fe - F0)/e and Jd = A0 (x_e - x0)/e, x_e as gathered.
      ierr = VecAXPY(Fe, -1.0, F0);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecScale(Fe, 1.0 / e);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(dx, -1.0, x0, ls.getSolution());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecScale(dx, 1.0 / e);
      assert(ierr == PETSC_SUCCESS);
      ierr = MatMult(A0, dx, Jd);
      assert(ierr == PETSC_SUCCESS);
      // By block: the rows of u, p, psi_D, psi_O, at the trial offsets (the
      // block vector holds field k at offset_k + its global field index).
      std::array<Real, 4> bJ{}, bE{};
      {
        const ::Vec fields[4] = {
          m_uIt.getData(), m_pIt.getData(), m_psiItD.getData(), m_psiItO.getData()};
        PetscInt off = 0;
        for (size_t b = 0; b < 4; ++b)
        {
          PetscInt rb = 0, re = 0, n = 0;
          ierr = VecGetOwnershipRange(fields[b], &rb, &re);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGetSize(fields[b], &n);
          assert(ierr == PETSC_SUCCESS);
          MPI_Comm comm;
          ierr = PetscObjectGetComm((PetscObject)Jd, &comm);
          assert(ierr == PETSC_SUCCESS);
          ::IS is = PETSC_NULLPTR;
          ierr = ISCreateStride(comm, re - rb, off + rb, 1, &is);
          assert(ierr == PETSC_SUCCESS);
          ::Vec sj = PETSC_NULLPTR, se = PETSC_NULLPTR, diff = PETSC_NULLPTR;
          ierr = VecGetSubVector(Jd, is, &sj);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGetSubVector(Fe, is, &se);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecDuplicate(se, &diff);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecWAXPY(diff, -1.0, sj, se);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecNorm(sj, NORM_2, &bJ[b]);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecNorm(diff, NORM_2, &bE[b]);
          assert(ierr == PETSC_SUCCESS);
          VecDestroy(&diff);
          VecRestoreSubVector(Fe, is, &se);
          VecRestoreSubVector(Jd, is, &sj);
          ISDestroy(&is);
          off += n;
        }
      }
      Real nJ = 0.0, nJmax = 0.0;
      ierr = VecNorm(Jd, NORM_2, &nJ);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(Jd, NORM_INFINITY, &nJmax);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecAXPY(Fe, -1.0, Jd);
      assert(ierr == PETSC_SUCCESS);
      Real nE = 0.0, nEmax = 0.0;
      ierr = VecNorm(Fe, NORM_2, &nE);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecNorm(Fe, NORM_INFINITY, &nEmax);
      assert(ierr == PETSC_SUCCESS);
      if (isRoot())
      {
        Alert::Info() << "[jacobian] e = " << e
                      << "  err_2 = " << nE / std::max<Real>(nJ, 1e-300)
                      << "  err_inf = " << nEmax / std::max<Real>(nJmax, 1e-300)
                      << "  | by block u, p, psi_D, psi_O: "
                      << bE[0] / std::max<Real>(bJ[0], 1e-300) << ", "
                      << bE[1] / std::max<Real>(bJ[1], 1e-300) << ", "
                      << bE[2] / std::max<Real>(bJ[2], 1e-300) << ", "
                      << bE[3] / std::max<Real>(bJ[3], 1e-300) << Alert::Raise;
      }
    }
    for (::Vec* v : {&u0, &u1, &d0, &d1, &o0, &o1, &p0, &p1, &F0, &Fe, &Jd, &dx})
      VecDestroy(v);
    (void)ierr;
  }

  void AtriumViscoelastic::computeWallShear()
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
    // Projected over the cells with the nodal wall indicator, so that every
    // rank writes every local node consistently (never Region::Boundary).
    const size_t dim = m_mesh.getSpaceDimension();
    m_wssRec.setData(m_wssTrial.getSolution().getData());
    m_wss.project(
      VectorFunction(dim, [this, dim](const Point& p) -> Math::SpatialVector<Real> {
        Math::SpatialVector<Real> out(dim);
        if (!(m_wallArea.getValue(p) > m_wallAreaThreshold))
        {
          out.setZero();
          return out;
        }
        const auto w = m_wssRec.getValue(p);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          out(c) = w(c);
        return out;
      }));

    axpy(m_dt, m_wss.getData(), m_netShear.getData());
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void AtriumViscoelastic::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // TAWSS and OSI on the wall only, projected over the cells with the nodal
    // wall indicator (see m_wallArea; never Region::Boundary).
    const auto onWallNode = [this](const Point& p) {
      return m_wallArea.getValue(p) > m_wallAreaThreshold;
    };
    m_tawss.project(RealFunction([this, elapsed, onWallNode](const Point& p) -> Real {
      return onWallNode(p) ? m_absShear.getValue(p) / elapsed : 0.0;
    }));

    m_osi.project(RealFunction([this, onWallNode](const Point& p) -> Real {
      if (!onWallNode(p))
        return 0.0;
      const Real abs = m_absShear.getValue(p);
      if (abs <= 0.0)
        return 0.0;
      const auto net = m_netShear.getValue(p);
      const Real mag = std::sqrt(Math::dot(net, net));
      return std::clamp<Real>(0.5 * (1.0 - mag / abs), 0.0, 0.5);
    }));

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

    // The same indices restricted to the wall of the appendage.
    const auto onAppendage = [this](const Point& p) {
      return m_appendageArea.getValue(p) > m_appendageAreaThreshold;
    };
    m_tawssAppendage.project(RealFunction([this, onAppendage](const Point& p) -> Real {
      return onAppendage(p) ? m_tawss.getValue(p) : 0.0;
    }));
    m_osiAppendage.project(RealFunction([this, onAppendage](const Point& p) -> Real {
      return onAppendage(p) ? m_osi.getValue(p) : 0.0;
    }));
  }

  void AtriumViscoelastic::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    // Volumetric flow rates (m^3/s). n is outward: the venous flux is negated
    // to read as the inflow, the mitral flux is the outflow.
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_qIn = -m_flux(m_one);

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);

    // Mean pressures of the openings. That of the mitral plane is the atrial
    // pressure read by the ventricle and by the valve at the next step.
    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_inletPressure = m_flux(m_one) / m_inletMeasure;

    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_outletMeanPressure = m_flux(m_one) / m_outletMeasure;
  }

  void AtriumViscoelastic::writeCSVHeader()
  {
    m_csv
      << "t,cycle,pVenous,pOut,valveOpen,qIn,qOut,pInletMean,pOutletMean,maxU,maxSigma,"
      << "newtonIts,dpsi,maxTAWSS,maxOSI,maxTAWSSAppendage,maxOSIAppendage,dt,retries,"
      << "nMin,nMax,lambdaMax,ageMean,ageAppendage,"
      << "kinetic,enstrophy,viscousDissipation,polymerPower,elasticEnergy,"
         "polymerDissipation,"
      << "ossConvective,ossGradDiv,subscaleKinetic,psiNorm,psiSubscale,"
      << "wiMax,wiVolumeFraction,wiAppendage,shearAppendage,kineticAppendage\n";
  }

  void AtriumViscoelastic::writeCSVRow(int cycle)
  {
    // Reductions over the communicator, taken on all ranks first.
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();
    const Real tawssAppendage = m_tawssAppendage.max();
    const Real osiAppendage = m_osiAppendage.max();

    if (!isRoot())
      return;

    const Metrics& m = m_metrics;
    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << m_pOut << ','
          << (m_valveOpen ? 1 : 0) << ',' << m_qIn << ',' << m_qOut << ','
          << m_inletPressure << ',' << m_outletMeanPressure << ',' << m_speed << ','
          << m_stress << ',' << m_conformationIts << ',' << m_psiIncrement << ',' << tawss
          << ',' << osi << ',' << tawssAppendage << ',' << osiAppendage << ',' << m_dt
          << ',' << m_timing.retries << ',' << m_nMin << ',' << m_nMax << ','
          << m_lambdaMax << ',' << m_ageMean << ',' << m_ageAppendage << ',' << m.kinetic
          << ',' << m.enstrophy << ',' << m.viscousDissipation << ',' << m.polymerPower
          << ',' << m.elasticEnergy << ',' << m.polymerDissipation << ','
          << m.convectiveSubscale << ',' << m.gradDivSubscale << ',' << m.velocitySubscale
          << ',' << m.psiNorm << ',' << m.psiSubscale << ',' << m.wiMax << ','
          << m.wiVolume << ',' << m.wiAppendage << ',' << m.shearAppendage << ','
          << m.kineticAppendage << '\n';
    m_csv.flush();
  }

  int AtriumViscoelastic::run()
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
                    << " steps (" << nominalSteps
                    << " total), opening pressures ramped over " << m_cfg.rampCycles
                    << " cycle(s); XDMF every " << m_cfg.outputEvery
                    << " steps and at every cycle boundary" << Alert::Raise;
      Alert::Info()
        << "[run] formulation " << m_cfg.formulation << ", "
        << (m_cfg.dynamicSubscales ? "dynamic" : "quasi-static") << " subscales"
        << (m_cfg.formulation == "psi"
               ? ""
               : (m_cfg.subscaleTimeResidual
                     ? ", d_t of the FE field in the stress subscale residual"
                     : ", stress subscale residual without d_t of the FE field"))
        << (twoField() ? "; two-field (u, p), generalised Newtonian, convection at u* "
                         "(BDF extrapolation): one linear solve per step;"
                       : "; Newton on (u, p, psi) with the momentum convection at u* "
                         "(BDF extrapolation);")
        << " the PC is an LU of an earlier Jacobian, GMRES solves the current one"
        << Alert::Raise;
      if (m_cfg.adaptiveTimeStep)
      {
        Alert::Info() << "[run] step control: retry with dt x " << m_cfg.dtShrink
                      << " (at most " << m_cfg.maxStepRetries << " times, dt >= " << dtMin
                      << " s), growth x " << m_cfg.dtGrow << " after "
                      << m_cfg.dtGrowAfter << " steps of <= " << m_cfg.newtonEasy
                      << " Newton iterations, up to dt0 = " << m_dtNominal
                      << " s; Newton prediction " << (m_cfg.newtonPredict ? "on" : "off")
                      << "; lagged preconditioner over " << m_cfg.pcReuseSolves
                      << " solves" << Alert::Raise;
      }
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
        if (!updateBoundaryData(m_t))
          m_failure = "the lumped ventricle did not converge";
        else
        {
          prepareStep();
          if (solveFlow())
            break;
        }
        if (m_jacobianChecked)
          return 0;

        rejectStep();
        const Real next = m_cfg.dtShrink * m_dt;
        if (!m_cfg.adaptiveTimeStep || retries >= m_cfg.maxStepRetries ||
          next < (1.0 - TimeRoundoff) * dtMin)
        {
          // The last accepted state (rejectStep restored it), for locate_blowup.py.
          m_t = tStart;
          m_xdmf.write(m_t).flush();
          Alert::Exception() << "[flow] diverged at step " << (step + 1)
                             << ", t = " << m_t << " s, dt = " << m_dt
                             << " s: " << m_failure << Alert::Raise;
          return 1;
        }
        if (isRoot())
        {
          Alert::Warning() << "[step] " << (step + 1) << " rejected at t = " << m_t
                           << " s, dt = " << m_dt << " s: " << m_failure
                           << "; retrying with dt = " << next << " s" << Alert::Raise;
        }
        ++m_timing.retries;
        dt = next;
      }
      if (isOwens())
      {
        solveStructure();
        accumulateKinetics();
      }
      if (m_cfg.residenceTime)
        solveResidence();
      m_dtPrev = m_dt;
      ++m_history;
      ++step;

      computeWallShear();
      cycleElapsed += m_dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      if (m_cfg.metrics)
        computeMetrics();
      writeProbes();
      m_timing.fluxes = secondsSince(fluxStart);

      const bool endOfCycle = cycleEnd - m_t <= TimeRoundoff * m_dtNominal;
      if (endOfCycle)
      {
        closeCycle(cycleElapsed);
        resetStructure(cycleElapsed);
        cycleElapsed = 0.0;

        const Real maxTawss = m_tawss.max();
        const Real maxOsi = m_osi.max();
        const Real appendageTawss = m_tawssAppendage.max();
        const Real appendageOsi = m_osiAppendage.max();
        if (isRoot())
        {
          Alert::Info() << "[cycle " << (cycle + 1) << "/" << m_cfg.cycles
                        << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                        << "  appendage: maxTAWSS=" << appendageTawss
                        << " Pa  maxOSI=" << appendageOsi
                        << "  mean age=" << m_ageAppendage << " s" << Alert::Raise;
        }
      }

      writeCSVRow(cycle);

      bool output = endOfCycle;
      // Diagnostic dumps: every accepted step whose max|sigma| exceeds
      // dumpStress, at most dumpMax of them (locate_blowup.py --last N).
      if (m_cfg.dumpStress > 0.0 && m_stress > m_cfg.dumpStress &&
        m_dumps < m_cfg.dumpMax)
      {
        output = true;
        ++m_dumps;
      }
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

      if (isRoot() &&
        (m_step < 3 || m_step % 20 == 0 || endOfCycle || m_timing.retries > 0))
      {
        const Real elapsed = secondsSince(runStart);
        const Real remainingTime = std::max<Real>(endTime - m_t, 0.0);

        Alert::Info() << "---- Step " << step << "/" << nominalSteps << "  cycle "
                      << (cycle + 1) << "/" << m_cfg.cycles << "  t = " << m_t << " s"
                      << "  dt = " << m_dt << " s" << Alert::Raise;
        Alert::Info() << "[LA] p_pv=" << m_pIn << "  p_out=" << m_pOut << " Pa  valve "
                      << (m_valveOpen ? "open" : "closed") << "  qIn=" << m_qIn
                      << "  qOut=" << m_qOut << " m^3/s  max|u|=" << m_speed
                      << " m/s  max|sigma|=" << m_stress
                      << " Pa  newton=" << m_conformationIts
                      << (m_newtonPredicted ? "p" : "") << " (|dpsi|=" << m_psiIncrement
                      << ")" << Alert::Raise;
        if (m_cfg.metrics)
        {
          Alert::Info()
            << "[metrics] K=" << m_metrics.kinetic
            << " J  viscous=" << m_metrics.viscousDissipation
            << " W  polymer=" << m_metrics.polymerDissipation
            << " W  OSS convective=" << m_metrics.convectiveSubscale << " W ("
            << (m_metrics.convectiveSubscale /
                 std::max<Real>(
                   m_metrics.viscousDissipation + m_metrics.polymerDissipation, 1.0e-300))
            << " of viscous + polymer)  max Wi=" << m_metrics.wiMax
            << "  V(Wi>1)/V=" << m_metrics.wiVolume << "  |psi~|/|psi|="
            << (m_metrics.psiSubscale / std::max<Real>(m_metrics.psiNorm, 1.0e-300))
            << Alert::Raise;
        }
        Alert::Info() << "[timing] vms=" << m_timing.vms << "  asm=" << m_timing.assembly
                      << "  ksp=" << m_timing.solve << " (" << m_timing.kspIterations
                      << " its, " << m_timing.pcBuilds
                      << " PC builds)  wss=" << m_timing.shear
                      << "  flux=" << m_timing.fluxes << "  out=" << m_timing.output
                      << "  retries=" << m_timing.retries << "  total=" << m_timing.total
                      << " s  |  ETA "
                      << (elapsed / std::max<Real>(m_t, 1e-30) * remainingTime / 60.0)
                      << " min" << Alert::Raise;
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
    {
      m_csv.close();
      if (m_probeCSV.is_open())
        m_probeCSV.close();
    }

    return 0;
  }

  // ==========================================================================
  // Openings: venous pressure, lumped ventricle and mitral valve
  // ==========================================================================
  AtriumViscoelastic::VentricleType::Input AtriumViscoelastic::makeVentricleInput(
    const Config& cfg, const Real* atrialPressure)
  {
    // The input of Atrium/Fluid, parameter by parameter.
    VentricleType::Input input;

    input.rho = cfg.lv.rho;
    input.R0 = cfg.lv.R0;
    input.d0 = cfg.lv.d0;
    input.Es = cfg.lv.Es;
    input.mu = cfg.lv.mu;
    input.eta = cfg.lv.eta;
    input.alpha = cfg.lv.alpha;
    input.alphaR = cfg.lv.alphaR;
    input.k0 = cfg.lv.k0;
    input.sigma0 = cfg.lv.sigma0;

    input.Rp = cfg.lv.Rp;
    input.Cp = cfg.lv.Cp;
    input.Rd = cfg.lv.Rd;
    input.Cd = cfg.lv.Cd;

    input.proximalRadius = cfg.lv.proximalRadius;
    input.proximalLength = cfg.lv.proximalLength;
    input.distalRadius = cfg.lv.distalRadius;
    input.distalLength = cfg.lv.distalLength;

    input.mu_0 = cfg.lv.mu_0;
    input.mu_Inf = cfg.lv.mu_Inf;
    input.lambda = cfg.lv.lambda;
    input.n = cfg.lv.n;
    input.yasuda = cfg.lv.yasuda;

    input.Kat = cfg.lv.Kat;
    input.Kp = cfg.lv.Kp;
    input.Kar = cfg.lv.Kar;

    input.cavityCapacity = cfg.lv.cavityCapacity;
    input.localTolerance = cfg.lv.localTolerance;
    input.localMaxIterations = static_cast<size_t>(cfg.lv.localMaxIterations);
    input.localDamping = cfg.lv.localDamping;
    input.absRegularization = cfg.lv.absRegularization;

    input.pSv = [p = cfg.lv.systemicVenousPressure](Real) { return p; };
    // With feedback the ventricle reads the mean pressure of the mitral plane
    // at the last accepted step; otherwise its own atrial waveform.
    if (cfg.ventricularFeedback)
      input.pAt = [atrialPressure](Real) { return *atrialPressure; };
    else
    {
      input.pAt = [p = cfg.atrialPressure, T = cfg.period](
                    Real t) { return atrialWave(p, T, t); };
    }
    input.u = [a = cfg.activation, T = cfg.period](
                Real t) { return activationWave(a, T, t); };
    input.m0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                 lowEc = cfg.lv.relaxationM0LowEc,
                 highEc = cfg.lv.relaxationM0HighEc](Real ec) {
      if (highEc <= lowEc || ec >= highEc)
        return high;
      if (ec <= lowEc)
        return low;
      const Real s = (ec - lowEc) / (highEc - lowEc);
      return (1.0 - s) * low + s * high;
    };
    input.dm0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                  lowEc = cfg.lv.relaxationM0LowEc,
                  highEc = cfg.lv.relaxationM0HighEc](Real ec) {
      if (highEc <= lowEc || ec <= lowEc || ec >= highEc)
        return 0.0;
      return (high - low) / (highEc - lowEc);
    };

    {
      using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
      typename PassiveEnergy::Parameters hp;
      hp.C0 = cfg.lv.passiveC0;
      hp.C1 = cfg.lv.passiveC1;
      hp.C2 = cfg.lv.passiveC2;
      hp.C3 = cfg.lv.passiveC3;
      input.passiveEnergy = PassiveEnergy(hp);
    }

    return input;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::activationWave(
    const Activation& cfg, Real period, Real t)
  {
    const Real T = period;
    const Real tau = t - T * std::floor(t / T);

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };

    if (tau < cfg.tRampStart)
      return 0.0;
    if (tau < cfg.tRampEnd)
    {
      return cfg.positiveValue *
        ss((tau - cfg.tRampStart) / (cfg.tRampEnd - cfg.tRampStart));
    }
    if (tau < cfg.tPlateauEnd)
      return cfg.positiveValue;
    if (tau < cfg.tRelaxEnd)
    {
      return cfg.positiveValue +
        (cfg.negativeValue - cfg.positiveValue) *
        ss((tau - cfg.tPlateauEnd) / (cfg.tRelaxEnd - cfg.tPlateauEnd));
    }
    if (tau < cfg.tNegativeEnd)
      return cfg.negativeValue;
    return cfg.negativeValue *
      (1.0 - ss((tau - cfg.tNegativeEnd) / (T - cfg.tNegativeEnd)));
  }

  AtriumViscoelastic::Real AtriumViscoelastic::atrialWave(
    const AtrialPressure& cfg, Real period, Real t)
  {
    const Real T = period;
    const Real tau = t - T * std::floor(t / T);

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };
    auto ramp = [&](Real a, Real b, Real s) { return a + (b - a) * ss(s); };

    if (tau < cfg.t1)
      return ramp(cfg.minValue, cfg.maxValue, tau / cfg.t1);
    if (tau < cfg.t2)
      return cfg.maxValue;
    if (tau < cfg.t3)
      return ramp(cfg.maxValue, cfg.minValue, (tau - cfg.t2) / (cfg.t3 - cfg.t2));
    if (tau < cfg.t4)
      return ramp(cfg.minValue, cfg.secondThreshold, (tau - cfg.t3) / (cfg.t4 - cfg.t3));
    if (tau < cfg.t5)
      return cfg.secondThreshold;
    if (tau < cfg.t6)
      return ramp(cfg.secondThreshold, cfg.minValue, (tau - cfg.t5) / (cfg.t6 - cfg.t5));
    return cfg.minValue;
  }

  AtriumViscoelastic::Real AtriumViscoelastic::mitralResistance(bool open) const
  {
    // The coupled valve applies the conductances of the ventricular model, so
    // that the volume leaving the atrium is the volume entering the ventricle.
    if (m_cfg.useVentricularModel && m_cfg.ventricularFeedback)
      return 1.0 / (open ? m_cfg.lv.Kat : m_cfg.lv.Kp);
    return open ? m_cfg.mitralResistanceOpen : m_cfg.mitralResistanceClosed;
  }

  bool AtriumViscoelastic::updateBoundaryData(Real t)
  {
    const Real r = ramp(t);
    m_pIn = r *
      (m_cfg.inletPressureMean +
        m_cfg.inletPressureAmplitude * std::sin(2.0 * M_PI * t / m_period));

    if (!m_cfg.useVentricularModel)
    {
      m_pOut = r * m_cfg.prescribedOutletPressure;
      m_valveOpen = true;
      m_mitralZ = m_cfg.prescribedOutletResistance * m_outletMeasure;
      return true;
    }

    // The ventricle is advanced, once, to any time the flow has not reached
    // before; it is never taken back. A step that is rejected and retried with
    // a smaller dt reads the ventricular pressure by linear interpolation in
    // the interval already covered, so that no state of the lumped model has
    // to be restored. Without rejections this is the scheme of Atrium/Fluid:
    // one ventricular step per flow step, taken before it.
    constexpr Real Roundoff = 1.0e-12;
    while (m_lvTime1 < t - Roundoff * m_period)
    {
      const auto report = m_model.step(t - m_lvTime1);
      if (!report.converged)
        return false;
      m_lvTime0 = m_lvTime1;
      m_lvPressure0 = m_lvPressure1;
      m_lvTime1 = t;
      m_lvPressure1 = m_model.getState().pv;
    }
    const Real w = m_lvTime1 > m_lvTime0
      ? std::clamp<Real>((t - m_lvTime0) / (m_lvTime1 - m_lvTime0), 0.0, 1.0)
      : 1.0;
    m_pOut = (1.0 - w) * m_lvPressure0 + w * m_lvPressure1;
    // A coupled ventricle keeps its own pressure during the ramp: the valve
    // stays closed while the ramped atrium is below it (Atrium/Fluid).
    if (!m_cfg.ventricularFeedback)
      m_pOut *= r;

    // Mitral valve, switched on the mean pressure of the mitral plane at the
    // last accepted step: the test of the ventricular model, p_v <= p_at.
    m_valveOpen = m_outletMeanPressure > m_pOut;
    m_mitralZ = mitralResistance(m_valveOpen) * m_outletMeasure;
    return true;
  }

  // ==========================================================================
  // Diagnostics
  // ==========================================================================
  bool AtriumViscoelastic::inAppendage(const Point& p) const
  {
    const auto& a = m_cfg.appendage;
    if (!(a.radius > 0.0))
      return false;
    const auto& x = p.getPhysicalCoordinates();
    const Real dx = x(0) - a.center[0], dy = x(1) - a.center[1], dz = x(2) - a.center[2];
    return dx * dx + dy * dy + dz * dz < a.radius * a.radius;
  }

  template <class F>
  AtriumViscoelastic::Real AtriumViscoelastic::integrate(F&& f)
  {
    // The quadrature of the log-conformation integrals, which is also that of
    // the subscales stored at the quadrature points.
    auto sweep = Integral(
      PointwiseScalar([&f](const Point& p) -> Real { return f(p); }, 2), m_qFlux);
    sweep.setOrder(static_cast<size_t>(std::max(1, m_cfg.psiQuadratureOrder)));
    m_flux = sweep;
    m_flux.assemble();
    return m_flux(m_one);
  }

  void AtriumViscoelastic::computeMetrics()
  {
    // All fields at the accepted step n+1 (m_uOld, m_psiOld hold it). The
    // projections Pi[.] and the stabilisation parameters are those of the step
    // just taken (step n), as in the forms.
    const Real rho = m_cfg.rho;
    const auto& ob = m_cfg.oldroydB;
    const bool polymer = !twoField();
    const Real s = ob.etaP / lambda0();
    Metrics m;

    const auto gradient = [this](const Point& p) -> Eigen::Matrix3d {
      const auto J = Jacobian(m_uOld).getValue(p);
      Eigen::Matrix3d L;
      for (int a = 0; a < 3; ++a)
      {
        for (int b = 0; b < 3; ++b)
          L(a, b) = J(a, b);
      }
      return L;
    };
    const auto shearRate = [](const Eigen::Matrix3d& L) -> Real {
      return std::sqrt(2.0 * (0.5 * (L + L.transpose())).squaredNorm());
    };
    const auto velocity = [this](const Point& p) -> Eigen::Vector3d {
      const auto u = m_uOld.getValue(p);
      return Eigen::Vector3d(u(0), u(1), u(2));
    };
    const auto relaxation = [this, &ob](const Point& p) -> Real {
      return isOwens() ? m_lambdaCell.getValue(p) : ob.lambda;
    };
    const auto expPsi = [this](const Point& p) -> Eigen::Matrix3d {
      return logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp();
    };

    m.kinetic = 0.5 * rho *
      integrate([&](const Point& p) -> Real { return velocity(p).squaredNorm(); });
    m.enstrophy = 0.5 * integrate([&](const Point& p) -> Real {
      const Eigen::Matrix3d L = gradient(p);
      return (L - L.transpose()).squaredNorm() / 2.0; // |curl u|^2
    });
    m.viscousDissipation = integrate([&](const Point& p) -> Real {
      const Real gd = shearRate(gradient(p));
      return solventViscosity(gd) * gd * gd;
    });
    m.convectiveSubscale = rho * rho * integrate([&](const Point& p) -> Real {
      const Eigen::Vector3d conv = gradient(p) * velocity(p);
      const auto pi = m_piConv.get().getValue(p);
      return m_tauK.getValue(p) *
        (conv - Eigen::Vector3d(pi(0), pi(1), pi(2))).squaredNorm();
    });
    m.gradDivSubscale = integrate([&](const Point& p) -> Real {
      const Real div = gradient(p).trace();
      return m_alpha2.getValue(p) * div * (div - m_piDiv.get().getValue(p));
    });
    m.velocitySubscale = 0.5 * rho * integrate([&](const Point& p) -> Real {
      const SubscaleStore::Entry& e = subscales().at(p);
      // The stored subscales are velocities (m/s): u~ = u1 + u2 + u3.
      return (e.u1 + e.u2 + e.u3).squaredNorm();
    });
    if (polymer)
    {
      m.polymerPower = s * integrate([&](const Point& p) -> Real {
        const Eigen::Matrix3d L = gradient(p);
        const Eigen::Matrix3d E = 0.5 * (L + L.transpose());
        return (expPsi(p) - Eigen::Matrix3d::Identity()).cwiseProduct(E).sum();
      });
      m.elasticEnergy = 0.5 * s *
        integrate([&](const Point& p) -> Real { return expPsi(p).trace() - 3.0; });
      m.polymerDissipation = 0.5 * s * integrate([&](const Point& p) -> Real {
        const Real tr = expPsi(p).trace();
        return pttFactor(tr) * (tr - 3.0) / relaxation(p);
      });
      m.psiNorm = std::sqrt(integrate([&](const Point& p) -> Real {
        return unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p)).squaredNorm();
      }) /
        m_volume);
      m.psiSubscale = std::sqrt(integrate([&](const Point& p) -> Real {
        return subscales().at(p).sig.squaredNorm();
      }) /
        m_volume);
      // Local Weissenberg number lambda gd, by cell.
      m_wiCell.project(RealFunction(
        [&](const Point& p) -> Real { return relaxation(p) * shearRate(gradient(p)); }));
      m.wiMax = m_wiCell.max();
      m.wiVolume = integrate([&](const Point& p) -> Real {
        return m_wiCell.getValue(p) > 1.0 ? 1.0 : 0.0;
      }) /
        m_volume;
    }
    if (m_appendageVolume > 0.0)
    {
      m.shearAppendage = integrate([&](const Point& p) -> Real {
        return inAppendage(p) ? shearRate(gradient(p)) : 0.0;
      }) /
        m_appendageVolume;
      m.kineticAppendage = 0.5 * rho * integrate([&](const Point& p) -> Real {
        return inAppendage(p) ? velocity(p).squaredNorm() : 0.0;
      });
      if (polymer)
      {
        m.wiAppendage = integrate([&](const Point& p) -> Real {
          return inAppendage(p) ? m_wiCell.getValue(p) : 0.0;
        }) /
          m_appendageVolume;
      }
    }
    m_metrics = m;
  }

  void AtriumViscoelastic::setupProbes()
  {
    // "x,y,z;x,y,z;...", metres.
    m_probes.clear();
    m_probeWeights.clear();
    std::stringstream list(m_cfg.probes);
    std::string item;
    while (std::getline(list, item, ';'))
    {
      if (item.find_first_not_of(" \t") == std::string::npos)
        continue;
      std::replace(item.begin(), item.end(), ',', ' ');
      std::stringstream coordinates(item);
      std::array<Real, 3> x{};
      if (!(coordinates >> x[0] >> x[1] >> x[2]))
      {
        throw std::runtime_error(
          "probes: expected \"x,y,z;x,y,z;...\", got \"" + m_cfg.probes + "\".");
      }
      m_probes.push_back(x);
    }
    if (m_probes.empty())
      return;
    if (!(m_cfg.probeRadius > 0.0))
      throw std::runtime_error("probe_radius must be positive.");
    const Real r2 = m_cfg.probeRadius * m_cfg.probeRadius;
    for (const auto& c : m_probes)
    {
      const Real w = integrate([&](const Point& p) -> Real {
        const auto& x = p.getPhysicalCoordinates();
        const Real d2 = (x(0) - c[0]) * (x(0) - c[0]) + (x(1) - c[1]) * (x(1) - c[1]) +
          (x(2) - c[2]) * (x(2) - c[2]);
        return d2 < 16.0 * r2 ? std::exp(-0.5 * d2 / r2) : 0.0;
      });
      if (!(w > 0.0))
      {
        throw std::runtime_error("probes: no fluid within 4 probe radii of (" +
          std::to_string(c[0]) + ", " + std::to_string(c[1]) + ", " +
          std::to_string(c[2]) + ") m.");
      }
      m_probeWeights.push_back(w);
    }
    if (isRoot())
    {
      const auto stem =
        std::filesystem::path(m_cfg.csvPath).replace_extension("").string();
      m_probeCSV.open(stem + "_probes.csv");
      if (!m_probeCSV)
        throw std::runtime_error("Failed to open " + stem + "_probes.csv");
      m_probeCSV << "t";
      for (size_t i = 0; i < m_probes.size(); ++i)
        m_probeCSV << ",ux" << i << ",uy" << i << ",uz" << i << ",trSigma" << i;
      m_probeCSV << '\n';
      Alert::Info() << "[probes] " << m_probes.size() << " Gaussian probes of radius "
                    << m_cfg.probeRadius << " m" << Alert::Raise;
    }
  }

  void AtriumViscoelastic::writeProbes()
  {
    if (m_probes.empty())
      return;
    const Real r2 = m_cfg.probeRadius * m_cfg.probeRadius;
    const Real s = m_cfg.oldroydB.etaP / lambda0();
    const bool polymer = !twoField();
    std::vector<Real> row;
    row.reserve(4 * m_probes.size());
    for (size_t i = 0; i < m_probes.size(); ++i)
    {
      const auto& c = m_probes[i];
      const auto weight = [&c, r2](const Point& p) -> Real {
        const auto& x = p.getPhysicalCoordinates();
        const Real d2 = (x(0) - c[0]) * (x(0) - c[0]) + (x(1) - c[1]) * (x(1) - c[1]) +
          (x(2) - c[2]) * (x(2) - c[2]);
        return d2 < 16.0 * r2 ? std::exp(-0.5 * d2 / r2) : 0.0;
      };
      for (Index k = 0; k < 3; ++k)
      {
        row.push_back(integrate([&](const Point& p) -> Real {
          const Real w = weight(p);
          return w > 0.0 ? w * m_uOld.getValue(p)(k) : 0.0;
        }) /
          m_probeWeights[i]);
      }
      row.push_back(polymer ? s * integrate([&](const Point& p) -> Real {
        const Real w = weight(p);
        if (!(w > 0.0))
          return 0.0;
        return w *
          (logEig(unpack(m_psiOldD.getValue(p), m_psiOldO.getValue(p))).exp().trace() -
            3.0);
      }) / m_probeWeights[i]
                            : 0.0);
    }
    if (!isRoot())
      return;
    m_probeCSV << m_t;
    for (const Real v : row)
      m_probeCSV << ',' << v;
    m_probeCSV << '\n';
    m_probeCSV.flush();
  }

  void AtriumViscoelastic::accumulateKinetics()
  {
    // int b dt and int b n_st (n_st - 1) dt at the nodes (the second is int a
    // dt wherever b is not floored), for Owens::resetCycle.
    if (m_cfg.owens.equilibrium || m_cfg.owens.resetCycle <= 0 ||
      m_cyclesDone >= m_cfg.owens.resetCycle)
      return;
    m_rateWork.project(RealFunction([this](const Point& p) -> Real {
      const Real gd = shearRateOf(m_uOld, p);
      const Real nst = owensNst(gd);
      return owensB(gd) * nst * (nst - 1.0);
    }));
    axpy(m_dt, m_rateWork.getData(), m_aSum.getData());
    m_rateWork.project(RealFunction(
      [this](const Point& p) -> Real { return owensB(shearRateOf(m_uOld, p)); }));
    axpy(m_dt, m_rateWork.getData(), m_bSum.getData());
  }

  void AtriumViscoelastic::resetStructure(Real elapsed)
  {
    ++m_cyclesDone;
    if (!isOwens() || m_cfg.owens.equilibrium || m_cfg.owens.resetCycle <= 0 ||
      elapsed <= 0.0)
      return;
    if (m_cyclesDone == m_cfg.owens.resetCycle)
    {
      // <b (n - n_st)(n + n_st - 1)> = 0 over the cycle: n (n - 1) = <b n_st
      // (n_st - 1)>/<b>, the positive root. Where b vanished all cycle the
      // structure has no rate and n is kept.
      const Real nMax = owensNst(0.0);
      m_rateWork.project(RealFunction([this, nMax](const Point& p) -> Real {
        const Real b = m_bSum.getValue(p);
        if (!(b > 0.0))
          return m_nOld.getValue(p);
        const Real q = std::max<Real>(m_aSum.getValue(p) / b, 0.0);
        return std::clamp<Real>(0.5 * (1.0 + std::sqrt(1.0 + 4.0 * q)), 1.0, nMax);
      }));
      m_nOld.setData(m_rateWork.getData());
      m_nNm1.setData(m_nOld.getData());
      m_nHist.setData(m_nOld.getData());
      m_n.getSolution().setData(m_nOld.getData());
      updateStructureFields();
      if (isRoot())
      {
        Alert::Info()
          << "[structure] n set to the periodic state of slow kinetics of cycle "
          << m_cyclesDone << ": n in [" << m_nMin << ", " << m_nMax << "]"
          << Alert::Raise;
      }
    }
    // The next cycle accumulates its own means.
    m_aSum = Real(0);
    m_bSum = Real(0);
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
      using Simulation = Rodin::Examples::ViscoelasticFluids::AtriumViscoelastic;
      Simulation::Config cfg;

      const auto getString = [](const char* name, std::string& out) {
        char buffer[1024];
        PetscBool set = PETSC_FALSE;
        PetscOptionsGetString(
          PETSC_NULLPTR, PETSC_NULLPTR, name, buffer, sizeof(buffer), &set);
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

      const auto getPoint = [](const char* name, std::array<Rodin::Real, 3>& out) {
        PetscReal values[3] = {0.0, 0.0, 0.0};
        PetscInt count = 3;
        PetscBool set = PETSC_FALSE;
        PetscOptionsGetRealArray(
          PETSC_NULLPTR, PETSC_NULLPTR, name, values, &count, &set);
        if (set && count != 3)
          throw std::runtime_error(std::string(name) + " expects x,y,z.");
        if (set)
          out = {values[0], values[1], values[2]};
      };

      getString("-av_mesh", cfg.meshPath);
      getString("-av_xdmf", cfg.xdmfBasename);
      getString("-av_csv", cfg.csvPath);
      getString("-av_fluid", cfg.fluid);
      getString("-av_formulation", cfg.formulation);
      getReal("-av_c1", cfg.c1);
      getReal("-av_c2", cfg.c2);
      getReal("-av_c3", cfg.c3);
      getReal("-av_c4", cfg.c4);
      getBool("-av_dynamic_subscales", cfg.dynamicSubscales);
      getBool("-av_adjoint_reaction", cfg.adjointReaction);
      getBool("-av_subscale_time_residual", cfg.subscaleTimeResidual);

      getReal("-av_mesh_scale", cfg.meshScale);
      getReal("-av_reference_length", cfg.referenceLength);
      getReal("-av_reference_velocity", cfg.referenceVelocity);
      getReal("-av_rho", cfg.rho);
      getReal("-av_eta_s", cfg.oldroydB.etaS);
      getReal("-av_eta_p", cfg.oldroydB.etaP);
      getReal("-av_lambda", cfg.oldroydB.lambda);
      getReal("-av_lambda0_factor", cfg.oldroydB.lambda0Factor);
      getReal("-av_lambda0_min", cfg.oldroydB.lambda0Min);
      getReal("-av_ptt_epsilon", cfg.oldroydB.pttEpsilon);
      getReal("-av_de", cfg.deborah);

      getReal("-av_period", cfg.period);
      getReal("-av_p_in", cfg.inletPressureMean);
      getReal("-av_p_in_amplitude", cfg.inletPressureAmplitude);
      getReal("-av_ramp_cycles", cfg.rampCycles);
      getReal("-av_inlet_impedance", cfg.inletImpedance);
      getReal("-av_inlet_tangential_damping", cfg.inletTangentialDamping);
      getReal("-av_inlet_backflow", cfg.inletBackflowStabilization);
      getReal("-av_outlet_backflow", cfg.outletBackflowStabilization);
      getReal("-av_inlet_shear_rate", cfg.inletShearRate);
      getBool("-av_lv_model", cfg.useVentricularModel);
      getBool("-av_lv_feedback", cfg.ventricularFeedback);
      getReal("-av_mitral_r_open", cfg.mitralResistanceOpen);
      getReal("-av_mitral_r_closed", cfg.mitralResistanceClosed);
      getReal("-av_p_out", cfg.prescribedOutletPressure);
      getReal("-av_r_out", cfg.prescribedOutletResistance);
      getReal("-av_lv_kat", cfg.lv.Kat);
      getReal("-av_lv_kp", cfg.lv.Kp);

      getInt("-av_steps_per_cycle", cfg.stepsPerCycle);
      getReal("-av_dt", cfg.timeStep);
      getInt("-av_cycles", cfg.cycles);
      getInt("-av_time_order", cfg.timeOrder);
      getReal("-av_eta_ref", cfg.etaRef);
      getReal("-av_cy_eta_zero", cfg.carreauYasuda.etaZero);
      getReal("-av_cy_eta_inf", cfg.carreauYasuda.etaInf);
      getReal("-av_cy_lambda", cfg.carreauYasuda.lambda);
      getReal("-av_cy_n", cfg.carreauYasuda.n);
      getReal("-av_cy_a", cfg.carreauYasuda.a);
      getInt("-av_bdf2_extrapolate", cfg.bdf2Extrapolate);
      getInt("-av_output_every", cfg.outputEvery);

      getReal("-av_vms_scale", cfg.vmsScale);
      getReal("-av_graddiv_scale", cfg.gradDivScale);
      getReal("-av_pressure_scale", cfg.pressureScale);
      getReal("-av_stress_div_scale", cfg.stressDivScale);
      getReal("-av_stress_scale", cfg.stressScale);
      getInt("-av_conformation_its", cfg.conformationIterations);
      getReal("-av_conformation_tol", cfg.conformationTolerance);
      getReal("-av_newton_max_step", cfg.newtonMaxStep);
      getReal("-av_max_velocity_factor", cfg.maxVelocityFactor);
      getReal("-av_dump_stress", cfg.dumpStress);
      getBool("-av_wall_kinematics", cfg.wallKinematics);
      getBool("-av_two_field", cfg.twoField);
      getBool("-av_residence_time", cfg.residenceTime);
      getInt("-av_jacobian_check", cfg.jacobianCheck);
      getInt("-av_dump_max", cfg.dumpMax);
      getBool("-av_vms", cfg.useVMS);
      getBool("-av_inlet_conformation", cfg.inletConformation);
      getBool("-av_trace_newton", cfg.traceNewton);
      getBool("-av_newton_d2", cfg.newtonD2);
      getInt("-av_psi_order", cfg.psiQuadratureOrder);
      getBool("-av_adaptive_dt", cfg.adaptiveTimeStep);
      getInt("-av_max_step_retries", cfg.maxStepRetries);
      getReal("-av_dt_shrink", cfg.dtShrink);
      getReal("-av_dt_grow", cfg.dtGrow);
      getInt("-av_dt_grow_after", cfg.dtGrowAfter);
      getInt("-av_newton_easy", cfg.newtonEasy);
      getInt("-av_newton_hard", cfg.newtonHard);
      getBool("-av_newton_predict", cfg.newtonPredict);
      getReal("-av_newton_predict_safety", cfg.newtonPredictSafety);
      getInt("-av_newton_calibrate_every", cfg.newtonCalibrateEvery);
      getInt("-av_pc_reuse_solves", cfg.pcReuseSolves);
      getReal("-av_pc_reuse_growth", cfg.pcReuseGrowth);
      getBool("-av_outlet_conformation_inflow", cfg.outletConformationInflow);
      getBool("-av_psi_skew", cfg.psiSkewAdvection);
      getBool("-av_owens_equilibrium", cfg.owens.equilibrium);
      getReal("-av_owens_lambda_h", cfg.owens.lambdaH);
      getReal("-av_owens_eta_plasma", cfg.owens.etaPlasma);
      getReal("-av_owens_a0", cfg.owens.a0);
      getReal("-av_owens_a_crit", cfg.owens.aCrit);
      getReal("-av_owens_gd_crit", cfg.owens.gdCrit);
      getReal("-av_owens_gd_max", cfg.owens.gdMax);
      getBool("-av_owens_b_floor", cfg.owens.bFloor);
      getReal("-av_owens_n_inlet", cfg.owens.nInlet);
      getReal("-av_owens_epsilon", cfg.owens.epsilon);
      getReal("-av_owens_n0", cfg.owens.nInitial);
      getInt("-av_owens_reset_cycle", cfg.owens.resetCycle);

      getPoint("-av_appendage_center", cfg.appendage.center);
      getReal("-av_appendage_radius", cfg.appendage.radius);
      getBool("-av_metrics", cfg.metrics);
      getString("-av_probes", cfg.probes);
      getReal("-av_probe_radius", cfg.probeRadius);

      // The XDMF and CSV writers do not create directories.
      if (world.rank() == 0)
      {
        for (const std::string& path : {cfg.xdmfBasename, cfg.csvPath})
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
