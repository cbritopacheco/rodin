/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RmuAxi.cpp
 * @brief Truncation benchmark of reduced outlet closures for a
 *        generalised-Newtonian fluid: pulsatile, pressure-driven flow through
 *        a stepped axisymmetric pipe, resolved in full or truncated and
 *        closed by a reduced resistance.
 *
 * Geometry. A series of coaxial tubes (x axial, y = r >= 0), radii R_k and
 * lengths L_k, joined by abrupt contractions. A symmetric vascular tree is a
 * series of generations, each made of N_k identical tubes in parallel; the
 * WRMS closure of a generation depends on its wall stress only, which the
 * identical tubes share, so a series of tubes is the axisymmetric surrogate
 * of a symmetric tree.
 *
 * Problem. The full geometry is resolved with the Carreau-Yasuda fluid,
 *
 *   rho (d_t u + u . grad u) - div(2 mu(gd) D(u)) + grad p = 0,  div u = 0,
 *
 * written on the meridian half-plane (every integral weighted by r, hoop
 * strain u_r/r included in D and gd). The inlet carries the traction
 * -p_in(t) n with p_in(t) = Pbar ramp(t) [1 + A sin(omega t)], A > 1 giving
 * flow reversal; the outlet of the full geometry carries p = 0.
 *
 * The truncated geometry stops at x_T inside one of the tubes; the
 * downstream tubes are replaced by a reduced resistance, the outlet carrying
 * the uniform traction -R(Q) Q n with Q the outlet flux:
 *
 *   tree     R(Q) = sum_k 8 mu_ap(4|Q|/(pi R_k^3)) L_k / (pi R_k^4)
 *            (exact WRMS series of the removed tubes);
 *   single   one compartment of radius R_s and an effective length fixed so
 *            that R(Q_0) equals the tree value: the reduction of the paper;
 *   rest     the constant R_tree(Q_0), the closure frozen at the a-priori
 *            resting operating point, times Config::constantScale;
 *   plateau  the constant sum_k 8 mu_inf L_k/(pi R_k^4).
 *
 * The resting operating point is fixed before solving, as in the paper: the
 * narrowest tube of the full geometry carries the mean apparent shear rate
 * gammaRef, Q_0 = pi R_min^3 gammaRef / 4, and Pbar is the WRMS pressure
 * budget of the whole series at Q_0. Both runs therefore see the same
 * forcing, and the full run is the reference the closures are measured
 * against.
 *
 * Time discretisation. BDF1. Each step is solved by Picard iterations on
 * the advecting velocity, the Carreau-Yasuda viscosity and the outlet
 * resistance, until the outlet flux and the velocity change converge; no
 * coefficient is lagged at convergence. The uniform outlet traction is
 * preconditioned by the local implicit term beta R A (u.n)(v.n) minus the same
 * term at the previous iterate, which cancels at convergence.
 */
#include <algorithm>
#include <chrono>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

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

namespace Rodin::Examples::Heart::RmuAxi
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  namespace
  {
    constexpr int RootRank = 0;

    constexpr Real Pi = std::numbers::pi_v<Real>;

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

    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }

    /// @brief e_r e_r: selects u_r in the axis penalty.
    MatrixFunction<Math::Matrix<Real>> radialProjector()
    {
      Math::Matrix<Real> P = Math::Matrix<Real>::Zero(2, 2);
      P(1, 1) = 1;
      return MatrixFunction(P);
    }

    /// @brief e_j. VectorFunction keeps a reference, hence the static storage.
    VectorFunction<Math::Vector<Real>> unit(size_t j)
    {
      static const Math::Vector<Real> e[2] = {
        Math::Vector<Real>::Unit(2, 0), Math::Vector<Real>::Unit(2, 1)};
      return VectorFunction(e[j]);
    }

    /// @brief A scalar field given by a callable, evaluated pointwise, with a
    ///        declared quadrature order.
    template <class F>
    class PointwiseScalar final : public RealFunctionBase<PointwiseScalar<F>>
    {
      public:
        using Parent = RealFunctionBase<PointwiseScalar<F>>;

        explicit PointwiseScalar(F f, size_t order)
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

    /// @brief One tube of the series: radius and length, in units of the
    ///        inlet radius.
    struct Tube
    {
      Real radius = 1.0;
      Real length = 1.0;
    };

    std::vector<Tube> parseTubes(const std::string& text)
    {
      std::vector<Tube> tubes;
      std::stringstream items(text);
      std::string item;
      while (std::getline(items, item, ','))
      {
        if (item.empty())
          continue;
        const auto colon = item.find(':');
        if (colon == std::string::npos)
          throw std::runtime_error("A tube is R:L (units of the inlet radius): " + item);
        tubes.push_back(
          Tube{std::stod(item.substr(0, colon)), std::stod(item.substr(colon + 1))});
      }
      return tubes;
    }
  } // namespace

  /**
   * @brief Hydraulic constitutive curve mu_ap(gd_a) of a Carreau-Yasuda fluid,
   *        tabulated once and log-log interpolated.
   *
   * @details mu_ap(tau_w) = tau_w^4/(4 I(tau_w)), gd_a = tau_w/mu_ap =
   *          4Q/(pi R^3), with I(tau_w) = int_0^{gd_w} gd^3 mu^2 (mu + gd mu')
   *          dgd. The same construction as the coronary outlet table.
   */
  class HydraulicCurve
  {
    public:
      struct CarreauYasuda
      {
        Real mu0 = 0.0712616;
        Real muInf = 0.0041;
        Real lambda = 9.2780;
        Real n = 0.3572;
        Real a = 1.8873;
      };

      /// Log-spaced table nodes in tau_w, Simpson panels per node and the
      /// bracket of the wall-stress table (Pa). The bounds cover every stress
      /// of the benchmark; outside them the Carreau-Yasuda plateaus are the
      /// correct continuation, hence clamping.
      static constexpr int TableNodes = 241;
      static constexpr int IntegralSteps = 2000;
      static constexpr Real TauMin = 1.0e-6;
      static constexpr Real TauMax = 1.0e4;
      static constexpr int ShearIterations = 200;
      static constexpr Real ShearTolerance = 1.0e-12;

      explicit HydraulicCurve(const CarreauYasuda& cy)
        : m_cy(cy)
      {
        const auto mu = [this](Real g) { return viscosity(g); };
        const auto dmu = [this](Real g) -> Real {
          if (g <= 0.0)
            return 0.0;
          const Real base = 1.0 + std::pow(m_cy.lambda * g, m_cy.a);
          return (m_cy.mu0 - m_cy.muInf) * (m_cy.n - 1.0) *
            std::pow(base, (m_cy.n - 1.0 - m_cy.a) / m_cy.a) *
            std::pow(m_cy.lambda, m_cy.a) * std::pow(g, m_cy.a - 1.0);
        };

        // gd_w(tau_w): tau = mu(gd) gd is increasing; bisection-safeguarded
        // Newton on log gd.
        const auto shearAt = [&](Real tauW) -> Real {
          Real lo = std::log(1e-14);
          Real hi = std::log(std::max<Real>(tauW / m_cy.muInf, 1e-12)) + 1.0;
          Real s = 0.5 * (lo + hi);
          for (int it = 0; it < ShearIterations; ++it)
          {
            const Real g = std::exp(s);
            const Real f = mu(g) * g - tauW;
            if (f < 0.0)
              lo = s;
            else
              hi = s;
            const Real df = (mu(g) + g * dmu(g)) * g;
            Real next = (std::abs(df) > 0.0) ? s - f / df : 0.5 * (lo + hi);
            if (!(next > lo && next < hi))
              next = 0.5 * (lo + hi);
            const Real step = std::abs(next - s);
            s = next;
            if (step < ShearTolerance)
              break;
          }
          return std::exp(s);
        };

        const auto integral = [&](Real gW) -> Real {
          const Real h = gW / static_cast<Real>(IntegralSteps);
          const auto f = [&](Real g) -> Real {
            if (g <= 0.0)
              return 0.0;
            const Real m = mu(g);
            return g * g * g * m * m * (m + g * dmu(g));
          };
          Real sum = f(0.0) + f(gW);
          for (int i = 1; i < IntegralSteps; ++i)
            sum += ((i % 2 == 1) ? 4.0 : 2.0) * f(static_cast<Real>(i) * h);
          return sum * h / 3.0;
        };

        for (int i = 0; i < TableNodes; ++i)
        {
          const Real tauW = std::exp(std::log(TauMin) +
            (std::log(TauMax) - std::log(TauMin)) * static_cast<Real>(i) /
              static_cast<Real>(TableNodes - 1));
          const Real I = integral(shearAt(tauW));
          if (!(I > 0.0) || !std::isfinite(I))
            continue;
          const Real muAp = std::pow(tauW, 4.0) / (4.0 * I);
          const Real gA = tauW / muAp;
          if (!m_logGamma.empty() && std::log(gA) <= m_logGamma.back())
            continue;
          m_logGamma.push_back(std::log(gA));
          m_logMu.push_back(std::log(muAp));
        }

        if (m_logGamma.size() < 2)
          throw std::runtime_error("Degenerate Carreau-Yasuda hydraulic curve.");
      }

      /// @brief Local Carreau-Yasuda viscosity mu(gd).
      Real viscosity(Real g) const
      {
        return m_cy.muInf +
          (m_cy.mu0 - m_cy.muInf) *
          std::pow(1.0 + std::pow(m_cy.lambda * std::max<Real>(g, 0.0), m_cy.a),
            (m_cy.n - 1.0) / m_cy.a);
      }

      /// @brief d mu / d gd of the local law.
      Real derivative(Real g) const
      {
        if (g <= 0.0)
          return 0.0;
        const Real base = 1.0 + std::pow(m_cy.lambda * g, m_cy.a);
        return (m_cy.mu0 - m_cy.muInf) * (m_cy.n - 1.0) *
          std::pow(base, (m_cy.n - 1.0 - m_cy.a) / m_cy.a) *
          std::pow(m_cy.lambda, m_cy.a) * std::pow(g, m_cy.a - 1.0);
      }

      /// @brief Hydraulic (effective) viscosity at the apparent shear rate.
      Real operator()(Real gA) const
      {
        const Real lg = std::log(std::max<Real>(gA, 1e-300));
        if (lg <= m_logGamma.front())
          return std::exp(m_logMu.front());
        if (lg >= m_logGamma.back())
          return std::exp(m_logMu.back());
        const auto it = std::upper_bound(m_logGamma.begin(), m_logGamma.end(), lg);
        const std::size_t i = static_cast<std::size_t>(it - m_logGamma.begin()) - 1;
        const Real w = (lg - m_logGamma[i]) / (m_logGamma[i + 1] - m_logGamma[i]);
        return std::exp(m_logMu[i] + w * (m_logMu[i + 1] - m_logMu[i]));
      }

      /// @brief Resistance of one tube (radius R, length L, SI) at flow Q.
      Real resistance(Real R, Real L, Real Q) const
      {
        const Real gA = 4.0 * std::abs(Q) / (Pi * R * R * R);
        return 8.0 * (*this)(gA)*L / (Pi * R * R * R * R);
      }

      const CarreauYasuda& getParameters() const
      {
        return m_cy;
      }

    private:
      CarreauYasuda m_cy;
      std::vector<Real> m_logGamma;
      std::vector<Real> m_logMu;
  };

  class Benchmark
  {
    public:
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;
      using Attribute = Rodin::Geometry::Attribute;
      using AttributeSet = Rodin::FlatSet<Attribute>;

      using VelocityFESType =
        Rodin::Variational::H1<2, Rodin::Math::SpatialVector<Real>, MeshType>;
      using PressureFESType = Rodin::Variational::H1<1, Real, MeshType>;

      using VelocityGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VelocityFESType>;
      using PressureGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<PressureFESType>;

      using VelocityTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VelocityGridFunctionType,
          VelocityFESType>;
      using VelocityTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VelocityFESType>;
      using PressureTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<PressureGridFunctionType,
          PressureFESType>;
      using PressureTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<PressureFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType, VelocityTrialFunctionType,
          PressureTrialFunctionType, VelocityTestFunctionType, PressureTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VelocityTrialFunctionType, VelocityTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<PressureFESType, ::Vec>;

      /// @brief Boundary labels written by make_rmu_axi_mesh.py.
      struct Labels
      {
        Attribute inlet = 1;
        Attribute outlet = 2;
        Attribute wall = 3;
        Attribute axis = 4;
        Attribute parentWall = 5;
        Attribute stubWall = 6;
      };

      enum class Closure
      {
        /// @brief Full geometry, p = 0 at its outlet: the reference.
        None,
        Tree,
        Single,
        Rest,
        Plateau
      };

      struct Config
      {
        std::string meshPath = "rmu_axi_full.mesh";
        std::string outputPrefix = "rmu_axi";

        /// @brief Inlet radius (m); the mesh is in units of it.
        Real scale = 1.0e-4;
        Real rho = 1050.0;
        HydraulicCurve::CarreauYasuda cy;

        /// @brief Shear-rate regularisation of the resolved viscosity (1/s).
        Real gammaRegularization = 1.0e-3;

        Closure closure = Closure::None;

        /// @brief Every tube of the full geometry (units of the inlet radius).
        std::vector<Tube> fullTubes;

        /// @brief The tubes the truncated geometry removes (downstream of x_T).
        std::vector<Tube> reducedTubes;

        /// @brief Axial position of the truncation (units of the inlet radius).
        Real truncation = 0.0;

        /// @brief Radius of the single compartment (units of the inlet
        ///        radius); <= 0 means the last reduced tube.
        Real singleRadius = 0.0;

        /// @brief Multiplier of the constant closures.
        Real constantScale = 1.0;

        /// @brief Mean apparent shear rate of the narrowest tube at the
        ///        a-priori resting operating point (1/s).
        Real gammaRef = 100.0;

        /// @brief p_in = Pbar [1 + A sin(omega t)].
        Real amplitude = 1.2;
        Real period = 0.85;
        int stepsPerCycle = 200;
        int cycles = 4;

        /// @brief Inlet pressure ramp (1 - cos)/2 over this many cycles.
        Real rampCycles = 1.0;

        /// @brief Picard iterations per step and their tolerance on the
        ///        relative velocity change and the outlet flux change
        ///        (relative to Q_0). The iteration starts from the linear
        ///        extrapolation 2 u^n - u^{n-1}.
        int picardIterations = 40;
        Real picardTolerance = 1.0e-6;

        /// @brief Newton tangent of the Carreau-Yasuda viscous term. The
        ///        step falls back to Picard if the residual increases.
        bool newton = true;

        Real backflowStabilization = 1.0;

        /// @brief u_r = 0 on the axis by the penalty gamma int u_r v_r, gamma =
        ///        axisPenalty mu_inf / R_1.
        Real axisPenalty = 1.0e8;
        Real pressurePenalty = 1.0e-12;

        /// @brief XDMF output period in steps (0 disables).
        int outputEvery = 0;

        Labels labels;
      };

      Benchmark(const Rodin::Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_curve(cfg.cy),
          m_mesh(makeMesh(context, cfg)),
          m_xdmf(context.getCommunicator(), cfg.outputPrefix + "/RmuAxi"),
          m_wallSet{cfg.labels.wall, cfg.labels.parentWall, cfg.labels.stubWall},
          m_uh(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_ph(std::integral_constant<size_t, 1>{}, m_mesh),
          m_rh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_radius(m_rh),
          m_u(m_uh),
          m_p(m_ph),
          m_v(m_uh),
          m_q(m_ph),
          m_uOld(m_uh),
          m_uNm1(m_uh),
          m_uIt(m_uh),
          m_du(m_uh),
          m_response(m_uh),
          m_s0(m_uh),
          m_s1(m_uh),
          m_wTrial(m_uh),
          m_wTest(m_uh),
          m_qFlux(m_ph),
          m_one(m_ph),
          m_flux(m_qFlux),
          m_flow(m_u, m_p, m_v, m_q),
          m_flowKSP(m_flow),
          m_projection(m_wTrial, m_wTest),
          m_projectionKSP(m_projection)
      {}

      int run();

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == RootRank)
        {
          Geometry::Mesh<Context::Local> mesh;
          mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
          if (mesh.getSpaceDimension() != 2 || mesh.getDimension() != 2)
          {
            throw std::runtime_error(
              "RmuAxi expects a triangular MEDIT mesh, Dimension 2.");
          }
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
        mesh.scale(cfg.scale);
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, 0);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        mesh.reconcile(1);
        return mesh;
      }

      bool isRoot() const
      {
        return m_mesh.getContext().getCommunicator().rank() == RootRank;
      }

      /// @brief sqrt(2 D:D) with the hoop strain, at p, for the field u.
      static Real shearRate(const VelocityGridFunctionType& u, const Point& p, Real reg)
      {
        const auto J = Jacobian(u).getValue(p);
        const Real exx = J(0, 0);
        const Real err = J(1, 1);
        const Real exr = 0.5 * (J(0, 1) + J(1, 0));
        const Real r = p.getPhysicalCoordinates()(1);
        const Real ett = (r > 1.0e-300) ? u.getValue(p)(1) / r : err;
        return std::sqrt(
          reg * reg + 2.0 * (exx * exx + err * err + ett * ett + 2.0 * exr * exr));
      }

      /// @brief Radius of the tube that contains the truncation (units of the
      ///        inlet radius): the outlet radius of the truncated geometry.
      Real outletRadius() const
      {
        Real x = 0.0;
        for (const auto& t : m_cfg.fullTubes)
        {
          x += t.length;
          if (m_cfg.truncation <= x * (1.0 + 1.0e-12))
            return t.radius;
        }
        return m_cfg.fullTubes.empty() ? 1.0 : m_cfg.fullTubes.back().radius;
      }

      /// @brief Viscous quantities of the iterate u^k at a quadrature point.
      struct ViscousPoint
      {
        Real shear = 0.0; ///< gd_k
        Real mu = 0.0; ///< mu(gd_k)
        Real c = 0.0; ///< mu'(gd_k) / gd_k
        Real hoop = 0.0; ///< u^k_r / r
      };

      /// @brief ViscousPoint of u^k at p, cached per quadrature point for the
      ///        current iterate: every coefficient of the viscous operator
      ///        reads the same velocity gradient, and the integrators of the
      ///        form visit the same points, so one evaluation serves all of
      ///        them. The cache is thread local and keyed on the revision of
      ///        u^k, which advances whenever u^k changes.
      const ViscousPoint& viscousAt(const Point& p) const
      {
        struct Key
        {
          std::uint64_t revision;
          size_t dimension;
          Index cell;
          Real r0;
          Real r1;

          bool operator==(const Key& o) const
          {
            return revision == o.revision && dimension == o.dimension &&
              cell == o.cell && r0 == o.r0 && r1 == o.r1;
          }
        };

        struct KeyHash
        {
          size_t operator()(const Key& k) const
          {
            size_t h = std::hash<Index>()(k.cell);
            h ^= std::hash<Real>()(k.r0) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<Real>()(k.r1) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<size_t>()(k.dimension) + (h << 6) + (h >> 2);
            return h;
          }
        };

        thread_local std::unordered_map<Key, ViscousPoint, KeyHash> cache;
        thread_local std::uint64_t cachedRevision = 0;
        thread_local const void* owner = nullptr;
        if (owner != this || cachedRevision != m_revision)
        {
          cache.clear();
          owner = this;
          cachedRevision = m_revision;
        }

        const auto& rc = p.getReferenceCoordinates();
        const Key key{m_revision, p.getPolytope().getDimension(),
          p.getPolytope().getIndex(), rc(0), rc.size() > 1 ? rc(1) : 0.0};
        const auto it = cache.find(key);
        if (it != cache.end())
          return it->second;

        ViscousPoint v;
        v.shear = shearRate(m_uIt, p, m_cfg.gammaRegularization);
        v.mu = m_curve.viscosity(v.shear);
        v.c = m_curve.derivative(v.shear) / v.shear;
        const Real y = std::max<Real>(p.getPhysicalCoordinates()(1), 1.0e-300);
        v.hoop = m_uIt.getValue(p)(1) / y;
        return cache.emplace(key, v).first->second;
      }

      Real tubeResistance(const Tube& t, Real Q) const
      {
        return m_curve.resistance(t.radius * m_cfg.scale, t.length * m_cfg.scale, Q);
      }

      /// @brief Resistance of a series of tubes at flow Q.
      Real seriesResistance(const std::vector<Tube>& tubes, Real Q) const
      {
        Real R = 0.0;
        for (const auto& t : tubes)
          R += tubeResistance(t, Q);
        return R;
      }

      /// @brief The reduced resistance of the configured closure at flow Q.
      Real closureResistance(Real Q) const
      {
        switch (m_cfg.closure)
        {
          case Closure::None:
            return 0.0;
          case Closure::Tree:
            return seriesResistance(m_cfg.reducedTubes, Q);
          case Closure::Single:
            return tubeResistance(Tube{m_singleRadius, m_singleLength}, Q);
          case Closure::Rest:
            return m_cfg.constantScale * m_restResistance;
          case Closure::Plateau:
            return m_cfg.constantScale * m_plateauResistance;
        }
        return 0.0;
      }

      /// @brief Outlet traction of the linearised step.
      /// @details By superposition the outlet flux is Q = q0 + q1 p_T, q0 the
      ///          flux with p_T = 0 and q1 < 0 the flux of a unit traction.
      ///          The closure p_T = R(Q) Q then gives the scalar equation
      ///          h(Q) = Q - q0 - q1 R(Q) Q = 0; R(Q) Q is increasing in Q (the
      ///          pressure-flow map is monotone), so h is increasing and its
      ///          root lies between 0 and q0. Bisection, then p_T = R(Q) Q.
      Real outletTraction(Real q0, Real q1) const
      {
        if (!(q1 < 0.0))
          return closureResistance(q0) * q0;
        const auto h = [&](Real Q) { return Q - q0 - q1 * closureResistance(Q) * Q; };
        Real lo = std::min<Real>(0.0, q0);
        Real hi = std::max<Real>(0.0, q0);
        constexpr int BisectionSteps = 200;
        for (int i = 0; i < BisectionSteps; ++i)
        {
          const Real mid = 0.5 * (lo + hi);
          if (h(mid) < 0.0)
            lo = mid;
          else
            hi = mid;
          if (hi - lo <= 1.0e-15 * std::max<Real>(std::abs(hi), 1e-300))
            break;
        }
        const Real Q = 0.5 * (lo + hi);
        return closureResistance(Q) * Q;
      }

      /// @brief Q = 2 pi int r u.n ds over a label (outward normal).
      Real flux(const VelocityGridFunctionType& u, Attribute tag)
      {
        const auto normal = BoundaryNormal(m_mesh);
        m_flux = BoundaryIntegral(m_radius * Dot(u, normal), m_qFlux).over(tag);
        m_flux.assemble();
        return 2.0 * Pi * m_flux(m_one);
      }

      /// @brief int r f ds / int r ds over a label.
      template <class F>
      Real boundaryMean(const F& f, Attribute tag)
      {
        m_flux = BoundaryIntegral(m_radius * f, m_qFlux).over(tag);
        m_flux.assemble();
        const Real num = m_flux(m_one);
        m_flux = BoundaryIntegral(m_radius, m_qFlux).over(tag);
        m_flux.assemble();
        return num / std::max<Real>(m_flux(m_one), 1e-300);
      }

      Real inletPressure(Real t) const
      {
        const Real T = m_cfg.period;
        const Real tr = m_cfg.rampCycles * T;
        const Real ramp =
          (tr > 0.0 && t < tr) ? 0.5 * (1.0 - std::cos(Pi * t / tr)) : 1.0;
        return m_pbar * ramp * (1.0 + m_cfg.amplitude * std::sin(2.0 * Pi * t / T));
      }

      void calibrate();
      void setupProblem();
      void diagnostics();
      void writeCSVHeader();
      void writeCSVRow(int cycle, int picard, Real residual);

      Config m_cfg;
      HydraulicCurve m_curve;
      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_wallSet;

      VelocityFESType m_uh;
      PressureFESType m_ph;
      PressureFESType m_rh;
      PressureGridFunctionType m_radius;

      VelocityTrialFunctionType m_u;
      PressureTrialFunctionType m_p;
      VelocityTestFunctionType m_v;
      PressureTestFunctionType m_q;

      /// @brief u^n, u^{n-1}, the Picard iterate u^k and its increment.
      VelocityGridFunctionType m_uOld;
      VelocityGridFunctionType m_uNm1;
      VelocityGridFunctionType m_uIt;
      VelocityGridFunctionType m_du;

      /// @brief Recovered rows of 2 D(u) (for the wall shear stress).
      VelocityGridFunctionType m_s0;
      VelocityGridFunctionType m_s1;
      VelocityTrialFunctionType m_wTrial;
      VelocityTestFunctionType m_wTest;

      PressureTestFunctionType m_qFlux;
      PressureGridFunctionType m_one;
      FluxFormType m_flux;

      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;
      VectorProblemType m_projection;
      Rodin::Solver::KSP m_projectionKSP;

      // Calibration (a-priori resting operating point).
      Real m_q0 = 0.0;
      Real m_pbar = 0.0;
      Real m_restResistance = 0.0;
      Real m_plateauResistance = 0.0;
      Real m_singleRadius = 0.0;
      Real m_singleLength = 0.0;

      // Step state.
      Real m_t = 0.0;
      Real m_dt = 0.0;
      Real m_pIn = 0.0;
      Real m_outletResistance = 0.0;

      /// @brief 1 for the Newton viscous tangent, 0 for Picard.
      Real m_newtonWeight = 1.0;

      /// @brief Revision of u^k, advanced whenever it changes (viscousAt).
      std::uint64_t m_revision = 1;

      /// @brief Outlet traction in the assembled form: 1 once, to extract the
      ///        unit-traction right-hand side, 0 otherwise.
      Real m_tractionAssembly = 0.0;

      /// @brief Right-hand side of a unit outlet traction, and the response
      ///        to it under the current linearisation.
      ::Vec m_unitTraction = PETSC_NULLPTR;
      ::Vec m_unitResponse = PETSC_NULLPTR;
      VelocityGridFunctionType m_response;

      /// @brief Accumulated wall-clock times (s).
      Real m_timeAssembly = 0.0;
      Real m_timeSolve = 0.0;
      Real m_timeDiagnostics = 0.0;
      Real m_outletTraction = 0.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_pInMean = 0.0;
      Real m_pOutMean = 0.0;
      Real m_dissipation = 0.0;
      Real m_wssParent = 0.0;
      Real m_wssStub = 0.0;

      std::ofstream m_csv;
  };

  void Benchmark::calibrate()
  {
    if (m_cfg.fullTubes.empty())
    {
      throw std::runtime_error(
        "-rmu_full_tubes is required (R:L,... units of the inlet radius).");
    }

    Real rMin = m_cfg.fullTubes.front().radius;
    for (const auto& t : m_cfg.fullTubes)
      rMin = std::min(rMin, t.radius);

    // A-priori resting operating point: gammaRef in the narrowest tube, and
    // the WRMS budget of the whole series at that flow (contraction losses
    // excluded; the full run includes them).
    const Real Rmin = rMin * m_cfg.scale;
    m_q0 = Pi * Rmin * Rmin * Rmin * m_cfg.gammaRef / 4.0;
    m_pbar = m_q0 * seriesResistance(m_cfg.fullTubes, m_q0);

    if (m_cfg.closure != Closure::None && m_cfg.reducedTubes.empty())
    {
      throw std::runtime_error(
        "-rmu_reduced_tubes is required by every closure but none.");
    }

    m_restResistance = seriesResistance(m_cfg.reducedTubes, m_q0);

    m_plateauResistance = 0.0;
    for (const auto& t : m_cfg.reducedTubes)
    {
      const Real R = t.radius * m_cfg.scale;
      m_plateauResistance +=
        8.0 * m_cfg.cy.muInf * t.length * m_cfg.scale / (Pi * std::pow(R, 4.0));
    }

    // Single compartment: radius R_s, and the effective length that gives the
    // tree resistance at the resting point (the calibration of the paper,
    // with the pressure budget of the removed tubes).
    if (!m_cfg.reducedTubes.empty())
    {
      m_singleRadius = (m_cfg.singleRadius > 0.0) ? m_cfg.singleRadius
                                                  : m_cfg.reducedTubes.back().radius;
      const Real Rs = m_singleRadius * m_cfg.scale;
      const Real muS = m_curve(4.0 * m_q0 / (Pi * Rs * Rs * Rs));
      m_singleLength =
        m_restResistance * Pi * std::pow(Rs, 4.0) / (8.0 * muS) / m_cfg.scale;
    }

    if (isRoot())
    {
      const Real omega = 2.0 * Pi / m_cfg.period;
      Alert::Info() << "[calib] Q0=" << m_q0 << " m^3/s  Pbar=" << m_pbar
                    << " Pa  gammaRef=" << m_cfg.gammaRef << " 1/s  A=" << m_cfg.amplitude
                    << "  T=" << m_cfg.period << " s" << Alert::Raise;
      for (const auto& t : m_cfg.fullTubes)
      {
        const Real R = t.radius * m_cfg.scale;
        const Real gA = 4.0 * m_q0 / (Pi * R * R * R);
        const Real wo = R * std::sqrt(omega * m_cfg.rho / m_curve(gA));
        const Real re = m_cfg.rho * (m_q0 / (Pi * R * R)) * 2.0 * R / m_curve(gA);
        Alert::Info() << "[calib] tube R=" << t.radius << " L=" << t.length
                      << "  gd_a(Q0)=" << gA << " 1/s  mu_ap=" << m_curve(gA)
                      << " Pa s  Wo=" << wo << "  Re=" << re
                      << "  share of Pbar=" << (m_q0 * tubeResistance(t, m_q0) / m_pbar)
                      << Alert::Raise;
      }
      if (m_cfg.closure != Closure::None)
      {
        Alert::Info() << "[calib] R_tree(Q0)=" << m_restResistance
                      << "  R_plateau=" << m_plateauResistance
                      << "  single: R_s=" << m_singleRadius << " L_s=" << m_singleLength
                      << "  constantScale=" << m_cfg.constantScale << Alert::Raise;
      }
    }
  }

  void Benchmark::setupProblem()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_dt;
    const auto& r = m_radius;
    const Real reg = m_cfg.gammaRegularization;

    // Order declared for the pointwise coefficients: they are not
    // polynomials, and the quadrature is chosen for the polynomial part of
    // each integrand (P2 x P2 and the weight r) plus this order. One keeps
    // the rule close to the Galerkin one; higher values multiply the
    // assembly cost without changing the converged solution appreciably.
    constexpr size_t CoefficientOrder = 1;

    const auto mu = PointwiseScalar(
      [this](const Point& p) -> Real { return viscousAt(p).mu; }, CoefficientOrder);
    const auto muOverR = PointwiseScalar(
      [this](const Point& p) -> Real {
        const Real y = std::max<Real>(p.getPhysicalCoordinates()(1), 1.0e-300);
        return viscousAt(p).mu / y;
      },
      CoefficientOrder);

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto ur = Component(m_u, 1);
    const auto vr = Component(m_v, 1);

    // Newton terms of the viscous operator. With DD_k = D(u^k) including the
    // hoop strain h_k = u^k_r/r, gd_k^2 = reg^2 + 2 DD_k:DD_k and
    // c_k = mu'(gd_k)/gd_k, the tangent of 2 mu(gd) DD(u):DD(v) adds
    //   4 int c_k (DD_k:DD(u)) (DD_k:DD(v)) r,
    // and the right-hand side removes its value at u^k,
    //   2 int c_k (gd_k^2 - reg^2) (DD_k:DD(v)) r,
    // so the two cancel at convergence and the converged state is that of the
    // Picard iteration. m_newtonWeight = 0 turns them off (Picard fallback).
    const auto newtonC = PointwiseScalar(
      [this](const Point& p) -> Real {
        return (m_newtonWeight == 0.0) ? 0.0 : m_newtonWeight * viscousAt(p).c;
      },
      CoefficientOrder);
    const auto newtonRhs = PointwiseScalar(
      [this, reg](const Point& p) -> Real {
        if (m_newtonWeight == 0.0)
          return 0.0;
        const auto& v = viscousAt(p);
        return m_newtonWeight * v.c * (v.shear * v.shear - reg * reg);
      },
      CoefficientOrder);
    const auto hoopK = PointwiseScalar(
      [this](const Point& p) -> Real { return viscousAt(p).hoop; }, CoefficientOrder);
    const auto hoopKOverR = PointwiseScalar(
      [this](const Point& p) -> Real {
        const Real y = std::max<Real>(p.getPhysicalCoordinates()(1), 1.0e-300);
        return viscousAt(p).hoop / y;
      },
      CoefficientOrder);
    const auto symUk = 0.5 * (Jacobian(m_uIt) + Transpose(Jacobian(m_uIt)));
    const auto sU = Dot(symUk, symU) + hoopKOverR * ur;

    const auto convU = Mult(Jacobian(m_u), m_uIt);
    const auto temam = (r * Div(m_uIt) + Dot(m_uIt, unit(1))) * Dot(m_u, m_v);
    const auto backflow =
      0.5 * rho * m_cfg.backflowStabilization * Max(-Dot(m_uIt, normal), 0.0);
    const Real axisGamma = m_cfg.axisPenalty * m_cfg.cy.muInf / m_cfg.scale;

    RealFunction pInFn = [this](const Point&) { return m_pIn; };

    // Outlet: the uniform traction p_T. It is assembled only once, with
    // p_T = 1, to obtain the right-hand side of a unit outlet traction; every
    // step is then solved by superposition (see run()), so the assembled value
    // is zero.
    RealFunction tractionFn = [this](const Point&) { return m_tractionAssembly; };

    const auto& L = m_cfg.labels;

    if (m_cfg.closure == Closure::None)
    {
      m_flow = (rho / dt) * Integral(r * m_u, m_v) -
        (rho / dt) * Integral(r * m_uOld, m_v) + rho * Integral(r * convU, m_v) +
        0.5 * rho * Integral(temam) + 2.0 * Integral(r * (mu * symU), symV) +
        2.0 * Integral(muOverR * ur, vr) +
        4.0 * Integral(r * (newtonC * sU) * symUk, symV) +
        4.0 * Integral((newtonC * hoopK) * sU, vr) -
        2.0 * Integral(r * newtonRhs * symUk, symV) -
        2.0 * Integral(newtonRhs * hoopK, vr) - Integral(r * m_p, Div(m_v)) -
        Integral(m_p, vr) + Integral(r * Div(m_u), m_q) + Integral(ur, m_q) +
        m_cfg.pressurePenalty * Integral(r * m_p, m_q) +
        BoundaryIntegral(r * pInFn * Dot(m_v, normal)).over(L.inlet) +
        BoundaryIntegral(r * backflow * Dot(m_u, m_v)).over(L.inlet) +
        BoundaryIntegral(r * backflow * Dot(m_u, m_v)).over(L.outlet) +
        BoundaryIntegral(RealFunction(axisGamma) * Dot(Mult(radialProjector(), m_u), m_v))
          .over(L.axis) +
        DirichletBC(m_u, Zero(dim)).on(m_wallSet);
    }
    else
    {
      m_flow = (rho / dt) * Integral(r * m_u, m_v) -
        (rho / dt) * Integral(r * m_uOld, m_v) + rho * Integral(r * convU, m_v) +
        0.5 * rho * Integral(temam) + 2.0 * Integral(r * (mu * symU), symV) +
        2.0 * Integral(muOverR * ur, vr) +
        4.0 * Integral(r * (newtonC * sU) * symUk, symV) +
        4.0 * Integral((newtonC * hoopK) * sU, vr) -
        2.0 * Integral(r * newtonRhs * symUk, symV) -
        2.0 * Integral(newtonRhs * hoopK, vr) - Integral(r * m_p, Div(m_v)) -
        Integral(m_p, vr) + Integral(r * Div(m_u), m_q) + Integral(ur, m_q) +
        m_cfg.pressurePenalty * Integral(r * m_p, m_q) +
        BoundaryIntegral(r * pInFn * Dot(m_v, normal)).over(L.inlet) +
        BoundaryIntegral(r * backflow * Dot(m_u, m_v)).over(L.inlet) +
        BoundaryIntegral(r * backflow * Dot(m_u, m_v)).over(L.outlet) +
        BoundaryIntegral(r * tractionFn * Dot(m_v, normal)).over(L.outlet) +
        BoundaryIntegral(RealFunction(axisGamma) * Dot(Mult(radialProjector(), m_u), m_v))
          .over(L.axis) +
        DirichletBC(m_u, Zero(dim)).on(m_wallSet);
    }
  }

  void Benchmark::diagnostics()
  {
    const auto& uSol = m_u.getSolution();
    const auto& pSol = m_p.getSolution();
    const Real reg = m_cfg.gammaRegularization;
    const Real xT = m_cfg.truncation * m_cfg.scale;
    const Real tol = 1.0e-9 * m_cfg.scale;

    m_qIn = -flux(uSol, m_cfg.labels.inlet);
    m_qOut = flux(uSol, m_cfg.labels.outlet);
    m_pInMean = boundaryMean(pSol, m_cfg.labels.inlet);
    m_pOutMean = boundaryMean(pSol, m_cfg.labels.outlet);

    // Viscous dissipation 2 pi int r 2 mu D:D over the region common to both
    // geometries, x <= x_T (a mesh column lies at x_T).
    constexpr size_t DissipationOrder = 2;
    const auto dissipation = PointwiseScalar(
      [&uSol, reg, xT, tol, this](const Point& p) -> Real {
        if (p.getPhysicalCoordinates()(0) > xT + tol)
          return 0.0;
        const Real g = shearRate(uSol, p, reg);
        return m_curve.viscosity(g) * g * g;
      },
      DissipationOrder);
    m_flux = Integral(m_radius * dissipation, m_qFlux);
    m_flux.assemble();
    m_dissipation = 2.0 * Pi * m_flux(m_one);

    // Wall shear stress on the parent and stub cylinders (n = e_r): the
    // rows of 2 D(u) are recovered by L2 projection and read on the wall;
    // tau_w = -mu (2 D)_xr, positive for forward flow.
    const auto jac = Jacobian(uSol);
    const auto offDiagonal = Component(jac, 0, 1) + Component(jac, 1, 0);
    m_projection = Integral(m_wTrial, m_wTest) -
      Integral(VectorFunction(2.0 * Component(jac, 0, 0), offDiagonal), m_wTest);
    m_projection.solve(m_projectionKSP);
    m_s0.setData(m_wTrial.getSolution().getData());
    m_projection = Integral(m_wTrial, m_wTest) -
      Integral(VectorFunction(offDiagonal, 2.0 * Component(jac, 1, 1)), m_wTest);
    m_projection.solve(m_projectionKSP);
    m_s1.setData(m_wTrial.getSolution().getData());

    const auto wss = PointwiseScalar(
      [this](const Point& p) -> Real {
        const auto s0 = m_s0.getValue(p);
        const auto s1 = m_s1.getValue(p);
        const Real g = std::sqrt(
          0.5 * (s0(0) * s0(0) + s0(1) * s0(1) + s1(0) * s1(0) + s1(1) * s1(1)));
        return -m_curve.viscosity(g) * s0(1);
      },
      DissipationOrder);
    m_wssParent = boundaryMean(wss, m_cfg.labels.parentWall);
    m_wssStub = boundaryMean(wss, m_cfg.labels.stubWall);
  }

  void Benchmark::writeCSVHeader()
  {
    m_csv << "t,cycle,pIn,qIn,qOut,pInMean,pOutMean,pT,R,picardIts,picardResidual,"
          << "dissipation,wssParent,wssStub\n";
  }

  void Benchmark::writeCSVRow(int cycle, int picard, Real residual)
  {
    if (!isRoot())
      return;
    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << m_qIn << ',' << m_qOut << ','
          << m_pInMean << ',' << m_pOutMean << ',' << m_outletTraction << ','
          << m_outletResistance << ',' << picard << ',' << residual << ','
          << m_dissipation << ',' << m_wssParent << ',' << m_wssStub << '\n';
    m_csv.flush();
  }

  int Benchmark::run()
  {
    m_dt = m_cfg.period / static_cast<Real>(m_cfg.stepsPerCycle);
    calibrate();

    m_radius.project(RealFunction([](const Point& p) -> Real {
      return std::max<Real>(p.getPhysicalCoordinates()(1), 0.0);
    }));
    m_one = Real(1);
    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0}};
    m_uNm1 = Math::SpatialVector<Real>{{0.0, 0.0}};
    m_uIt = Math::SpatialVector<Real>{{0.0, 0.0}};
    m_du = Math::SpatialVector<Real>{{0.0, 0.0}};
    configureMassSolver(m_projectionKSP, "rmu_proj_");

    m_u.setName("velocity");
    m_p.setName("pressure");
    if (m_cfg.outputEvery > 0)
    {
      m_xdmf.setMesh(m_mesh);
      m_xdmf.add("velocity", m_u.getSolution());
      m_xdmf.add("pressure", m_p.getSolution());
    }

    setupProblem();

    // Right-hand side of a unit outlet traction: the difference of two
    // assemblies, so that its sign and its Dirichlet rows are those of the
    // assembled system.
    if (m_cfg.closure != Closure::None)
    {
      auto& axb = m_flow.getLinearSystem();
      m_tractionAssembly = 1.0;
      m_flow.assemble();
      PetscErrorCode ierr = VecDuplicate(axb.getVector(), &m_unitTraction);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecCopy(axb.getVector(), m_unitTraction);
      assert(ierr == PETSC_SUCCESS);
      m_tractionAssembly = 0.0;
      m_flow.assemble();
      ierr = VecAXPY(m_unitTraction, -1.0, axb.getVector());
      assert(ierr == PETSC_SUCCESS);
      ierr = VecDuplicate(axb.getSolution(), &m_unitResponse);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }

    if (isRoot())
    {
      m_csv.open(m_cfg.outputPrefix + "/RmuAxi.csv");
      if (!m_csv)
        throw std::runtime_error("Cannot open " + m_cfg.outputPrefix + "/RmuAxi.csv");
      m_csv.precision(12);
      writeCSVHeader();
    }

    const int total = m_cfg.cycles * m_cfg.stepsPerCycle;
    bool fieldSplits = false;

    for (int step = 0; step < total; ++step)
    {
      m_t += m_dt;
      const int cycle = step / m_cfg.stepsPerCycle;
      m_pIn = inletPressure(m_t);

      // Picard on advection, viscosity and outlet resistance, from the
      // extrapolation 2 u^n - u^{n-1} (u^0 = u^{-1} = 0 at rest).
      {
        ::Vec& z = m_uIt.getData();
        PetscErrorCode ierr = VecCopy(m_uOld.getData(), z);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecScale(z, 2.0);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAXPY(z, -1.0, m_uNm1.getData());
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateBegin(z, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecGhostUpdateEnd(z, INSERT_VALUES, SCATTER_FORWARD);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        ++m_revision;
      }
      Real qk = (m_cfg.closure == Closure::None) ? 0.0 : flux(m_uIt, m_cfg.labels.outlet);
      int k = 0;
      Real residual = 0.0;
      Real previous = std::numeric_limits<Real>::infinity();
      m_newtonWeight = m_cfg.newton ? 1.0 : 0.0;
      for (; k < std::max(1, m_cfg.picardIterations); ++k)
      {
        const auto t0 = std::chrono::steady_clock::now();
        m_flow.assemble();
        const auto t1 = std::chrono::steady_clock::now();
        m_timeAssembly += std::chrono::duration<Real>(t1 - t0).count();
        if (!fieldSplits)
        {
          m_flow.setFieldSplits();
          fieldSplits = true;
        }
        m_flow.solve(m_flowKSP);
        m_timeSolve +=
          std::chrono::duration<Real>(std::chrono::steady_clock::now() - t1).count();

        ::KSPConvergedReason reason = KSP_CONVERGED_ITS;
        PetscErrorCode ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
        assert(ierr == PETSC_SUCCESS);
        if (reason < 0)
        {
          Alert::Exception() << "[flow] linear solve failed at step " << (step + 1)
                             << " (KSP reason " << static_cast<int>(reason) << ")."
                             << Alert::Raise;
        }

        // Outlet by superposition: x = x0 + p_T x1, x1 the response to the
        // unit traction under the same factorised operator. The closure is
        // thus imposed exactly in the linearised step, whatever R(Q).
        if (m_cfg.closure != Closure::None)
        {
          auto& axb = m_flow.getLinearSystem();
          const auto ts = std::chrono::steady_clock::now();
          ierr = KSPSolve(m_flowKSP.getHandle(), m_unitTraction, m_unitResponse);
          assert(ierr == PETSC_SUCCESS);
          m_timeSolve +=
            std::chrono::duration<Real>(std::chrono::steady_clock::now() - ts).count();
          m_response.setData(m_unitResponse, 0);
          const Real q0 = flux(m_u.getSolution(), m_cfg.labels.outlet);
          const Real q1 = flux(m_response, m_cfg.labels.outlet);
          m_outletTraction = outletTraction(q0, q1);
          ierr = VecAXPY(axb.getSolution(), m_outletTraction, m_unitResponse);
          assert(ierr == PETSC_SUCCESS);
          m_u.getSolution().setData(axb.getSolution(), 0);
          m_p.getSolution().setData(axb.getSolution(), m_uh.getSize());
        }

        // Convergence: velocity change relative to the velocity, and outlet
        // flux change relative to the resting flow.
        ierr =
          VecWAXPY(m_du.getData(), -1.0, m_uIt.getData(), m_u.getSolution().getData());
        assert(ierr == PETSC_SUCCESS);
        PetscReal duNorm = 0.0;
        PetscReal uNorm = 0.0;
        ierr = VecNorm(m_du.getData(), NORM_INFINITY, &duNorm);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecNorm(m_u.getSolution().getData(), NORM_INFINITY, &uNorm);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;

        m_uIt.setData(m_u.getSolution().getData());
        ++m_revision;
        const Real qNext =
          (m_cfg.closure == Closure::None) ? 0.0 : flux(m_uIt, m_cfg.labels.outlet);
        const Real dq = std::abs(qNext - qk) / m_q0;
        qk = qNext;

        residual = std::max<Real>(duNorm / std::max<Real>(uNorm, 1e-300), dq);
        if (residual < m_cfg.picardTolerance)
        {
          ++k;
          break;
        }
        if (k >= 1 && residual > previous && m_newtonWeight > 0.0)
        {
          m_newtonWeight = 0.0;
          if (isRoot())
          {
            Alert::Warning() << "[step " << (step + 1) << "] Newton residual increased ("
                             << previous << " -> " << residual
                             << "): Picard for this step." << Alert::Raise;
          }
        }
        previous = residual;
      }

      // The resistance of the converged state.
      m_outletResistance = closureResistance(qk);

      m_uNm1.setData(m_uOld.getData());
      m_uOld.setData(m_uIt.getData());
      const auto td = std::chrono::steady_clock::now();
      diagnostics();
      m_timeDiagnostics +=
        std::chrono::duration<Real>(std::chrono::steady_clock::now() - td).count();
      writeCSVRow(cycle, k, residual);

      if (m_cfg.outputEvery > 0 && (step + 1) % m_cfg.outputEvery == 0)
        m_xdmf.write(m_t).flush();

      if (isRoot() && ((step + 1) % 20 == 0 || step < 3))
      {
        Alert::Info() << "[step " << (step + 1) << "/" << total << "] t=" << m_t
                      << "  pIn=" << m_pIn << "  qIn=" << m_qIn << "  qOut=" << m_qOut
                      << "  pT=" << m_outletTraction << "  picard=" << k << " (res "
                      << residual << ")  time: assembly " << m_timeAssembly
                      << " s, solve " << m_timeSolve << " s, diagnostics "
                      << m_timeDiagnostics << " s" << Alert::Raise;
      }
    }

    if (m_cfg.outputEvery > 0)
      m_xdmf.close();
    if (isRoot())
      m_csv.close();
    if (m_unitTraction)
      VecDestroy(&m_unitTraction);
    if (m_unitResponse)
      VecDestroy(&m_unitResponse);
    return 0;
  }
} // namespace Rodin::Examples::Heart::RmuAxi

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

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;
    {
      using Benchmark = Rodin::Examples::Heart::RmuAxi::Benchmark;
      Benchmark::Config cfg;

      const auto getString = [](const char* name, std::string& out) {
        char buffer[4096] = {};
        PetscBool set = PETSC_FALSE;
        PetscOptionsGetString(
          PETSC_NULLPTR, PETSC_NULLPTR, name, buffer, sizeof(buffer), &set);
        if (set)
          out = buffer;
        return set == PETSC_TRUE;
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

      getString("-rmu_mesh", cfg.meshPath);
      getString("-rmu_output_prefix", cfg.outputPrefix);
      getReal("-rmu_scale", cfg.scale);
      getReal("-rmu_rho", cfg.rho);
      getReal("-rmu_cy_mu0", cfg.cy.mu0);
      getReal("-rmu_cy_muinf", cfg.cy.muInf);
      getReal("-rmu_cy_lambda", cfg.cy.lambda);
      getReal("-rmu_cy_n", cfg.cy.n);
      getReal("-rmu_cy_a", cfg.cy.a);
      getReal("-rmu_gamma_regularization", cfg.gammaRegularization);

      std::string closure;
      if (getString("-rmu_closure", closure))
      {
        if (closure == "none")
          cfg.closure = Benchmark::Closure::None;
        else if (closure == "tree")
          cfg.closure = Benchmark::Closure::Tree;
        else if (closure == "single")
          cfg.closure = Benchmark::Closure::Single;
        else if (closure == "rest")
          cfg.closure = Benchmark::Closure::Rest;
        else if (closure == "plateau")
          cfg.closure = Benchmark::Closure::Plateau;
        else
          throw std::runtime_error("-rmu_closure: none, tree, single, rest or plateau.");
      }

      std::string tubes;
      if (getString("-rmu_full_tubes", tubes))
        cfg.fullTubes = Rodin::Examples::Heart::RmuAxi::parseTubes(tubes);
      if (getString("-rmu_reduced_tubes", tubes))
        cfg.reducedTubes = Rodin::Examples::Heart::RmuAxi::parseTubes(tubes);
      getReal("-rmu_truncation", cfg.truncation);
      getReal("-rmu_single_radius", cfg.singleRadius);
      getReal("-rmu_constant_scale", cfg.constantScale);
      getReal("-rmu_gamma_ref", cfg.gammaRef);
      getReal("-rmu_amplitude", cfg.amplitude);
      getReal("-rmu_period", cfg.period);
      getInt("-rmu_steps_per_cycle", cfg.stepsPerCycle);
      getInt("-rmu_cycles", cfg.cycles);
      getReal("-rmu_ramp_cycles", cfg.rampCycles);
      getInt("-rmu_picard_its", cfg.picardIterations);
      {
        PetscBool set = PETSC_FALSE;
        PetscBool value = PETSC_TRUE;
        PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-rmu_newton", &value, &set);
        if (set)
          cfg.newton = (value == PETSC_TRUE);
      }
      getReal("-rmu_picard_tol", cfg.picardTolerance);
      getReal("-rmu_backflow", cfg.backflowStabilization);
      getReal("-rmu_axis_penalty", cfg.axisPenalty);
      getInt("-rmu_output_every", cfg.outputEvery);

      Benchmark benchmark(context, cfg);
      status = benchmark.run();
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
