/*
 * Uses Rodin's RmuAxi pressure-driven Taylor–Hood formulation and WRMS table.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

namespace
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;
  constexpr Real Pi = std::numbers::pi_v<Real>;
  // Coordinates and velocity are scaled so the mixed system is well conditioned.
  constexpr Real Length = 1e-4;
  constexpr Real Speed = 1e-2;
  constexpr Real MuInf = 0.0041;
  constexpr Real Density = 1050;
  constexpr Real Period = 0.85;
  constexpr Real Pressure = MuInf * Speed / Length;
  constexpr Real Shear = Speed / Length;
  constexpr Real StrouhalCoefficient = Density * Length * Length / MuInf;
  constexpr Real Reynolds = Density * Speed * Length / MuInf;
  constexpr Real Regularization = 1e-3;
  constexpr Real PicardTolerance = 2e-7;
  constexpr int MaxIterations = 80;
  constexpr Real WaveAmplitude = 0.45;

  constexpr Attribute Inlet = 1;
  constexpr Attribute Wall = 2;
  constexpr std::array<Attribute, 2> Cut{3, 4};
  constexpr std::array<Attribute, 2> Distal{5, 6};
  constexpr std::array<Real, 2> Radius{0.72, 0.46};

#include "HydraulicCurve.inc"

  class SlitCurve
  {
    public:
      explicit SlitCurve(const HydraulicCurve& blood)
      {
        // WRMS for a slit: gamma_a = 3/tau_w^2 times the integral of tau gamma(tau).
        constexpr int Nodes = 16385;
        Real previousStress = 0;
        Real previousIntegrand = 0;
        Real integral = 0;
        for (int node = 0; node < Nodes; ++node)
        {
          const Real gamma = std::pow(10.0, -10.0 + 18.0 * node / (Nodes - 1));
          const Real stress = blood.viscosity(gamma) * gamma;
          const Real integrand = stress * gamma;
          integral += node == 0
            ? stress * integrand / 3
            : (stress - previousStress) * (previousIntegrand + integrand) / 2;
          const Real apparent = 3 * integral / (stress * stress);
          m_rate.push_back(std::log(apparent));
          m_viscosity.push_back(std::log(stress / apparent));
          previousStress = stress;
          previousIntegrand = integrand;
        }
      }

      Real operator()(Real apparent) const
      {
        const Real rate = std::log(std::max(apparent, 1e-300));
        const auto last = std::upper_bound(m_rate.begin(), m_rate.end(), rate);
        if (last == m_rate.begin())
          return std::exp(m_viscosity.front());
        if (last == m_rate.end())
          return std::exp(m_viscosity.back());
        const auto index = static_cast<size_t>(last - m_rate.begin());
        const Real weight =
          (rate - m_rate[index - 1]) / (m_rate[index] - m_rate[index - 1]);
        return std::exp(std::lerp(m_viscosity[index - 1], m_viscosity[index], weight));
      }

    private:
      std::vector<Real> m_rate, m_viscosity;
  };

  template <class F>
  class QuadratureCoefficient final : public RealFunctionBase<QuadratureCoefficient<F>>
  {
    public:
      explicit QuadratureCoefficient(F function)
        : m_function(std::move(function))
      {}
      Real getValue(const Point& p) const
      {
        return m_function(p);
      }
      Optional<size_t> getOrder(const Polytope&) const noexcept
      {
        return 1;
      }
      QuadratureCoefficient* copy() const noexcept override
      {
        return new QuadratureCoefficient(*this);
      }

    private:
      F m_function;
  };

  void check(PetscErrorCode value)
  {
    if (value != PETSC_SUCCESS)
      throw std::runtime_error("PETSc error " + std::to_string(value));
  }

  using MeshType = Mesh<Context::MPI>;
  MeshType loadMesh(const Context::MPI& context, const std::string& path)
  {
    Rodin::MPI::Sharder sharder(context);
    Mesh<Context::Local> mesh;
    mesh.load(path, IO::FileFormat::MEDIT);
    const auto dim = mesh.getDimension();
    mesh.getConnectivity().compute(dim, dim);
    mesh.getConnectivity().compute(dim, 0);
    mesh.getConnectivity().compute(dim, dim - 1);
    mesh.getConnectivity().compute(dim - 1, dim);
    mesh.getConnectivity().compute(dim - 1, 0);
    if (dim == 3)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
    }
    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(1);
    sharder.shard(partitioner);
    sharder.scatter(0);
    auto result = sharder.gather(0);
    auto& connectivity = result.getConnectivity();
    connectivity.compute(dim, dim);
    connectivity.compute(dim, 0);
    connectivity.compute(dim, dim - 1);
    connectivity.compute(dim - 1, dim);
    connectivity.compute(dim - 1, 0);
    if (dim == 3)
    {
      connectivity.compute(2, 1);
      connectivity.compute(1, 0);
      connectivity.compute(3, 1);
      result.reconcile(2);
    }
    result.reconcile(1);
    return result;
  }

  using VelocitySpace = H1<2, Math::SpatialVector<Real>, MeshType>;
  using PressureSpace = H1<1, Real, MeshType>;
  template <class Space>
  using Field = PETSc::Variational::GridFunction<Space>;
  template <class Space>
  using Trial = PETSc::Variational::TrialFunction<Field<Space>, Space>;
  template <class Space>
  using Test = PETSc::Variational::TestFunction<Space>;
  using FlowProblem = Problem<PETSc::Math::LinearSystem, Trial<VelocitySpace>,
    Trial<PressureSpace>, Test<VelocitySpace>, Test<PressureSpace>>;

  enum class Closure
  {
    Full,
    Constant,
    Variable
  };

  class Benchmark
  {
    public:
      Benchmark(
        const Context::MPI& context, const std::string& path, Closure closure, Real dt)
        : m_mesh(loadMesh(context, path)),
          m_uh(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_ph(std::integral_constant<size_t, 1>{}, m_mesh),
          m_u(m_uh),
          m_p(m_ph),
          m_v(m_uh),
          m_q(m_ph),
          m_old(m_uh),
          m_iterate(m_uh),
          m_difference(m_uh),
          m_one(m_ph),
          m_section(m_q),
          m_flow(m_u, m_p, m_v, m_q),
          m_ksp(m_flow),
          m_curve(
            HydraulicCurve::CarreauYasuda{0.0712616, MuInf, 9.2780, 0.3572, 1.8773}),
          m_slitCurve(m_curve),
          m_closure(closure),
          m_dt(dt),
          m_caps(closure == Closure::Full ? Distal : Cut)
      {
        setup();
      }

      ~Benchmark()
      {
        for (auto* vector :
          {&m_inletLoad, &m_loads[0], &m_loads[1], &m_response[0], &m_response[1]})
        {
          if (*vector)
            VecDestroy(vector);
        }
      }

      void run(const std::string& path, Real head, int steps, int cycles,
        const std::string& calibration)
      {
        PetscInt snapshotStep = 0;
        check(PetscOptionsGetInt(
          nullptr, nullptr, "-y_snapshot_step", &snapshotStep, nullptr));
        if (snapshotStep < 0 || snapshotStep > steps * cycles ||
          (snapshotStep > 0 && m_mesh.getDimension() != 3))
        {
          throw std::runtime_error(
            "A velocity snapshot requires a valid step on a 3D mesh");
        }
        if (m_closure != Closure::Full)
        {
          std::ifstream input(calibration);
          if (!(input >> m_referenceFlow[0] >> m_referenceFlow[1] >> m_referenceDrop[0] >>
                m_referenceDrop[1]) ||
            m_referenceFlow.minCoeff() <= 0 || m_referenceDrop.minCoeff() <= 0)
          {
            throw std::runtime_error(
              "A positive steady calibration is required from " + calibration);
          }
        }
        if (steps == 0)
        {
          if (m_closure != Closure::Full)
            throw std::runtime_error("The full geometry supplies the steady calibration");
          advance(head);
          const auto flows = outletFlux(m_flow.getLinearSystem().getSolution());
          std::array<Real, 2> drops{sectionPressure(Cut[0]), sectionPressure(Cut[1])};
          std::ofstream out(calibration);
          out.precision(16);
          out << flows.transpose() << ' ' << drops[0] << ' ' << drops[1] << '\n';
          if (!out || flows.minCoeff() <= 0 ||
            *std::min_element(drops.begin(), drops.end()) <= 0)
            throw std::runtime_error("Invalid calibration at the two cuts");
          std::cout << "Steady calibration: q = " << flows.transpose()
                    << ", cut pressure = " << drops[0] << ' ' << drops[1]
                    << " (scaled units)\n";
          return;
        }

        std::ofstream out(path);
        if (!out)
          throw std::runtime_error("Cannot write " + path);
        out.precision(13);
        const auto dimensions = m_mesh.getDimension();
        const std::string units = dimensions == 3 ? "m3_s" : "m2_s";
        const Real flowScale = Speed * std::pow(Length, dimensions - 1);
        out << "step,time_s,inlet_pressure_pa,inlet_" << units << ",wide_" << units
            << ",narrow_" << units << ','
            << "wide_cut_pressure_pa,narrow_cut_pressure_pa,iterations,residual\n";
        for (int step = 1; step <= steps * cycles; ++step)
        {
          const Real time = step * m_dt;
          const Real forcing =
            head * (1 + WaveAmplitude * std::sin(2 * Pi * time / Period));
          const auto [iterations, residual] = advance(forcing);
          if (step == snapshotStep)
            writeVelocity(path + ".vtk", time, forcing * Pressure);
          const auto solution = m_flow.getLinearSystem().getSolution();
          const auto fluxes = outletFlux(solution);
          Real inlet = 0;
          check(VecDot(m_inletLoad, solution, &inlet));
          out << step << ',' << time << ',' << forcing * Pressure << ','
              << inlet * flowScale << ',' << fluxes[0] * flowScale << ','
              << fluxes[1] * flowScale << ',' << sectionPressure(Cut[0]) * Pressure << ','
              << sectionPressure(Cut[1]) * Pressure << ',' << iterations << ','
              << residual << '\n';
          out.flush();
          if (step == 1 || step % 10 == 0)
            std::cout << path << " step " << step << '/' << steps * cycles << std::endl;
        }
      }

    private:
      void writeVelocity(const std::string& path, Real time, Real inletPressure)
      {
        // VTK quadratic tetrahedron order: four vertices then six edge midpoints.
        const Polytope::Traits traits(Polytope::Type::Tetrahedron);
        std::array<Math::SpatialPoint, 10> samples;
        for (int i = 0; i < 4; ++i)
          samples[i] = traits.getVertex(i);
        constexpr int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
        for (int i = 0; i < 6; ++i)
          samples[4 + i] = 0.5 * (samples[edges[i][0]] + samples[edges[i][1]]);

        const size_t cells = m_mesh.getCellCount();
        std::ofstream out(path);
        out.precision(16);
        out << "# vtk DataFile Version 3.0\nRefined Y velocity; coordinates in m, "
               "velocity in m/s\n"
            << "ASCII\nDATASET UNSTRUCTURED_GRID\nFIELD FieldData 2\n"
            << "time_s 1 1 double\n"
            << time << "\ninlet_pressure_pa 1 1 double\n"
            << inletPressure << "\nPOINTS " << 10 * cells << " double\n";
        std::vector<Eigen::Vector3d> velocities;
        velocities.reserve(10 * cells);
        for (auto cell = m_mesh.getCell(); cell; ++cell)
        {
          for (const auto& rc : samples)
          {
            const Point point(*cell, rc);
            out << (Length * point.getPhysicalCoordinates()).transpose() << '\n';
            const auto velocity = m_u.getSolution().getValue(point);
            velocities.emplace_back(
              Speed * velocity(0), Speed * velocity(1), Speed * velocity(2));
          }
        }
        out << "CELLS " << cells << ' ' << 11 * cells << '\n';
        for (size_t i = 0; i < cells; ++i)
        {
          out << 10;
          for (int j = 0; j < 10; ++j)
            out << ' ' << 10 * i + j;
          out << '\n';
        }
        out << "CELL_TYPES " << cells << '\n';
        for (size_t i = 0; i < cells; ++i)
          out << "24\n";
        out << "POINT_DATA " << velocities.size() << "\nVECTORS velocity_m_s double\n";
        for (const auto& v : velocities)
          out << v.transpose() << '\n';
        if (!out)
        {
          throw std::runtime_error(
            "Could not write the quadratic velocity field to " + path);
        }
      }

      void setup()
      {
        const auto dim = m_mesh.getSpaceDimension();
        const auto normal = BoundaryNormal(m_mesh);
        const auto velocityGradient = Jacobian(m_iterate);
        const auto mu =
          QuadratureCoefficient([&, velocityGradient](const Point& point) -> Real {
            const auto jac = velocityGradient.getValue(point);
            const Real gamma = std::hypot(Regularization,
              Shear * (jac + jac.transpose()).norm() / std::sqrt(Real(2)));
            return m_curve.viscosity(gamma) / MuInf;
          });
        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto transport = Mult(Jacobian(m_u), m_iterate);
        const Real inertia = m_dt > 0 ? StrouhalCoefficient / m_dt : 0;
        const auto in = RealFunction([this](const Point&) { return m_applied[0]; });
        const auto a = RealFunction([this](const Point&) { return m_applied[1]; });
        const auto b = RealFunction([this](const Point&) { return m_applied[2]; });
        const FlatSet<Attribute> open{Inlet, m_caps[0], m_caps[1]};
        const auto backflow = (0.5 * Reynolds) * Max(-Dot(m_iterate, normal), 0.0);
        const auto zero = dim == 3 ? Math::SpatialVector<Real>{{0, 0, 0}}
                                   : Math::SpatialVector<Real>{{0, 0}};
        m_old = zero;
        m_iterate = zero;
        m_one = 1.0;
        m_flow = inertia * Integral(m_u, m_v) - inertia * Integral(m_old, m_v) +
          Reynolds * Integral(transport, m_v) +
          (0.5 * Reynolds) * Integral(Div(m_iterate) * Dot(m_u, m_v)) +
          2.0 * Integral(mu * symU, symV) - Integral(m_p, Div(m_v)) +
          Integral(Div(m_u), m_q) +
          BoundaryIntegral(backflow * Dot(m_u, m_v)).over(open) +
          BoundaryIntegral(in * Dot(m_v, normal)).over(Inlet) +
          BoundaryIntegral(a * Dot(m_v, normal)).over(m_caps[0]) +
          BoundaryIntegral(b * Dot(m_v, normal)).over(m_caps[1]) +
          DirichletBC(m_u, Zero(dim)).on(Wall);

        std::array<::Vec*, 3> loads{&m_inletLoad, &m_loads[0], &m_loads[1]};
        for (size_t j = 0; j < loads.size(); ++j)
        {
          m_applied[j] = 1;
          m_flow.assemble();
          auto& system = m_flow.getLinearSystem();
          check(VecDuplicate(system.getVector(), loads[j]));
          check(VecCopy(system.getVector(), *loads[j]));
          m_applied[j] = 0;
        }
        for (auto& response : m_response)
          check(VecDuplicate(m_flow.getLinearSystem().getSolution(), &response));
      }

      Real sectionPressure(Attribute section)
      {
        if (m_closure == Closure::Full)
        {
          m_section = FaceIntegral(m_p.getSolution(), m_q).over(section);
          m_section.assemble();
          const Real integrated = m_section(m_one);
          m_section = FaceIntegral(RealFunction(1), m_q).over(section);
          m_section.assemble();
          const Real area = m_section(m_one);
          if (!(area > 0))
            throw std::runtime_error("Missing interior cut section");
          return integrated / area;
        }
        return m_outletPressure[section == Cut[0] ? 0 : 1];
      }

      Eigen::Vector2d outletFlux(::Vec state) const
      {
        Eigen::Vector2d result;
        for (int j = 0; j < 2; ++j)
        {
          check(VecDot(m_loads[j], state, &result[j]));
          result[j] *= -1;
        }
        return result;
      }

      Real drop(int branch, Real flow) const
      {
        const Real base = m_referenceDrop[branch] * flow / m_referenceFlow[branch];
        if (m_closure == Closure::Constant)
          return base;
        const Real radius = Radius[branch];
        if (m_mesh.getDimension() == 3)
        {
          const Real ratePerFlow = 4 * Shear / (Pi * radius * radius * radius);
          return base * m_curve(ratePerFlow * std::abs(flow)) /
            m_curve(ratePerFlow * m_referenceFlow[branch]);
        }
        const Real ratePerFlow = 3 * Shear / (2 * radius * radius);
        return base * m_slitCurve(ratePerFlow * std::abs(flow)) /
          m_slitCurve(ratePerFlow * m_referenceFlow[branch]);
      }

      void solveOutlets()
      {
        auto& system = m_flow.getLinearSystem();
        const Eigen::Vector2d free = outletFlux(system.getSolution());
        Eigen::Matrix2d responseFlux;
        for (int j = 0; j < 2; ++j)
        {
          check(KSPSolve(m_ksp.getHandle(), m_loads[j], m_response[j]));
          responseFlux.col(j) = outletFlux(m_response[j]);
        }
        const auto drops = [this](const Eigen::Vector2d& q) -> Eigen::Vector2d {
          return {drop(0, q[0]), drop(1, q[1])};
        };
        Eigen::Vector2d flow = m_outletFlow;
        constexpr Real Tolerance = 1e-11;
        for (int iter = 0; iter < 30; ++iter)
        {
          const Eigen::Vector2d residual = flow - free - responseFlux * drops(flow);
          if (residual.lpNorm<Eigen::Infinity>() < Tolerance)
          {
            m_outletPressure = drops(flow);
            m_outletFlow = flow;
            for (int j = 0; j < 2; ++j)
              check(VecAXPY(system.getSolution(), m_outletPressure[j], m_response[j]));
            m_u.getSolution().setData(system.getSolution(), 0);
            m_p.getSolution().setData(system.getSolution(), m_uh.getSize());
            return;
          }
          Eigen::Vector2d derivative;
          for (int j = 0; j < 2; ++j)
          {
            const Real h = 1e-6 * m_referenceFlow[j];
            derivative[j] = (drop(j, flow[j] + h) - drop(j, flow[j] - h)) / (2 * h);
          }
          const Eigen::Matrix2d jacobian =
            Eigen::Matrix2d::Identity() - responseFlux * derivative.asDiagonal();
          const Eigen::Vector2d update = jacobian.fullPivLu().solve(residual);
          Real scale = 1;
          while (scale > 1.0 / 256)
          {
            const Eigen::Vector2d trial = flow - scale * update;
            if ((trial - free - responseFlux * drops(trial)).norm() < residual.norm())
              break;
            scale *= 0.5;
          }
          flow -= scale * update;
        }
        throw std::runtime_error("Two-outlet traction solve did not converge");
      }

      std::pair<int, Real> advance(Real head)
      {
        m_applied[0] = head;
        Real residual = 0;
        for (int iteration = 1; iteration <= MaxIterations; ++iteration)
        {
          m_flow.assemble();
          m_flow.solve(m_ksp);
          KSPConvergedReason reason;
          check(KSPGetConvergedReason(m_ksp.getHandle(), &reason));
          if (reason < 0)
            throw std::runtime_error("Flow linear solver failed");
          if (m_closure != Closure::Full)
            solveOutlets();
          const auto solution = m_u.getSolution().getData();
          check(VecWAXPY(m_difference.getData(), -1, m_iterate.getData(), solution));
          Real difference = 0;
          Real magnitude = 0;
          check(VecNorm(m_difference.getData(), NORM_INFINITY, &difference));
          check(VecNorm(solution, NORM_INFINITY, &magnitude));
          residual = difference / std::max(magnitude, 1e-12);
          m_iterate.setData(solution);
          if (residual < PicardTolerance)
          {
            m_old.setData(solution);
            return {iteration, residual};
          }
        }
        throw std::runtime_error(
          "Flow Picard iteration did not converge: " + std::to_string(residual));
      }

      MeshType m_mesh;
      VelocitySpace m_uh;
      PressureSpace m_ph;
      Trial<VelocitySpace> m_u;
      Trial<PressureSpace> m_p;
      Test<VelocitySpace> m_v;
      Test<PressureSpace> m_q;
      Field<VelocitySpace> m_old, m_iterate, m_difference;
      Field<PressureSpace> m_one;
      LinearForm<PressureSpace, ::Vec> m_section;
      FlowProblem m_flow;
      Rodin::Solver::KSP m_ksp;
      HydraulicCurve m_curve;
      SlitCurve m_slitCurve;
      Closure m_closure;
      Real m_dt;
      std::array<Attribute, 2> m_caps;
      std::array<Real, 3> m_applied{};
      ::Vec m_inletLoad = nullptr;
      std::array<::Vec, 2> m_loads{}, m_response{};
      Eigen::Vector2d m_referenceFlow, m_referenceDrop;
      Eigen::Vector2d m_outletFlow = Eigen::Vector2d::Zero();
      Eigen::Vector2d m_outletPressure = Eigen::Vector2d::Zero();
  };
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int result = 0;
  try
  {
    if (argc < 8)
    {
      throw std::runtime_error("Usage: TruncationY mesh full|constant|variable output "
                               "head_pa steps cycles calibration [PETSc options]");
    }
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm(PETSC_COMM_WORLD, boost::mpi::comm_attach);
    if (comm.size() != 1)
      throw std::runtime_error("Use one process per benchmark case");
    const Context::MPI context(env, comm);
    const std::string mode(argv[2]);
    if (mode != "full" && mode != "constant" && mode != "variable")
      throw std::runtime_error("Unknown closure: " + mode);
    const Closure closure = mode == "full"
      ? Closure::Full
      : (mode == "constant" ? Closure::Constant : Closure::Variable);
    const Real head = std::stod(argv[4]);
    const int steps = std::stoi(argv[5]);
    const int cycles = std::stoi(argv[6]);
    if (!std::isfinite(head) || head <= 0 || steps < 0 || cycles < 1)
    {
      throw std::runtime_error(
        "Head and cycles must be positive, steps must be nonnegative");
    }
    Benchmark benchmark(context, argv[1], closure, steps > 0 ? Period / steps : 0);
    benchmark.run(argv[3], head / Pressure, steps, cycles, argv[7]);
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    result = 1;
  }
  PetscFinalize();
  return result;
}
