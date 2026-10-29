/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <utility>
#include <cmath>
#include <limits>
#include <filesystem>
#include <iomanip>
#include <Rodin/Alert.h>
#include "Atrium.h"
#include "Log.h"
#include <boost/mpi/collectives.hpp>
#include <Rodin/Scotch/MeshPartitioner.h>

namespace Rodin::Examples::Heart::Atrium
{
  Atrium::Atrium(const Context::MPI& context, const Parameters& cfg)
    : m_parameters(cfg.validated()),
      m_mesh(makeMesh(context, m_parameters)),
      m_cavity(m_mesh, m_parameters.fluidMaterial),
      m_csvPath(m_parameters.outputDirectory + "/Diagnostics.csv")
  {
    if (m_parameters.mode != Mode::Fluid)
    {
      Stage stage("Construct solid", context.getCommunicator().rank() == 0);
      m_solid.emplace(
        m_mesh, m_parameters.solid, m_parameters, m_parameters.solidMaterial);
    }
    if (m_parameters.mode != Mode::Solid)
    {
      Stage stage("Construct fluid", context.getCommunicator().rank() == 0);
      m_fluid.emplace(
        m_mesh, m_parameters.fluid, m_parameters, m_parameters.fluidMaterial);
    }
    if (m_parameters.mode == Mode::FSI)
      m_fsi.emplace(*m_solid, *m_fluid, m_parameters.fsi);
  }

  Atrium::MeshType Atrium::makeMesh(const Context::MPI& context, const Parameters& cfg)
  {
    const auto& comm = context.getCommunicator();
    Stage stage("Mesh loading and Scotch partitioning", comm.rank() == 0);
    constexpr int Root = 0;
    Geometry::LocalMesh local;
    std::string failure;
    if (comm.rank() == Root)
    {
      try
      {
        local.load(cfg.meshPath, IO::FileFormat::MEDIT);
        if (local.getDimension() != 3 || local.getSpaceDimension() != 3 ||
          local.getCellCount() == 0)
          Alert::Exception() << "Expected a nonempty combined 3D mesh" << Alert::Raise;
        for (size_t d = 1; d <= 3; ++d)
        {
          local.getConnectivity().compute(d, 0);
          local.getConnectivity().compute(d, d);
          local.getConnectivity().compute(d, d - 1);
          local.getConnectivity().compute(d - 1, d);
        }
      }
      catch (const std::exception& error)
      {
        failure = error.what();
      }
    }
    boost::mpi::broadcast(comm, failure, Root);
    if (!failure.empty())
      Alert::Exception() << failure << Alert::Raise;
    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == Root)
    {
      try
      {
        Scotch::Partitioner partitioner(local);
        partitioner.partition(static_cast<size_t>(comm.size()));
        sharder.shard(partitioner);
      }
      catch (const std::exception& error)
      {
        failure = error.what();
      }
    }
    boost::mpi::broadcast(comm, failure, Root);
    if (!failure.empty())
      Alert::Exception() << failure << Alert::Raise;
    if (comm.rank() == Root)
      sharder.scatter(Root);
    auto mesh = sharder.gather(Root);
    for (size_t d = 1; d <= 3; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(d, d);
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
    }
    mesh.reconcile(2);
    mesh.reconcile(1);
    return mesh;
  }

  Atrium::~Atrium() = default;

  Atrium& Atrium::initialize()
  {
    if (m_finished || m_failed)
    {
      Alert::Exception() << "Cannot initialize a finished or failed atrium"
                         << Alert::Raise;
    }
    if (m_initialized)
      return *this;
    Stage stage("Initialization", m_mesh.getContext().getCommunicator().rank() == 0);
    m_failed = true;
    if (m_solid)
      m_solid->initialize();
    if (m_fluid)
      m_fluid->initialize();
    if (m_fsi)
      m_fsi->initialize();
    std::string failure;
    const auto& comm = m_mesh.getContext().getCommunicator();
    if (comm.rank() == 0)
    {
      try
      {
        const auto parent = std::filesystem::path(m_csvPath).parent_path();
        if (!parent.empty())
          std::filesystem::create_directories(parent);
        m_csvOutput.open(m_csvPath);
        if (!m_csvOutput)
          failure = "Cannot open run diagnostics " + m_csvPath;
      }
      catch (const std::exception& error)
      {
        failure = error.what();
      }
    }
    boost::mpi::broadcast(comm, failure, 0);
    if (!failure.empty())
      Alert::Exception() << failure << Alert::Raise;
    writeDiagnostics();
    m_initialized = true;
    m_failed = false;
    return *this;
  }

  void Atrium::step()
  {
    initialize();
    const bool root = m_mesh.getContext().getCommunicator().rank() == 0;
    Stage stage("Timestep " + std::to_string(m_steps + 1) + "/" +
        std::to_string(m_parameters.steps),
      root);
    if (root)
    {
      Log("Time") << "time=" << getTime() + m_parameters.timeStep
                  << " s dt=" << m_parameters.timeStep << " s" << Alert::Raise;
    }
    m_failed = true;
    if (m_fsi)
      m_fsi->step();
    else
    {
      if (m_solid)
        m_solid->step();
      if (m_fluid)
        m_fluid->step();
    }
    writeDiagnostics();
    ++m_steps;
    m_failed = false;
  }

  void Atrium::finish()
  {
    if (m_finished)
      return;
    if (m_solid)
      m_solid->finish();
    if (m_fluid)
      m_fluid->finish();
    if (m_csvOutput.is_open())
      m_csvOutput.close();
    m_finished = true;
  }

  void Atrium::writeDiagnostics()
  {
    Stage stage("Cavity volume and CSV diagnostics",
      m_mesh.getContext().getCommunicator().rank() == 0);
    const Real reference = m_cavity.getReferenceVolume() * CubicMillimetresToMillilitres;
    if (m_solid)
    {
      m_solidCavityVolume =
        m_cavity.evaluate(m_solid->getDisplacement()) * CubicMillimetresToMillilitres;
    }
    if (m_fluid)
    {
      m_fluidDomainVolume =
        CavityVolume::measure(m_fluid->getMesh()) * CubicMetresToMillilitres;
    }
    const auto& comm = m_mesh.getContext().getCommunicator();
    std::string failure;
    if (comm.rank() == 0)
    {
      const auto solid =
        m_solid ? m_solid->getDiagnostics() : Optional<Solid::Diagnostics>{};
      const auto fluid =
        m_fluid ? m_fluid->getDiagnostics() : Optional<Fluid::Diagnostics>{};
      const auto fields = std::to_array<std::pair<StringView, Optional<Real>>>(
        {{"time_s", getTime()}, {"reference_volume_ml", reference},
          {"solid_cavity_volume_ml", m_solidCavityVolume},
          {"fluid_domain_volume_ml", m_fluidDomainVolume},
          {"solid_volume_change_percent",
            m_solidCavityVolume
              ? Optional<Real>(100 * (*m_solidCavityVolume / reference - 1))
              : Optional<Real>{}},
          {"max_displacement_mm",
            solid ? Optional<Real>(solid->maxDisplacement) : Optional<Real>{}},
          {"activation_max",
            solid ? Optional<Real>(solid->activationMax) : Optional<Real>{}},
          {"newton_iterations",
            solid ? Optional<Real>(solid->newtonIterations) : Optional<Real>{}},
          {"residual", solid ? Optional<Real>(solid->residual) : Optional<Real>{}},
          {"min_scaled_jacobian_mm",
            solid ? Optional<Real>(solid->minimumJacobian) : Optional<Real>{}},
          {"max_velocity_mm_per_s",
            solid ? Optional<Real>(solid->maxVelocity) : Optional<Real>{}},
          {"max_acceleration_mm_per_s2",
            solid ? Optional<Real>(solid->maxAcceleration) : Optional<Real>{}},
          {"min_active_extension",
            solid ? Optional<Real>(solid->minimumExtension) : Optional<Real>{}},
          {"max_active_extension",
            solid ? Optional<Real>(solid->maximumExtension) : Optional<Real>{}},
          {"max_contractile_stress_kPa",
            solid ? Optional<Real>(solid->maximumStress) : Optional<Real>{}},
          {"fluid_cycle", fluid ? Optional<Real>(fluid->cycle) : Optional<Real>{}},
          {"fluid_inlet_pressure_pa",
            fluid ? Optional<Real>(fluid->inletPressure) : Optional<Real>{}},
          {"fluid_outlet_pressure_pa",
            fluid ? Optional<Real>(fluid->outletPressure) : Optional<Real>{}},
          {"fluid_arterial_pressure_pa",
            fluid ? fluid->arterialPressure : Optional<Real>{}},
          {"fluid_distal_pressure_pa", fluid ? fluid->distalPressure : Optional<Real>{}},
          {"fluid_valve_open",
            fluid ? Optional<Real>(fluid->valveOpen) : Optional<Real>{}},
          {"fluid_inlet_flux_m3_per_s",
            fluid ? Optional<Real>(fluid->inletFlux) : Optional<Real>{}},
          {"fluid_outlet_flux_m3_per_s",
            fluid ? Optional<Real>(fluid->outletFlux) : Optional<Real>{}},
          {"fluid_max_velocity_m_per_s",
            fluid ? Optional<Real>(fluid->maximumVelocity) : Optional<Real>{}},
          {"fluid_linear_iterations",
            fluid ? Optional<Real>(fluid->linearIterations) : Optional<Real>{}},
          {"fluid_pc_age_timesteps",
            fluid ? Optional<Real>(fluid->preconditionerAge) : Optional<Real>{}},
          {"fluid_pc_refreshes",
            fluid ? Optional<Real>(fluid->preconditionerRefreshes) : Optional<Real>{}},
          {"fluid_pc_refreshed",
            fluid ? Optional<Real>(fluid->preconditionerRefreshed) : Optional<Real>{}},
          {"fluid_pc_retried",
            fluid ? Optional<Real>(fluid->preconditionerRetried) : Optional<Real>{}},
          {"fsi_iterations",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().iterations)
                  : Optional<Real>{}},
          {"fsi_traction_data_update_pa",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().tractionDataUpdate)
                  : Optional<Real>{}},
          {"fsi_displacement_change",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().displacementChange)
                  : Optional<Real>{}},
          {"fsi_displacement_tolerance",
            m_fsi ? Optional<Real>(m_parameters.fsi.relativeTolerance)
                  : Optional<Real>{}},
          {"fsi_converged",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().converged) : Optional<Real>{}},
          {"fsi_budget_exhausted",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().budgetExhausted)
                  : Optional<Real>{}},
          {"fsi_velocity_mismatch_m_per_s",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().velocityMismatch)
                  : Optional<Real>{}},
          {"fsi_newmark_velocity_mismatch_m_per_s",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().newmarkVelocityMismatch)
                  : Optional<Real>{}},
          {"fsi_robin_coefficient_pa_s_per_m",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().robinCoefficient)
                  : Optional<Real>{}},
          {"fsi_displacement_mismatch_m",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().displacementMismatch)
                  : Optional<Real>{}},
          {"fsi_min_mesh_jacobian",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().minimumMeshJacobian)
                  : Optional<Real>{}},
          {"fluid_wall_flux_m3_per_s",
            m_fsi ? Optional<Real>(m_fsi->getDiagnostics().wallFlux) : Optional<Real>{}},
          {"fluid_total_flux_m3_per_s",
            m_fsi && fluid ? Optional<Real>(fluid->inletFlux + fluid->outletFlux +
                               m_fsi->getDiagnostics().wallFlux)
                           : Optional<Real>{}},
          {"fsi_cavity_volume_difference_ml",
            m_fsi && m_solidCavityVolume && m_fluidDomainVolume
              ? Optional<Real>(*m_fluidDomainVolume - *m_solidCavityVolume)
              : Optional<Real>{}},
          {"fluid_max_TAWSS_pa",
            fluid ? Optional<Real>(fluid->maximumTAWSS) : Optional<Real>{}},
          {"fluid_max_OSI", fluid ? Optional<Real>(fluid->maximumOSI) : Optional<Real>{}},
          {"fluid_max_activation",
            fluid ? Optional<Real>(fluid->maximumActivation) : Optional<Real>{}},
          {"fluid_max_thrombin",
            fluid ? Optional<Real>(fluid->maximumThrombin) : Optional<Real>{}},
          {"fluid_min_fibrinogen",
            fluid ? Optional<Real>(fluid->minimumFibrinogen) : Optional<Real>{}},
          {"fluid_max_fibrin",
            fluid ? Optional<Real>(fluid->maximumFibrin) : Optional<Real>{}}});
      if (m_csvOutput.tellp() == 0)
      {
        for (size_t i = 0; i < fields.size(); ++i)
          m_csvOutput << (i == 0 ? "" : ",") << fields[i].first;
        m_csvOutput << '\n';
      }
      m_csvOutput << std::setprecision(17);
      for (size_t i = 0; i < fields.size(); ++i)
      {
        if (i != 0)
          m_csvOutput << ',';
        if (fields[i].second)
          m_csvOutput << *fields[i].second;
      }
      m_csvOutput << '\n';
      m_csvOutput.flush();
      if (!m_csvOutput)
        failure = "Cannot write run diagnostics " + m_csvPath;
    }
    boost::mpi::broadcast(comm, failure, 0);
    if (!failure.empty())
      Alert::Exception() << failure << Alert::Raise;
  }

  Optional<Real> Atrium::getSolidCavityVolume() const noexcept
  {
    return m_solidCavityVolume;
  }

  Optional<Real> Atrium::getFluidDomainVolume() const noexcept
  {
    return m_fluidDomainVolume;
  }

  size_t Atrium::getStepCount() const noexcept
  {
    return m_steps;
  }
  Real Atrium::getTime() const noexcept
  {
    if (m_solid)
      return m_solid->getTime();
    if (m_fluid)
      return m_fluid->getTime();
    return 0;
  }
  const Atrium::MeshType& Atrium::getMesh() const noexcept
  {
    return m_mesh;
  }

  const Atrium::Solid* Atrium::getSolid() const noexcept
  {
    return m_solid ? &*m_solid : nullptr;
  }
  const Atrium::Fluid* Atrium::getFluid() const noexcept
  {
    return m_fluid ? &*m_fluid : nullptr;
  }
}
