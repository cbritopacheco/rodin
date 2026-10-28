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
#include <boost/mpi/collectives.hpp>
#include <Rodin/Scotch/MeshPartitioner.h>

namespace Rodin::Examples::Heart::Atrium
{
  Atrium::Atrium(const Context::MPI& context, const Parameters& cfg)
    : m_mesh(makeMesh(context, cfg)),
      m_cavity(m_mesh, cfg.fluidMaterial),
      m_csvPath(cfg.mode == Mode::Fluid ? cfg.fluid.csvPath : cfg.solid.output + ".csv")
  {
    if (cfg.mode != Mode::Solid && cfg.mode != Mode::Fluid && cfg.mode != Mode::FSI)
      Alert::Exception() << "Invalid atrium mode" << Alert::Raise;
    const bool hasSolid = cfg.mode != Mode::Fluid;
    const bool hasFluid = cfg.mode != Mode::Solid;
    if (cfg.mode == Mode::FSI)
    {
      const auto dt = cfg.solid.timeStep;
      const auto fluidDt = cfg.fluid.dt;
      const auto tolerance = 64 * std::numeric_limits<Real>::epsilon() *
        std::max(std::abs(dt), std::abs(fluidDt));
      if (cfg.solid.fibersOnly || !std::isfinite(dt) || !std::isfinite(fluidDt) ||
        dt <= 0 || fluidDt <= 0 || std::abs(dt - fluidDt) > tolerance)
      {
        Alert::Exception() << "FSI components require matching positive "
                              "time steps and an active solid mechanics problem"
                           << Alert::Raise;
      }
    }
    if (hasSolid)
    {
      auto parameters = cfg.solid;
      parameters.writeCSV = false;
      m_solid.emplace(m_mesh, parameters, cfg.solidMaterial);
    }
    if (hasFluid)
    {
      auto parameters = cfg.fluid;
      parameters.writeCSV = false;
      m_fluid.emplace(m_mesh, parameters, cfg.fluidMaterial);
    }
  }

  Atrium::MeshType Atrium::makeMesh(const Context::MPI& context, const Parameters& cfg)
  {
    const auto& comm = context.getCommunicator();
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
    m_failed = true;
    if (m_solid)
      m_solid->initialize();
    if (m_fluid)
      m_fluid->initialize();
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
    m_failed = true;
    if (m_solid)
      m_solid->step();
    if (m_fluid)
      m_fluid->step();
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
