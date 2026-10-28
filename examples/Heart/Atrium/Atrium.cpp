/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cmath>
#include <limits>
#include <filesystem>
#include <iomanip>
#include <Rodin/Alert.h>
#include "Atrium.h"
#include <boost/mpi/collectives.hpp>
#include <Rodin/Scotch/MeshPartitioner.h>

namespace Rodin::Examples::Heart
{
  Atrium::Atrium(const Context::MPI& context, const Parameters& cfg)
    : m_mesh(makeMesh(context, cfg)),
      m_cavity(m_mesh, cfg.fluidMaterial),
      m_volumePath((cfg.mode == Mode::Fluid ? cfg.fluid.xdmfBasename : cfg.solid.output) +
        ".volume.csv")
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
      m_solid.emplace(m_mesh, cfg.solid, cfg.solidMaterial);
    if (hasFluid)
      m_fluid.emplace(m_mesh, cfg.fluid, cfg.fluidMaterial);
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
        const auto parent = std::filesystem::path(m_volumePath).parent_path();
        if (!parent.empty())
          std::filesystem::create_directories(parent);
        m_volumeOutput.open(m_volumePath);
        if (!m_volumeOutput)
          failure = "Cannot open cavity-volume output " + m_volumePath;
        else
        {
          m_volumeOutput << "time_s,reference_volume_ml,solid_cavity_volume_ml,fluid_"
                            "domain_volume_ml,solid_volume_change_percent\n";
        }
      }
      catch (const std::exception& error)
      {
        failure = error.what();
      }
    }
    boost::mpi::broadcast(comm, failure, 0);
    if (!failure.empty())
      Alert::Exception() << failure << Alert::Raise;
    writeVolume();
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
    writeVolume();
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
    if (m_volumeOutput.is_open())
      m_volumeOutput.close();
    m_finished = true;
  }

  void Atrium::writeVolume()
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
      m_volumeOutput << std::setprecision(17) << getTime() << ',' << reference << ',';
      if (m_solidCavityVolume)
        m_volumeOutput << *m_solidCavityVolume;
      m_volumeOutput << ',';
      if (m_fluidDomainVolume)
        m_volumeOutput << *m_fluidDomainVolume;
      m_volumeOutput << ',';
      if (m_solidCavityVolume)
        m_volumeOutput << 100 * (*m_solidCavityVolume / reference - 1);
      m_volumeOutput << '\n';
      m_volumeOutput.flush();
      if (!m_volumeOutput)
        failure = "Cannot write cavity-volume output " + m_volumePath;
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
