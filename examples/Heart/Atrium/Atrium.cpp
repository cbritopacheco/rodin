/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cmath>
#include <limits>
#include <Rodin/Alert.h>
#include "Atrium.h"

namespace Rodin::Examples::Heart
{
  Atrium::Atrium(const Parameters& cfg)
  {
    if (cfg.mode != Mode::SolidOnly)
    {
      Alert::Exception() << "Fluid modes require an MPI context and PETSc"
                         << Alert::Raise;
    }
    m_solid = std::make_unique<Solid>(cfg.solid);
  }

#ifdef RODIN_ATRIUM_HAS_FLUID
  Atrium::Atrium(const Context::MPI& context, const Parameters& cfg)
  {
    if (cfg.mode != Mode::SolidOnly && cfg.mode != Mode::FluidOnly &&
      cfg.mode != Mode::Independent)
      Alert::Exception() << "Invalid atrium mode" << Alert::Raise;
    const bool hasSolid = cfg.mode != Mode::FluidOnly;
    const bool hasFluid = cfg.mode != Mode::SolidOnly;
    if (hasSolid && context.getCommunicator().size() != 1)
    {
      Alert::Exception() << "The local atrial solid requires one MPI rank"
                         << Alert::Raise;
    }
    if (cfg.mode == Mode::Independent)
    {
      const auto dt = cfg.solid.timeStep;
      const auto fluidDt = cfg.fluid.dt;
      const auto tolerance = 64 * std::numeric_limits<Real>::epsilon() *
        std::max(std::abs(dt), std::abs(fluidDt));
      if (cfg.solid.fibersOnly || !std::isfinite(dt) || !std::isfinite(fluidDt) ||
        dt <= 0 || fluidDt <= 0 || std::abs(dt - fluidDt) > tolerance)
      {
        Alert::Exception() << "Independent components require matching positive "
                              "time steps and an active solid mechanics problem"
                           << Alert::Raise;
      }
    }
    if (hasSolid)
      m_solid = std::make_unique<Solid>(cfg.solid);
    if (hasFluid)
      m_fluid = std::make_unique<Fluid>(context, cfg.fluid);
  }

#endif

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
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_fluid)
      m_fluid->initialize();
#endif
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
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_fluid)
      m_fluid->step();
#endif
    ++m_steps;
    m_failed = false;
  }

  void Atrium::finish()
  {
    if (m_finished)
      return;
    if (m_solid)
      m_solid->finish();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_fluid)
      m_fluid->finish();
#endif
    m_finished = true;
  }

  size_t Atrium::getStepCount() const noexcept
  {
    return m_steps;
  }
  Real Atrium::getTime() const noexcept
  {
    if (m_solid)
      return m_solid->getTime();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_fluid)
      return m_fluid->getTime();
#endif
    return 0;
  }
  const Atrium::Solid* Atrium::getSolid() const noexcept
  {
    return m_solid.get();
  }
#ifdef RODIN_ATRIUM_HAS_FLUID
  const Atrium::Fluid* Atrium::getFluid() const noexcept
  {
    return m_fluid.get();
  }
#endif
}
