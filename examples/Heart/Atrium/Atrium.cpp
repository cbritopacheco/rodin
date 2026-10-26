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
  struct Atrium::Impl
  {
    explicit Impl(const Parameters& cfg)
    {
      if (cfg.mode != Mode::SolidOnly)
      {
        Alert::Exception() << "Fluid modes require an MPI context and PETSc"
                           << Alert::Raise;
      }
      solid = std::make_unique<Solid>(cfg.solid);
    }

#ifdef RODIN_ATRIUM_HAS_FLUID
    Impl(const Context::MPI& context, const Parameters& cfg)
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
        solid = std::make_unique<Solid>(cfg.solid);
      if (hasFluid)
        fluid = std::make_unique<Fluid>(context, cfg.fluid);
    }
    std::unique_ptr<Fluid> fluid;
#endif
    std::unique_ptr<Solid> solid;
    size_t steps = 0;
    bool initialized = false;
    bool finished = false;
    bool failed = false;
  };

  Atrium::Atrium(const Parameters& parameters)
    : m_impl(std::make_unique<Impl>(parameters))
  {}
#ifdef RODIN_ATRIUM_HAS_FLUID
  Atrium::Atrium(const Context::MPI& context, const Parameters& parameters)
    : m_impl(std::make_unique<Impl>(context, parameters))
  {}
#endif
  Atrium::~Atrium() = default;

  Atrium& Atrium::initialize()
  {
    if (m_impl->finished || m_impl->failed)
    {
      Alert::Exception() << "Cannot initialize a finished or failed atrium"
                         << Alert::Raise;
    }
    if (m_impl->initialized)
      return *this;
    m_impl->failed = true;
    if (m_impl->solid)
      m_impl->solid->initialize();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_impl->fluid)
      m_impl->fluid->initialize();
#endif
    m_impl->initialized = true;
    m_impl->failed = false;
    return *this;
  }

  void Atrium::step()
  {
    initialize();
    m_impl->failed = true;
    if (m_impl->solid)
      m_impl->solid->step();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_impl->fluid)
      m_impl->fluid->step();
#endif
    ++m_impl->steps;
    m_impl->failed = false;
  }

  void Atrium::finish()
  {
    if (m_impl->finished)
      return;
    if (m_impl->solid)
      m_impl->solid->finish();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_impl->fluid)
      m_impl->fluid->finish();
#endif
    m_impl->finished = true;
  }

  size_t Atrium::getStepCount() const noexcept
  {
    return m_impl->steps;
  }
  Real Atrium::getTime() const noexcept
  {
    if (m_impl->solid)
      return m_impl->solid->getTime();
#ifdef RODIN_ATRIUM_HAS_FLUID
    if (m_impl->fluid)
      return m_impl->fluid->getTime();
#endif
    return 0;
  }
  const Atrium::Solid* Atrium::getSolid() const noexcept
  {
    return m_impl->solid.get();
  }
#ifdef RODIN_ATRIUM_HAS_FLUID
  const Atrium::Fluid* Atrium::getFluid() const noexcept
  {
    return m_impl->fluid.get();
  }
#endif
}
