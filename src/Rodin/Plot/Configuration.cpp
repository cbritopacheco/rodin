/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Configuration.h"

namespace Rodin::Plot
{
  Configuration::Configuration()
  {
    setSampleCount(8);
  }

  Configuration& Configuration::setSampleCount(int sampleCount)
  {
    m_sampleCount = sampleCount;
    return *this;
  }

  int Configuration::getSampleCount() const
  {
    return m_sampleCount;
  }
}
