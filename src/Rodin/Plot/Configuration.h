/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_CONFIGURATION_H
#define RODIN_PLOT_CONFIGURATION_H

namespace Rodin::Plot
{
  class Configuration
  {
    public:
      Configuration();

      /**
       * Sets the sample count to be used when rendering anti-aliased objects.
       * @param sampleCount Number of samples to use.
       * @returns Reference to self for method chaining.
       */
      Configuration& setSampleCount(int sampleCount);

      /**
       * @returns Number of samples used when rendering anti-aliased objects.
       */
      int getSampleCount() const;

    private:
      int m_sampleCount;
  };
}

#endif
