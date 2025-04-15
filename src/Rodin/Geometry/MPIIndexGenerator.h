#ifndef RODIN_GEOMETRY_MPIINDEXGENERATOR_H
#define RODIN_GEOMETRY_MPIINDEXGENERATOR_H

#include "IndexGenerator.h"

namespace Rodin::Geometry
{
  class DistributedIndexGenerator : public IndexGeneratorBase
  {
    public:
      DistributedIndexGenerator(size_t globalOffset, size_t localCount)
        : m_globalOffset(globalOffset),
          m_localCount(localCount),
          m_current(globalOffset)
      {}

      bool end() const override
      {
        return m_current >= m_globalOffset + m_localCount;
      }

      DistributedIndexGenerator& operator++() override
      {
        ++m_current;
        return *this;
      }

      Index operator*() const noexcept override
      {
        return m_current - m_globalOffset;
      }

      DistributedIndexGenerator* copy() const noexcept override
      {
        return new DistributedIndexGenerator(*this);
      }

      DistributedIndexGenerator* move() noexcept override
      {
        return new DistributedIndexGenerator(std::move(*this));
      }

    private:
      size_t m_globalOffset;  // Global offset for this rank.
      size_t m_localCount;    // How many indices belong to this rank.
      size_t m_current;       // Current global index.
  };
}

#endif
