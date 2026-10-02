/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Distributed matrix finite element ownership and shard maps.
 */
#ifndef RODIN_MPI_VARIATIONAL_MATRIXRANGE_H
#define RODIN_MPI_VARIATIONAL_MATRIXRANGE_H

#include "Rodin/Variational/MatrixRange.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational::Detail
{
  /** Distributed component expansion of an existing scalar space. */
  template <class Derived, class ScalarSpace, class LocalSpace, bool CellsOnly = false>
  class DistributedMatrixSpace : public MatrixSpace<Derived, ScalarSpace,
                                   typename LocalSpace::ElementType, CellsOnly>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent =
        MatrixSpace<Derived, ScalarSpace, typename LocalSpace::ElementType, CellsOnly>;
      /// @brief Finite element space type.
      using FESType = LocalSpace;
      using Parent::getGlobalIndex;

      /// @brief Expands scalar ownership and ghost maps for a local matrix shard.
      DistributedMatrixSpace(
        ScalarSpace scalar, LocalSpace shard, size_t rows, size_t cols)
        : Parent(std::move(scalar), rows, cols),
          m_shard(std::move(shard))
      {
        for (size_t i = 0; i < m_shard.getSize(); ++i)
          m_globalToLocal.emplace(getGlobalIndex(i), i);
      }

      /// @brief Returns the matrix space on the local mesh shard.
      const LocalSpace& getShard() const
      {
        return m_shard;
      }

      /// @brief Returns the half-open range of owned global component DOFs.
      void getOwnershipRange(Index& begin, Index& end) const
      {
        this->getScalarSpace().getOwnershipRange(begin, end);
        begin *= this->getVectorDimension();
        end *= this->getVectorDimension();
      }

      /// @brief Maps a local component DOF to its global coefficient index.
      Index getGlobalIndex(Index local) const
      {
        assert(static_cast<size_t>(local) < m_shard.getSize());
        const Index components = this->getVectorDimension();
        return this->getScalarSpace().getGlobalIndex(local / components) * components +
          local % components;
      }

      /// @brief Returns the shard-local index of a global DOF, if present.
      Optional<Index> getLocalIndex(Index global) const
      {
        const auto it = m_globalToLocal.find(global);
        if (it == m_globalToLocal.end())
          return std::nullopt;
        return it->second;
      }

    private:
      LocalSpace m_shard;
      std::map<Index, Index> m_globalToLocal;
  };
}

#endif
