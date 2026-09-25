/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_ASSEMBLY_BOUNDARYDOFS_H
#define RODIN_MPI_ASSEMBLY_BOUNDARYDOFS_H

#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include "Rodin/Types.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @brief Selects boundary DOF functionals using the existing mesh overlap.
   *
   * Architecture: for locally supported spaces, the required set is the union
   * of owned DOFs and DOFs of owned cells. Vertex-overlap shards contain every
   * incident face of these DOFs. Select the smallest distributed face index
   * among eligible faces, including ghosts. Artificial boundaries at the outer
   * edge of the overlap cannot touch this required set. No communication is
   * performed for these spaces.
   *
   * P0g is globally supported and is deliberately different: select the lowest
   * eligible owned face globally and broadcast its evaluated payload. Ownership
   * and numbering are never changed. No numerical comparison selects a source.
   */
  template <class FES>
  class MPIBoundaryDOFs
  {
    public:
      using Scalar = typename FES::ScalarType;
      static constexpr bool Global = std::is_same_v<FES,
        Variational::P0g<typename FES::RangeType, typename FES::MeshType>>;

      MPIBoundaryDOFs(const FES& fes, const FlatSet<Geometry::Attribute>& attributes)
        : m_fes(fes)
      {
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const size_t dim = mesh.getDimension();
        IndexSet required;
        if constexpr (!Global)
        {
          Index begin, end;
          fes.getOwnershipRange(begin, end);
          for (Index i = begin; i < end; ++i)
            required.insert(i);
          for (auto cell = mesh.getCell(); cell; ++cell)
            if (shard.isOwned(dim, cell->getIndex()))
              for (Index dof : fes.getDOFs(dim, cell->getIndex()))
                required.insert(dof);
        }

        IndexMap<Index> sources;
        Index first = std::numeric_limits<Index>::max();
        for (auto face = mesh.getFace(); face; ++face)
        {
          const Index i = face->getIndex();
          if constexpr (Global)
            if (!shard.isOwned(dim - 1, i))
              continue;
          if (attributes.empty())
          {
            if (!shard.isBoundary(i))
              continue;
          }
          else if (!face->getAttribute() || !attributes.contains(*face->getAttribute()))
            continue;
          const Index source = shard.getPolytopeMap(dim - 1).left.at(i);
          first = std::min(first, source);
          const auto dofs = fes.getDOFs(dim - 1, i);
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
          {
            const Index global = dofs[local];
            if constexpr (!Global)
              if (!required.contains(global))
                continue;
            const auto found = sources.find(global);
            if (found == sources.end() || source < found->second)
            {
              sources[global] = source;
              m_dofs[global] = {i, local};
            }
          }
        }
        if constexpr (Global)
        {
          const auto& comm = mesh.getContext().getCommunicator();
          const Index selected = boost::mpi::all_reduce(comm, first, boost::mpi::minimum<Index>());
          if (selected == std::numeric_limits<Index>::max())
            return;
          m_sourceRank = boost::mpi::all_reduce(comm,
            first == selected ? comm.rank() : comm.size(), boost::mpi::minimum<int>());
          if (comm.rank() != m_sourceRank)
            m_dofs.clear();
        }
      }

      /// Global DOF -> (shard-local face, face-local functional ordinal).
      const auto& getDOFs() const { return m_dofs; }

      /** Broadcasts only the globally supported P0g payload; otherwise a no-op. */
      template <class Payload>
      void synchronize(IndexMap<Payload>& values) const
      {
        if constexpr (Global)
        {
          if (m_sourceRank < 0)
            return;
          std::vector<std::pair<Index, Payload>> entries(values.begin(), values.end());
          boost::mpi::broadcast(m_fes.getMesh().getContext().getCommunicator(), entries, m_sourceRank);
          values.clear();
          for (auto& entry : entries)
            values.emplace(std::move(entry));
        }
      }

      /** Evaluates prescribed values or affine offsets at the selected functionals. */
      template <class Function>
      void assemble(IndexMap<Scalar>& values, const Function& function) const
      {
        values.clear();
        const size_t dim = m_fes.getMesh().getDimension() - 1;
        for (const auto& [global, indices] : m_dofs)
        {
          const auto [face, local] = indices;
          const auto& element = m_fes.getFiniteElement(dim, face);
          const auto mapping = m_fes.getPullback({dim, face}, function);
          values.emplace(global, element.getLinearForm(local)(mapping));
        }
        synchronize(values);
      }

    private:
      const FES& m_fes;
      IndexMap<std::pair<Index, Index>> m_dofs;
      int m_sourceRank = -1;
  };
}

#endif
