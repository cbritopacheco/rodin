/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_VARIATIONAL_INTERPOLATION_H
#define RODIN_MPI_VARIATIONAL_INTERPOLATION_H

#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/Variational/Interpolation.h"
#include "Rodin/Variational/P0/ForwardDecls.h"
#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Selects and evaluates distributed interpolation functionals.
   *
   * Architecture: for each owned DOF @f$ i @f$, select the eligible incident
   * entity @f$ K_i @f$ with the smallest distributed index and evaluate
   * @f$ c_i = \ell_{i,K_i}(f) @f$ using the space's pullback and functional.
   * The vertex overlap contains the complete incident star of owned DOFs
   * for P0, P1, and H1. Selection and evaluation for these spaces require no
   * communication. Only owned coefficients are returned; the storage backend
   * commits them and refreshes ghost copies from their owners.
   *
   * P0g is globally supported. Its source is the smallest eligible owned
   * entity globally. The selected rank broadcasts the evaluated coefficients,
   * which are retained only by their DOF owner, even when its shard is empty.
   * An empty eligible set returns no coefficients and therefore preserves
   * the destination. This operation is DOF interpolation, not an
   * @f$ L^2 @f$ projection.
   *
   * @pre The predicate and source function describe the same mathematical
   * entity on every replica. Rank-dependent predicates do not define a
   * distributed region. Required incidences must have been computed before
   * constructing the space, as for ordinary region traversal.
   *
   * @note Boundary and interface selection examines overlap faces incident
   * to owned DOFs. Their complete stars exclude artificial overlap boundaries.
   * Globally supported spaces select only owned faces.
   *
   * @par Communication and storage
   * Locally supported spaces store only source records for owned DOFs and
   * perform no communication during selection or evaluation. The storage
   * backend refreshes its existing ghost layer. P0g uses two scalar reductions
   * to select its source and broadcasts only its coefficient payload; no
   * mesh entities or mesh-sized coefficient arrays are gathered.
   */
  template <class FES>
    requires std::is_same_v<typename FES::ContextType, Context::MPI>
  class Interpolation<FES> final
  {
    public:
      /// Coefficient scalar type, independent of the storage backend.
      using Scalar = typename FES::ScalarType;

      /// Rank identifier in the mesh communicator.
      using Rank = int;

      /// Whether the space has globally supported constant DOFs.
      static constexpr bool Global =
        std::is_same_v<FES, P0g<typename FES::RangeType, typename FES::MeshType>>;

      /** Selects the functional source for each coefficient to be updated. */
      template <class Pred>
      Interpolation(const FES& fes, Geometry::Region region, const Pred& pred)
        : m_fes(fes),
          m_dimension(0)
      {
        if constexpr (std::is_same_v<FES,
                        P0<typename FES::RangeType, typename FES::MeshType>>)
        {
          if (region != Geometry::Region::Cells)
            Alert::Exception() << "P0 interpolation requires the cell region."
                               << Alert::Raise;
        }
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const size_t dim = mesh.getDimension();
        if (region != Geometry::Region::Cells && dim == 0)
          return;
        m_dimension = region == Geometry::Region::Cells ? dim : dim - 1;
        Index begin, end;
        fes.getOwnershipRange(begin, end);

        using Candidate = std::pair<Index, std::pair<Index, Index>>;
        UnorderedMap<Index, Candidate> candidates;
        Index first = std::numeric_limits<Index>::max();
        for (auto entity = mesh.getPolytope(m_dimension); entity; ++entity)
        {
          const Index local = entity->getIndex();
          if constexpr (Global)
            if (!shard.isOwned(m_dimension, local))
              continue;
          if (region == Geometry::Region::Boundary && !shard.isBoundary(local))
            continue;
          if (region == Geometry::Region::Interface && !shard.isInterface(local))
            continue;
          if (!pred(*entity))
            continue;
          const Index source = mesh.getGlobalIndex(m_dimension, local);
          first = std::min(first, source);
          const auto dofs = fes.getDOFs(m_dimension, local);
          for (Index ordinal = 0; ordinal < static_cast<Index>(dofs.size()); ++ordinal)
          {
            const Index global = dofs[ordinal];
            if constexpr (!Global)
              if (global < begin || global >= end)
                continue;
            const auto found = candidates.find(global);
            if (found == candidates.end() || source < found->second.first)
              candidates[global] = {source, {local, ordinal}};
          }
        }
        std::vector<std::pair<Index, std::pair<Index, Index>>> ordered;
        ordered.reserve(candidates.size());
        for (const auto& [global, candidate] : candidates)
          ordered.emplace_back(global, candidate.second);
        std::sort(ordered.begin(), ordered.end());
        m_dofs.reserve(ordered.size());
        for (const auto& [global, indices] : ordered)
          m_dofs.emplace_hint(m_dofs.end(), global, indices);

        if constexpr (Global)
        {
          const auto& comm = mesh.getContext().getCommunicator();
          const Index selected =
            boost::mpi::all_reduce(comm, first, boost::mpi::minimum<Index>());
          if (selected == std::numeric_limits<Index>::max())
            return;
          m_source = boost::mpi::all_reduce(comm,
            first == selected ? comm.rank() : comm.size(), boost::mpi::minimum<Rank>());
          if (comm.rank() != *m_source)
            m_dofs.clear();
        }
      }

      /// Global DOF -> (shard-local entity, entity-local functional ordinal).
      const auto& getDOFs() const
      {
        return m_dofs;
      }

      /** Evaluates selected coefficients and returns only locally owned entries. */
      template <class Function>
      void assemble(IndexMap<Scalar>& values, const Function& function) const
      {
        values.clear();
        for (const auto& [global, indices] : m_dofs)
        {
          const auto [entity, ordinal] = indices;
          const auto& element = m_fes.getFiniteElement(m_dimension, entity);
          const auto mapping = m_fes.getPullback({m_dimension, entity}, function);
          values.emplace(global, element.getLinearForm(ordinal)(mapping));
        }
        if constexpr (Global)
        {
          if (!m_source)
            return;
          std::vector<std::pair<Index, Scalar>> entries(values.begin(), values.end());
          boost::mpi::broadcast(
            m_fes.getMesh().getContext().getCommunicator(), entries, *m_source);
          Index begin, end;
          m_fes.getOwnershipRange(begin, end);
          values.clear();
          for (const auto& [global, value] : entries)
            if (begin <= global && global < end)
              values.emplace(global, value);
        }
      }

    private:
      const FES& m_fes;
      size_t m_dimension;
      IndexMap<std::pair<Index, Index>> m_dofs;
      Optional<Rank> m_source;
  };

  /** Deduces the space that defines the interpolation functionals. */
  template <class FES, class Pred>
    requires std::is_same_v<typename FES::ContextType, Context::MPI>
  Interpolation(const FES&, Geometry::Region, const Pred&) -> Interpolation<FES>;
}

#endif
