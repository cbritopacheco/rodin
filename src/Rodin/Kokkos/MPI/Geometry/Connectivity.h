/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MPI_GEOMETRY_CONNECTIVITY_H
#define RODIN_KOKKOS_MPI_GEOMETRY_CONNECTIVITY_H

/**
 * @file
 * @brief Incidence preparation for resident MPI shards with fixed numbering.
 */

#include "Rodin/Kokkos/Geometry/Connectivity.h"
#include "Rodin/MPI/Context/MPI.h"

namespace Rodin::Geometry
{
  /**
   * @brief Borrows a distributed shard's incidence relations.
   *
   * Only the owning mesh can replace the numbered entities and their MPI
   * metadata. Incidence preparation can use existing entities on this rank.
   *
   * @tparam Backend Kokkos execution backend.
   */
  template <class Backend>
  class Connectivity<Context::Kokkos<Context::MPI, Backend>> final
    : public ConnectivityBase
  {
    public:
      /// @brief Owner of the resident shard's connectivity arrays.
      using LocalType = Connectivity<Context::Kokkos<Context::Local, Backend>>;

      /// @brief Complete or restricted incidence computation.
      using Mode = typename LocalType::Mode;

      /// @param local Shard connectivity whose owner must outlive this object.
      explicit Connectivity(LocalType& local)
        : m_local(local)
      {}

      /// @param other Borrow tied to another mesh.
      Connectivity(const Connectivity& other) = delete;

      /// @param other Borrow tied to another mesh. @returns This object.
      Connectivity& operator=(const Connectivity& other) = delete;

      /**
       * @param d Source dimension.
       * @param dp Target dimension.
       * @param mode Discovery policy when rebuilding a prerequisite relation.
       * @returns This connectivity.
       * @throws Alert::Exception If discovery would add new entities.
       */
      Connectivity& compute(size_t d, size_t dp, Mode mode = Mode::Discover)
      {
        if (mode == Mode::Restrict)
          m_local.restrict(d, dp);
        else
          m_local.computeExisting(d, dp);
        return *this;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This object after checked discovery.
      Connectivity& discover(size_t d, size_t dp)
      {
        return compute(d, dp, Mode::Discover);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This object after computing incidence among existing entities.
      Connectivity& restrict(size_t d, size_t dp)
      {
        return compute(d, dp, Mode::Restrict);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This object after clearing the relation.
      Connectivity& clear(size_t d, size_t dp)
      {
        m_local.clear(d, dp);
        return *this;
      }

      /// @param d Source dimension. @param dp Target dimension. @param dirty New dirty state. @returns This object.
      Connectivity& setDirty(size_t d, size_t dp, bool dirty = true)
      {
        m_local.setDirty(d, dp, dirty);
        return *this;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Whether relation storage is present.
      bool hasIncidence(size_t d, size_t dp) const
      {
        return m_local.hasIncidence(d, dp);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Whether the relation is present and current.
      bool isPrepared(size_t d, size_t dp) const
      {
        return m_local.isPrepared(d, dp);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Whether recomputation is required.
      bool isDirty(size_t d, size_t dp) const
      {
        return m_local.isDirty(d, dp);
      }

      /// @returns The resident connectivity revision.
      std::uint64_t getRevision() const
      {
        return m_local.getRevision();
      }

      /// @returns Borrowed device tables, invalidated by incidence mutation.
      typename LocalType::View getView() const
      {
        return m_local.getView();
      }

      Polytope::Type getGeometry(size_t d, Index i) const override
      {
        return m_local.getGeometry(d, i);
      }

      const Polytope::Key& getPolytope(size_t d, Index i) const override
      {
        return m_local.getPolytope(d, i);
      }

      size_t getCount(size_t d) const override
      {
        return m_local.getCount(d);
      }

      size_t getCount(Polytope::Type g) const override
      {
        return m_local.getCount(g);
      }

      size_t getDimension() const override
      {
        return m_local.getDimension();
      }

      const Incidence& getIncidence(size_t d, size_t dp) const override
      {
        return m_local.getIncidence(d, dp);
      }

      const IndexVector& getIncidence(
        std::pair<size_t, size_t> pair, Index i) const override
      {
        return m_local.getIncidence(pair, i);
      }

    private:
      LocalType& m_local;
  };
}

#endif
