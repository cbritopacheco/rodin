/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_CONNECTIVITY_H
#define RODIN_GEOMETRY_CONNECTIVITY_H

#include <vector>

#include "Rodin/Array.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents the set of incidence relations.
   *
   * Stores the set of incidence relations:
   * @f[
   *  d \rightarrow d'
   * @f]
   * for a fixed pair of topological dimensions @f$ (d, d') @f$.
   */
  class Connectivity
  {
    public:
      using Incidence = std::vector<Array<Index>>;

      Connectivity() = default;

      Connectivity(size_t meshDim)
        : m_meshDimension(meshDim)
      {
        initialize(meshDim);
      }

      Connectivity(const Connectivity&) = default;

      Connectivity(Connectivity&&) = default;

      Connectivity& operator=(const Connectivity&) = default;

      Connectivity& operator=(Connectivity&&) = default;

      inline
      Connectivity& initialize(size_t meshDim)
      {
        m_connectivity.resize(meshDim + 1);
        for (auto& v : m_connectivity)
          v.resize(meshDim + 1);

        m_dirty.resize(meshDim + 1);
        for (auto& v : m_dirty)
          v.resize(meshDim + 1, true);

        return *this;
      }

      inline
      Connectivity& reserve(size_t d, size_t count)
      {
        assert(d < m_connectivity.size());
        for (auto& v : m_connectivity[d])
          v.reserve(count);
        return *this;
      }

      inline
      Connectivity& reserve(const std::pair<size_t, size_t> p, size_t count)
      {
        const auto& [d, dp] = p;
        assert(d < m_connectivity.size());
        assert(dp < m_connectivity[d].size());
        m_connectivity[d][dp].reserve(count);
        return *this;
      }

      inline
      Connectivity& connect(const std::pair<size_t, size_t> p, std::initializer_list<Index> in)
      {
        Array<Index> a(in.size());
        std::copy(in.begin(), in.end(), a.begin());
        return connect(p, std::move(a));
      }

      inline
      Connectivity& connect(const std::pair<size_t, size_t> p, const Array<Index>& in)
      {
        const auto& [d, dp] = p;
        assert(d < m_connectivity.size());
        assert(dp < m_connectivity[d].size());
        m_dirty[d][dp] = true;
        m_connectivity[d][dp].push_back(in);
        return *this;
      }

      inline
      Connectivity& connect(const std::pair<size_t, size_t> p, Array<Index>&& in)
      {
        const auto& [d, dp] = p;
        assert(d < m_connectivity.size());
        assert(dp < m_connectivity[d].size());
        m_dirty[d][dp] = true;
        m_connectivity[d][dp].push_back(in);
        return *this;
      }

      inline
      const Incidence& getIncidence(size_t d, size_t dp) const
      {
        assert(d < m_connectivity.size());
        assert(dp < m_connectivity[d].size());
        return m_connectivity[d][dp];
      }

      inline
      const Array<Index>& getIncidence(const std::pair<size_t, size_t> p, Index idx) const
      {
        const auto& [d, dp] = p;
        assert(d < m_connectivity.size());
        assert(dp < m_connectivity[d].size());
        assert(idx < m_connectivity[d][dp].size());
        return m_connectivity[d][dp][idx];
      }

      inline
      size_t getMeshDimension() const
      {
        return m_meshDimension;
      }

      void compute(size_t d, size_t dp)
      {
      }

      void build(size_t d)
      {
        assert(false);
      }

      void transpose(size_t d, size_t dp)
      {
        assert(false);
      }

      void intersection(size_t d, size_t dp)
      {
        assert(false);
      }

    private:
      size_t m_meshDimension;
      std::vector<std::vector<Incidence>> m_connectivity;
      std::vector<std::vector<bool>> m_dirty;
  };
}

#endif
