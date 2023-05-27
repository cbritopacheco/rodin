/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Connectivity.h"

namespace Rodin::Geometry
{
  MeshConnectivity& MeshConnectivity::initialize(size_t maximalDimension)
  {
    m_maximalDimension = maximalDimension;

    assert(m_count.size() == 0);
    m_count.resize(maximalDimension + 1, 0);

    assert(m_connectivity.size() == 0);
    m_connectivity.resize(maximalDimension + 1);
    for (auto& v : m_connectivity)
      v.resize(maximalDimension + 1);

    assert(m_dirty.size() == 0);
    m_dirty.resize(maximalDimension + 1);
    for (auto& v : m_dirty)
      v.resize(maximalDimension + 1, true);

    assert(m_index.size() == 0);
    m_index.resize(maximalDimension + 1);

    assert(m_geometry.size() == 0);
    m_geometry.resize(maximalDimension + 1);

    assert(m_polytopes.size() == 0);
    m_polytopes.resize(maximalDimension + 1);

    m_gcount.reserve(Polytope::Geometries.size());
    for (const auto& g : Polytope::Geometries)
      m_gcount.insert({g, 0});

    return *this;
  }

  MeshConnectivity& MeshConnectivity::reserve(size_t d, size_t count)
  {
    assert(d < m_connectivity.size());
    for (auto& v : m_connectivity[d])
      v.reserve(count);
    m_index[d].reserve(count);
    m_geometry[d].reserve(count);
    m_polytopes[d].reserve(count);
    return *this;
  }

  MeshConnectivity& MeshConnectivity::nodes(size_t count)
  {
    m_count[0] = count;
    m_gcount[Geometry::Polytope::Geometry::Point] = count;
    return *this;
  }

  MeshConnectivity&
  MeshConnectivity::polytope(Polytope::Geometry t, const Array<Index>& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);
    m_connectivity[d][0].emplace_back(in.begin(), in.end());
    m_polytopes[d].push_back(in);
    m_geometry[d].push_back(t);
    m_dirty[d][0] = true;
    m_count[d] += 1;
    m_gcount[t] += 1;
    return *this;
  }

  MeshConnectivity&
  MeshConnectivity::polytope(Polytope::Geometry t, Array<Index>&& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);
    m_connectivity[d][0].emplace_back(in.begin(), in.end());
    m_polytopes[d].push_back(std::move(in));
    m_geometry[d].push_back(t);
    m_dirty[d][0] = true;
    m_count[d] += 1;
    m_gcount[t] += 1;
    return *this;
  }

  const MeshConnectivity::IndexMap& MeshConnectivity::getIndexMap(size_t dim) const
  {
    return m_index[dim];
  }

  const std::optional<Index> MeshConnectivity::getIndex(size_t dim, const IndexSet& key) const
  {
    auto it = m_index[dim].find(key);
    if (it == m_index[dim].end())
      return {};
    else
      return it->second;
  }

  const Incidence& MeshConnectivity::getIncidence(size_t d, size_t dp) const
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    return m_connectivity[d][dp];
  }

  const FlatSet<Index>& MeshConnectivity::getIncidence(
      const std::pair<size_t, size_t> p, Index idx) const
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    assert(idx < m_connectivity[d][dp].size());
    return m_connectivity[d][dp][idx];
  }

  size_t MeshConnectivity::getCount(size_t dim) const
  {
    return m_count[dim];
  }

  size_t MeshConnectivity::getCount(Polytope::Geometry g) const
  {
    return m_gcount.at(g);
  }

  size_t MeshConnectivity::getMeshDimension() const
  {
    for (size_t i = m_count.size() - 1; i >= 0; i--)
    {
      if (m_count[i] > 0)
        return i;
    }
    return 0;
  }

  Polytope::Geometry MeshConnectivity::getGeometry(size_t d, Index idx) const
  {
    if (d == 0)
      return Polytope::Geometry::Point;
    else
      return m_geometry[d][idx];
  }

  const Array<Index>& MeshConnectivity::getPolytope(size_t d, Index idx) const
  {
    return m_polytopes[d][idx];
  }

  MeshConnectivity& MeshConnectivity::compute(size_t d, size_t dp)
  {
    const size_t D = getMeshDimension();
    if (m_dirty[D][D])
      transpose(0, D).intersection(D, D, 0);
    assert(!m_dirty[D][D]);
    if (!(d == D && dp == 0))
    {
      if (m_dirty[d][0])
        build(d);
      if (m_dirty[dp][0])
        build(dp);
      if (m_dirty[d][dp])
      {
        if (d < dp)
        {
          compute(dp, d).transpose(d, dp);
        }
        else
        {
          size_t dpp;
          if (d == 0 && dp == 0)
            dpp = D;
          else
            dpp = 0;
          compute(d, dpp).compute(dpp, dp).intersection(d, dp, dpp);
        }
      }
    }
    m_dirty[d][dp] = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::build(size_t d)
  {
    for (size_t i = 0; i < m_polytopes[d].size(); i++)
    {
      const auto& vertices = m_polytopes[d][i];
      m_index[d].insert({ IndexSet(vertices.begin(), vertices.end()), i });
    }

    const size_t D = getMeshDimension();
    for (Index i = 0; i < m_count[D]; i++)
    {
      IndexSet s;
      for (auto [geometry, vertices] : local(d, i))
      {
        auto [it, inserted] =
          m_index[d].insert({ IndexSet(vertices.begin(), vertices.end()), m_count[d] });
        const auto& [v, idx] = *it;
        if (inserted)
        {
          m_polytopes[d].push_back(std::move(vertices));
          m_geometry[d].push_back(geometry);
          m_connectivity[d][0].push_back(v);
        }
        m_count[d] += inserted && !(d == D || d == 0);
        m_gcount[geometry] += inserted && !(d == D || d == 0);
        s.insert(idx);
      }
      m_connectivity[D][d].push_back(std::move(s));
    }
    m_dirty[D][d] = false;
    m_dirty[d][0] = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::transpose(size_t d, size_t dp)
  {
    assert(d < dp);
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_connectivity[d][dp].resize(m_count[d]);
    for (Index j = 0; j < m_count[dp]; j++)
    {
      for (Index i : m_connectivity[dp][d][j])
        m_connectivity[d][dp][i].insert(j);
    }
    m_dirty[d][dp] = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::intersection(size_t d, size_t dp, size_t dpp)
  {
    assert(d >= dp);
    m_connectivity[d][dp].resize(m_count[d]);
    for (Index i = 0; i < m_count[d]; i++)
    {
      assert(i < m_connectivity[d][dpp].size());
      for (Index k : m_connectivity[d][dpp][i])
      {
        assert(k < m_connectivity[dpp][dp].size());
        for (Index j : m_connectivity[dpp][dp][k])
        {
          assert(i < m_connectivity[dp][0].size());
          assert(j < m_connectivity[dp][0].size());
          const auto& d0i = m_connectivity[d][0][i];
          const auto& d0j = m_connectivity[dp][0][j];
          if ((d == dp && i != j) ||
              (d > dp && std::includes(d0i.begin(), d0i.end(), d0j.begin(), d0j.end())))
          {
            m_connectivity[d][dp][i].insert(j);
          }
        }
      }
    }
    m_dirty[d][dp] = false;
    return *this;
  }

  std::vector<MeshConnectivity::SubPolytope> MeshConnectivity::local(size_t dim, Index i)
  {
    const size_t D = getMeshDimension();
    const auto& p = m_polytopes[D][i];
    switch (m_geometry[D][i])
    {
      case Polytope::Geometry::Point:
      {
        assert(dim == 0);
        assert(p.size()  == 0);
        return { {Polytope::Geometry::Point, { { p(0) } } } };
      }
      case Polytope::Geometry::Segment:
      {
        assert(dim <= 1);
        assert(p.size() == 2);
        if (dim == 0)
        {
          return {
            { Polytope::Geometry::Point, { { p(0) } } },
            { Polytope::Geometry::Point, { { p(1) } } }
          };
        }
        else if (dim == 1)
        {
          return { { Polytope::Geometry::Segment, { { p(0), p(1) } } } };
        }
        else
        {
          assert(false);
          return {};
        }
        break;
      }
      case Polytope::Geometry::Triangle:
      {
        assert(dim <= 2);
        assert(p.size() == 3);
        if (dim == 0)
        {
          return {
            { Polytope::Geometry::Point, { { p(0) } } },
            { Polytope::Geometry::Point, { { p(1) } } },
            { Polytope::Geometry::Point, { { p(2) } } }
          };
        }
        else if (dim == 1)
        {
          return {
            { Polytope::Geometry::Segment, { { p(0), p(1) } } },
            { Polytope::Geometry::Segment, { { p(1), p(2) } } },
            { Polytope::Geometry::Segment, { { p(0), p(2) } } }
          };
        }
        else if (dim == 2)
        {
          return { { Polytope::Geometry::Triangle, { { p(0), p(1), p(2) } } } };
        }
        else
        {
          assert(false);
          return {};
        }
        break;
      }
      case Polytope::Geometry::Quadrilateral:
      {
        assert(dim <= 2);
        assert(p.size() == 4);
        if (dim == 0)
        {
          return {
            { Polytope::Geometry::Point, { { p(0) } } },
            { Polytope::Geometry::Point, { { p(1) } } },
            { Polytope::Geometry::Point, { { p(2) } } },
            { Polytope::Geometry::Point, { { p(3) } } }
          };
        }
        else if (dim == 1)
        {
          return {
            { Polytope::Geometry::Segment, { { p(0), p(1) } } },
            { Polytope::Geometry::Segment, { { p(1), p(3) } } },
            { Polytope::Geometry::Segment, { { p(3), p(2) } } },
            { Polytope::Geometry::Segment, { { p(2), p(0) } } }
          };
        }
        else if (dim == 2)
        {
          return { { Polytope::Geometry::Quadrilateral, { { p(0), p(1), p(2), p(3) } } } };
        }
        else
        {
          assert(false);
          return {};
        }
        break;
      }
      case Polytope::Geometry::Tetrahedron:
      {
        assert(dim <= 3);
        assert(p.size() == 4);
        if (dim == 0)
        {
          return {
            { Polytope::Geometry::Point, { { p(0) } } },
            { Polytope::Geometry::Point, { { p(1) } } },
            { Polytope::Geometry::Point, { { p(2) } } },
            { Polytope::Geometry::Point, { { p(3) } } }
          };
        }
        else if (dim == 1)
        {
          return {
            { Polytope::Geometry::Segment, { { p(0), p(1) } } },
            { Polytope::Geometry::Segment, { { p(1), p(2) } } },
            { Polytope::Geometry::Segment, { { p(0), p(2) } } },
            { Polytope::Geometry::Segment, { { p(0), p(3) } } },
            { Polytope::Geometry::Segment, { { p(1), p(3) } } },
            { Polytope::Geometry::Segment, { { p(2), p(3) } } }
          };
        }
        else if (dim == 2)
        {
          return {
            { Polytope::Geometry::Triangle, { { p(0), p(1), p(2) } } },
            { Polytope::Geometry::Triangle, { { p(0), p(1), p(3) } } },
            { Polytope::Geometry::Triangle, { { p(0), p(2), p(3) } } },
            { Polytope::Geometry::Triangle, { { p(1), p(2), p(3) } } }
          };
        }
        else if (dim == 3)
        {
          return { { Polytope::Geometry::Tetrahedron, { { p(0), p(1), p(2), p(3) } } } };
        }
        else
        {
          assert(false);
          return {};
        }
        break;
      }
    }
    assert(false); // We should not reach here. There is an unhandled case.
    return {};
  }
}
