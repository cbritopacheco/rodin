/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Connectivity.h"

namespace Rodin::Geometry
{
  Connectivity<Context::Sequential>::Connectivity()
  {
    m_count.resize(1, 0);
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::initialize(size_t maximalDimension)
  {
    m_maximalDimension = maximalDimension;

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

    return *this;
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::reserve(size_t d, size_t count)
  {
    assert(d < m_connectivity.size());
    m_index[d].left.rehash(count); m_index[d].right.rehash(count);
    m_geometry[d].reserve(count);
    m_connectivity[d][0].reserve(count);
    return *this;
  }

  Connectivity<Context::Sequential>& Connectivity<Context::Sequential>::nodes(size_t count)
  {
    m_count[0] = count;
    m_gcount[Geometry::Polytope::Type::Point] = count;
    for (size_t i = 0; i < count; i++)
    {
      auto p = m_index[0].left.insert({ IndexArray{{ i }}, i });
      assert(p.second);
    }
    return *this;
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::polytope(Polytope::Type t, const Array<Index>& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);
    auto [it, inserted] = m_index[d].left.insert({ in, m_count[d] });
    if (inserted)
    {
      m_connectivity[d][0].emplace_back().insert_unique(it->first.begin(), it->first.end());
      m_geometry[d].push_back(t);
      m_count[d] += 1;
      m_gcount[t] += 1;
      m_dirty[d][0] = false;
    }
    return *this;
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::polytope(Polytope::Type t, Array<Index>&& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);
    auto [it, inserted] = m_index[d].left.insert({ std::move(in), m_count[d] });
    if (inserted)
    {
      m_connectivity[d][0].emplace_back().insert_unique(it->first.begin(), it->first.end());
      m_geometry[d].push_back(t);
      m_count[d] += 1;
      m_gcount[t] += 1;
      m_dirty[d][0] = false;
    }
    return *this;
  }

  const Connectivity<Context::Sequential>::PolytopeIndex&
  Connectivity<Context::Sequential>::getIndexMap(size_t dim) const
  {
    return m_index[dim];
  }

  const std::optional<Index>
  Connectivity<Context::Sequential>::getIndex(size_t dim, const IndexArray& key) const
  {
    auto it = m_index[dim].left.find(key);
    if (it == m_index[dim].left.end())
      return {};
    else
      return it->second;
  }

  const Incidence& Connectivity<Context::Sequential>::getIncidence(size_t d, size_t dp) const
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    return m_connectivity[d][dp];
  }

  const IndexSet& Connectivity<Context::Sequential>::getIncidence(
      const std::pair<size_t, size_t> p, Index idx) const
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    assert(idx < m_connectivity[d][dp].size());
    return m_connectivity[d][dp][idx];
  }

  size_t Connectivity<Context::Sequential>::getCount(size_t dim) const
  {
    return m_count[dim];
  }

  size_t Connectivity<Context::Sequential>::getCount(Polytope::Type g) const
  {
    return m_gcount[g];
  }

  size_t Connectivity<Context::Sequential>::getMeshDimension() const
  {
    for (int i = m_count.size() - 1; i >= 0; i--)
    {
      if (m_count[i] > 0)
        return i;
    }
    return 0;
  }

  Polytope::Type Connectivity<Context::Sequential>::getGeometry(size_t d, Index idx) const
  {
    if (d == 0)
      return Polytope::Type::Point;
    else
      return m_geometry[d][idx];
  }

  const Array<Index>& Connectivity<Context::Sequential>::getPolytope(size_t d, Index idx) const
  {
    return m_index[d].right.at(idx);
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::setIncidence(const std::pair<size_t, size_t>& p, Incidence&& inc)
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_connectivity[d][dp] = std::move(inc);
    return *this;
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::compute(size_t d, size_t dp)
  {
    const size_t D = getMeshDimension();
    if (d == D && dp == 0)
      return *this;
    if (m_dirty[D][D])
      transpose(0, D).intersection(D, D, 0);
    assert(!m_dirty[D][D]);
    if (d != D && d != 0 && (m_dirty[D][d] || m_dirty[d][0]))
      build(d);
    assert(!m_dirty[D][d]);
    assert(!m_dirty[d][0] || d == D || d == 0);
    if (dp != D && dp != 0 && (m_dirty[D][dp] || m_dirty[dp][0]))
      build(dp);
    assert(!m_dirty[D][dp]);
    assert(!m_dirty[dp][0] || dp == D || dp == 0);
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
    m_dirty[d][dp] = false;
    return *this;
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::build(size_t d)
  {
    const size_t D = getMeshDimension();
    assert(d > 0);
    assert(d < D);
    assert(!m_dirty[D][0]);
    assert(!m_dirty[D][D]);
    for (Index i = 0; i < m_count[D]; i++)
    {
      IndexSet s;
      std::vector<SubPolytope> subpolytopes;
      local(subpolytopes, d, i);
      for (auto& [geometry, vertices] : subpolytopes)
      {
        auto insert = m_index[d].left.insert({ std::move(vertices), m_count[d] });
        const PolytopeIndex::left_iterator it = insert.first;
        const bool inserted = insert.second;
        const auto& [v, idx] = *it;
        if (inserted)
        {
          m_geometry[d].push_back(geometry);
          m_connectivity[d][0].emplace_back().insert_unique(v.begin(), v.end());
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

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::transpose(size_t d, size_t dp)
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

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::intersection(size_t d, size_t dp, size_t dpp)
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
          assert(i < m_connectivity[d][0].size());
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

  void Connectivity<Context::Sequential>::local(
      std::vector<Connectivity<Context::Sequential>::SubPolytope>& out, size_t dim, Index i)
  {
    const size_t D = getMeshDimension();

    const auto& p = m_index[D].right.at(i);
    switch (m_geometry[D][i])
    {
      case Polytope::Type::Point:
      {
        assert(dim == 0);
        assert(p.size()  == 0);
        out.resize(1);
        out[0].geometry = Polytope::Type::Point;
        out[0].vertices.resize(1);
        out[0].vertices.coeffRef(0) = i;
        return;
      }
      case Polytope::Type::Segment:
      {
        assert(dim <= 1);
        assert(p.size() == 2);
        if (dim == 0)
        {
          out.resize(2);
          out[0].geometry = Polytope::Type::Point;
          out[0].vertices.resize(1);
          out[0].vertices.coeffRef(0) = p.coeff(0);
          out[1].geometry = Polytope::Type::Point;
          out[1].vertices.resize(1);
          out[1].vertices.coeffRef(0) = p.coeff(1);
        }
        else if (dim == 1)
        {
          out.resize(1);
          out[0].geometry = Polytope::Type::Segment;
          out[0].vertices = p;
        }
        else
        {
          assert(false);
          out = {};
        }
        return;
      }
      case Polytope::Type::Triangle:
      {
        assert(dim <= 2);
        assert(p.size() == 3);
        if (dim == 0)
        {
          out.resize(3);
          out[0].geometry = Polytope::Type::Point;
          out[0].vertices.resize(1);
          out[0].vertices.coeffRef(0) = p.coeff(0);

          out[1].geometry = Polytope::Type::Point;
          out[1].vertices.resize(1);
          out[1].vertices.coeffRef(0) = p.coeff(1);

          out[2].geometry = Polytope::Type::Point;
          out[2].vertices.resize(1);
          out[2].vertices.coeffRef(0) = p.coeff(2);
        }
        else if (dim == 1)
        {
          out.resize(3);

          out[0].geometry = Polytope::Type::Segment;
          out[0].vertices.resize(2);
          out[0].vertices.coeffRef(0) = p.coeff(0);
          out[0].vertices.coeffRef(1) = p.coeff(1);

          out[1].geometry = Polytope::Type::Segment;
          out[1].vertices.resize(2);
          out[1].vertices.coeffRef(0) = p.coeff(1);
          out[1].vertices.coeffRef(1) = p.coeff(2);

          out[2].geometry = Polytope::Type::Segment;
          out[2].vertices.resize(2);
          out[2].vertices.coeffRef(0) = p.coeff(0);
          out[2].vertices.coeffRef(1) = p.coeff(2);
        }
        else if (dim == 2)
        {
          out.resize(1);
          out[0].geometry = Polytope::Type::Triangle;
          out[0].vertices = p;
        }
        else
        {
          assert(false);
          out = {};
        }
        return;
      }
      case Polytope::Type::Quadrilateral:
      {
        assert(dim <= 2);
        assert(p.size() == 4);
        if (dim == 0)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Point, { { p(0) } } };
          out[1] = { Polytope::Type::Point, { { p(1) } } };
          out[2] = { Polytope::Type::Point, { { p(2) } } };
          out[3] = { Polytope::Type::Point, { { p(3) } } };
        }
        else if (dim == 1)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Segment, { { p(0), p(1) } } };
          out[1] = { Polytope::Type::Segment, { { p(1), p(3) } } };
          out[2] = { Polytope::Type::Segment, { { p(3), p(2) } } };
          out[3] = { Polytope::Type::Segment, { { p(2), p(0) } } };
        }
        else if (dim == 2)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Quadrilateral, p };
        }
        else
        {
          assert(false);
          out = {};
        }
        return;
      }
      case Polytope::Type::Tetrahedron:
      {
        assert(dim <= 3);
        assert(p.size() == 4);
        if (dim == 0)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Point, { { p(0) } } };
          out[1] = { Polytope::Type::Point, { { p(1) } } };
          out[2] = { Polytope::Type::Point, { { p(2) } } };
          out[3] = { Polytope::Type::Point, { { p(3) } } };
        }
        else if (dim == 1)
        {
          out.resize(6);
          out[0] = { Polytope::Type::Segment, { { p(0), p(1) } } };
          out[1] = { Polytope::Type::Segment, { { p(1), p(2) } } };
          out[2] = { Polytope::Type::Segment, { { p(0), p(2) } } };
          out[3] = { Polytope::Type::Segment, { { p(0), p(3) } } };
          out[4] = { Polytope::Type::Segment, { { p(1), p(3) } } };
          out[5] = { Polytope::Type::Segment, { { p(2), p(3) } } };
        }
        else if (dim == 2)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Triangle, { { p(0), p(1), p(2) } } };
          out[1] = { Polytope::Type::Triangle, { { p(0), p(1), p(3) } } };
          out[2] = { Polytope::Type::Triangle, { { p(0), p(2), p(3) } } };
          out[3] = { Polytope::Type::Triangle, { { p(1), p(2), p(3) } } };
        }
        else if (dim == 3)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Tetrahedron, p };
        }
        else
        {
          assert(false);
          out = {};
        }
        return;
      }
    }
    assert(false); // We should not reach here. There is an unhandled case.
    out = {};
  }

  Connectivity<Context::Sequential>&
  Connectivity<Context::Sequential>::clear(size_t d, size_t dp)
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_dirty[d][dp] = true;
    m_connectivity[d][dp].clear();
    return *this;
  }
}
