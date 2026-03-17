/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Connectivity.h"
#include "Rodin/Configure.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Geometry
{
  Connectivity<Context::Local>::Connectivity()
    : m_maximalDimension(0), m_dimension(0)
  {
    m_count.resize(1, 0);
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::initialize(size_t maximalDimension)
  {
    m_maximalDimension = maximalDimension;
    m_dimension = 0;

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

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::reserve(size_t d, size_t count)
  {
    assert(d < m_connectivity.size());
    m_index[d].left.reserve(count);
    m_index[d].right.reserve(count);
    m_geometry[d].reserve(count);
    m_connectivity[d][0].reserve(count);
    return *this;
  }

  Connectivity<Context::Local>& Connectivity<Context::Local>::nodes(size_t count)
  {
    m_count[0] = count;
    m_index[0].left.reserve(count);
    m_index[0].right.reserve(count);
    m_gcount[Geometry::Polytope::Type::Point] = count;
    for (size_t i = 0; i < count; i++)
    {
      const auto p = m_index[0].right.insert({ Polytope::Key( {i }), i });
      assert(p.second);
      m_index[0].left.push_back(&p.first->first);
    }
    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::polytope(Polytope::Type t, const Polytope::Key& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::Traits(t).getDimension();
    assert(d > 0);
    assert(d <= m_maximalDimension);
    const auto [it, inserted] = m_index[d].right.insert({ in, m_count[d]});
    if (inserted)
    {
      m_index[d].left.push_back(&it->first);
      auto& v = m_connectivity[d][0].emplace_back();
      for (const Index& j : it->first)
        v.push_back(j);
      m_geometry[d].push_back(t);
      m_count[d] += 1;
      m_gcount[t] += 1;
      m_dirty[d][0] = false;
      m_dimension = std::max(m_dimension, d);
    }
    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::polytope(Polytope::Type t, Polytope::Key&& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::Traits(t).getDimension();
    assert(d > 0);
    assert(d <= m_maximalDimension);
    const auto [it, inserted] = m_index[d].right.insert({ std::move(in), m_count[d] });
    if (inserted)
    {
      m_index[d].left.push_back(&it->first);
      auto& v = m_connectivity[d][0].emplace_back();
      for (const Index& j : it->first)
        v.push_back(j);
      m_geometry[d].push_back(t);
      m_count[d] += 1;
      m_gcount[t] += 1;
      m_dirty[d][0] = false;
      m_dimension = std::max(m_dimension, d);
    }
    return *this;
  }

  const Connectivity<Context::Local>::PolytopeIndex&
  Connectivity<Context::Local>::getIndexMap(size_t dim) const
  {
    return m_index[dim];
  }

  const Optional<Index>
  Connectivity<Context::Local>::getIndex(size_t dim, const Polytope::Key& key) const
  {
    const auto it = m_index[dim].right.find(key);
    if (it == m_index[dim].right.end())
      return {};
    else
      return it->second;
  }

  const Incidence& Connectivity<Context::Local>::getIncidence(size_t d, size_t dp) const
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    return m_connectivity[d][dp];
  }

  const IndexVector& Connectivity<Context::Local>::getIncidence(
      const std::pair<size_t, size_t> p, Index idx) const
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    assert(idx < m_connectivity[d][dp].size());
    return m_connectivity[d][dp][idx];
  }

  size_t Connectivity<Context::Local>::getCount(size_t dim) const
  {
    return m_count[dim];
  }

  size_t Connectivity<Context::Local>::getCount(Polytope::Type g) const
  {
    return m_gcount[g];
  }

  size_t Connectivity<Context::Local>::getDimension() const
  {
    return m_dimension;
  }

  Polytope::Type Connectivity<Context::Local>::getGeometry(size_t d, Index idx) const
  {
    if (d == 0)
      return Polytope::Type::Point;
    else
      return m_geometry[d][idx];
  }

  const Polytope::Key& Connectivity<Context::Local>::getPolytope(size_t d, Index idx) const
  {
    const Polytope::Key* p = m_index[d].left[idx];
    assert(p);
    return *p;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::setIncidence(const std::pair<size_t, size_t>& p, Incidence&& inc)
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_connectivity[d][dp] = std::move(inc);
    // This incidence relation is now explicitly provided (e.g. from serialized IO),
    // so it must not be considered stale and recomputed on next access.
    m_dirty[d][dp] = false;
    return *this;
  }

  Connectivity<Context::Local>& Connectivity<Context::Local>::discover(size_t d, size_t dp)
  {
    return this->compute(d, dp, Mode::Discover);
  }

  Connectivity<Context::Local>& Connectivity<Context::Local>::restrict(size_t d, size_t dp)
  {
    return this->compute(d, dp, Mode::Restrict);
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::compute(size_t d, size_t dp, Mode mode)
  {
    if (getCount(0) == 0)
      return *this;
    const size_t D = getDimension();
    if (d == D && dp == 0)
      return *this;
    if (m_dirty[D][D])
      transpose(0, D).intersection(D, D, 0);
    assert(!m_dirty[D][D]);
    if (d != D && d != 0 && (m_dirty[D][d] || m_dirty[d][0]))
      build(d, mode);
    assert(!m_dirty[D][d]);
    assert(!m_dirty[d][0] || d == D || d == 0);
    if (dp != D && dp != 0 && (m_dirty[D][dp] || m_dirty[dp][0]))
      build(dp, mode);
    assert(!m_dirty[D][dp]);
    assert(!m_dirty[dp][0] || dp == D || dp == 0);
    if (m_dirty[d][dp])
    {
      if (d < dp)
      {
        compute(dp, d, mode).transpose(d, dp);
      }
      else
      {
        size_t dpp;
        if (d == 0 && dp == 0)
          dpp = D;
        else
          dpp = 0;
        compute(d, dpp).compute(dpp, dp, mode).intersection(d, dp, dpp);
      }
    }
    m_dirty[d][dp] = false;
    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::build(size_t d, Mode mode)
  {
    const size_t D = getDimension();
    assert(d > 0);
    assert(d < D);
    assert(!m_dirty[D][0]);
    assert(!m_dirty[D][D]);

    m_connectivity[D][d].clear();
    m_connectivity[D][d].reserve(m_count[D]);

    for (Index i = 0; i < m_count[D]; ++i)
      local(i, d, mode);

    m_dirty[D][d] = false;

    switch (mode)
    {
      case Mode::Discover:
      {
        m_dirty[d][0] = false;
        break;
      }
      case Mode::Restrict:
      {
        assert(!m_dirty[d][0]); // precondition for restricted mode
        break;
      }
    }

    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::local(size_t i, size_t d, Mode mode)
  {
    static thread_local std::vector<SubPolytope> subpolytopes;
    std::vector<Index> s;

    const size_t D = getDimension();
    assert(d > 0);
    assert(d < D);

    this->getSubPolytopes(subpolytopes, i, d);

    s.reserve(subpolytopes.size());

    for (auto& [geometry, vertices] : subpolytopes)
    {
      Index idx;

      switch (mode)
      {
        case Mode::Discover:
        {
          auto insert = m_index[d].right.insert({ std::move(vertices), m_count[d] });
          const auto it = insert.first;
          const bool inserted = insert.second;
          const auto& arr = it->first;
          idx = it->second;

          if (inserted)
          {
            auto& v = m_connectivity[d][0].emplace_back();
            v.reserve(arr.size());
            for (Index j : arr)
              v.push_back(j);
            m_index[d].left.push_back(&it->first);
            m_geometry[d].push_back(geometry);

            if (d != 0 && d != D)
            {
              ++m_count[d];
              ++m_gcount[geometry];
            }
          }
          break;
        }
        case Mode::Restrict:
        {
          // LOOKUP ONLY
          auto it = m_index[d].right.find(vertices);
          if (it == m_index[d].right.end())
            continue;       // this subpolytope is not in the tracked set
          idx = it->second; // existing entity index
          break;
        }
      }
      s.push_back(idx);
    }

    m_connectivity[D][d].push_back(std::move(s));
    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::transpose(size_t d, size_t dp)
  {
    static thread_local std::vector<size_t> s_mark;
    static thread_local size_t s_epoch = 1;

    assert(d < dp);
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity.size());
    assert(d < m_connectivity[d].size());
    assert(dp < m_connectivity[dp].size());

    const Index n_dst = m_count[d];
    const Index n_src = m_count[dp];

    m_connectivity[d][dp].assign(n_dst, {});

    std::vector<size_t> counts(n_dst, 0);
    for (Index j = 0; j < n_src; ++j)
    {
      for (Index i : m_connectivity[dp][d][j])
        ++counts[i];
    }

    for (Index i = 0; i < n_dst; ++i)
      m_connectivity[d][dp][i].reserve(counts[i]);

    for (Index j = 0; j < n_src; ++j)
    {
      for (Index i : m_connectivity[dp][d][j])
        m_connectivity[d][dp][i].push_back(j);
    }

    if (s_mark.size() < static_cast<size_t>(n_src))
      s_mark.assign(n_src, 0);

    for (Index i = 0; i < n_dst; ++i)
    {
      if (++s_epoch == 0)
      {
        std::fill(s_mark.begin(), s_mark.end(), 0);
        s_epoch = 1;
      }

      auto& v = m_connectivity[d][dp][i];
      size_t write = 0;
      for (size_t r = 0; r < v.size(); ++r)
      {
        const Index j = v[r];
        if (s_mark[j] != s_epoch)
        {
            s_mark[j] = s_epoch;
            v[write++] = j;
        }
      }
      v.resize(write);
    }

    m_dirty[d][dp] = false;
    return *this;
  }

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::intersection(size_t d, size_t dp, size_t dpp)
  {
    static thread_local std::vector<size_t> s_markJ;
    static thread_local size_t s_epochJ = 1;

    static thread_local std::vector<size_t> s_markV;
    static thread_local size_t s_epochV = 1;

    assert(d >= dp);
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity.size());
    assert(dpp < m_connectivity.size());
    assert(d < m_connectivity[d].size());
    assert(dp < m_connectivity[dp].size());
    assert(dpp < m_connectivity[dpp].size());
    assert(0 < m_connectivity[d].size());
    assert(0 < m_connectivity[dp].size());

    const Index n_d  = m_count[d];
    const Index n_dp = m_count[dp];
    const Index n_0  = m_count[0];

    m_connectivity[d][dp].assign(n_d, {});

    if (s_markJ.size() < static_cast<size_t>(n_dp))
      s_markJ.assign(n_dp, 0);

    if (s_markV.size() < static_cast<size_t>(n_0))
      s_markV.assign(n_0, 0);

    for (Index i = 0; i < n_d; ++i)
    {
      auto& out = m_connectivity[d][dp][i];

      ++s_epochJ;
      if (s_epochJ == 0)
      {
        std::fill(s_markJ.begin(), s_markJ.end(), 0);
        s_epochJ = 1;
      }

      ++s_epochV;
      if (s_epochV == 0)
      {
        std::fill(s_markV.begin(), s_markV.end(), 0);
        s_epochV = 1;
      }

      // mark vertices of i
      const auto& d0i = m_connectivity[d][0][i];
      for (Index v : d0i)
        s_markV[v] = s_epochV;

      // traverse via dpp
      const auto& idpp = m_connectivity[d][dpp][i];
      for (Index k : idpp)
      {
        const auto& kdps = m_connectivity[dpp][dp][k];
        for (Index j : kdps)
        {
          if (s_markJ[j] == s_epochJ)
            continue;  // already added for this i
          if (d == dp && i == j)
            continue;  // exclude self

          bool subset = true;
          if (d > dp)
          {
            const auto& d0j = m_connectivity[dp][0][j];
            for (Index v : d0j)
            {
              if (s_markV[v] != s_epochV)
              {
                subset = false;
                break;
              }
            }
          }
          if (d == dp || subset)
          {
            s_markJ[j] = s_epochJ;
            out.push_back(j); // preserve discovery order
          }
        }
      }
    }

    m_dirty[d][dp] = false;
    return *this;
  }

  void Connectivity<Context::Local>::getSubPolytopes(
      std::vector<SubPolytope>& out, Index i, size_t dim) const
  {
    const size_t D = getDimension();
    const auto& p = *m_index[D].left[i];
    switch (m_geometry[D][i])
    {
      case Polytope::Type::Point:
      {
        assert(dim == 0);
        assert(p.size() == 1);
        out.resize(1);
        out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
        return;
      }
      case Polytope::Type::Segment:
      {
        assert(dim <= 1);
        assert(p.size() == 2);
        if (dim == 0)
        {
          out.resize(2);
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
        }
        else if (dim == 1)
        {
          // Local segment orientation: (0 -> 1)
          out.resize(1);
          out[0] = { Polytope::Type::Segment, p };
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
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
          out[2] = { Polytope::Type::Point, Polytope::Key({ p(2) }) };
        }
        else if (dim == 1)
        {
          // Triangle edges in CCW loop: (0->1), (1->2), (2->0)
          out.resize(3);
          out[0] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(1) }) };
          out[1] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(2) }) };
          out[2] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(0) }) };
        }
        else if (dim == 2)
        {
          // Face orientation: (0,1,2) = reference CCW triangle
          out.resize(1);
          out[0] = { Polytope::Type::Triangle, p };
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
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
          out[2] = { Polytope::Type::Point, Polytope::Key({ p(2) }) };
          out[3] = { Polytope::Type::Point, Polytope::Key({ p(3) }) };
        }
        else if (dim == 1)
        {
          // Quad edges in CCW loop: (0->1->2->3->0)
          out.resize(4);
          out[0] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(1) }) };
          out[1] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(2) }) };
          out[2] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(3) }) };
          out[3] = { Polytope::Type::Segment, Polytope::Key({ p(3), p(0) }) };
        }
        else if (dim == 2)
        {
          // Face orientation: (0,1,2,3) = CCW quad
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
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
          out[2] = { Polytope::Type::Point, Polytope::Key({ p(2) }) };
          out[3] = { Polytope::Type::Point, Polytope::Key({ p(3) }) };
        }
        else if (dim == 1)
        {
          // Edge list: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
          out.resize(6);
          out[0] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(1) }) };
          out[1] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(2) }) };
          out[2] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(3) }) };
          out[3] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(2) }) };
          out[4] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(3) }) };
          out[5] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(3) }) };
        }
        else if (dim == 2)
        {
          // Faces oriented to match ∂[0,1,2,3] = [1,2,3] - [0,2,3] + [0,1,3] - [0,1,2]
          out.resize(4);
          out[0] = { Polytope::Type::Triangle, Polytope::Key({ p(1), p(2), p(3) }) }; // +[1,2,3]
          out[1] = { Polytope::Type::Triangle, Polytope::Key({ p(0), p(3), p(2) }) }; // -[0,2,3]
          out[2] = { Polytope::Type::Triangle, Polytope::Key({ p(0), p(1), p(3) }) }; // +[0,1,3]
          out[3] = { Polytope::Type::Triangle, Polytope::Key({ p(0), p(2), p(1) }) }; // -[0,1,2]
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
      case Polytope::Type::Hexahedron:
      {
        assert(dim <= 3);
        assert(p.size() == 8);

        if (dim == 0)
        {
          // 8 vertices
          out.resize(8);
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
          out[2] = { Polytope::Type::Point, Polytope::Key({ p(2) }) };
          out[3] = { Polytope::Type::Point, Polytope::Key({ p(3) }) };
          out[4] = { Polytope::Type::Point, Polytope::Key({ p(4) }) };
          out[5] = { Polytope::Type::Point, Polytope::Key({ p(5) }) };
          out[6] = { Polytope::Type::Point, Polytope::Key({ p(6) }) };
          out[7] = { Polytope::Type::Point, Polytope::Key({ p(7) }) };
        }
        else if (dim == 1)
        {
          // 12 edges:
          // bottom: (0->1),(1->2),(2->3),(3->0)  (CCW in (x,y))
          // top:    (4->5),(5->6),(6->7),(7->4)  (CCW in (x,y))
          // verts:  (0->4),(1->5),(2->6),(3->7)
          out.resize(12);
          out[0]  = { Polytope::Type::Segment, Polytope::Key({ p(0), p(1) }) };
          out[1]  = { Polytope::Type::Segment, Polytope::Key({ p(1), p(2) }) };
          out[2]  = { Polytope::Type::Segment, Polytope::Key({ p(2), p(3) }) };
          out[3]  = { Polytope::Type::Segment, Polytope::Key({ p(3), p(0) }) };

          out[4]  = { Polytope::Type::Segment, Polytope::Key({ p(4), p(5) }) };
          out[5]  = { Polytope::Type::Segment, Polytope::Key({ p(5), p(6) }) };
          out[6]  = { Polytope::Type::Segment, Polytope::Key({ p(6), p(7) }) };
          out[7]  = { Polytope::Type::Segment, Polytope::Key({ p(7), p(4) }) };

          out[8]  = { Polytope::Type::Segment, Polytope::Key({ p(0), p(4) }) };
          out[9]  = { Polytope::Type::Segment, Polytope::Key({ p(1), p(5) }) };
          out[10] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(6) }) };
          out[11] = { Polytope::Type::Segment, Polytope::Key({ p(3), p(7) }) };
        }
        else if (dim == 2)
        {
          // 6 quad faces, analogous to the wedge:
          //  - start with bottom,
          //  - then wrap the 4 side faces following the bottom edges
          //  - end with top.
          //
          // bottom: (0,1,2,3)
          // side 0: (0,1,5,4) (attached to edge 0->1)
          // side 1: (1,2,6,5) (attached to edge 1->2)
          // side 2: (2,3,7,6) (attached to edge 2->3)
          // side 3: (3,0,4,7) (attached to edge 3->0)
          // top:    (4,5,6,7)
          out.resize(6);
          out[0] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(0), p(1), p(2), p(3) }) };
          out[1] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(0), p(1), p(5), p(4) }) };
          out[2] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(1), p(2), p(6), p(5) }) };
          out[3] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(2), p(3), p(7), p(6) }) };
          out[4] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(3), p(0), p(4), p(7) }) };
          out[5] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(4), p(5), p(6), p(7) }) };
        }
        else if (dim == 3)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Hexahedron, p };
        }
        else
        {
          assert(false);
          out = {};
        }
        return;
      }
      case Polytope::Type::Wedge:
      {
        assert(dim <= 3);
        assert(p.size() == 6);
        if (dim == 0)
        {
          out.resize(6);
          out[0] = { Polytope::Type::Point, Polytope::Key({ p(0) }) };
          out[1] = { Polytope::Type::Point, Polytope::Key({ p(1) }) };
          out[2] = { Polytope::Type::Point, Polytope::Key({ p(2) }) };
          out[3] = { Polytope::Type::Point, Polytope::Key({ p(3) }) };
          out[4] = { Polytope::Type::Point, Polytope::Key({ p(4) }) };
          out[5] = { Polytope::Type::Point, Polytope::Key({ p(5) }) };
        }
        else if (dim == 1)
        {
          // Edges: bottom tri, top tri, verticals
          // bottom: (0->1),(1->2),(2->0)
          // top:    (3->4),(4->5),(5->3)
          // verts:  (0->3),(1->4),(2->5)
          out.resize(9);
          out[0] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(1) }) };
          out[1] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(2) }) };
          out[2] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(0) }) };
          out[3] = { Polytope::Type::Segment, Polytope::Key({ p(3), p(4) }) };
          out[4] = { Polytope::Type::Segment, Polytope::Key({ p(4), p(5) }) };
          out[5] = { Polytope::Type::Segment, Polytope::Key({ p(5), p(3) }) };
          out[6] = { Polytope::Type::Segment, Polytope::Key({ p(0), p(3) }) };
          out[7] = { Polytope::Type::Segment, Polytope::Key({ p(1), p(4) }) };
          out[8] = { Polytope::Type::Segment, Polytope::Key({ p(2), p(5) }) };
        }
        else if (dim == 2)
        {
          // Faces: two triangles + three quads, oriented consistently
          // bottom tri: (0,1,2)
          // top    tri: (3,5,4) (flipped to match prism boundary orientation)
          // quads: (0,1,4,3), (1,2,5,4), (2,0,3,5)
          out.resize(5);
          out[0] = { Polytope::Type::Triangle,      Polytope::Key({ p(0), p(1), p(2) }) };
          out[1] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(0), p(1), p(4), p(3) }) };
          out[2] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(1), p(2), p(5), p(4) }) };
          out[3] = { Polytope::Type::Quadrilateral, Polytope::Key({ p(2), p(0), p(3), p(5) }) };
          out[4] = { Polytope::Type::Triangle,      Polytope::Key({ p(3), p(5), p(4) }) };
        }
        else if (dim == 3)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Wedge, p };
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

  Connectivity<Context::Local>&
  Connectivity<Context::Local>::clear(size_t d, size_t dp)
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_dirty[d][dp] = true;
    m_connectivity[d][dp].clear();
    return *this;
  }

  Connectivity<Context::Local>& Connectivity<Context::Local>::setDirty(size_t d, size_t dp, bool dirty)
  {
    if (d >= m_dirty.size() || dp >= m_dirty[d].size())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Invalid dimensions for setting dirty flag: " << d << " -> " << dp
        << Alert::Raise;
    }
    m_dirty[d][dp] = dirty;
    return *this;
  }

  bool Connectivity<Context::Local>::isDirty(size_t d, size_t dp) const
  {
    if (d >= m_dirty.size() || dp >= m_dirty[d].size())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Invalid dimensions for checking dirty flag: " << d << " -> " << dp
        << Alert::Raise;
    }
    return m_dirty[d][dp];
  }

  size_t Connectivity<Context::Local>::getMaximalDimension() const
  {
    return m_maximalDimension;
  }
}
