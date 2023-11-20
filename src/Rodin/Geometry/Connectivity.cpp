/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Connectivity.h"

namespace Rodin::Geometry
{
  MeshConnectivity::MeshConnectivity()
  {
    m_count.resize(1, 0);
  }

  MeshConnectivity& MeshConnectivity::initialize(size_t maximalDimension)
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

  MeshConnectivity& MeshConnectivity::reserve(size_t d, size_t count)
  {
    assert(d < m_connectivity.size());
    m_index[d].write(
        [&](auto& obj)
        {
          obj.left.rehash(count);
          obj.right.rehash(count);
        });
    m_geometry[d].write(
        [&](auto& obj)
        {
          obj.reserve(count);
        });
    m_connectivity[d][0].write(
        [&](auto& obj)
        {
          obj.reserve(count);
        });
    return *this;
  }

  MeshConnectivity& MeshConnectivity::nodes(size_t count)
  {
    m_count[0].write() = count;
    m_gcount[Geometry::Polytope::Type::Point].write() = count;
    for (size_t i = 0; i < count; i++)
    {
      m_index[0].write(
          [&](auto& obj)
          {
            auto p = obj.left.insert({ IndexArray{{ i }}, i });
            assert(p.second);
          });
    }
    return *this;
  }

  MeshConnectivity&
  MeshConnectivity::polytope(Polytope::Type t, const Array<Index>& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);
    std::pair<PolytopeIndex::left_iterator, bool> res;
    m_index[d].write(
        [&](PolytopeIndex& obj)
        {
          size_t countd;
          m_count[d].read([&](const auto& obj) { countd = obj; });
          res = obj.left.insert({ in, countd });
        });
    const auto& it = res.first;
    const bool& inserted = res.second;
    if (inserted)
    {
      m_connectivity[d][0].write(
          [&](auto& obj)
          {
            obj.emplace_back().insert_unique(it->first.begin(), it->first.end());
          });
      m_geometry[d].write(
          [&](auto& obj)
          {
            obj.push_back(t);
          });
      m_count[d].write([](auto& obj){ obj += 1; });
      m_gcount[t].write([](auto& obj){ obj += 1; });
      m_dirty[d][0].write() = false;
    }
    return *this;
  }

  MeshConnectivity&
  MeshConnectivity::polytope(Polytope::Type t, Array<Index>&& in)
  {
    assert(in.size() > 0);
    const size_t d = Polytope::getGeometryDimension(t);
    assert(d > 0);
    assert(d <= m_maximalDimension);

    size_t countd;
    m_count[d].read([&](const auto& obj) { countd = obj; });

    std::pair<PolytopeIndex::left_iterator, bool> res;
    m_index[d].write(
        [&](auto& obj)
        {
          res = obj.left.insert({ std::move(in), countd });
        });
    const auto& it = res.first;
    const bool& inserted = res.second;
    if (inserted)
    {
      m_connectivity[d][0].write(
          [&](auto& obj)
          {
            obj.emplace_back().insert_unique(it->first.begin(), it->first.end());
          });
      m_geometry[d].write(
          [&](auto& obj)
          {
            obj.push_back(t);
          });
      m_count[d].write([](auto& obj){ obj += 1; });
      m_gcount[t].write([](auto& obj){ obj += 1; });
      m_dirty[d][0].write() = false;
    }
    return *this;
  }

  const MeshConnectivity::PolytopeIndex& MeshConnectivity::getIndexMap(size_t dim) const
  {
    return m_index[dim].read();
  }

  const std::optional<Index> MeshConnectivity::getIndex(size_t dim, const IndexArray& key) const
  {
    bool found;
    PolytopeIndex::left_const_iterator it;
    m_index[dim].read(
        [&](const auto& obj)
        {
          it = obj.left.find(key);
          found = (it == obj.left.end());
        });
    if (found)
      return it->second;
    else
      return {};
  }

  const Incidence& MeshConnectivity::getIncidence(size_t d, size_t dp) const
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    return m_connectivity[d][dp].read();
  }

  const IndexSet& MeshConnectivity::getIncidence(
      const std::pair<size_t, size_t> p, Index idx) const
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    assert(idx < m_connectivity[d][dp].read().size());
    return m_connectivity[d][dp].read()[idx];
  }

  size_t MeshConnectivity::getCount(size_t dim) const
  {
    return m_count[dim].read();
  }

  size_t MeshConnectivity::getCount(Polytope::Type g) const
  {
    return m_gcount[g].read();
  }

  size_t MeshConnectivity::getMeshDimension() const
  {
    for (int i = m_count.size() - 1; i >= 0; i--)
    {
      bool b;
      m_count[i].read([&](const auto& obj) { b = (obj > 0); });
      if (b)
        return i;
    }
    return 0;
  }

  Polytope::Type MeshConnectivity::getGeometry(size_t d, Index idx) const
  {
    if (d == 0)
      return Polytope::Type::Point;
    else
    {
      Polytope::Type t;
      m_geometry[d].read([&](const auto& obj) { t = obj[idx]; });
      return t;
    }
  }

  const Array<Index>& MeshConnectivity::getPolytope(size_t d, Index idx) const
  {
    return m_index[d].read().right.at(idx);
  }

  MeshConnectivity& MeshConnectivity::setIncidence(const std::pair<size_t, size_t>& p, Incidence&& inc)
  {
    const auto& [d, dp] = p;
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_connectivity[d][dp].write() = std::move(inc);
    return *this;
  }

  MeshConnectivity& MeshConnectivity::compute(size_t d, size_t dp)
  {
    const size_t D = getMeshDimension();
    if (d == D && dp == 0)
      return *this;
    if (m_dirty[D][D].read())
      transpose(0, D).intersection(D, D, 0);
    assert(!m_dirty[D][D].read());
    if (d != D && d != 0 && (m_dirty[D][d].read() || m_dirty[d][0].read()))
      build(d);
    assert(!m_dirty[D][d].read());
    assert(!m_dirty[d][0].read() || d == D || d == 0);
    if (dp != D && dp != 0 && (m_dirty[D][dp].read() || m_dirty[dp][0].read()))
      build(dp);
    assert(!m_dirty[D][dp].read());
    assert(!m_dirty[dp][0].read() || dp == D || dp == 0);
    if (m_dirty[d][dp].read())
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
    m_dirty[d][dp].write() = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::build(size_t d)
  {
    const size_t D = getMeshDimension();
    assert(d > 0);
    assert(d < D);
    assert(!m_dirty[D][0].read());
    assert(!m_dirty[D][D].read());
    for (Index i = 0; i < m_count[D].read(); i++)
    {
      IndexSet s;
      std::vector<SubPolytope> subpolytopes;
      local(subpolytopes, d, i);
      for (auto& gv : subpolytopes)
      {
        const auto& geometry = gv.geometry;
        auto& vertices = gv.vertices;
        size_t countd;
        m_count[d].read([&](const auto& obj) { countd = obj; });
        std::pair<PolytopeIndex::left_const_iterator, bool> insert;
        m_index[d].write([&](auto& obj) { insert = obj.left.insert({ std::move(vertices), countd }); });
        const PolytopeIndex::left_const_iterator it = insert.first;
        const bool inserted = insert.second;
        const auto& v = it->first;
        const auto& idx = it->second;
        if (inserted)
        {
          m_geometry[d].write(
              [&](auto& obj)
              {
                obj.push_back(geometry);
              });
          m_connectivity[d][0].write(
              [&](auto& obj)
              {
                obj.emplace_back().insert_unique(v.begin(), v.end());
              });
        }
        m_count[d].write([&](auto& obj) { obj += inserted && !(d == D || d == 0); } );
        m_gcount[geometry].write([&](auto& obj) { obj += inserted && !(d == D || d == 0); });
        s.insert(idx);
      }
      m_connectivity[D][d].write(
          [&](auto& obj)
          {
            obj.push_back(std::move(s));
          });
    }
    m_dirty[D][d].write() = false;
    m_dirty[d][0].write() = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::transpose(size_t d, size_t dp)
  {
    assert(d < dp);
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    size_t countd;
    m_count[d].read([&](const auto& obj) { countd = obj; });
    m_connectivity[d][dp].write([&](auto& obj) { obj.resize(countd); });
    size_t countdp;
    m_count[dp].read([&](const auto& obj) { countdp = obj; });
    for (Index j = 0; j < countdp; j++)
    {
      const IndexSet* isj = nullptr;
      m_connectivity[dp][d].read([&](const auto& obj) { isj = &obj[j]; });
      for (Index i : *isj)
        m_connectivity[d][dp].write([&](auto& obj) { obj[i].insert(j); });
    }
    m_dirty[d][dp].write() = false;
    return *this;
  }

  MeshConnectivity& MeshConnectivity::intersection(size_t d, size_t dp, size_t dpp)
  {
    assert(d >= dp);
    size_t countd;
    m_count[d].read([&](const auto& obj) { countd = obj; });
    m_connectivity[d][dp].write([&](auto& obj) { obj.resize(countd); });
    for (Index i = 0; i < countd; i++)
    {
      m_connectivity[d][dpp].read(
          [&](const auto& obj) { assert(i < obj.size()); });
      const IndexSet* isi = nullptr;
      m_connectivity[d][dpp].read([&](const auto& obj) { isi = &obj[i]; });
      for (Index k : *isi)
      {
        m_connectivity[dpp][dp].read(
            [&](const auto& obj) { assert(k < obj.size()); });
        const IndexSet* isk = nullptr;
        m_connectivity[dpp][dp].read([&](const auto& obj) { isk = &obj[k]; });
        for (Index j : *isk)
        {
          m_connectivity[d][0].read(
              [&](const auto& obj) { assert(i < obj.size()); });
          m_connectivity[dp][0].read(
              [&](const auto& obj) { assert(j < obj.size()); });
          const IndexSet* isd0i = nullptr;
          m_connectivity[d][0].read([&](const auto& obj) { isd0i = &obj[i]; });
          const IndexSet* isdp0j = nullptr;
          m_connectivity[dp][0].read([&](const auto& obj) { isdp0j = &obj[j]; });
          if ((d == dp && i != j) ||
              (d > dp && std::includes(isd0i->begin(), isd0i->end(), isdp0j->begin(), isdp0j->end())))
          {
            m_connectivity[d][dp].write([&](auto& obj) { obj[i].insert(j); });
          }
        }
      }
    }
    m_dirty[d][dp].write() = false;
    return *this;
  }

  void MeshConnectivity::local(std::vector<MeshConnectivity::SubPolytope>& out, size_t dim, Index i)
  {
    const size_t D = getMeshDimension();
    const IndexArray* p = nullptr;
    m_index[D].read([&](const auto& obj) { p = &obj.right.at(i); });
    Polytope::Type geometryD;
    m_geometry[D].read([&](const auto& obj) { geometryD = obj[i]; });
    switch (geometryD)
    {
      case Polytope::Type::Point:
      {
        assert(dim == 0);
        assert(p->size()  == 0);
        out.resize(1);
        out[0].geometry = Polytope::Type::Point;
        out[0].vertices.resize(1);
        out[0].vertices.coeffRef(0) = i;
        return;
      }
      case Polytope::Type::Segment:
      {
        assert(dim <= 1);
        assert(p->size() == 2);
        if (dim == 0)
        {
          out.resize(2);
          out[0].geometry = Polytope::Type::Point;
          out[0].vertices.resize(1);
          out[0].vertices.coeffRef(0) = p->coeff(0);
          out[1].geometry = Polytope::Type::Point;
          out[1].vertices.resize(1);
          out[1].vertices.coeffRef(0) = p->coeff(1);
        }
        else if (dim == 1)
        {
          out.resize(1);
          out[0].geometry = Polytope::Type::Segment;
          out[0].vertices = *p;
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
        assert(p->size() == 3);
        if (dim == 0)
        {
          out.resize(3);
          out[0].geometry = Polytope::Type::Point;
          out[0].vertices.resize(1);
          out[0].vertices.coeffRef(0) = p->coeff(0);

          out[1].geometry = Polytope::Type::Point;
          out[1].vertices.resize(1);
          out[1].vertices.coeffRef(0) = p->coeff(1);

          out[2].geometry = Polytope::Type::Point;
          out[2].vertices.resize(1);
          out[2].vertices.coeffRef(0) = p->coeff(2);
        }
        else if (dim == 1)
        {
          out.resize(3);

          out[0].geometry = Polytope::Type::Segment;
          out[0].vertices.resize(2);
          out[0].vertices.coeffRef(0) = p->coeff(0);
          out[0].vertices.coeffRef(1) = p->coeff(1);

          out[1].geometry = Polytope::Type::Segment;
          out[1].vertices.resize(2);
          out[1].vertices.coeffRef(0) = p->coeff(1);
          out[1].vertices.coeffRef(1) = p->coeff(2);

          out[2].geometry = Polytope::Type::Segment;
          out[2].vertices.resize(2);
          out[2].vertices.coeffRef(0) = p->coeff(0);
          out[2].vertices.coeffRef(1) = p->coeff(2);
        }
        else if (dim == 2)
        {
          out.resize(1);
          out[0].geometry = Polytope::Type::Triangle;
          out[0].vertices = *p;
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
        assert(p->size() == 4);
        if (dim == 0)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Point, { { (*p)(0) } } };
          out[1] = { Polytope::Type::Point, { { (*p)(1) } } };
          out[2] = { Polytope::Type::Point, { { (*p)(2) } } };
          out[3] = { Polytope::Type::Point, { { (*p)(3) } } };
        }
        else if (dim == 1)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Segment, { { (*p)(0), (*p)(1) } } };
          out[1] = { Polytope::Type::Segment, { { (*p)(1), (*p)(3) } } };
          out[2] = { Polytope::Type::Segment, { { (*p)(3), (*p)(2) } } };
          out[3] = { Polytope::Type::Segment, { { (*p)(2), (*p)(0) } } };
        }
        else if (dim == 2)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Quadrilateral, *p };
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
        assert(p->size() == 4);
        if (dim == 0)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Point, { { (*p)(0) } } };
          out[1] = { Polytope::Type::Point, { { (*p)(1) } } };
          out[2] = { Polytope::Type::Point, { { (*p)(2) } } };
          out[3] = { Polytope::Type::Point, { { (*p)(3) } } };
        }
        else if (dim == 1)
        {
          out.resize(6);
          out[0] = { Polytope::Type::Segment, { { (*p)(0), (*p)(1) } } };
          out[1] = { Polytope::Type::Segment, { { (*p)(1), (*p)(2) } } };
          out[2] = { Polytope::Type::Segment, { { (*p)(0), (*p)(2) } } };
          out[3] = { Polytope::Type::Segment, { { (*p)(0), (*p)(3) } } };
          out[4] = { Polytope::Type::Segment, { { (*p)(1), (*p)(3) } } };
          out[5] = { Polytope::Type::Segment, { { (*p)(2), (*p)(3) } } };
        }
        else if (dim == 2)
        {
          out.resize(4);
          out[0] = { Polytope::Type::Triangle, { { (*p)(0), (*p)(1), (*p)(2) } } };
          out[1] = { Polytope::Type::Triangle, { { (*p)(0), (*p)(1), (*p)(3) } } };
          out[2] = { Polytope::Type::Triangle, { { (*p)(0), (*p)(2), (*p)(3) } } };
          out[3] = { Polytope::Type::Triangle, { { (*p)(1), (*p)(2), (*p)(3) } } };
        }
        else if (dim == 3)
        {
          out.resize(1);
          out[0] = { Polytope::Type::Tetrahedron, *p };
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

  MeshConnectivity& MeshConnectivity::clear(size_t d, size_t dp)
  {
    assert(d < m_connectivity.size());
    assert(dp < m_connectivity[d].size());
    m_dirty[d][dp].write() = true;
    m_connectivity[d][dp].write().clear();
    return *this;
  }
}
