/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESH_H
#define RODIN_EXTERNAL_MMG_MESH_H

#include "Rodin/Geometry.h"

namespace Rodin::External::MMG
{
  class Mesh : public Geometry::Mesh<Context::Serial>
  {
    public:
      struct Edge
      {
        std::pair<int, int> endpoints;
        int ref;
      };

      using Parent = Geometry::Mesh<Context::Serial>;

      Mesh() : Parent()
      {}

      Mesh(const Mesh& other)
        : Parent(other),
          m_corners(other.m_corners),
          m_ridges(other.m_ridges),
          m_edges(other.m_edges)
      {}

      Mesh(Mesh&& other)
        : Parent(std::move(other)),
          m_corners(std::move(other.m_corners)),
          m_ridges(std::move(other.m_ridges)),
          m_edges(std::move(other.m_edges))
      {}

      Mesh& operator=(Mesh&& other)
      {
        Parent::operator=(std::move(other));
        m_edges = std::move(other.m_edges);
        m_ridges = std::move(other.m_ridges);
        m_corners = std::move(other.m_corners);
        return *this;
      }

      void save(
          const boost::filesystem::path& filename,
          IO::FileFormat fmt = IO::FileFormat::MFEM,
          size_t precison = 16) const override;

      Mesh& load(
          const boost::filesystem::path& filename,
          IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      Mesh& corner(int vertexIdx);

      Mesh& ridge(int edgeIdx);

      Mesh& edge(const std::pair<int, int>& endpoints, int ref);

      const std::set<int>& getCorners() const
      {
        return m_corners;
      }

      const std::set<int>& getRidges() const
      {
        return m_ridges;
      }

      const std::vector<Edge>& getEdges() const
      {
        return m_edges;
      }

    private:
       std::set<int> m_corners;

       std::set<int> m_ridges;
       std::vector<Edge> m_edges;
  };
}

#endif
