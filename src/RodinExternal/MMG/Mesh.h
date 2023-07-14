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
  /**
   * @brief Mesh class which has support for MMG functionalities.
   */
  class Mesh : public Geometry::Mesh<Context::Serial>
  {
    public:
      /// Parent class
      using Parent = Geometry::Mesh<Context::Serial>;

      /// Index of corners in the mesh
      using CornerIndex = IndexSet;

      /// Index of ridges in the mesh
      using RidgeIndex = IndexSet;

      /**
      * @brief Constructs an empty mesh with no elements.
      */
      Mesh() = default;

      /**
       * @brief Copy constructor.
       */
      Mesh(const Mesh& other)
        : Parent(other),
          m_cornerIndex(other.m_cornerIndex),
          m_ridgeIndex(other.m_ridgeIndex)
      {}

      /**
       * @brief Move constructor.
       */
      Mesh(Mesh&& other)
        : Parent(std::move(other)),
          m_cornerIndex(std::move(other.m_cornerIndex)),
          m_ridgeIndex(std::move(other.m_ridgeIndex))
      {}

      /**
       * @brief Move assignment.
       */
      Mesh& operator=(Mesh&& other)
      {
        Parent::operator=(std::move(other));
        m_ridgeIndex = std::move(other.m_ridgeIndex);
        m_cornerIndex = std::move(other.m_cornerIndex);
        return *this;
      }

      /**
       * @brief Move assignment from a Mesh<Context::Serial> object.
       */
      Mesh& operator=(Parent&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      void save(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt = IO::FileFormat::MFEM,
         size_t precison = 16) const override;

      Mesh& load(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      /**
       * @brief Adds the vertex to the corner index.
       */
      Mesh& setCorner(Index vertexIdx);

      /**
       * @brief Adds the edge to the corner index.
       */
      Mesh& setRidge(Index edgeIdx);

      /**
       * @brief Gets the index of corners.
       */
      const CornerIndex& getCorners() const
      {
        return m_cornerIndex;
      }

      /**
       * @brief Gets the index of ridges.
       */
      const RidgeIndex& getRidges() const
      {
        return m_ridgeIndex;
      }

    private:
      CornerIndex m_cornerIndex;
      RidgeIndex  m_ridgeIndex;
  };
}

#endif
