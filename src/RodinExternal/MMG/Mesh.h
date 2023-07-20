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
       * @brief Class used to build MMG::Mesh instances.
       */
      class Builder : public Parent::Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          /**
           * @brief Deleted copy constructor.
           */
          Builder(const Builder&) = delete;

          /**
           * @brief Move constructor.
           */
          Builder(Builder&& other)
            : Parent::Builder(std::move(other)),
              m_cornerIndex(std::move(other.m_cornerIndex)),
              m_ridgeIndex(std::move(other.m_ridgeIndex))
          {}

          /**
           * @brief Move assignment.
           */
          Builder& operator=(Builder&& other);

          /**
           * @brief Adds the vertex to the corner index.
           */
          Builder& corner(Index vertexIdx);

          /**
           * @brief Adds the edge to the ridge index.
           */
          Builder& ridge(Index edgeIdx);

          /**
           * @brief Finishes construction of the MMG::Mesh and returns it.
           */
          Mesh finalize();

        private:
          CornerIndex m_cornerIndex;
          RidgeIndex  m_ridgeIndex;
      };

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

      Mesh& operator=(Parent&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      inline
      MMG::Mesh::Builder build() const
      {
        return MMG::Mesh::Builder();
      }

      /**
       * @brief Adds the vertex to the corner index.
       */
      Mesh& setCorner(Index vertexIdx);

      /**
       * @brief Adds the edge to the ridge index.
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
       * @brief Gets the index of corners.
       */
      CornerIndex& getCorners()
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

      /**
       * @brief Gets the index of ridges.
       */
      RidgeIndex& getRidges()
      {
        return m_ridgeIndex;
      }

      void save(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt = IO::FileFormat::MFEM,
         size_t precison = 16) const override;

      Mesh& load(
         const boost::filesystem::path& filename,
         IO::FileFormat fmt = IO::FileFormat::MFEM) override;

    private:
      CornerIndex m_cornerIndex;
      RidgeIndex  m_ridgeIndex;
  };
}

#endif
