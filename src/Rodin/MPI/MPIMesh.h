/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MPI_MESH_H
#define RODIN_GEOMETRY_MPI_MESH_H

#include "Rodin/Configure.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Shard.h"

#include "MPIContext.h"

#include "MPIConnectivity.h"

namespace Rodin::Geometry
{
  using MPIMesh = Mesh<Context::MPI>;

  template <>
  class Mesh<Context::MPI> : public MeshBase
  {
    public:
      class Builder
      {
        public:
          Builder(const Context::MPI& context);

          Builder& initialize(Shard&& shard);

          Mesh finalize();

        private:
          Shard m_shard;
          Context::MPI m_context;
      };

      Mesh(const Context::MPI& context)
        : m_context(context)
      {}

      Shard& getShard();

      const Shard& getShard() const;

      Mesh& scale(Real c) override;

      void flush() override;

      bool isSubMesh() const override;

      size_t getDimension() const override;

      size_t getSpaceDimension() const override;

      const Context::MPI& getContext() const override;

      size_t getPolytopeCount(size_t d) const override;

      size_t getPolytopeCount(Polytope::Type g) const override;

      /**
       * @brief Returns a PolytopeIterator for the given global polytope index.
       * @param dimension The topological dimension of the polytopes to iterate over.
       * @param globalIdx The global index (consistent across processes) of the
       * polytope to start iteration from.
       *
       * - If the provided global index is within the range that belongs to
       *   this process, the function constructs a PolytopeIterator that starts
       *   at the corresponding local index.
       * - Otherwise, the function returns an "empty" PolytopeIterator,
       *   indicating that the requested global polytope is not available on
       *   the calling process.
       *
       * @return A PolytopeIterator starting at the appropriate local index if
       * this process owns that global index; otherwise, an empty iterator.
       */
      PolytopeIterator getPolytope(size_t dimension, Index globalIdx) const override;

      CellIterator getCell(Index globalIdx) const override;

      FaceIterator getFace(Index globalIdx) const override;

      VertexIterator getVertex(Index globalIdx) const override;

      PolytopeIterator getPolytope(size_t dimension) const override;

      CellIterator getCell() const override;

      FaceIterator getFace() const override;

      VertexIterator getVertex() const override;

      FaceIterator getBoundary() const override;

      FaceIterator getInterface() const override;

      std::optional<Index> getLocalIndex(size_t dimension, Index globalIdx) const;

      bool isInterface(Index globalFaceIdx) const override;

      bool isBoundary(Index globalFaceIdx) const override;

      Real getVolume() const override;

      Real getVolume(Attribute attr) const override;

      Real getVolume(const FlatSet<Attribute>& attr) const override;

      Real getPerimeter() const override;

      Real getPerimeter(Attribute attr) const override;

      Real getPerimeter(const FlatSet<Attribute>& attr) const override;

      Real getArea() const override;

      Real getArea(Attribute attr) const override;

      Real getArea(const FlatSet<Attribute>& attr) const override;

      Real getMeasure(size_t d) const override;

      Real getMeasure(size_t d, Attribute attr) const override;

      Real getMeasure(size_t d, const FlatSet<Attribute>& attr) const override;

      Polytope::Type getGeometry(size_t dimension, Index globalIdx) const override;

      Attribute getAttribute(size_t dimension, Index globalIdx) const override;

      Mesh& setAttribute(const std::pair<size_t, Index>& p, Attribute attr) override;

      Mesh& setVertexCoordinates(Index globalIdx, Real s, size_t i) override;

      Mesh& setVertexCoordinates(Index globalIdx, const Math::PointVector& coords) override;

      virtual const PolytopeTransformation& getPolytopeTransformation(size_t dimension, Index globalIdx) const override;

      virtual Mesh& setPolytopeTransformation(const std::pair<size_t, Index> p, PolytopeTransformation* trans) override;

      MPIMesh& load(
        std::function<boost::filesystem::path(size_t)> filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM);

      MPIMesh& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      void save(
        std::function<boost::filesystem::path(size_t)> filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM);

      void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const override;

      Eigen::Map<const Math::PointVector> getVertexCoordinates(Index globalIdx) const override;

      SubMeshBase& asSubMesh() override
      {
        throw std::runtime_error("asSubMesh() not implemented");
      }

      const SubMeshBase& asSubMesh() const override
      {
        throw std::runtime_error("asSubMesh() not implemented");
      }

      Connectivity<Context::MPI>& getConnectivity() override
      {
        throw std::runtime_error("getConnectivity() not implemented");
        // return m_connectivity;
      }

      const Connectivity<Context::MPI>& getConnectivity() const override
      {
        throw std::runtime_error("getConnectivity() not implemented");
        // return m_connectivity;
      }

    private:
      Context::MPI m_context;
      Shard m_shard;

      mutable FlatMap<Index, Math::PointVector> m_vertices;
      mutable TransformationIndex m_transformationIndex;
      // Connectivity<Context::MPI> m_connectivity;
  };
}

#endif
