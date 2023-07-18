/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MESH_H
#define RODIN_GEOMETRY_MESH_H

#include <set>
#include <string>
#include <deque>

#include <boost/filesystem.hpp>

#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Utility/IsSpecialization.h"
#include "Rodin/Variational/Traits.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1/ForwardDecls.h"

#include "ForwardDecls.h"
#include "Connectivity.h"
#include "Simplex.h"
#include "SimplexCount.h"
#include "SimplexIndexed.h"
#include "SimplexIterator.h"
#include "SimplexTransformation.h"

namespace Rodin::Geometry
{
  /**
  * @brief Abstract base class for Mesh objects.
  */
  class MeshBase
  {
    public:
      virtual ~MeshBase() = default;

      virtual MeshBase& scale(Scalar c) = 0;

      /**
       * @brief Gets the maximum number @f$ t @f$ by which the mesh will
       * remain valid, when displacing by @f$ u @f$.
       * @param[in] u Displacement at each node
       *
       * This function will calculate the maximum number @f$ t @f$ so that
       * the displacement
       * @f[
       *   x \mapsto x + t u(x)
       * @f]
       * gives a valid mesh without actually displacing the mesh.
       *
       * @note The vector dimension of @f$ u @f$ must be equal to the
       * space dimension.
       *
       * @returns Maximum time so that the mesh remains valid.
       */
      template <class FES>
      Scalar getMaximumDisplacement(const Variational::GridFunction<FES>& u)
      {
        Scalar res;
        assert(false);
        // getHandle().CheckDisplacements(u.getHandle(), res);
        return res;
      }

      /**
       * @brief Indicates whether the mesh is a surface or not.
       *
       * A mesh is considered a surface mesh if has codimension of 1, meaning
       * the difference between its space dimension and its topological
       * dimension is 1.
       *
       * @returns True if mesh is a surface, false otherwise.
       */
      bool isSurface() const;

      /**
       * @brief Gets the number of vertices in the mesh.
       */
      inline
      size_t getVertexCount() const
      {
        return getCount(0);
      }

      /**
       * @brief Gets the number of faces in the mesh.
       */
      inline
      size_t getFaceCount() const
      {
        return getCount(getDimension() - 1);
      }

      /**
       * @brief Gets the number of elements in the mesh.
       */
      inline
      size_t getElementCount() const
      {
        return getCount(getDimension());
      }

      inline
      Attribute getFaceAttribute(Index index) const
      {
        return getAttribute(getDimension() - 1, index);
      }

      inline
      Attribute getElementAttribute(Index index) const
      {
        return getAttribute(getDimension(), index);
      }

      /**
       * @brief Gets the total volume of the mesh.
       * @returns Sum of all element volumes.
       */
      Scalar getVolume();

      /**
       * @brief Gets the sum of the volumes of the elements given by the
       * specified attribute.
       * @param[in] attr Attribute of elements
       * @returns Sum of element volumes with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the volume.
       */
      Scalar getVolume(Attribute attr);

      /**
       * @brief Gets the total perimeter of the mesh.
       * @returns Sum of all element perimeters.
       */
      Scalar getPerimeter();

      /**
       * @brief Gets the sum of the perimeters of the elements given by the
       * specified attribute.
       * @param[in] attr Attribute of elements
       * @returns Sum of element perimeters with given attribute
       * @note If the element attribute does not exist then this function
       * will return 0 as the perimeter.
       */
      Scalar getPerimeter(Attribute attr);

      /**
       * @brief Gets the labels of the domain elements in the mesh.
       * @returns Set of all the attributes in the mesh object.
       * @see getBoundaryAttributes() const
       */
      virtual const FlatSet<Attribute>& getAttributes(size_t d) const = 0;

      bool operator==(const MeshBase& other) const
      {
        return this == &other;
      }

      bool operator!=(const MeshBase& other) const
      {
        return this != &other;
      }

      /**
       * @brief Gets the dimension of the elements.
       * @returns Dimension of the elements.
       * @see getSpaceDimension() const
       */
      virtual size_t getDimension() const = 0;

      /**
       * @brief Gets the dimension of the ambient space
       * @returns Dimension of the space which the mesh is embedded in
       * @see getDimension() const
       */
      virtual size_t getSpaceDimension() const = 0;

      virtual MeshBase& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) = 0;

      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const = 0;

      /**
       * @brief Indicates whether the mesh is a submesh or not.
       * @returns True if mesh is a submesh, false otherwise.
       *
       * A Mesh which is also a SubMesh may be casted into down to access
       * the SubMesh functionality. For example:
       * @code{.cpp}
       * if (mesh.isSubMesh())
       * {
       *   // Cast is well defined
       *   auto& submesh = static_cast<SubMesh&>(mesh);
       * }
       * @endcode
       *
       */
      virtual bool isSubMesh() const = 0;

      virtual bool isInterface(Index faceIdx) const = 0;

      virtual bool isBoundary(Index faceIdx) const = 0;

      virtual FaceIterator getBoundary() const = 0;

      virtual FaceIterator getInterface() const = 0;

      virtual size_t getCount(size_t dim) const = 0;

      virtual size_t getCount(Polytope::Geometry g) const = 0;

      virtual ElementIterator getElement(Index idx = 0) const = 0;

      virtual FaceIterator getFace(Index idx = 0) const = 0;

      virtual VertexIterator getVertex(Index idx = 0) const = 0;

      virtual PolytopeIterator getPolytope(size_t dimension, Index idx) const = 0;

      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const = 0;

      virtual Polytope::Geometry getGeometry(size_t dimension, Index idx) const = 0;

      virtual Attribute getAttribute(size_t dimension, Index index) const = 0;

      virtual MeshBase& setAttribute(const std::pair<size_t, Index>&, Attribute attr) = 0;

      virtual MeshConnectivity& getConnectivity() = 0;

      virtual const MeshConnectivity& getConnectivity() const = 0;

      virtual Eigen::Map<const Math::SpatialVector> getVertexCoordinates(Index idx) const = 0;

      virtual void flush() = 0;
  };

  using BoundaryIndex = IndexSet;

  using AttributeIndex = PolytopeIndexed<Geometry::Attribute>;

  using TransformationIndex = PolytopeIndexed<std::unique_ptr<PolytopeTransformation>>;

  /**
  * @brief Represents the subdivision of some domain into faces of (possibly)
  * different geometries.
  */
  template <>
  class Mesh<Context::Serial> : public MeshBase
  {
    public:
      /**
       * @brief Class used to build Mesh<Context::Serial> instances.
       */
      class Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          virtual ~Builder() = default;

          /**
           * @brief Deleted copy constructor.
           */
          Builder(const Builder&) = delete;

          /**
           * @brief Move constructor.
           */
          Builder(Builder&&) = default;

          /**
           * @brief Deleted copy assignment.
           */
          Builder& operator=(const Builder&) = delete;

          /**
           * @brief Move assignment.
           */
          Builder& operator=(Builder&&) = default;

          /**
           * @brief Reserves memory accross the data structure for the polytopes of
           * the given dimension.
           */
          Builder& reserve(size_t d, size_t count);

          /**
           * @brief Initializes construction of the mesh object.
           */
          Builder& initialize(size_t sdim);

          /**
           * @brief Sets the number of nodes in the mesh.
           */
          Builder& nodes(size_t n);

          /**
           * @brief Adds vertex with coordinates given by the fixed size array.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          template <size_t Size>
          inline
          Builder& vertex(const Scalar (&data)[Size])
          {
            assert(Size == m_sdim);
            m_vertices.col(m_nodes++) = data;
            return *this;
          }

          /**
           * @brief Adds vertex with coordinates given by the initializer list.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(std::initializer_list<Scalar> l);

          /**
           * @brief Adds vertex with coordinates given by the array pointer.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Scalar* data);

          /**
           * @brief Adds vertex with coordinates given by the mapped memory.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Eigen::Map<const Math::Vector>& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(Math::Vector&& x);

          /**
           * @brief Adds vertex with coordinates given by the vector.
           *
           * This method requires nodes(size_t) to be called beforehand.
           */
          Builder& vertex(const Math::Vector& x);

          /**
           * @brief Sets the attribute of the given polytope.
           */
          Builder& attribute(const std::pair<size_t, Index>& p, Attribute attr);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Geometry t, std::initializer_list<Index> vs)
          {
            return polytope(t, IndexArray({ vs }));
          }

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Geometry t, const IndexArray& vs);

          /**
           * @brief Adds polytope defined by the given vertices.
           */
          Builder& polytope(Polytope::Geometry t, IndexArray&& vs);

          /**
           * @brief Finalizes construction of the Mesh<Context::Serial> object.
           */
          Mesh finalize();

          Builder& setVertices(Math::Matrix&& connectivity);

          Builder& setConnectivity(MeshConnectivity&& connectivity);

          Builder& setAttributeIndex(AttributeIndex&& connectivity);

          Builder& setTransformationIndex(TransformationIndex&& connectivity);

        private:
          size_t m_sdim;
          size_t m_nodes;

          Math::Matrix m_vertices;
          MeshConnectivity m_connectivity;

          AttributeIndex m_attributeIndex;
          TransformationIndex m_transformationIndex;

          std::vector<FlatSet<Attribute>> m_attributes;
      };

      /**
       * @brief Generates a uniform grid for a given geometry.
       */
      static Mesh UniformGrid(Polytope::Geometry g, size_t h, size_t w);

      /**
      * @brief Constructs an empty mesh with no elements.
      */
      Mesh()
        : m_sdim(0)
      {}

      Mesh(const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        load(filename, fmt);
      }

      /**
      * @brief Performs a copy of another mesh.
      */
      Mesh(const Mesh& other)
        : m_sdim(other.m_sdim),
          m_vertices(other.m_vertices),
          m_connectivity(other.m_connectivity),
          m_attributeIndex(other.m_attributeIndex),
          m_attributes(other.m_attributes)
      {}

      /**
      * @brief Move constructs the mesh from another mesh.
      */
      Mesh(Mesh&& other) = default;

      Mesh& operator=(const Mesh& other) = delete;

      /**
      * @brief Move assigns the mesh from another mesh.
      */
      Mesh& operator=(Mesh&&) = default;

      /**
       * @brief Displaces the mesh nodes by the displacement @f$ u @f$.
       * @param[in] u Displacement at each node
       *
       * Given a vector valued function @f$ \vec{u} @f$, the method will perform the
       * displacement
       * @f[
       *   x \mapsto x + \vec{u}(x)
       * @f]
       * at each node @f$ x @f$ of the mesh.
       *
       * @note The vector dimension of @f$ \vec{u} @f$ must be equal to the
       * space dimension.
       *
       * @returns Reference to this (for method chaining)
       */
      template <class FunctionDerived>
      Mesh& displace(const Variational::FunctionBase<FunctionDerived>& u)
      {
        for (auto it = getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(*it, it->getTransformation(),
              Polytope::getVertices(Polytope::Geometry::Point).col(0), it->getCoordinates());
          m_vertices.col(it->getIndex()) += u(p);
        }
        return *this;
      }

      /**
       * @internal
       */
      Mesh& displace(const
          Variational::GridFunction<Variational::P1<Math::Vector,
          Context::Serial, Geometry::Mesh<Context::Serial>>>& u);

      const PolytopeIndexed<Attribute>& getAttributeIndex() const
      {
        return m_attributeIndex;
      }

      const PolytopeIndexed<std::unique_ptr<PolytopeTransformation>>& getTransformationIndex() const
      {
        return m_transformationIndex;
      }

      inline
      const Math::Matrix& getVertices() const
      {
        return m_vertices;
      }

      /**
      * @brief Skins the mesh to obtain its boundary mesh
      * @returns SubMesh object to the boundary region of the mesh
      *
      * This function "skins" the mesh to return the boundary as a new SubMesh
      * object. The resulting mesh will be embedded in the original space
      * dimension.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 2
      * @f$ are included if the connectivity @f$ (D - 1) \longrightarrow d @f$
      * is already computed in the mesh.
      */
      virtual SubMesh<Context::Serial> skin() const;

      /**
      * @brief Trims the elements with the given attribute.
      * @param[in] attr Attribute to trim
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call trim(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context::Serial> trim(Attribute attr) const;

      /**
      * @brief Trims the elements with the given attribute.
      * @param[in] attrs Attributes to trim
      * @returns SubMesh object to the remaining region of the mesh
      *
      * This function will trim discard all the elements that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of elements that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context::Serial> trim(const FlatSet<Attribute>& attrs) const;

      /**
      * @brief Keeps the elements with the given attribute.
      * @param[in] attr Attribute to keep
      * @returns SubMesh of the remaining region mesh
      *
      * Convenience function to call keep(const std::FlatSet<Attribute>&) with
      * only one attribute.
      */
      virtual SubMesh<Context::Serial> keep(Attribute attr) const;

      /**
      * @brief Trims the elements with the given attributes.
      * @param[in] attrs Attributes to trim
      * @returns SubMesh object to the remaining region of the mesh
      *
      * This function will trim keep only the elements that have an attribute
      * in the given set of attributes.
      *
      * The lower dimensional polytopes of dimension @f$ 1 \leq d \leq D - 1
      * @f$ are included if the connectivity @f$ D \longrightarrow d @f$ is
      * already computed in the mesh.
      *
      * @returns A SubMesh object consisting of elements that have attributes
      * not in the given set.
      */
      virtual SubMesh<Context::Serial> keep(const FlatSet<Attribute>& attrs) const;

      /**
      * @brief Loads a mesh from file in the given format.
      * @param[in] filename Name of file to read
      * @param[in] fmt Mesh file format
      * @returns Reference to this (for method chaining)
      */
      virtual Mesh& load(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM) override;

      /**
      * @brief Saves a mesh to file in the given format.
      * @param[in] filename Name of file to write
      * @param[in] fmt Mesh file format
      * @returns Reference to this (for method chaining)
      */
      virtual void save(
        const boost::filesystem::path& filename,
        IO::FileFormat fmt = IO::FileFormat::MFEM, size_t precison = 16) const override;

      virtual Mesh& scale(Scalar c) override;

      virtual Mesh& setAttribute(const std::pair<size_t, Index>&, Attribute attr) override;

      virtual size_t getCount(size_t dim) const override;

      virtual size_t getCount(Polytope::Geometry g) const override;

      virtual FaceIterator getBoundary() const override;

      virtual FaceIterator getInterface() const override;

      virtual ElementIterator getElement(Index idx = 0) const override;

      virtual FaceIterator getFace(Index idx = 0) const override;

      virtual VertexIterator getVertex(Index idx = 0) const override;

      virtual PolytopeIterator getPolytope(size_t dimension, Index idx = 0) const override;

      virtual bool isSubMesh() const override
      {
        return false;
      }

      virtual bool isInterface(Index faceIdx) const override;

      virtual bool isBoundary(Index faceIdx) const override;

      virtual size_t getDimension() const override;

      virtual size_t getSpaceDimension() const override;

      virtual const PolytopeTransformation& getPolytopeTransformation(
          size_t dimension, Index idx) const override;

      virtual Polytope::Geometry getGeometry(size_t dimension, Index idx) const override;

      virtual Attribute getAttribute(size_t dimension, Index index) const override;

      virtual MeshConnectivity& getConnectivity() override
      {
        return m_connectivity;
      }

      virtual const MeshConnectivity& getConnectivity() const override
      {
        return m_connectivity;
      }

      virtual void flush() override
      {
        m_transformationIndex.clear();
      }

      virtual Eigen::Map<const Math::SpatialVector> getVertexCoordinates(Index idx) const override;

      virtual const FlatSet<Attribute>& getAttributes(size_t d) const override;

    private:
      static const GeometryIndexed<Math::Matrix> s_vertices;

      size_t m_sdim;

      Math::Matrix m_vertices;
      MeshConnectivity m_connectivity;

      PolytopeIndexed<Geometry::Attribute> m_attributeIndex;
      mutable PolytopeIndexed<std::unique_ptr<PolytopeTransformation>> m_transformationIndex;

      std::vector<FlatSet<Attribute>> m_attributes;
  };

  using SerialMesh = Mesh<Context::Serial>;
}

#include "Mesh.hpp"

#endif
