/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_GEOMETRY_MESH_H
#define RODIN_KOKKOS_GEOMETRY_MESH_H

/**
 * @file
 * @brief Resident Kokkos specialization of MeshBase.
 */

#include <cmath>
#include <memory>
#include <utility>
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/IdentityTransformation.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Connectivity.h"

namespace Rodin::Geometry
{
  /**
   * @brief Local mesh whose authoritative data resides in Kokkos memory.
   * @tparam Backend Kokkos execution backend, including the CUDA tag.
   * @ingroup MeshTypes
   *
   * ## Architecture
   *
   * This specialization implements MeshBase directly. Connectivity, vertices,
   * attributes and geometric maps are owned Kokkos arrays in the context's
   * memory space. Host meshes are temporary construction and file-format
   * adapters only; no LocalMesh is retained. There is no upload lifecycle.
   * Mutators execute on the selected execution instance and finish before
   * returning. Device kernels borrow getDeviceView(), which contains no host
   * objects. Coordinate changes invalidate host evaluation caches using the
   * same flush semantics as the local mesh. Copies deep-copy resident arrays.
   *
   * Legacy host iterators and reference-returning queries use the resident
   * data through explicit transfers and the existing evaluation caches.
   * Connectivity computation discovers entities on the device. Attribute and
   * map tables grow lazily in resident memory, preserving existing local IDs.
   * Identity and real P1 maps on all Rodin reference
   * families are supported. Higher-order and custom maps are rejected.
   * Geometry must remain unchanged during evaluation. Borrowed views expire
   * when the owner is destroyed, moved from, or changed by mutation, load or assignment.
   * The application initializes Kokkos before construction and finalizes it
   * after all meshes and kernels have finished.
   */
  template <class Backend>
  class Mesh<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>> final
    : public MeshBase
  {
    public:
      /// Host polymorphic interface.
      using Parent = MeshBase;

      /// Distribution and execution context.
      using ContextType = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

      /// Context alias used by geometry traits.
      using Context = ContextType;

      /// Execution-space type.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// Resident memory space.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// Resident connectivity specialization.
      using ConnectivityType = Connectivity<ContextType>;

      /// Borrowed device mesh handle.
      using DeviceView = Rodin::Kokkos::Geometry::MeshView;

      /**
       * @brief Constructs resident meshes with the existing local builder vocabulary.
       *
       * Construction inputs are host data, discarded by finalize(). The returned
       * mesh owns only resident arrays and small metadata, independently of the builder.
       */
      class Builder
      {
        public:
          /// @param context Context of the resulting mesh.
          explicit Builder(ContextType context = ContextType())
            : m_context(std::move(context))
          {}

          /// @param dimension Embedding dimension. @returns This builder.
          Builder& initialize(size_t dimension)
          {
            m_builder.initialize(dimension);
            return *this;
          }

          /// @param count Vertex count. @returns This builder.
          Builder& nodes(size_t count)
          {
            m_builder.nodes(count);
            return *this;
          }

          /// @param d Entity dimension. @param count Expected count. @returns This builder.
          Builder& reserve(size_t d, size_t count)
          {
            m_builder.reserve(d, count);
            return *this;
          }

          /// @param point Vertex coordinates. @returns This builder.
          Builder& vertex(std::initializer_list<Real> point)
          {
            m_builder.vertex(point);
            return *this;
          }

          /// @param point Vertex coordinates. @returns This builder.
          Builder& vertex(const Math::SpatialPoint& point)
          {
            m_builder.vertex(point);
            return *this;
          }

          /// @param geometry Reference family. @param vertices Ordered vertices. @returns This builder.
          Builder& polytope(
            Polytope::Type geometry, std::initializer_list<Index> vertices)
          {
            m_builder.polytope(geometry, vertices);
            return *this;
          }

          /// @param geometry Reference family. @param vertices Ordered vertices. @returns This builder.
          Builder& polytope(Polytope::Type geometry, const IndexArray& vertices)
          {
            m_builder.polytope(geometry, vertices);
            return *this;
          }

          /// @param entity Dimension and index. @param attribute Optional label. @returns This builder.
          Builder& attribute(
            const std::pair<size_t, Index>& entity, Optional<Attribute> attribute)
          {
            m_builder.attribute(entity, attribute);
            return *this;
          }

          /// @returns Independently owned resident mesh.
          Mesh finalize()
          {
            return Mesh(m_context, m_builder.finalize());
          }

        private:
          ContextType m_context;
          typename Mesh<Rodin::Context::Local>::Builder m_builder;
      };

      /// @param context Execution association; Kokkos must already be initialized.
      explicit Mesh(ContextType context = ContextType())
        : m_context(std::move(context)),
          m_connectivity(m_context)
      {
        initializeCaches();
      }

      /**
       * @param context Execution association.
       * @param source Construction input; its lifetime is independent of this mesh.
       */
      Mesh(ContextType context, const MeshBase& source)
        : MeshBase(source),
          m_context(std::move(context)),
          m_connectivity(m_context, source.getConnectivity()),
          m_sdim(source.getSpaceDimension())
      {
        if (m_sdim > 3 || getDimension() > m_sdim)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid Kokkos mesh embedding dimension." << Alert::Raise;
        }
        std::vector<Real> coordinates;
        std::vector<Attribute> attributes;
        std::vector<unsigned char> present;
        std::vector<GeometryKind> kinds;
        IndexVector controlOffsets(1, 0);
        std::vector<Real> controls;
        for (Index vertex = 0; vertex < source.getVertexCount(); ++vertex)
        {
          const auto point = source.getVertexCoordinates(vertex);
          for (size_t s = 0; s < m_sdim; ++s)
          {
            requireFinite(point(s));
            coordinates.push_back(point(s));
          }
        }
        for (size_t d = 0; d <= getDimension(); ++d)
        {
          for (Index i = 0; i < getPolytopeCount(d); ++i)
          {
            const auto attribute = source.getAttribute(d, i);
            attributes.push_back(attribute.value_or(Attribute{}));
            present.push_back(attribute.has_value());
            const auto& map = source.getPolytopeTransformation(d, i);
            validateTransformation(d, source.getGeometry(d, i), map);
            const auto* parametric =
              dynamic_cast<const ParametricTransformation<Variational::RealP1Element>*>(
                &map);
            const auto& key = source.getConnectivity().getPolytope(d, i);
            bool vertexMap = parametric != nullptr;
            std::vector<Real> entityControls;
            for (Index a = 0; a < key.size(); ++a)
            {
              const auto point = source.getVertexCoordinates(key[a]);
              for (size_t s = 0; s < m_sdim; ++s)
              {
                const Real value =
                  parametric ? parametric->getPointMatrix()(s, a) : Real(0);
                requireFinite(value);
                entityControls.push_back(value);
                vertexMap = vertexMap && value == point(s);
              }
            }
            kinds.push_back(parametric
                ? (vertexMap ? GeometryKind::VertexP1 : GeometryKind::ParametricP1)
                : GeometryKind::Identity);
            if (parametric && !vertexMap)
            {
              controls.insert(
                controls.end(), entityControls.begin(), entityControls.end());
            }
            controlOffsets.push_back(
              controlOffsets.back() + (parametric && !vertexMap ? key.size() : 0));
          }
        }
        m_coordinates = upload(coordinates, "mesh coordinates");
        m_attributes = upload(attributes, "mesh attributes");
        m_attributePresent = upload(present, "mesh attribute presence");
        m_geometryKinds = upload(kinds, "mesh geometry kinds");
        m_controlOffsets = upload(controlOffsets, "mesh control offsets");
        m_controls = upload(controls, "mesh controls");
        initializeCaches();
      }

      /// @param other Mesh to deep-copy, without retaining evaluation caches.
      Mesh(const Mesh& other)
        : MeshBase(other),
          m_context(other.m_context),
          m_connectivity(other.m_connectivity),
          m_sdim(other.m_sdim)
      {
        other.synchronizeEntities();
        const auto& execution = m_context.getExecutionSpace();
        m_coordinates = Rodin::Kokkos::Internal::clone(execution, other.m_coordinates);
        m_attributes = Rodin::Kokkos::Internal::clone(execution, other.m_attributes);
        m_attributePresent =
          Rodin::Kokkos::Internal::clone(execution, other.m_attributePresent);
        m_geometryKinds =
          Rodin::Kokkos::Internal::clone(execution, other.m_geometryKinds);
        m_controlOffsets =
          Rodin::Kokkos::Internal::clone(execution, other.m_controlOffsets);
        m_controls = Rodin::Kokkos::Internal::clone(execution, other.m_controls);
        initializeCaches();
      }

      /// @param other Mesh whose arrays are transferred; the source becomes empty.
      Mesh(Mesh&& other)
        : Mesh(other.m_context)
      {
        *this = std::move(other);
      }

      /// @param other Mesh to deep-copy. @returns This mesh.
      Mesh& operator=(const Mesh& other)
      {
        if (this != &other)
        {
          Mesh next(other);
          *this = std::move(next);
        }
        return *this;
      }

      /// @param other Mesh whose arrays are transferred; the source becomes empty. @returns This mesh.
      Mesh& operator=(Mesh&& other)
      {
        if (this != &other)
        {
          other.synchronizeEntities();
          MeshBase::operator=(other);
          m_context.getExecutionSpace().fence("Rodin mesh replacement");
          m_context = other.m_context;
          m_connectivity = std::move(other.m_connectivity);
          m_sdim = other.m_sdim;
          m_coordinates = std::move(other.m_coordinates);
          m_attributes = std::move(other.m_attributes);
          m_attributePresent = std::move(other.m_attributePresent);
          m_geometryKinds = std::move(other.m_geometryKinds);
          m_controlOffsets = std::move(other.m_controlOffsets);
          m_controls = std::move(other.m_controls);
          m_transformations.clear();
          m_quadratures.clear();
          initializeCaches();
          other.m_connectivity = ConnectivityType(other.m_context);
          other.m_sdim = 0;
          other.m_transformations.clear();
          other.m_quadratures.clear();
          other.initializeCaches();
        }
        return *this;
      }

      /// @returns Associated execution context.
      const ContextType& getContext() const override
      {
        return m_context;
      }

      /// @returns Topological dimension.
      size_t getDimension() const override
      {
        return m_connectivity.getDimension();
      }

      /// @returns Embedding dimension.
      size_t getSpaceDimension() const override
      {
        return m_sdim;
      }

      /// @returns Resident connectivity owner.
      ConnectivityType& getConnectivity() override
      {
        return m_connectivity;
      }

      /// @returns Resident connectivity owner.
      const ConnectivityType& getConnectivity() const override
      {
        return m_connectivity;
      }

      /// @returns Borrowed resident data; kernels must finish before owner replacement.
      DeviceView getDeviceView() const
      {
        synchronizeEntities();
        const auto topology = m_connectivity.getDeviceView();
        return {getDimension(), m_sdim, topology.entityOffsets, topology.geometries,
          topology.vertexOffsets, topology.vertexIndices, m_coordinates.data(),
          m_attributes.data(), m_attributePresent.data(), m_geometryKinds.data(),
          m_controlOffsets.data(), m_controls.data(), topology.incidenceDescriptors};
      }

      /// @param d Entity dimension. @returns Entity count.
      size_t getPolytopeCount(size_t d) const override
      {
        return m_connectivity.getCount(d);
      }

      /// @param g Reference family. @returns Entity count.
      size_t getPolytopeCount(Polytope::Type g) const override
      {
        return m_connectivity.getCount(g);
      }

      /// @param d Entity dimension. @param i Entity index. @returns Reference family.
      Polytope::Type getGeometry(size_t d, Index i) const override
      {
        return m_connectivity.getGeometry(d, i);
      }

      /// @returns False; this specialization owns a complete local mesh.
      bool isSubMesh() const override
      {
        return false;
      }

      /// @returns No value; a complete mesh cannot be cast to a submesh. @throws Alert::Exception Always.
      SubMeshBase& asSubMesh() override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      /// @returns No value; a complete mesh cannot be cast to a submesh. @throws Alert::Exception Always.
      const SubMeshBase& asSubMesh() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      /// @param factor Coordinate scale. @returns This mesh after resident execution finishes.
      Mesh& scale(Real factor) override
      {
        requireFinite(factor);
        const auto coordinates = m_coordinates;
        ::Kokkos::parallel_for(
          "Rodin mesh scale", policy(coordinates.extent(0)),
          KOKKOS_LAMBDA(Index i) { coordinates(i) *= factor; });
        flush();
        return *this;
      }

      /**
       * @brief Invalidates evaluation caches and restores vertex-based default maps.
       * Resident updates finish before returning, as for the local mesh's geometry edits.
       */
      void flush() override
      {
        synchronizeEntities();
        ::Kokkos::deep_copy(
          m_context.getExecutionSpace(), m_geometryKinds, GeometryKind::VertexP1);
        ::Kokkos::deep_copy(m_context.getExecutionSpace(), m_controlOffsets, Index(0));
        m_context.getExecutionSpace().fence("Rodin geometry update");
        m_controls = DeviceArray<Real>();
        m_transformations.clear();
        m_quadratures.clear();
      }

      /// @param vertex Vertex index. @returns Coordinates copied from resident memory.
      Math::SpatialPoint getVertexCoordinates(Index vertex) const override
      {
        assert(vertex < getVertexCount());
        const auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
          ::Kokkos::subview(
            m_coordinates, std::make_pair(vertex * m_sdim, (vertex + 1) * m_sdim)));
        Math::SpatialPoint point(m_sdim);
        for (size_t s = 0; s < m_sdim; ++s)
          point(s) = host(s);
        return point;
      }

      /// @param vertex Vertex index. @param value Coordinate. @param component Coordinate index. @returns This mesh.
      Mesh& setVertexCoordinates(Index vertex, Real value, size_t component) override
      {
        assert(vertex < getVertexCount() && component < m_sdim);
        requireFinite(value);
        ::Kokkos::deep_copy(m_context.getExecutionSpace(),
          ::Kokkos::subview(m_coordinates, vertex * m_sdim + component), value);
        flush();
        return *this;
      }

      /// @param vertex Vertex index. @param point New coordinates. @returns This mesh.
      Mesh& setVertexCoordinates(Index vertex, const Math::SpatialPoint& point) override
      {
        assert(vertex < getVertexCount() && point.size() == m_sdim);
        for (size_t s = 0; s < m_sdim; ++s)
          requireFinite(point(s));
        Rodin::Kokkos::Geometry::Vector value(m_sdim);
        for (size_t s = 0; s < m_sdim; ++s)
          value(s) = point(s);
        const auto coordinates = m_coordinates;
        const size_t sdim = m_sdim;
        ::Kokkos::parallel_for(
          "Rodin vertex update", policy(sdim),
          KOKKOS_LAMBDA(Index s) { coordinates(vertex * sdim + s) = value(s); });
        flush();
        return *this;
      }

      /// @param d Entity dimension. @param i Entity index. @returns Optional resident label.
      Optional<Attribute> getAttribute(size_t d, Index i) const override
      {
        synchronizeEntities();
        assert(i < getPolytopeCount(d));
        const auto entity = m_connectivity.getOffset(d) + i;
        if (!Rodin::Kokkos::Internal::read(
              m_context.getExecutionSpace(), m_attributePresent, entity))
          return {};
        return Rodin::Kokkos::Internal::read(
          m_context.getExecutionSpace(), m_attributes, entity);
      }

      /// @param entity Dimension and index. @param attribute Optional label. @returns This mesh.
      Mesh& setAttribute(const std::pair<size_t, Index>& entity,
        const Optional<Attribute>& attribute) override
      {
        synchronizeEntities();
        assert(entity.second < getPolytopeCount(entity.first));
        const Index index = m_connectivity.getOffset(entity.first) + entity.second;
        const auto& execution = m_context.getExecutionSpace();
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(m_attributes, index),
          attribute.value_or(Attribute{}));
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(m_attributePresent, index),
          static_cast<unsigned char>(attribute.has_value()));
        execution.fence("Rodin attribute update");
        return *this;
      }

      /// @param entity Dimension and index. @param transformation Owned map, or null for the default map. @returns This mesh.
      Mesh& setPolytopeTransformation(const std::pair<size_t, Index> entity,
        PolytopeTransformation* transformation) override
      {
        synchronizeEntities();
        std::unique_ptr<PolytopeTransformation> map(transformation);
        assert(entity.second < getPolytopeCount(entity.first));
        const Index index = m_connectivity.getOffset(entity.first) + entity.second;
        GeometryKind kind = GeometryKind::VertexP1;
        std::vector<Real> values;
        Index count = 0;
        if (map)
        {
          validateTransformation(
            entity.first, getGeometry(entity.first, entity.second), *map);
          if (const auto* parametric =
                dynamic_cast<const ParametricTransformation<Variational::RealP1Element>*>(
                  map.get()))
          {
            const auto& points = parametric->getPointMatrix();
            count = points.getCount();
            for (Index a = 0; a < count; ++a)
            {
              for (size_t s = 0; s < m_sdim; ++s)
              {
                requireFinite(points(s, a));
                values.push_back(points(s, a));
              }
            }
            kind = GeometryKind::ParametricP1;
          }
          else
            kind = GeometryKind::Identity;
        }
        const auto& execution = m_context.getExecutionSpace();
        const Index first =
          Rodin::Kokkos::Internal::read(execution, m_controlOffsets, index);
        const Index end =
          Rodin::Kokkos::Internal::read(execution, m_controlOffsets, index + 1);
        const Index oldCount = end - first;
        const Index total = Rodin::Kokkos::Internal::read(
          execution, m_controlOffsets, m_controlOffsets.extent(0) - 1);
        const auto temporary = upload(values, "new P1 controls");
        // The compact table stores independent maps only. Default maps borrow vertices.
        DeviceArray<Index> offsets(
          ::Kokkos::view_alloc(execution, std::string("mesh control offsets")),
          m_controlOffsets.extent(0));
        DeviceArray<Real> controls(
          ::Kokkos::view_alloc(execution, std::string("mesh controls")),
          (total - oldCount + count) * m_sdim);
        const auto previous = m_controlOffsets;
        ::Kokkos::parallel_for(
          "Rodin control offsets", policy(offsets.extent(0)), KOKKOS_LAMBDA(Index i) {
            offsets(i) = i > index ? previous(i) - oldCount + count : previous(i);
          });
        ::Kokkos::deep_copy(execution,
          ::Kokkos::subview(controls, std::make_pair(Index(0), first * m_sdim)),
          ::Kokkos::subview(m_controls, std::make_pair(Index(0), first * m_sdim)));
        ::Kokkos::deep_copy(execution,
          ::Kokkos::subview(
            controls, std::make_pair(first * m_sdim, (first + count) * m_sdim)),
          temporary);
        ::Kokkos::deep_copy(execution,
          ::Kokkos::subview(
            controls, std::make_pair((first + count) * m_sdim, controls.extent(0))),
          ::Kokkos::subview(
            m_controls, std::make_pair(end * m_sdim, m_controls.extent(0))));
        execution.fence("Rodin control table replacement");
        m_controlOffsets = std::move(offsets);
        m_controls = std::move(controls);
        ::Kokkos::deep_copy(execution, ::Kokkos::subview(m_geometryKinds, index), kind);
        execution.fence("Rodin transformation update");
        m_transformations.clear();
        m_quadratures.clear();
        return *this;
      }

      /// @param d Entity dimension. @param i Entity index. @returns Cached host adapter reconstructed from resident geometry.
      const PolytopeTransformation& getPolytopeTransformation(
        size_t d, Index i) const override
      {
        synchronizeEntities();
        return m_transformations.get({d, i}, getPolytopeCount(d),
          [this](
            size_t dimension, Index index) -> std::unique_ptr<PolytopeTransformation> {
            const Index entity = m_connectivity.getOffset(dimension) + index;
            const auto kind = Rodin::Kokkos::Internal::read(
              m_context.getExecutionSpace(), m_geometryKinds, entity);
            if (kind == GeometryKind::Identity)
              return std::make_unique<IdentityTransformation>(dimension);
            const auto& vertices = m_connectivity.getPolytope(dimension, index);
            PointCloud points(m_sdim, vertices.size());
            if (kind == GeometryKind::VertexP1)
            {
              for (Index a = 0; a < vertices.size(); ++a)
              {
                for (size_t s = 0; s < m_sdim; ++s)
                {
                  points(s, a) =
                    Rodin::Kokkos::Internal::read(m_context.getExecutionSpace(),
                      m_coordinates, vertices[a] * m_sdim + s);
                }
              }
            }
            else
            {
              const auto first =
                Rodin::Kokkos::Internal::read(
                  m_context.getExecutionSpace(), m_controlOffsets, entity) *
                m_sdim;
              const auto host =
                ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
                  ::Kokkos::subview(
                    m_controls, std::make_pair(first, first + vertices.size() * m_sdim)));
              for (Index a = 0; a < vertices.size(); ++a)
              {
                for (size_t s = 0; s < m_sdim; ++s)
                  points(s, a) = host(a * m_sdim + s);
              }
            }
            return std::make_unique<ParametricTransformation<Variational::RealP1Element>>(
              std::move(points),
              Variational::RealP1Element(getGeometry(dimension, index)));
          });
      }

      /// @param d Entity dimension. @param i Entity index. @param formula Persistent reference quadrature. @returns Host point cache attached to this mesh.
      const PolytopeQuadrature& getQuadrature(
        size_t d, Index i, const QF::QuadratureFormula& formula) const override
      {
        return m_quadratures.get(
          {d, i}, getPolytopeCount(d), formula, [this, d, i, &formula]() {
            return std::make_unique<PolytopeQuadrature>(Polytope(d, i, *this), formula);
          });
      }

      /// @param path Input path. @param format Existing Rodin file format. @returns This mesh, replaced only after successful construction.
      Mesh& load(const boost::filesystem::path& path, IO::FileFormat format) override
      {
        Mesh<Rodin::Context::Local> input;
        input.load(path, format);
        Mesh next(m_context, input);
        *this = std::move(next);
        return *this;
      }

      /// @param path Output path. @param format Existing Rodin file format.
      void save(const boost::filesystem::path& path, IO::FileFormat format) const override
      {
        typename Mesh<Rodin::Context::Local>::Builder builder;
        builder.initialize(m_sdim).nodes(getVertexCount());
        for (Index v = 0; v < getVertexCount(); ++v)
          builder.vertex(getVertexCoordinates(v));
        for (size_t d = 1; d <= getDimension(); ++d)
        {
          for (Index i = 0; i < getPolytopeCount(d); ++i)
          {
            const auto& key = m_connectivity.getPolytope(d, i);
            IndexArray vertices(key.size());
            for (Index a = 0; a < key.size(); ++a)
              vertices(a) = key[a];
            builder.polytope(getGeometry(d, i), vertices);
          }
        }
        auto output = builder.finalize();
        for (size_t d = 0; d <= getDimension(); ++d)
        {
          for (Index i = 0; i < getPolytopeCount(d); ++i)
          {
            output.setAttribute({d, i}, getAttribute(d, i));
            output.setPolytopeTransformation(
              {d, i}, getPolytopeTransformation(d, i).copy());
          }
        }
        for (size_t d = 0; d <= getDimension(); ++d)
        {
          for (size_t dp = 0; dp <= getDimension(); ++dp)
          {
            const auto& relation = m_connectivity.getIncidence(d, dp);
            if (!relation.empty())
              output.getConnectivity().setIncidence({d, dp}, Incidence(relation));
          }
        }
        if (const auto name = getName())
          output.setName(std::string(*name));
        output.save(path, format);
      }

      /**
       * @brief Traverses boundary faces.
       * @returns Host iterator attached to this mesh.
       */
      FaceIterator getBoundary() const override
      {
        std::vector<Index> indices;
        const size_t count = getDimension() ? getFaceCount() : 0;
        for (Index i = 0; i < count; i++)
        {
          if (isBoundary(i))
            indices.push_back(i);
        }
        if (indices.size() == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Mesh has an empty boundary." << Alert::Raise;
        }
        return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
      }

      /**
       * @brief Traverses interface faces.
       * @returns Host iterator attached to this mesh.
       */
      FaceIterator getInterface() const override
      {
        std::vector<Index> indices;
        const size_t count = getDimension() ? getFaceCount() : 0;
        for (Index i = 0; i < count; i++)
        {
          if (isInterface(i))
            indices.push_back(i);
        }
        if (indices.size() == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Mesh has an empty interface." << Alert::Raise;
        }
        return FaceIterator(*this, VectorIndexGenerator(std::move(indices)));
      }

      /**
       * @brief Traverses cells.
       * @returns Host iterator attached to this mesh.
       */
      CellIterator getCell() const override
      {
        return getCell(0);
      }

      /**
       * @brief Traverses faces.
       * @returns Host iterator attached to this mesh.
       */
      FaceIterator getFace() const override
      {
        return getFace(0);
      }

      /**
       * @brief Traverses vertices.
       * @returns Host iterator attached to this mesh.
       */
      VertexIterator getVertex() const override
      {
        return getVertex(0);
      }

      /**
       * @brief Traverses polytopes.
       * @param dimension Entity dimension.
       * @returns Host iterator attached to this mesh.
       */
      PolytopeIterator getPolytope(size_t dimension) const override
      {
        return getPolytope(dimension, 0);
      }

      /**
       * @brief Traverses cells.
       * @param idx First entity index.
       * @returns Host iterator attached to this mesh.
       */
      CellIterator getCell(Index idx) const override
      {
        return CellIterator(*this, BoundedIndexGenerator(idx, this->getCellCount()));
      }

      /**
       * @brief Traverses faces.
       * @param idx First entity index.
       * @returns Host iterator attached to this mesh.
       */
      FaceIterator getFace(Index idx) const override
      {
        return FaceIterator(*this, BoundedIndexGenerator(idx, this->getFaceCount()));
      }

      /**
       * @brief Traverses vertices.
       * @param idx First entity index.
       * @returns Host iterator attached to this mesh.
       */
      VertexIterator getVertex(Index idx) const override
      {
        return VertexIterator(*this, BoundedIndexGenerator(idx, this->getVertexCount()));
      }

      /**
       * @brief Traverses polytopes.
       * @param dimension Entity dimension.
       * @param idx First entity index.
       * @returns Host iterator attached to this mesh.
       */
      PolytopeIterator getPolytope(size_t dimension, Index idx) const override
      {
        return PolytopeIterator(dimension, *this,
          BoundedIndexGenerator(idx, this->getPolytopeCount(dimension)));
      }

      /// @param face Face index. @returns Whether it meets more than one cell.
      bool isInterface(Index face) const override
      {
        return getIncidentCellCount(face) > 1;
      }

      /// @param face Face index. @returns Whether it meets exactly one cell.
      bool isBoundary(Index face) const override
      {
        return getIncidentCellCount(face) == 1;
      }

      /**
       * @brief Computes the volume by resident quadrature reduction.
       * @returns Integrated geometric measure.
       */
      Real getVolume() const override
      {
        return measure(3, nullptr, false);
      }

      /**
       * @brief Computes the volume by resident quadrature reduction.
       * @param attribute Selected material label.
       * @returns Integrated geometric measure.
       */
      Real getVolume(Attribute attribute) const override
      {
        const FlatSet<Attribute> selected{attribute};
        return measure(3, &selected, false);
      }

      /**
       * @brief Computes the volume by resident quadrature reduction.
       * @param attributes Selected material labels.
       * @returns Integrated geometric measure.
       */
      Real getVolume(const FlatSet<Attribute>& attributes) const override
      {
        return measure(3, &attributes, false);
      }

      /**
       * @brief Computes the area by resident quadrature reduction.
       * @returns Integrated geometric measure.
       */
      Real getArea() const override
      {
        return measure(2, nullptr, false);
      }

      /**
       * @brief Computes the area by resident quadrature reduction.
       * @param attribute Selected material label.
       * @returns Integrated geometric measure.
       */
      Real getArea(Attribute attribute) const override
      {
        const FlatSet<Attribute> selected{attribute};
        return measure(2, &selected, false);
      }

      /**
       * @brief Computes the area by resident quadrature reduction.
       * @param attributes Selected material labels.
       * @returns Integrated geometric measure.
       */
      Real getArea(const FlatSet<Attribute>& attributes) const override
      {
        return measure(2, &attributes, false);
      }

      /**
       * @brief Computes the measure by resident quadrature reduction.
       * @param d Entity dimension.
       * @returns Integrated geometric measure.
       */
      Real getMeasure(size_t d) const override
      {
        return measure(d, nullptr, false);
      }

      /**
       * @brief Computes the measure by resident quadrature reduction.
       * @param d Entity dimension.
       * @param attribute Selected material label.
       * @returns Integrated geometric measure.
       */
      Real getMeasure(size_t d, Attribute attribute) const override
      {
        const FlatSet<Attribute> selected{attribute};
        return measure(d, &selected, false);
      }

      /**
       * @brief Computes the measure by resident quadrature reduction.
       * @param d Entity dimension.
       * @param attributes Selected material labels.
       * @returns Integrated geometric measure.
       */
      Real getMeasure(size_t d, const FlatSet<Attribute>& attributes) const override
      {
        return measure(d, &attributes, false);
      }

      /**
       * @brief Computes the perimeter by resident quadrature reduction.
       * @returns Integrated geometric measure.
       */
      Real getPerimeter() const override
      {
        return measure(getDimension() ? getDimension() - 1 : 0, nullptr, true);
      }

      /**
       * @brief Computes the perimeter by resident quadrature reduction.
       * @param attribute Selected material label.
       * @returns Integrated geometric measure.
       */
      Real getPerimeter(Attribute attribute) const override
      {
        const FlatSet<Attribute> selected{attribute};
        return measure(getDimension() ? getDimension() - 1 : 0, &selected, true);
      }

      /**
       * @brief Computes the perimeter by resident quadrature reduction.
       * @param attributes Selected material labels.
       * @returns Integrated geometric measure.
       */
      Real getPerimeter(const FlatSet<Attribute>& attributes) const override
      {
        return measure(getDimension() ? getDimension() - 1 : 0, &attributes, true);
      }

    private:
      using GeometryKind = Rodin::Kokkos::Geometry::GeometryKind;
      template <class T>
      using DeviceArray = ::Kokkos::View<T*, MemorySpace>;

      /**
       * @brief One resident quadrature reduction over entities of a reference family.
       * Named functors permit use from private host methods with NVCC.
       */
      struct Measure
      {
        DeviceView view; ///< Resident mesh data.
        Polytope::Type geometry; ///< Selected reference family.
        DeviceArray<Attribute> labels; ///< Selected labels, possibly empty.
        DeviceArray<Real> references; ///< Packed reference samples.
        DeviceArray<Real> quadratureWeights; ///< Reference weights.
        size_t dimension; ///< Entity dimension.
        bool filter; ///< Whether a label selection was requested.
        bool boundary; ///< Whether only boundary entities contribute.

        /// @param i Entity index. @param sum Reduction accumulator.
        KOKKOS_INLINE_FUNCTION
        void operator()(Index i, Real& sum) const
        {
          if (view.getConnectivity().getGeometry(dimension, i) != geometry)
            return;
          if (boundary &&
            view.getConnectivity().getIncidence(dimension, dimension + 1, i).size() !=
              1)
            return;
          const Index entity = view.entityOffsets[dimension] + i;
          if (filter)
          {
            bool selected = false;
            if (view.attributePresent[entity])
            {
              for (Index a = 0; a < labels.extent(0); ++a)
                selected = selected || labels(a) == view.attributes[entity];
            }
            if (!selected)
              return;
          }
          for (Index q = 0; q < quadratureWeights.extent(0); ++q)
          {
            Rodin::Kokkos::Geometry::Vector reference(dimension);
            for (size_t s = 0; s < dimension; ++s)
              reference(s) = references(q * dimension + s);
            sum += quadratureWeights(q) *
              view.getPolytope(dimension, i).getPoint(reference).getDistortion();
          }
        }
      };

      /// @brief Grows entity data after device discovery without changing existing maps.
      struct GrowEntities
      {
        ::Kokkos::Array<Index, 5> oldOffsets,
          newOffsets; ///< Dimension offsets before and after discovery.
        DeviceArray<Attribute> oldAttributes,
          attributes; ///< Existing and replacement labels.
        DeviceArray<unsigned char> oldPresent, present; ///< Explicit label presence.
        DeviceArray<GeometryKind> oldKinds, kinds; ///< Existing and replacement maps.
        DeviceArray<Index> oldControls,
          controls; ///< Existing and replacement control-point offsets.
        size_t dimension; ///< Topological dimension.

        /// @param i New flattened entity index, including the terminal control offset.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          if (i == newOffsets[dimension + 1])
          {
            controls(i) =
              oldControls.extent(0) ? oldControls(oldOffsets[dimension + 1]) : 0;
            return;
          }
          size_t d = 0;
          while (i >= newOffsets[d + 1])
            ++d;
          const Index local = i - newOffsets[d];
          const bool existing = local < oldOffsets[d + 1] - oldOffsets[d];
          if (existing)
          {
            const Index old = oldOffsets[d] + local;
            attributes(i) = oldAttributes(old);
            present(i) = oldPresent(old);
            kinds(i) = oldKinds(old);
            controls(i) = oldControls(old);
          }
          else
          {
            present(i) = 0;
            kinds(i) = GeometryKind::VertexP1;
            controls(i) = oldControls.extent(0) ? oldControls(oldOffsets[d + 1]) : 0;
          }
        }
      };

      /**
       * @brief Extends resident entity slots on first use after connectivity discovery.
       * Like the core attribute and transformation indices, const queries may
       * materialize newly required slots. Actual attributes and geometry stay
       * resident; only the small dimension-count metadata is examined on the host.
       */
      void synchronizeEntities() const
      {
        std::lock_guard<std::mutex> lock(m_entityMutex);
        ::Kokkos::Array<Index, 5> offsets{};
        bool changed = false;
        for (size_t d = 0; d <= getDimension() + 1; ++d)
        {
          offsets[d] = m_connectivity.getOffset(d);
          changed = changed || offsets[d] != m_layoutOffsets[d];
        }
        if (!changed)
          return;
        const auto& execution = m_context.getExecutionSpace();
        const Index total = offsets[getDimension() + 1];
        GrowEntities work;
        work.oldOffsets = m_layoutOffsets;
        work.newOffsets = offsets;
        work.oldAttributes = m_attributes;
        work.oldPresent = m_attributePresent;
        work.oldKinds = m_geometryKinds;
        work.oldControls = m_controlOffsets;
        work.dimension = getDimension();
        work.attributes = DeviceArray<Attribute>(
          ::Kokkos::view_alloc(execution, std::string("mesh attributes")), total);
        work.present = DeviceArray<unsigned char>(
          ::Kokkos::view_alloc(execution, std::string("mesh attribute presence")), total);
        work.kinds = DeviceArray<GeometryKind>(
          ::Kokkos::view_alloc(execution, std::string("mesh geometry kinds")), total);
        work.controls = DeviceArray<Index>(
          ::Kokkos::view_alloc(execution, std::string("mesh control offsets")),
          total + 1);
        ::Kokkos::parallel_for("Rodin discovered entity slots", policy(total + 1), work);
        execution.fence("Rodin discovered mesh data");
        m_attributes = std::move(work.attributes);
        m_attributePresent = std::move(work.present);
        m_geometryKinds = std::move(work.kinds);
        m_controlOffsets = std::move(work.controls);
        m_layoutOffsets = offsets;
      }

      /// @param count Iteration count. @returns Policy on this mesh's execution instance.
      auto policy(Index count) const
      {
        return ::Kokkos::RangePolicy<ExecutionSpace>(
          m_context.getExecutionSpace(), 0, count);
      }

      /**
       * @tparam T Entry type.
       * @param values Temporary construction values.
       * @param label Allocation label.
       * @returns Independent resident array.
       */
      template <class T>
      auto upload(const std::vector<T>& values, const char* label) const
      {
        return Rodin::Kokkos::Internal::upload(
          m_context.getExecutionSpace(), values, label);
      }

      /// @brief Initializes host compatibility cache buckets.
      void initializeCaches()
      {
        m_transformations.initialize(getDimension());
        m_quadratures.initialize(getDimension());
        m_layoutOffsets = {};
        for (size_t d = 0; d <= getDimension() + 1; ++d)
          m_layoutOffsets[d] = m_connectivity.getOffset(d);
      }

      /// @param value Coordinate or scale to validate.
      void requireFinite(Real value) const
      {
        if (!std::isfinite(value))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "A mesh coordinate or scale is nonfinite." << Alert::Raise;
        }
      }

      /**
       * @param d Reference dimension.
       * @param geometry Reference family.
       * @param map Transformation to validate before any device update.
       */
      void validateTransformation(
        size_t d, Polytope::Type geometry, const PolytopeTransformation& map) const
      {
        if (map.getReferenceDimension() != d || map.getPhysicalDimension() != m_sdim)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Inconsistent geometry-map dimensions." << Alert::Raise;
        }
        if (!dynamic_cast<const IdentityTransformation*>(&map) &&
          !dynamic_cast<const ParametricTransformation<Variational::RealP1Element>*>(
            &map))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos geometry supports identity and real P1 transformations."
            << Alert::Raise;
        }
        if (const auto* parametric =
              dynamic_cast<const ParametricTransformation<Variational::RealP1Element>*>(
                &map))
        {
          if (static_cast<size_t>(parametric->getPointMatrix().rows()) != m_sdim ||
            parametric->getPointMatrix().getCount() !=
              Polytope::Traits(geometry).getVertexCount())
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid P1 geometry control points." << Alert::Raise;
          }
        }
      }

      /// @param face Face index. @returns Incident cell count for classification.
      size_t getIncidentCellCount(Index face) const
      {
        const size_t d = getDimension();
        if (!d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "A zero-dimensional mesh has no faces." << Alert::Raise;
        }
        RODIN_GEOMETRY_MESH_REQUIRE_INCIDENCE(d - 1, d);
        return m_connectivity.getIncidence({d - 1, d}, face).size();
      }

      /**
       * @param dimension Entity dimension.
       * @param attributes Optional label filter; null selects every label.
       * @param boundary Whether only boundary entities contribute.
       * @returns Resident quadrature sum copied to the host.
       */
      Real measure(
        size_t dimension, const FlatSet<Attribute>* attributes, bool boundary) const
      {
        if (dimension > getDimension() || (boundary && !getDimension()))
          return 0;
        if (boundary)
        {
          const size_t d = getDimension();
          if (!m_connectivity.hasIncidence(d - 1, d))
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Compute face-to-cell incidence before measuring the boundary."
              << Alert::Raise;
          }
        }
        const auto view = getDeviceView();
        const auto labels = upload(attributes
            ? std::vector<Attribute>(attributes->begin(), attributes->end())
            : std::vector<Attribute>{},
          "measure labels");
        const bool filter = attributes != nullptr;
        Real result = 0;
        for (size_t g = 0; g < 8; ++g)
        {
          const auto geometry = static_cast<Polytope::Type>(g);
          if (Polytope::Traits(geometry).getDimension() != dimension ||
            !getPolytopeCount(geometry))
            continue;
          const Polytope::Traits traits(geometry);
          const Variational::RealP1Element element(geometry);
          // Match Polytope::getMeasure's degree policy, including its curved-map floor.
          constexpr size_t MinimumCurvedMeasureOrder = 8;
          size_t order = std::max(MinimumCurvedMeasureOrder, 2 * element.getOrder() + 2);
          if (traits.isSimplex())
            order = 0;
          else if (dimension == m_sdim && geometry != Polytope::Type::Pyramid)
            order = dimension * element.getFactorOrder() - 1;
          const auto& formula = QF::PolytopeQuadratureFormula::get(order, geometry);
          std::vector<Real> points, weights;
          for (size_t q = 0; q < formula.getSize(); ++q)
          {
            const auto& point = formula.getPoint(q);
            for (size_t s = 0; s < dimension; ++s)
              points.push_back(point(s));
            weights.push_back(formula.getWeight(q));
          }
          const auto references = upload(points, "measure reference points");
          const auto quadratureWeights = upload(weights, "measure reference weights");
          Real subtotal = 0;
          ::Kokkos::parallel_reduce("Rodin mesh measure",
            policy(getPolytopeCount(dimension)),
            Measure{view, geometry, labels, references, quadratureWeights, dimension,
              filter, boundary},
            subtotal);
          result += subtotal;
        }
        return result;
      }

      ContextType m_context;
      ConnectivityType m_connectivity;
      size_t m_sdim = 0;
      DeviceArray<Real> m_coordinates, m_controls;
      mutable DeviceArray<Attribute> m_attributes;
      mutable DeviceArray<unsigned char> m_attributePresent;
      mutable DeviceArray<GeometryKind> m_geometryKinds;
      mutable DeviceArray<Index> m_controlOffsets;
      mutable ::Kokkos::Array<Index, 5> m_layoutOffsets{};
      mutable std::mutex m_entityMutex;
      mutable PolytopeTransformationIndex m_transformations;
      mutable PolytopeQuadratureIndex m_quadratures;
  };
}

#endif
