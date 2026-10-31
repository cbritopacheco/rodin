/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_H1_H
#define RODIN_KOKKOS_VARIATIONAL_H1_H

/**
 * @file
 * @brief Quadratic H1 Lagrange spaces on a resident Kokkos mesh.
 */

#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <string>

#include "Rodin/Variational/H1/H1.h"
#include "Rodin/Variational/H1/Div.h"
#include "Rodin/Variational/H1/Grad.h"
#include "Rodin/Variational/H1/Jacobian.h"
#include "Rodin/Variational/H1/ShapeFunction.h"

#include "Rodin/Kokkos/Geometry/Mesh.h"

#include "DOFs.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Quadratic Lagrange basis of the simplices.
   *
   * Node @f$ q @f$ of a simplex follows the node order of
   * Rodin::Variational::H1Element of degree two. It is either the vertex
   * @f$ i @f$ or the midpoint of the edge @f$ (i, j) @f$, and its basis
   * function is a product of barycentric coordinates
   * @f[
   *   \phi_q = \lambda_i (2 \lambda_i - 1) \quad \text{or} \quad
   *   \phi_q = 4 \lambda_i \lambda_j .
   * @f]
   * The barycentric coordinates of a simplex are its P1 basis, so the
   * values and physical gradients are combinations of the P1 basis cached
   * by IntegrationPointView.
   */
  class QuadraticBasis
  {
    public:
      /// @brief Maximum number of nodes, reached by the tetrahedron.
      static constexpr Index Capacity = 10;

      /**
       * @brief Gets the basis of every simplex.
       *
       * The table is built once on the host from the reference nodes of the
       * host elements.
       *
       * @return Quadratic basis of the point, segment, triangle and
       * tetrahedron.
       */
      static const QuadraticBasis& get()
      {
        static const QuadraticBasis s_basis = [] {
          using Type = ::Rodin::Geometry::Polytope::Type;
          constexpr std::array<Type, 4> simplices = {
            Type::Point, Type::Segment, Type::Triangle, Type::Tetrahedron};
          QuadraticBasis basis;
          for (size_t d = 0; d < simplices.size(); ++d)
          {
            const auto& nodes =
              ::Rodin::Variational::H1Element<2, Real>::getNodes(simplices[d]);
            assert(nodes.size() == getCount(d));
            for (size_t q = 0; q < nodes.size(); ++q)
            {
              std::array<Real, 4> lambda = {1, 0, 0, 0};
              for (size_t k = 0; k < d; ++k)
              {
                lambda[k + 1] = nodes[q](k);
                lambda[0] -= nodes[q](k);
              }
              std::array<std::uint8_t, 2> pair = {0, 0};
              size_t found = 0;
              for (size_t a = 0; a <= d; ++a)
              {
                if (lambda[a] > 0.25)
                  pair[found++] = static_cast<std::uint8_t>(a);
              }
              assert(found == 1 || found == 2);
              basis.m_first[d][q] = pair[0];
              basis.m_second[d][q] = found == 1 ? pair[0] : pair[1];
            }
          }
          return basis;
        }();
        return s_basis;
      }

      /**
       * @brief Gets the number of nodes of a simplex.
       * @param d Dimension of the simplex.
       * @return @f$ (d + 1)(d + 2) / 2 @f$.
       */
      KOKKOS_INLINE_FUNCTION
      static constexpr Index getCount(Index d)
      {
        return (d + 1) * (d + 2) / 2;
      }

      /**
       * @brief Evaluates a basis function.
       * @param ip Integration point on a simplex.
       * @param q Local node index.
       * @return Basis value @f$ \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
      {
        const Index d = ip.getPoint().getPolytope().getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        const Index i = m_first[d][q];
        const Index j = m_second[d][q];
        const Real li = ip.getBasis(i);
        if (i == j)
          return li * (2 * li - 1);
        return 4 * li * ip.getBasis(j);
      }

      /**
       * @brief Evaluates the physical gradient of a basis function.
       * @param ip Integration point on a simplex.
       * @param q Local node index.
       * @return Gradient @f$ \nabla \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(const Point& ip, Index q) const
      {
        const Index d = ip.getPoint().getPolytope().getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        const Index i = m_first[d][q];
        const Index j = m_second[d][q];
        const Real li = ip.getBasis(i);
        const auto& gi = ip.getGradient(i);
        Geometry::Vector result(gi.size());
        if (i == j)
        {
          for (Index c = 0; c < gi.size(); ++c)
            result(c) = (4 * li - 1) * gi(c);
          return result;
        }
        const Real lj = ip.getBasis(j);
        const auto& gj = ip.getGradient(j);
        for (Index c = 0; c < gi.size(); ++c)
          result(c) = 4 * (lj * gi(c) + li * gj(c));
        return result;
      }

      /**
       * @brief Gets the reference coordinates of a node.
       * @param polytope Simplex of the mesh.
       * @param q Local node index.
       * @return Vertex @f$ i @f$ or midpoint of the edge @f$ (i, j) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
      {
        const Index d = polytope.getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        Geometry::Vector result = getReferenceVertex(polytope, m_first[d][q]);
        const Geometry::Vector other = getReferenceVertex(polytope, m_second[d][q]);
        for (Index k = 0; k < result.size(); ++k)
          result(k) = (result(k) + other(k)) / 2;
        return result;
      }

    private:
      ::Kokkos::Array<::Kokkos::Array<std::uint8_t, Capacity>, 4> m_first;
      ::Kokkos::Array<::Kokkos::Array<std::uint8_t, Capacity>, 4> m_second;
  };

  /**
   * @brief Degrees of freedom of a quadratic H1 space on the device.
   *
   * Holds the local to global numbering of every polytope of a simplicial
   * mesh, copied from the host space, and the quadratic basis. Component
   * @f$ r @f$ of the scalar degree of freedom @f$ q @f$ has the global index
   * @f$ q + r n @f$, where @f$ n @f$ is the stride, the size of the scalar
   * space.
   *
   * @tparam MemorySpace Kokkos memory space of the tables.
   */
  template <class MemorySpace>
  class QuadraticDOFs
  {
    public:
      /// @brief Empty tables.
      QuadraticDOFs() = default;

      /**
       * @brief Copies the numbering of a scalar host space to the device.
       * @param fes Scalar quadratic host space.
       * @param execution Execution space ordering the transfers.
       */
      template <class FES, class ExecutionSpace>
      QuadraticDOFs(const FES& fes, const ExecutionSpace& execution)
        : m_stride(fes.getSize()),
          m_basis(QuadraticBasis::get())
      {
        const auto& mesh = fes.getMesh();
        assert(mesh.getDimension() < m_dofs.size());
        for (size_t d = 0; d <= mesh.getDimension(); ++d)
        {
          const Index count = mesh.getPolytopeCount(d);
          const Index n = QuadraticBasis::getCount(d);
          m_dofs[d] = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, std::string("Rodin H1 degrees of freedom")),
            count * n);
          auto host = ::Kokkos::create_mirror_view(m_dofs[d]);
          for (Index i = 0; i < count; ++i)
          {
            const auto& dofs = fes.getDOFs(d, i);
            assert(static_cast<Index>(dofs.size()) == n);
            for (Index q = 0; q < n; ++q)
              host(i * n + q) = dofs(q);
          }
          ::Kokkos::deep_copy(execution, m_dofs[d], host);
        }
        execution.fence("Rodin H1 degrees of freedom transfer");
      }

      /**
       * @brief Gets the number of local degrees of freedom of a polytope.
       * @param vertices Number of vertices of a simplex.
       * @return @f$ n (n + 1) / 2 @f$ for @f$ n @f$ vertices.
       */
      KOKKOS_INLINE_FUNCTION
      static constexpr Index getCapacity(Index vertices)
      {
        return vertices * (vertices + 1) / 2;
      }

      /// @brief Gets the offset between components. @return Scalar space size.
      KOKKOS_INLINE_FUNCTION
      Index getStride() const
      {
        return m_stride;
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param polytope Simplex of the mesh.
       * @return Scalar global indices in local order.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::IndexView bind(const Geometry::PolytopeView& polytope) const
      {
        const Index d = polytope.getDimension();
        const Index n = QuadraticBasis::getCount(d);
        return {m_dofs[d].data() + polytope.getIndex() * n, n};
      }

      /// @param ip Integration point. @param q Local index. @return Basis value.
      template <class Point>
      KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
      {
        return m_basis.getBasis(ip, q);
      }

      /// @param ip Integration point. @param q Local index. @return Physical gradient.
      template <class Point>
      KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(const Point& ip, Index q) const
      {
        return m_basis.getGradient(ip, q);
      }

      /// @param polytope Simplex. @param q Local index. @return Reference node.
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
      {
        return m_basis.getNode(polytope, q);
      }

    private:
      Index m_stride = 0;
      QuadraticBasis m_basis;
      ::Kokkos::Array<::Kokkos::View<Index*, MemorySpace>, 4> m_dofs;
  };
}

namespace Rodin::Kokkos::Variational::Internal
{
  /**
   * @brief Host copy of a resident mesh with the quadratic spaces built on
   * it.
   *
   * The host spaces number the degrees of freedom and answer the form
   * language queries; the scalar numbering is copied to the device.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class QuadraticMirror
  {
    public:
      /// @brief Resident mesh type.
      using DeviceMeshType = ::Rodin::Geometry::Mesh<
        ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>;

      /// @brief Host mesh type.
      using HostMeshType = ::Rodin::Geometry::Mesh<::Rodin::Context::Local>;

      /// @brief Memory space of the device tables.
      using MemorySpace = typename DeviceMeshType::MemorySpace;

      /// @brief Scalar host space.
      using ScalarSpace = ::Rodin::Variational::H1<2, Real, HostMeshType>;

      /// @brief Vector host space.
      using VectorSpace =
        ::Rodin::Variational::H1<2, Math::SpatialVector<Real>, HostMeshType>;

      /**
       * @brief Builds the host spaces and copies the numbering to the device.
       * @param mesh Resident simplicial mesh. The incidences
       * @f$ d \to d - 1 @f$ it lacks are computed on the host copy.
       * @param vdim Vector dimension, zero for a scalar space only.
       */
      QuadraticMirror(const DeviceMeshType& mesh, size_t vdim)
        : m_mesh(simplicial(mesh)),
          m_scalar(std::integral_constant<size_t, 2>{}, m_mesh),
          m_dofs(m_scalar, mesh.getContext().getExecutionSpace())
      {
        if (vdim > 0)
          m_vector.emplace(std::integral_constant<size_t, 2>{}, m_mesh, vdim);
      }

      QuadraticMirror(const QuadraticMirror&) = delete;

      QuadraticMirror& operator=(const QuadraticMirror&) = delete;

      /// @brief Gets the scalar host space. @return Scalar space.
      const ScalarSpace& getScalar() const
      {
        return m_scalar;
      }

      /// @brief Gets the vector host space. @return Vector space.
      const VectorSpace& getVector() const
      {
        assert(m_vector);
        return *m_vector;
      }

      /// @brief Gets the device numbering. @return Scalar degrees of freedom.
      const QuadraticDOFs<MemorySpace>& getDOFs() const
      {
        return m_dofs;
      }

    private:
      /**
       * @brief Copies a mesh to the host after checking its geometries.
       * @param mesh Resident mesh.
       * @return Host copy of @p mesh, with the incidences @f$ d \to d - 1 @f$.
       */
      static HostMeshType simplicial(const DeviceMeshType& mesh)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        for (size_t d = 1; d <= mesh.getDimension(); ++d)
        {
          for (Index i = 0; i < mesh.getPolytopeCount(d); ++i)
          {
            const Type g = mesh.getGeometry(d, i);
            if (g != Type::Segment && g != Type::Triangle && g != Type::Tetrahedron)
            {
              Alert::Exception()
                << "Kokkos H1 spaces of degree 2 require a simplicial mesh, found " << g
                << " at (" << d << ", " << i << ")." << Alert::Raise;
            }
          }
        }
        // The incidences missing on the device only number the host
        // closure, the kernels never bind the polytopes they create.
        HostMeshType host = mesh.toHost();
        for (size_t d = host.getDimension(); d > 1; --d)
        {
          if (host.getConnectivity().getIncidence(d, d - 1).size() == 0)
            host.getConnectivity().compute(d, d - 1);
        }
        return host;
      }

      HostMeshType m_mesh;
      ScalarSpace m_scalar;
      std::optional<VectorSpace> m_vector;
      QuadraticDOFs<MemorySpace> m_dofs;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup H1Specializations
   * @brief Real valued quadratic H1 Lagrange space on a resident simplicial
   * Kokkos mesh.
   *
   * The space numbers its degrees of freedom as the host space of the same
   * degree built on Geometry::Mesh::toHost(), so the two share their global
   * numbering. The host space is kept for the form language queries, and
   * its numbering is copied to the device, where the basis is evaluated by
   * Rodin::Kokkos::Variational::QuadraticBasis.
   *
   * Only the degree two is available on the device, on meshes made of
   * segments, triangles and tetrahedra.
   *
   * @tparam K Polynomial degree, two.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Backend>
  class H1<K, Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        H1<K, Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    static_assert(K == 2, "Kokkos H1 spaces are available for K = 2.");

    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = ScalarType;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = H1Element<K, RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType =
        Rodin::Kokkos::Variational::QuadraticDOFs<typename MeshType::MemorySpace>;

      /**
       * @brief Constructs the space on a resident mesh.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Simplicial mesh supporting the space, which must
       * outlive it.
       */
      H1([[maybe_unused]] std::integral_constant<size_t, K> degree, const MeshType& mesh)
        : m_mesh(mesh),
          m_mirror(std::make_shared<const MirrorType>(mesh, 0))
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_mirror(other.m_mirror)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_mirror(std::move(other.m_mirror))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_mirror = other.m_mirror;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_mirror = std::move(other.m_mirror);
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_mirror->getScalar().getFiniteElement(d, i);
      }

      /// @brief Gets the number of degrees of freedom. @return Size of the space.
      size_t getSize() const override
      {
        return m_mirror->getScalar().getSize();
      }

      /// @brief Gets the vector dimension. @return One, the space being scalar.
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Global indices in local order.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_mirror->getScalar().getDOFs(d, i);
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local index within the polytope.
       * @return Global index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        return m_mirror->getScalar().getGlobalIndex(idx, local);
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Scalar H1 space on @p mesh.
       */
      H1<K, RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return H1<K, RangeType, Geometry::Mesh<Context::Local>>(
          std::integral_constant<size_t, K>{}, mesh);
      }

      /// @brief Gets the device degrees of freedom. @return Device numbering and basis.
      const DOFsType& getDeviceDOFs() const
      {
        return m_mirror->getDOFs();
      }

    private:
      using MirrorType = Rodin::Kokkos::Variational::Internal::QuadraticMirror<Backend>;

      std::reference_wrapper<const MeshType> m_mesh;
      std::shared_ptr<const MirrorType> m_mirror;
  };

  /**
   * @ingroup H1Specializations
   * @brief Vector valued quadratic H1 Lagrange space on a resident
   * simplicial Kokkos mesh.
   *
   * The space numbers its degrees of freedom as the host vector space:
   * component @f$ r @f$ of the scalar degree of freedom @f$ q @f$ has the
   * global index @f$ q + r n @f$, where @f$ n @f$ is the size of the scalar
   * space, and the local index @f$ a\,d + r @f$ when @f$ q @f$ is the local
   * node @f$ a @f$ of a polytope and @f$ d @f$ the vector dimension.
   *
   * @tparam K Polynomial degree, two.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Backend>
  class H1<K, Math::SpatialVector<Real>,
    Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>
    final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        H1<K, Math::SpatialVector<Real>,
          Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    static_assert(K == 2, "Kokkos H1 spaces are available for K = 2.");

    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = Math::SpatialVector<ScalarType>;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = H1Element<K, RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType =
        Rodin::Kokkos::Variational::QuadraticDOFs<typename MeshType::MemorySpace>;

      /**
       * @brief Constructs the vector space on a resident mesh.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Simplicial mesh supporting the space, which must
       * outlive it.
       * @param[in] vdim Vector dimension, at most three.
       */
      H1([[maybe_unused]] std::integral_constant<size_t, K> degree, const MeshType& mesh,
        size_t vdim)
        : m_mesh(mesh),
          m_vdim(vdim),
          m_mirror(std::make_shared<const MirrorType>(mesh, vdim))
      {
        assert(vdim > 0 && vdim <= RODIN_MAXIMAL_SPACE_DIMENSION);
      }

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_mirror(other.m_mirror)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(other.m_vdim),
          m_mirror(std::move(other.m_mirror))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
          m_mirror = other.m_mirror;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = other.m_vdim;
          m_mirror = std::move(other.m_mirror);
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Vector element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_mirror->getVector().getFiniteElement(d, i);
      }

      /// @brief Gets the number of degrees of freedom. @return Size of the space.
      size_t getSize() const override
      {
        return m_mirror->getVector().getSize();
      }

      /// @brief Gets the vector dimension. @return Number of components.
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Global indices in local order.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_mirror->getVector().getDOFs(d, i);
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local index within the polytope.
       * @return Global index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        return m_mirror->getVector().getGlobalIndex(idx, local);
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Vector H1 space on @p mesh.
       */
      H1<K, RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return H1<K, RangeType, Geometry::Mesh<Context::Local>>(
          std::integral_constant<size_t, K>{}, mesh, m_vdim);
      }

      /// @brief Gets the device degrees of freedom. @return Scalar device numbering and basis.
      const DOFsType& getDeviceDOFs() const
      {
        return m_mirror->getDOFs();
      }

    private:
      using MirrorType = Rodin::Kokkos::Variational::Internal::QuadraticMirror<Backend>;

      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
      std::shared_ptr<const MirrorType> m_mirror;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief Quadratic H1 space on a resident mesh.
  template <class Range, class Backend>
  using DeviceH1Space = ::Rodin::Variational::H1<2, Range,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;

  /// @brief Real valued quadratic H1 space on a resident mesh.
  template <class Backend>
  using DeviceH1 = DeviceH1Space<Real, Backend>;

  /// @brief Vector valued quadratic H1 space on a resident mesh.
  template <class Backend>
  using DeviceVectorH1 = DeviceH1Space<Math::SpatialVector<Real>, Backend>;
}

#endif
