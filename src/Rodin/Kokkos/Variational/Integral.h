/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTEGRAL_H
#define RODIN_KOKKOS_VARIATIONAL_INTEGRAL_H

/**
 * @file
 * @brief Device integration of bilinear and linear forms on resident Kokkos
 * meshes.
 */

#include <algorithm>
#include <array>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/Dot.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/BoundaryIntegral.h"
#include "Rodin/Variational/FaceIntegral.h"
#include "Rodin/Variational/InterfaceIntegral.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/H1/ForwardDecls.h"

#include "Rodin/Kokkos/Utility.h"

#include "DeviceView.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Coordinate list of a sparse matrix in device memory.
   *
   * Duplicate entries are summed by the consumer.
   *
   * @tparam MemorySpace Kokkos memory space of the arrays.
   */
  template <class MemorySpace>
  struct DeviceTriplets
  {
    /// @brief Row indices, test degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> rows;

    /// @brief Column indices, trial degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> cols;

    /// @brief Entry values.
    ::Kokkos::View<Real*, MemorySpace> values;
  };

  /**
   * @brief Arrangement of the entries of a device integrator.
   *
   * Integrators with equal layouts list the same coordinates in the same
   * order, so their values can be summed entry by entry.
   */
  struct DeviceLayout
  {
    /// @brief Integrated mesh.
    const void* mesh = nullptr;

    /// @brief Vector dimension of the trial space, zero for linear forms.
    Index trialVDim = 0;

    /// @brief Vector dimension of the test space.
    Index testVDim = 0;

    /// @brief Selected attributes; empty selects every polytope of the region.
    std::vector<::Rodin::Geometry::Attribute> attributes;

    /// @brief Integration region.
    ::Rodin::Geometry::Region region = ::Rodin::Geometry::Region::Cells;

    /// @brief Offset between the components of the trial space, zero for linear forms.
    Index trialStride = 0;

    /// @brief Offset between the components of the test space.
    Index testStride = 0;

    /// @brief Maximum scalar local count of the trial space, zero for linear forms.
    Index trialCapacity = 0;

    /// @brief Maximum scalar local count of the test space.
    Index testCapacity = 0;

    /// @brief Whether the entries couple pairs of polytopes, as global integrators do.
    bool global = false;

    /// @brief Selected trial attributes of a global integrator, the
    /// attributes and the region above selecting its test polytopes.
    std::vector<::Rodin::Geometry::Attribute> trialAttributes = {};

    bool operator==(const DeviceLayout&) const = default;
  };

  /**
   * @brief Finds where an entry of a device integrator is stored.
   * @param rank Rank of the entry, the entries of each polytope following
   * the entries of the previous polytope.
   * @param count Number of entries.
   * @param spacing Spacing of the entries, as given by the integrator.
   * @return Position of the entry past the offset of the integrator.
   *
   * When every integrated polytope owns the same number n of entries and
   * there are M > 1 polytopes, entry j of the k-th polytope is stored at
   * k + j * M and the spacing is M, so that the threads of consecutive
   * polytopes write contiguous memory. Otherwise the spacing is one and
   * the entries are stored in rank order.
   */
  KOKKOS_INLINE_FUNCTION size_t getEntryPosition(
    size_t rank, size_t count, size_t spacing)
  {
    const size_t n = count / spacing;
    return rank / n + (rank % n) * spacing;
  }

  /**
   * @brief Finds the rank of the entry stored at a position, the inverse
   * of getEntryPosition().
   * @param position Position of the entry past the offset of the integrator.
   * @param count Number of entries.
   * @param spacing Spacing of the entries, as given by the integrator.
   * @return Rank of the entry.
   */
  KOKKOS_INLINE_FUNCTION size_t getEntryRank(
    size_t position, size_t count, size_t spacing)
  {
    const size_t n = count / spacing;
    return (position % spacing) * n + position / spacing;
  }

  /**
   * @brief Interface of bilinear form integrators evaluated on the device.
   *
   * Assembly::Kokkos discovers this interface on the local integrators of a
   * bilinear form and asks each one for its contributions as device triplets.
   * An integration is a call to prepare(), which locates the entries of the
   * cells, followed by fill(), which writes them.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class DeviceBilinearFormIntegrator
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Memory space of the triplets.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device triplet type.
      using TripletsType = DeviceTriplets<MemorySpace>;

      virtual ~DeviceBilinearFormIntegrator() = default;

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      virtual DeviceLayout getLayout() const = 0;

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries.
       */
      virtual size_t prepare() const = 0;

      /**
       * @brief Gets the spacing of the entries located by the last call to
       * prepare().
       * @return Spacing given to getEntryPosition().
       */
      virtual size_t getSpacing() const = 0;

      /**
       * @brief Integrates the form over its region.
       * @param values Entry values, the entry of rank r being written at
       * offset + getEntryPosition(r, count, getSpacing()).
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values
       * instead of replacing them.
       * @param rows Row indices, test degrees of freedom, written as the
       * values unless empty.
       * @param cols Column indices, trial degrees of freedom, written as the
       * values unless empty.
       *
       * The entries are the ones located by the last call to prepare().
       */
      virtual void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const = 0;

      /**
       * @brief Integrates the form over its region into new triplets.
       * @return Local contributions, with rows indexing test and columns
       * indexing trial degrees of freedom, arranged as by fill() at offset
       * zero.
       */
      TripletsType assemble() const
      {
        const size_t count = prepare();
        TripletsType result;
        result.rows = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet rows")),
          count);
        result.cols = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet columns")),
          count);
        result.values = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin triplet values")),
          count);
        fill(result.values, 0, Real(1), false, result.rows, result.cols);
        ::Kokkos::fence("Rodin bilinear form integration");
        return result;
      }
  };

  /**
   * @brief Coordinate list of a vector in device memory.
   *
   * Duplicate entries are summed by the consumer.
   *
   * @tparam MemorySpace Kokkos memory space of the arrays.
   */
  template <class MemorySpace>
  struct DeviceEntries
  {
    /// @brief Entry indices, test degrees of freedom.
    ::Kokkos::View<Index*, MemorySpace> indices;

    /// @brief Entry values.
    ::Kokkos::View<Real*, MemorySpace> values;
  };

  /**
   * @brief Interface of linear form integrators evaluated on the device.
   *
   * Assembly::Kokkos discovers this interface on the integrators of a
   * linear form and asks each one for its contributions as device entries.
   * An integration is a call to prepare(), which locates the entries of the
   * cells, followed by fill(), which writes them.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class DeviceLinearFormIntegrator
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Memory space of the entries.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device entry type.
      using EntriesType = DeviceEntries<MemorySpace>;

      virtual ~DeviceLinearFormIntegrator() = default;

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      virtual DeviceLayout getLayout() const = 0;

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries.
       */
      virtual size_t prepare() const = 0;

      /**
       * @brief Gets the spacing of the entries located by the last call to
       * prepare().
       * @return Spacing given to getEntryPosition().
       */
      virtual size_t getSpacing() const = 0;

      /**
       * @brief Integrates the form over its region.
       * @param values Entry values, the entry of rank r being written at
       * offset + getEntryPosition(r, count, getSpacing()).
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values
       * instead of replacing them.
       * @param indices Entry indices, test degrees of freedom, written as the
       * values unless empty.
       *
       * The entries are the ones located by the last call to prepare().
       */
      virtual void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const = 0;

      /**
       * @brief Integrates the form over its region into new entries.
       * @return Local contributions, indexed by test degrees of freedom and
       * arranged as by fill() at offset zero.
       */
      EntriesType assemble() const
      {
        const size_t count = prepare();
        EntriesType result;
        result.indices = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin entry indices")),
          count);
        result.values = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(
            ::Kokkos::WithoutInitializing, std::string("Rodin entry values")),
          count);
        fill(result.values, 0, Real(1), false, result.indices);
        ::Kokkos::fence("Rodin linear form integration");
        return result;
      }
  };

  /**
   * @brief Tests whether a cell belongs to an attribute selection.
   * @param cell Device cell.
   * @param attributes Selected attributes; empty selects every cell.
   * @return Whether the cell is integrated.
   */
  template <class AttributeView>
  KOKKOS_INLINE_FUNCTION bool isSelected(
    const Geometry::PolytopeView& cell, const AttributeView& attributes)
  {
    if (attributes.extent(0) == 0)
      return true;
    if (!cell.hasAttribute())
      return false;
    const auto attribute = cell.getAttribute();
    for (Index i = 0; i < attributes.extent(0); ++i)
    {
      if (attributes(i) == attribute)
        return true;
    }
    return false;
  }

  /**
   * @brief Tests whether a polytope is integrated over a region.
   * @param view Device mesh.
   * @param region Integration region.
   * @param polytope Device polytope of the dimension of the region.
   * @param attributes Selected attributes; empty selects every polytope of
   * the region.
   * @return Whether the polytope is integrated.
   *
   * Boundary faces are incident to one cell and interface faces to more
   * than one, as on the host. These regions need the face to cell
   * incidence.
   */
  template <class AttributeView>
  KOKKOS_INLINE_FUNCTION bool isIntegrated(const Geometry::MeshView& view,
    ::Rodin::Geometry::Region region, const Geometry::PolytopeView& polytope,
    const AttributeView& attributes)
  {
    if (region == ::Rodin::Geometry::Region::Boundary ||
      region == ::Rodin::Geometry::Region::Interface)
    {
      const Index d = view.getDimension();
      const Index n =
        view.getConnectivity().getIncidence(d - 1, d, polytope.getIndex()).size();
      if (region == ::Rodin::Geometry::Region::Boundary ? n != 1 : n <= 1)
        return false;
    }
    return isSelected(polytope, attributes);
  }
}

namespace Rodin::Kokkos::Variational::Internal
{
  /// @brief Kokkos backend of a device space.
  template <class Mesh>
  struct MeshBackend;

  /// @brief Backend of a resident mesh.
  template <class Backend>
  struct MeshBackend<
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>
  {
    using Type = Backend;
  };

  /// @brief Kokkos backend of the mesh of a device space.
  template <class FES>
  using BackendOf = typename MeshBackend<typename FES::MeshType>::Type;

  /// @brief Whether a space is a P1 space.
  template <class FES>
  struct IsP1Space
  {
    static constexpr bool Value = false;
  };

  /// @brief P1 space.
  template <class Range, class Mesh>
  struct IsP1Space<::Rodin::Variational::P1<Range, Mesh>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Whether a space is an H1 space.
  template <class FES>
  struct IsH1Space
  {
    static constexpr bool Value = false;
  };

  /// @brief H1 space.
  template <size_t K, class Range, class Mesh>
  struct IsH1Space<::Rodin::Variational::H1<K, Range, Mesh>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Derivatives of a trial/test pair handled by a dedicated host handler.
  enum class DerivativeKind
  {
    /// @brief No dedicated handler.
    None,

    /// @brief Gradients or Jacobians on both sides.
    Gradients,

    /// @brief Divergence of the trial function against a scalar test function.
    TrialDivergence,

    /// @brief Scalar trial function against the divergence of the test function.
    TestDivergence,

    /// @brief Jacobian of the trial function times a function against a test function.
    JacobianValue
  };

  /**
   * @brief Derivatives of a trial/test pair, selecting the fallback
   * quadrature order of the host handlers.
   *
   * When both spaces are P1, the host handlers of forms with derivatives,
   * such as @f$ \nabla u \cdot \nabla v @f$, @f$ \nabla u : \nabla v @f$
   * or @f$ (\nabla \cdot u) q @f$, fall back to the quadrature order
   * @f$ k_u + k_v - 2 @f$. When both are H1, the order loses one degree per
   * derivative. The device reproduces both.
   */
  template <class LHSDerived, class RHSDerived>
  struct DerivativeForm
  {
    static constexpr DerivativeKind Kind = DerivativeKind::None;
  };

  /// @brief Gradient trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Scaled gradient trial operator and gradient test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Divergence trial operator and scalar test function.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    static constexpr DerivativeKind Kind =
      std::is_same_v<typename TestFES::RangeType, Real>
      ? DerivativeKind::TrialDivergence
      : DerivativeKind::None;
  };

  /// @brief Scalar trial function and divergence test operator.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
                              ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr DerivativeKind Kind =
      std::is_same_v<typename TrialFES::RangeType, Real>
      ? DerivativeKind::TestDivergence
      : DerivativeKind::None;
  };

  /// @brief Jacobian trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
      TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Scaled Jacobian trial operator and Jacobian test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr DerivativeKind Kind = DerivativeKind::Gradients;
  };

  /// @brief Jacobian trial operator times a function, against a test function.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct DerivativeForm<
    ::Rodin::Variational::Mult<
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
      ::Rodin::Variational::FunctionBase<CoefficientDerived>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    static constexpr DerivativeKind Kind = DerivativeKind::JacobianValue;
  };

  /**
   * @brief Quadrature tables of the polytopes of one dimension of a resident
   * mesh, one formula per reference geometry.
   *
   * The quadrature order is selected on the host, once per reference
   * geometry, using the first polytope of each geometry. It must therefore
   * depend only on the geometry, which holds for meshes whose polytopes of a
   * given geometry share one transformation degree. The reference
   * coordinates are flattened and padded to three components.
   *
   * @tparam Mesh Resident Kokkos mesh type.
   */
  template <class Mesh>
  class DeviceQuadrature
  {
    public:
      /// @brief Kokkos execution space.
      using ExecutionSpace = typename Mesh::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename Mesh::MemorySpace;

      /// @brief Number of reference geometries.
      static constexpr size_t GeometryCount = ::Rodin::Geometry::Polytope::Types.size();

      /// @brief Stride of the flattened reference coordinates.
      static constexpr Index ReferenceStride = 3;

      /**
       * @brief Builds and uploads the tables of every geometry present.
       * @param mesh Resident mesh.
       * @param dimension Dimension of the integrated polytopes.
       * @param order Callable giving the quadrature order on a host polytope.
       */
      template <class Order>
      DeviceQuadrature(const Mesh& mesh, size_t dimension, const Order& order)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        m_representatives = getRepresentatives(mesh, dimension);
        const auto& representatives = m_representatives;
        std::vector<Index> offsets(GeometryCount + 1, 0);
        std::vector<Real> weights;
        std::vector<Real> points;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          offsets[g] = weights.size();
          if (!representatives[g])
            continue;
          const auto it = mesh.getPolytope(dimension, *representatives[g]);
          const ::Rodin::Geometry::Polytope& polytope = *it;
          m_vertices = std::max<Index>(m_vertices, polytope.getVertices().size());
          m_minimum = m_minimum
            ? std::min<Index>(m_minimum, polytope.getVertices().size())
            : polytope.getVertices().size();
          m_orders[g] = order(polytope);
          const auto& qf =
            QF::QuadratureFormula::get(m_orders[g], polytope.getGeometry());
          for (size_t qp = 0; qp < qf.getSize(); ++qp)
          {
            weights.push_back(qf.getWeight(qp));
            const auto& point = qf.getPoint(qp);
            for (size_t k = 0; k < ReferenceStride; ++k)
              points.push_back(
                k < static_cast<size_t>(point.size()) ? point(k) : Real(0));
          }
        }
        offsets[GeometryCount] = weights.size();
        m_offsets =
          Rodin::Kokkos::Internal::upload(execution, offsets, "Rodin quadrature offsets");
        m_weights =
          Rodin::Kokkos::Internal::upload(execution, weights, "Rodin quadrature weights");
        m_points =
          Rodin::Kokkos::Internal::upload(execution, points, "Rodin quadrature points");
      }

      /**
       * @brief Tests whether the tables were built with the given orders.
       * @param mesh Resident mesh the tables were built on.
       * @param dimension Dimension of the integrated polytopes.
       * @param order Callable giving the quadrature order on a host polytope.
       * @return Whether every geometry present has the same order.
       */
      template <class Order>
      bool hasOrders(const Mesh& mesh, size_t dimension, const Order& order) const
      {
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (!m_representatives[g])
            continue;
          const auto it = mesh.getPolytope(dimension, *m_representatives[g]);
          if (order(*it) != m_orders[g])
            return false;
        }
        return true;
      }

      /// @brief Gets the largest vertex count of the polytopes. @return Vertex count.
      Index getMaximumVertexCount() const
      {
        return m_vertices;
      }

      /**
       * @brief Returns the smallest vertex count of the geometries present.
       * @return Minimum vertex count, zero on an empty mesh.
       */
      Index getMinimumVertexCount() const
      {
        return m_minimum;
      }

      /// @param g Reference geometry. @return First quadrature point of the geometry.
      KOKKOS_INLINE_FUNCTION
      Index begin(Index g) const
      {
        return m_offsets(g);
      }

      /// @param g Reference geometry. @return Past the last quadrature point of the geometry.
      KOKKOS_INLINE_FUNCTION
      Index end(Index g) const
      {
        return m_offsets(g + 1);
      }

      /// @param q Quadrature point. @return Reference weight.
      KOKKOS_INLINE_FUNCTION
      Real getWeight(Index q) const
      {
        return m_weights(q);
      }

      /**
       * @brief Reads the reference coordinates of a quadrature point.
       * @param cell Device polytope.
       * @param q Quadrature point of the polytope geometry.
       * @return Reference coordinates, of the polytope dimension.
       */
      KOKKOS_INLINE_FUNCTION
      Rodin::Kokkos::Geometry::Vector getReference(
        const Rodin::Kokkos::Geometry::PolytopeView& cell, Index q) const
      {
        const Index rdim = cell.getDimension();
        Rodin::Kokkos::Geometry::Vector reference(rdim);
        for (Index k = 0; k < rdim; ++k)
          reference(k) = m_points(q * ReferenceStride + k);
        return reference;
      }

      /**
       * @brief Maps a quadrature point onto a polytope.
       * @param cell Device polytope.
       * @param q Quadrature point of the polytope geometry.
       * @return Point on the polytope.
       */
      KOKKOS_INLINE_FUNCTION
      Rodin::Kokkos::Geometry::PointView getPoint(
        const Rodin::Kokkos::Geometry::PolytopeView& cell, Index q) const
      {
        return cell.getPoint(getReference(cell, q));
      }

      /**
       * @brief Calls a functor with the static geometry tag of the polytopes.
       *
       * The tag applies when every polytope has one geometry with a static
       * basis, of dimension the spatial dimension.
       *
       * @param spaceDimension Spatial dimension of the mesh.
       * @param f Callable taking a StaticGeometry tag and returning whether it
       * integrated the polytopes.
       * @return Result of @p f, false when no tag applies.
       */
      template <class F>
      bool visitStaticGeometry(size_t spaceDimension, F&& f) const
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        Optional<Type> geometry;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (!m_representatives[g])
            continue;
          if (geometry)
            return false;
          geometry = static_cast<Type>(g);
        }
        if (!geometry ||
          ::Rodin::Geometry::Polytope::Traits(*geometry).getDimension() != spaceDimension)
          return false;
        switch (*geometry)
        {
          case Type::Triangle:
            return f(StaticGeometry<Type::Triangle>{});
          case Type::Quadrilateral:
            return f(StaticGeometry<Type::Quadrilateral>{});
          case Type::Tetrahedron:
            return f(StaticGeometry<Type::Tetrahedron>{});
          case Type::Hexahedron:
            return f(StaticGeometry<Type::Hexahedron>{});
          default:
            return false;
        }
      }

      /**
       * @brief Builds a range policy on the mesh execution space.
       * @param execution Execution space instance.
       * @param count Iteration count.
       * @return Range policy over @f$ [0, count) @f$.
       */
      static auto policy(const ExecutionSpace& execution, Index count)
      {
        return ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
          execution, 0, count);
      }

      /**
       * @brief Finds the first polytope of every reference geometry.
       * @param mesh Resident mesh.
       * @param dimension Dimension of the polytopes.
       * @return Index of the first polytope per geometry, empty when absent.
       */
      static std::array<Optional<Index>, GeometryCount> getRepresentatives(
        const Mesh& mesh, size_t dimension)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        constexpr Index None = std::numeric_limits<Index>::max();
        ::Kokkos::View<Index[GeometryCount], MemorySpace> first(
          ::Kokkos::view_alloc(execution, std::string("Rodin geometry representatives")));
        ::Kokkos::deep_copy(execution, first, None);
        const auto view = mesh.getDeviceView();
        findRepresentatives(
          execution, view, first, dimension, mesh.getPolytopeCount(dimension));
        auto host = ::Kokkos::create_mirror_view(first);
        ::Kokkos::deep_copy(execution, host, first);
        execution.fence("Rodin geometry representatives");
        std::array<Optional<Index>, GeometryCount> result;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (host(g) != None)
            result[g] = host(g);
        }
        return result;
      }

      /**
       * @brief Records the minimum polytope index of every geometry.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param first Per-geometry minimum, initialized to the maximum index.
       * @param dimension Dimension of the polytopes.
       * @param count Polytope count.
       *
       * Public as the enclosing function of a device lambda.
       */
      static void findRepresentatives(const ExecutionSpace& execution,
        const Rodin::Kokkos::Geometry::MeshView& view,
        const ::Kokkos::View<Index[GeometryCount], MemorySpace>& first, Index dimension,
        Index count)
      {
        ::Kokkos::parallel_for(
          "Rodin geometry representatives", policy(execution, count),
          KOKKOS_LAMBDA(Index p) {
            const Index g =
              static_cast<Index>(view.getPolytope(dimension, p).getGeometry());
            ::Kokkos::atomic_min(&first(g), p);
          });
      }

    private:
      ::Kokkos::View<Index*, MemorySpace> m_offsets;
      ::Kokkos::View<Real*, MemorySpace> m_weights;
      ::Kokkos::View<Real*, MemorySpace> m_points;
      std::array<Optional<Index>, GeometryCount> m_representatives;
      std::array<size_t, GeometryCount> m_orders{};
      Index m_vertices = 0;
      Index m_minimum = 0;
  };

  /**
   * @brief Uploads the attribute selection of an integrator.
   * @param execution Execution space instance.
   * @param attributes Selected attributes; empty selects every cell.
   * @return Device copy of the selection.
   */
  template <class Execution>
  auto uploadAttributes(
    const Execution& execution, const FlatSet<::Rodin::Geometry::Attribute>& attributes)
  {
    const std::vector<::Rodin::Geometry::Attribute> list(
      attributes.begin(), attributes.end());
    return Rodin::Kokkos::Internal::upload(
      execution, list, "Rodin integration attributes");
  }

  /**
   * @brief Dimension of the polytopes of an integration region.
   * @param dimension Mesh dimension.
   * @param region Integration region.
   * @return The mesh dimension for cells, the face dimension otherwise.
   */
  inline size_t getRegionDimension(size_t dimension, ::Rodin::Geometry::Region region)
  {
    return region == ::Rodin::Geometry::Region::Cells ? dimension : dimension - 1;
  }

  /**
   * @brief Polytopes integrated by a device integrator and the position of
   * their entries.
   *
   * The polytopes are the cells of the mesh or its faces, as selected by the
   * integration region. The quadrature tables and the attribute selection
   * are uploaded on the first integration and again only when the mesh, the
   * region, the quadrature orders or the selection change, so that repeated
   * integrations do not allocate device memory.
   *
   * @tparam Mesh Resident Kokkos mesh type.
   */
  template <class Mesh>
  class DevicePolytopes
  {
    public:
      /// @brief Kokkos execution space.
      using ExecutionSpace = typename Mesh::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename Mesh::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = DeviceQuadrature<Mesh>;

      /// @brief Device attribute selection.
      using AttributesType = ::Kokkos::View<::Rodin::Geometry::Attribute*, MemorySpace>;

      /**
       * @brief Locates the entries of every selected polytope.
       * @param mesh Resident mesh.
       * @param region Integration region.
       * @param attributes Selected attributes; empty selects every polytope
       * of the region.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param order Callable giving the quadrature order on a host polytope.
       * @return Number of entries.
       * @throws Alert::Exception If the region is the boundary or the
       * interface and the face to cell incidence is missing.
       */
      template <class TrialDOFs, class TestDOFs, class Order>
      size_t prepare(const Mesh& mesh, ::Rodin::Geometry::Region region,
        const FlatSet<::Rodin::Geometry::Attribute>& attributes,
        const TrialDOFs& trialDOFs, Index trialVDim, const TestDOFs& testDOFs,
        Index testVDim, const Order& order)
      {
        const size_t meshDimension = mesh.getDimension();
        if ((region == ::Rodin::Geometry::Region::Boundary ||
              region == ::Rodin::Geometry::Region::Interface) &&
          (meshDimension == 0 ||
            !mesh.getConnectivity().hasIncidence(meshDimension - 1, meshDimension)))
        {
          Alert::Exception()
            << "Compute face-to-cell incidence before integrating over the "
            << "boundary or the interface of a Kokkos mesh." << Alert::Raise;
        }
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t dimension = getRegionDimension(meshDimension, region);
        const Index count = mesh.getPolytopeCount(dimension);
        if (m_mesh != &mesh || m_dimension != dimension || m_polytopeCount != count ||
          !m_quadrature || !m_quadrature->hasOrders(mesh, dimension, order))
        {
          m_quadrature.emplace(mesh, dimension, order);
          m_offsets = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin entry offsets")),
            count);
          m_mesh = &mesh;
          m_dimension = dimension;
          m_polytopeCount = count;
        }
        if (!m_selection || *m_selection != attributes)
        {
          m_attributes = uploadAttributes(execution, attributes);
          m_selection = attributes;
        }
        m_region = region;
        m_count = locate(execution, mesh.getDeviceView(), region, dimension, m_attributes,
          m_offsets, count, trialDOFs, trialVDim, testDOFs, testVDim, m_spacing);
        return m_count;
      }

      /// @brief Gets the quadrature tables. @return Tables of the last integration.
      const QuadratureType& getQuadrature() const
      {
        assert(m_quadrature);
        return *m_quadrature;
      }

      /// @brief Gets the selection. @return Device attributes.
      const AttributesType& getAttributes() const
      {
        return m_attributes;
      }

      /// @brief Gets the position of the first entry of every polytope. @return Offsets.
      const ::Kokkos::View<Index*, MemorySpace>& getOffsets() const
      {
        return m_offsets;
      }

      /**
       * @brief Gets the spacing of the entries of the last integration.
       * @return Distance between consecutive entries of one polytope.
       */
      size_t getSpacing() const
      {
        return m_spacing;
      }

      /// @brief Gets the integration region. @return Region of the last integration.
      ::Rodin::Geometry::Region getRegion() const
      {
        return m_region;
      }

      /// @brief Gets the dimension of the polytopes. @return Polytope dimension.
      Index getDimension() const
      {
        return m_dimension;
      }

      /// @brief Gets the number of polytopes of the dimension. @return Polytope count.
      Index getPolytopeCount() const
      {
        return m_polytopeCount;
      }

      /// @brief Gets the entry count of the last integration. @return Count.
      size_t getCount() const
      {
        return m_count;
      }

      /**
       * @brief Computes the entry offsets of the polytopes.
       *
       * A polytope owns the product of its trial and test local counts. When
       * every integrated polytope owns the same count and there are more
       * than one, the offset of a polytope is its rank among them and the
       * spacing is their number, as described by getEntryPosition().
       * Otherwise the offset of a polytope is the count of the previous
       * ones and the spacing is one.
       *
       * Public as the enclosing function of a device lambda.
       */
      template <class TrialDOFs, class TestDOFs>
      static size_t locate(const ExecutionSpace& execution,
        const Rodin::Kokkos::Geometry::MeshView& view, ::Rodin::Geometry::Region region,
        Index dimension, const AttributesType& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& offsets, Index count,
        const TrialDOFs& trialDOFs, Index trialVDim, const TestDOFs& testDOFs,
        Index testVDim, size_t& spacing)
      {
        const auto local = KOKKOS_LAMBDA(Index p)->Index
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return 0;
          return trialDOFs.bind(polytope).size() * trialVDim *
            testDOFs.bind(polytope).size() * testVDim;
        };
        Index selected = 0;
        Index smallest = 0;
        Index largest = 0;
        ::Kokkos::parallel_reduce(
          "Rodin entry counts", QuadratureType::policy(execution, count),
          KOKKOS_LAMBDA(Index p, Index & m, Index & lo, Index & hi) {
            const Index n = local(p);
            if (n == 0)
              return;
            m += 1;
            lo = n < lo ? n : lo;
            hi = n > hi ? n : hi;
          },
          ::Kokkos::Sum<Index>(selected), ::Kokkos::Min<Index>(smallest),
          ::Kokkos::Max<Index>(largest));
        const bool uniform = selected > 1 && smallest == largest;
        Index total = 0;
        ::Kokkos::parallel_scan(
          "Rodin entry offsets", QuadratureType::policy(execution, count),
          KOKKOS_LAMBDA(Index p, Index & update, bool final) {
            if (final)
              offsets(p) = update;
            const Index n = local(p);
            update += uniform && n > 0 ? 1 : n;
          },
          total);
        spacing = uniform ? selected : 1;
        return uniform ? size_t(selected) * smallest : total;
      }

    private:
      const Mesh* m_mesh = nullptr;
      ::Rodin::Geometry::Region m_region = ::Rodin::Geometry::Region::Cells;
      Index m_dimension = 0;
      Index m_polytopeCount = 0;
      Optional<QuadratureType> m_quadrature;
      Optional<FlatSet<::Rodin::Geometry::Attribute>> m_selection;
      AttributesType m_attributes;
      ::Kokkos::View<Index*, MemorySpace> m_offsets;
      size_t m_count = 0;
      size_t m_spacing = 1;
  };
}

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Device integration of a trial/test dot product over a region of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{R} A(u) \cdot B(v) \ dx
   * @f]
   * for scalar or vector device trial and test spaces on a resident Kokkos
   * mesh, possibly different, such as the velocity and pressure spaces of a
   * Stokes problem. The region @f$ R @f$ is given by getRegion(): the cells,
   * the faces, the boundary or the interface of the mesh. The operators
   * @f$ A @f$ and @f$ B @f$ are lowered once to DeviceView, then a kernel
   * integrates every polytope of the region with one thread and writes its
   * local matrix as device triplets.
   *
   * A face is integrated with the basis of the face, whose degrees of
   * freedom are the ones of the face, as on the host.
   *
   * The quadrature order is selected on the host as by the QuadratureRule
   * the host would use, once per reference geometry, as described in
   * Internal::DeviceQuadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   *
   * Rodin::Variational::Integral, BoundaryIntegral, FaceIntegral and
   * InterfaceIntegral specialize this class for device spaces.
   */
  template <class LHSDerived, class RHSDerived, DeviceSpace TrialFES, DeviceSpace TestFES>
  class DeviceBilinearIntegral
    : public ::Rodin::Variational::LocalBilinearFormIntegratorBase<Real>,
      public DeviceBilinearFormIntegrator<Internal::BackendOf<TrialFES>>
  {
    public:
      /// @brief Finite element space of the trial functions.
      using TrialFESType = TrialFES;

      /// @brief Finite element space of the test functions.
      using TestFESType = TestFES;

      /// @brief Range type of the trial space.
      using TrialRange = typename TrialFESType::RangeType;

      /// @brief Range type of the test space.
      using TestRange = typename TestFESType::RangeType;

      /// @brief Device degrees of freedom of the trial space.
      using TrialDOFsType = typename TrialFESType::DOFsType;

      /// @brief Device degrees of freedom of the test space.
      using TestDOFsType = typename TestFESType::DOFsType;

      /// @brief Type of the left operand of the dot product.
      using LHSType = ::Rodin::Variational::ShapeFunctionBase<LHSDerived, TrialFESType,
        ::Rodin::Variational::ShapeFunctionSpaceType::Trial>;

      /// @brief Type of the right operand of the dot product.
      using RHSType = ::Rodin::Variational::ShapeFunctionBase<RHSDerived, TestFESType,
        ::Rodin::Variational::ShapeFunctionSpaceType::Test>;

      /// @brief Type of the integrand.
      using IntegrandType = ::Rodin::Variational::Dot<LHSType, RHSType>;

      /// @brief Host integrator interface.
      using Parent = ::Rodin::Variational::LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = DeviceBilinearFormIntegrator<Internal::BackendOf<TrialFES>>;

      /// @brief Device triplet type.
      using TripletsType = typename DeviceParent::TripletsType;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos execution space.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Internal::DevicePolytopes<MeshType>;

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      DeviceBilinearIntegral(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      DeviceBilinearIntegral(const DeviceBilinearIntegral& other)
        : Parent(other),
          DeviceParent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      DeviceBilinearIntegral(DeviceBilinearIntegral&& other)
        : Parent(std::move(other)),
          DeviceParent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      /// @brief Gets the integrand. @return Dot product being integrated.
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      const ::Rodin::Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      DeviceBilinearIntegral& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      DeviceLayout getLayout() const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&trialfes.getMesh(), static_cast<Index>(trialfes.getVectorDimension()),
          static_cast<Index>(testfes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(),
          trialfes.getDeviceDOFs().getStride(), testfes.getDeviceDOFs().getStride(),
          TrialDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES),
          TestDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /**
       * @brief Locates the entries of every selected polytope.
       * @return Number of entries, each selected polytope with m trial and
       * n test scalar degrees of freedom owning m * trialVDim * n * testVDim
       * entries.
       */
      size_t prepare() const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const MeshType& mesh = trialfes.getMesh();
        if (&testfes.getMesh() != &mesh)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos integrals require trial and test spaces on the same mesh."
            << Alert::Raise;
        }
        return m_polytopes.prepare(mesh, this->getRegion(), this->getAttributes(),
          trialfes.getDeviceDOFs(), trialfes.getVectorDimension(),
          testfes.getDeviceDOFs(), testfes.getVectorDimension(),
          [&](const ::Rodin::Geometry::Polytope& polytope) {
            return getOrder(integrand, polytope);
          });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /**
       * @brief Integrates every selected polytope on the device.
       *
       * Each selected polytope writes its local matrix, the entry of local
       * test index te and local trial index tr having rank te * ntr + tr
       * among the entries of the polytope.
       */
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const MeshType& mesh = trialfes.getMesh();
        const DeviceView<LHSType> trial(integrand.getLHS());
        const DeviceView<RHSType> test(integrand.getRHS());
        const TrialDOFsType trialDOFs = trialfes.getDeviceDOFs();
        const TestDOFsType testDOFs = testfes.getDeviceDOFs();
        const Index trialVDim = trialfes.getVectorDimension();
        const Index testVDim = testfes.getVectorDimension();
        const auto& quadrature = m_polytopes.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        constexpr bool p1 = std::is_same_v<TrialDOFsType, VertexDOFs> &&
          std::is_same_v<TestDOFsType, VertexDOFs>;
        const bool maximal = trialVDim == VectorCapacity<TrialRange> &&
          testVDim == VectorCapacity<TestRange>;
        // Polytopes of one geometry filling the space take kernels of static size.
        if ((!p1 || maximal) &&
          quadrature.visitStaticGeometry(mesh.getSpaceDimension(), [&](auto shape) {
            using Shape = decltype(shape);
            constexpr Index N = Rodin::Kokkos::Geometry::StaticBasis<Shape::Value>::Count;
            constexpr Index Capacity = IntegrationPointView::Capacity;
            if constexpr (TrialDOFsType::getCapacity(N) <= Capacity &&
              TestDOFsType::getCapacity(N) <= Capacity)
            {
              integratePolytopes<N, p1, Shape>(mesh, m_polytopes, trial, test, trialDOFs,
                trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, rows,
                cols);
              return true;
            }
            return false;
          }))
          return;
        if constexpr (p1)
        {
          if (vertices == 4 && quadrature.getMinimumVertexCount() == 4 && maximal)
          {
            integratePolytopes<4, true>(mesh, m_polytopes, trial, test, trialDOFs,
              trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, rows,
              cols);
            return;
          }
        }
        constexpr Index Capacity = IntegrationPointView::Capacity;
        if (vertices <= 4)
        {
          if constexpr (getLocalSize<4>() <= LocalCapacity)
          {
            integratePolytopes<4>(mesh, m_polytopes, trial, test, trialDOFs, trialVDim,
              testDOFs, testVDim, values, offset, sign, accumulate, rows, cols);
            return;
          }
        }
        else if constexpr (getLocalSize<Capacity>() <= LocalCapacity)
        {
          integratePolytopes<Capacity>(mesh, m_polytopes, trial, test, trialDOFs,
            trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, rows, cols);
          return;
        }
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals of these spaces on polytopes of " << vertices
          << " vertices need local matrices of more than " << LocalCapacity
          << " entries, lower the polynomial degree." << Alert::Raise;
      }

      /**
       * @brief Gets the number of entries of the local matrix of a polytope.
       * @tparam Vertices Maximum vertex count of the polytopes.
       * @return Entry count of the local matrix in the kernel.
       */
      template <Index Vertices>
      static constexpr Index getLocalSize()
      {
        return TrialDOFsType::getCapacity(Vertices) * VectorCapacity<TrialRange> *
          TestDOFsType::getCapacity(Vertices) * VectorCapacity<TestRange>;
      }

      /**
       * @brief Integrates the polytopes of at most Vertices vertices.
       *
       * When Exact, both spaces are P1, every polytope has Vertices vertices
       * and both vector dimensions are maximal, so the local loops have
       * static bounds. Shape is the geometry tag of the integration points.
       * Public as the enclosing function of a device lambda.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Trial, class Test>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const Trial& trial, const Test& test, const TrialDOFsType& trialDOFs,
        Index trialVDim, const TestDOFsType& testDOFs, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        using TestValue = std::decay_t<decltype(test.getBasis(
          std::declval<const IntegrationPointView&>(), Index(0), Index(0)))>;
        constexpr Index TrialCapacity =
          TrialDOFsType::getCapacity(Vertices) * VectorCapacity<TrialRange>;
        constexpr Index TestCapacity =
          TestDOFsType::getCapacity(Vertices) * VectorCapacity<TestRange>;
        const auto& execution = mesh.getContext().getExecutionSpace();
        const auto view = mesh.getDeviceView();
        const auto quadrature = polytopes.getQuadrature();
        const auto attributes = polytopes.getAttributes();
        const auto offsets = polytopes.getOffsets();
        const size_t spacing = polytopes.getSpacing();
        const auto region = polytopes.getRegion();
        const Index dimension = polytopes.getDimension();
        const bool index = rows.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin bilinear form integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()),
          KOKKOS_LAMBDA(Index p) {
            const auto polytope = view.getPolytope(dimension, p);
            if (!Rodin::Kokkos::Variational::isIntegrated(
                  view, region, polytope, attributes))
              return;
            KOKKOS_ASSERT(Exact ? polytope.getVertices().size() == Vertices
                                : polytope.getVertices().size() <= Vertices);
            const auto trialIndices = trialDOFs.bind(polytope);
            const auto testIndices = testDOFs.bind(polytope);
            const Index ntr = Exact ? TrialCapacity : trialIndices.size() * trialVDim;
            const Index nte = Exact ? TestCapacity : testIndices.size() * testVDim;
            const Index trialStride = Exact ? VectorCapacity<TrialRange> : trialVDim;
            const Index testStride = Exact ? VectorCapacity<TestRange> : testVDim;
            const Index g = static_cast<Index>(polytope.getGeometry());
            ::Kokkos::Array<Real, TrialCapacity * TestCapacity> local{};
            ::Kokkos::Array<TestValue, TestCapacity> phi;
            const Index begin = quadrature.begin(g);
            const Index end = quadrature.end(g);
            if (begin < end)
            {
              IntegrationPointView ip(
                polytope, quadrature.getReference(polytope, begin), Shape{});
              for (Index q = begin; q < end; ++q)
              {
                if (q > begin)
                  ip.setPoint(quadrature.getReference(polytope, q), Shape{});
                const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
                const auto testBasis = Rodin::Kokkos::Variational::bindBasis(test, ip);
                const auto trialBasis = Rodin::Kokkos::Variational::bindBasis(trial, ip);
                RODIN_KOKKOS_UNROLL
                for (Index te = 0; te < nte; ++te)
                  phi[te] = testBasis.getBasis(te / testStride, te % testStride);
                RODIN_KOKKOS_UNROLL
                for (Index tr = 0; tr < ntr; ++tr)
                {
                  const auto psi =
                    trialBasis.getBasis(tr / trialStride, tr % trialStride);
                  RODIN_KOKKOS_UNROLL
                  for (Index te = 0; te < nte; ++te)
                  {
                    if (Rodin::Kokkos::Variational::overlaps(psi, phi[te]))
                      local[te * ntr + tr] +=
                        wdet * Rodin::Kokkos::Variational::dot(psi, phi[te]);
                  }
                }
              }
            }
            // Local index a * vdim + r is component r of the scalar degree of
            // freedom a, numbered globally as dof + r * stride.
            const size_t base = offset + offsets(p);
            RODIN_KOKKOS_UNROLL
            for (Index te = 0; te < nte; ++te)
            {
              const Index row =
                testIndices[te / testStride] + (te % testStride) * testDOFs.getStride();
              RODIN_KOKKOS_UNROLL
              for (Index tr = 0; tr < ntr; ++tr)
              {
                const size_t e = base + size_t(te * ntr + tr) * spacing;
                const Real v = sign * local[te * ntr + tr];
                if (accumulate)
                  values(e) += v;
                else
                  values(e) = v;
                if (index)
                {
                  rows(e) = row;
                  cols(e) = trialIndices[tr / trialStride] +
                    (tr % trialStride) * trialDOFs.getStride();
                }
              }
            }
          });
      }

      DeviceBilinearIntegral* copy() const noexcept override = 0;

    private:
      /// @brief Maximum number of components of a value of the given range.
      template <class Range>
      static constexpr Index VectorCapacity =
        std::is_same_v<Range, Real> ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION;

      /**
       * @brief Maximum number of entries of the local matrix of a kernel.
       *
       * Each thread holds its local matrix, which spills beyond the
       * registers, so higher degrees on large polytopes are rejected.
       */
      static constexpr Index LocalCapacity = 8192;

      /**
       * @brief Selects the quadrature order on a polytope.
       * @param integrand Integrated dot product.
       * @param polytope Host polytope.
       * @return Order chosen as by the host QuadratureRule handling the same
       * integrand on a local mesh.
       *
       * The P1 derivative handlers fall back to @f$ k_u + k_v - 2 @f$, the H1
       * ones lose one degree per derivative, and the generic handler falls
       * back to the sum of the leaf orders.
       */
      size_t getOrder(
        const IntegrandType& integrand, const ::Rodin::Geometry::Polytope& polytope) const
      {
        auto order = this->Parent::getOrder(polytope);
        if (!order)
          order = integrand.getOrder(polytope);
        constexpr auto kind = Internal::DerivativeForm<LHSDerived, RHSDerived>::Kind;
        constexpr bool p1 =
          Internal::IsP1Space<TrialFES>::Value && Internal::IsP1Space<TestFES>::Value;
        constexpr bool h1 =
          Internal::IsH1Space<TrialFES>::Value && Internal::IsH1Space<TestFES>::Value;
        if (!order && kind != Internal::DerivativeKind::None && (p1 || h1))
        {
          const size_t d = polytope.getDimension();
          const Index i = polytope.getIndex();
          const auto& trialfe =
            integrand.getLHS().getFiniteElementSpace().getFiniteElement(d, i);
          const auto& testfe =
            integrand.getRHS().getFiniteElementSpace().getFiniteElement(d, i);
          const size_t trialDegree = trialfe.getDegree().value_or(2 * trialfe.getOrder());
          const size_t testDegree = testfe.getDegree().value_or(2 * testfe.getOrder());
          if (p1 || kind == Internal::DerivativeKind::Gradients)
            order =
              (trialDegree == 0 || testDegree == 0) ? 0 : trialDegree + testDegree - 2;
          else if (kind == Internal::DerivativeKind::TrialDivergence)
            order = trialDegree == 0 ? 0 : trialDegree + testDegree - 1;
          else if (kind == Internal::DerivativeKind::TestDivergence)
            order = testDegree == 0 ? 0 : trialDegree + testDegree - 1;
          else
            order = trialDegree == 0 ? testDegree : trialDegree + testDegree - 1;
        }
        if (!order)
        {
          const size_t d = polytope.getDimension();
          const Index i = polytope.getIndex();
          auto trialOrder = integrand.getLHS().getLeaf().getOrder(polytope);
          auto testOrder = integrand.getRHS().getLeaf().getOrder(polytope);
          if (!trialOrder)
            trialOrder = integrand.getLHS()
                           .getFiniteElementSpace()
                           .getFiniteElement(d, i)
                           .getDegree();
          if (!testOrder)
            testOrder = integrand.getRHS()
                          .getFiniteElementSpace()
                          .getFiniteElement(d, i)
                          .getDegree();
          if (trialOrder && testOrder)
            order = *trialOrder + *testOrder;
        }
        if (!order)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic quadrature order is established. "
            << "Specify an integration order with setOrder()." << Alert::Raise;
        }
        return *order;
      }

      std::unique_ptr<IntegrandType> m_integrand;
      mutable PolytopesType m_polytopes;
  };

  /**
   * @brief Device integration of a test operator over a region of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{R} A(v) \ dx
   * @f]
   * for a scalar valued operator @f$ A @f$ of a scalar or vector device test
   * function on a resident Kokkos mesh, such as @f$ f \cdot v @f$ or
   * @f$ g \nabla \cdot v @f$. The region @f$ R @f$ is given by
   * getRegion(). The operator is lowered once to DeviceView, then a kernel
   * integrates every polytope of the region with one thread and writes its
   * local vector as device entries.
   *
   * A face is integrated with the basis of the face, whose degrees of
   * freedom are the ones of the face, as on the host.
   *
   * The quadrature order is selected on the host as by the generic
   * QuadratureRule, once per reference geometry, as described in
   * Internal::DeviceQuadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   *
   * Rodin::Variational::Integral, BoundaryIntegral, FaceIntegral and
   * InterfaceIntegral specialize this class for device spaces.
   */
  template <class NestedDerived, DeviceSpace FES>
  class DeviceLinearIntegral
    : public ::Rodin::Variational::LinearFormIntegratorBase<Real>,
      public DeviceLinearFormIntegrator<Internal::BackendOf<FES>>
  {
    public:
      /// @brief Finite element space of the test functions.
      using FESType = FES;

      /// @brief Range type of the test space.
      using Range = typename FESType::RangeType;

      /// @brief Device degrees of freedom of the test space.
      using DOFsType = typename FESType::DOFsType;

      /// @brief Type of the integrand.
      using IntegrandType = ::Rodin::Variational::ShapeFunctionBase<NestedDerived,
        FESType, ::Rodin::Variational::ShapeFunctionSpaceType::Test>;

      /// @brief Host integrator interface.
      using Parent = ::Rodin::Variational::LinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = DeviceLinearFormIntegrator<Internal::BackendOf<FES>>;

      /// @brief Device entry type.
      using EntriesType = typename DeviceParent::EntriesType;

      /// @brief Mesh type.
      using MeshType = typename FESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Internal::DevicePolytopes<MeshType>;

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      DeviceLinearIntegral(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      DeviceLinearIntegral(const DeviceLinearIntegral& other)
        : Parent(other),
          DeviceParent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      DeviceLinearIntegral(DeviceLinearIntegral&& other)
        : Parent(std::move(other)),
          DeviceParent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      /// @brief Gets the integrand. @return Test operator being integrated.
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      const ::Rodin::Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      DeviceLinearIntegral& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t local) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      DeviceLayout getLayout() const override
      {
        const auto& fes = getIntegrand().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&fes.getMesh(), 0, static_cast<Index>(fes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(), 0,
          fes.getDeviceDOFs().getStride(), 0,
          DOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /**
       * @brief Locates the entries of every selected polytope.
       * @return Number of entries, each selected polytope with n scalar
       * degrees of freedom owning n * vdim entries.
       */
      size_t prepare() const override
      {
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        return m_polytopes.prepare(fes.getMesh(), this->getRegion(),
          this->getAttributes(), ConstantDOFs{}, 1, fes.getDeviceDOFs(),
          fes.getVectorDimension(), [&](const ::Rodin::Geometry::Polytope& polytope) {
            return getOrder(integrand, polytope);
          });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /// @brief Integrates every selected polytope on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const override
      {
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const DeviceView<IntegrandType> test(integrand);
        const DOFsType dofs = fes.getDeviceDOFs();
        const Index vdim = fes.getVectorDimension();
        const auto& quadrature = m_polytopes.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        constexpr bool p1 = std::is_same_v<DOFsType, VertexDOFs>;
        // Polytopes of one geometry filling the space take kernels of static size.
        if ((!p1 || vdim == VectorCapacity) &&
          quadrature.visitStaticGeometry(
            fes.getMesh().getSpaceDimension(), [&](auto shape) {
              using Shape = decltype(shape);
              constexpr Index N =
                Rodin::Kokkos::Geometry::StaticBasis<Shape::Value>::Count;
              if constexpr (DOFsType::getCapacity(N) <= IntegrationPointView::Capacity)
              {
                integratePolytopes<N, p1, Shape>(fes.getMesh(), m_polytopes, test, dofs,
                  vdim, values, offset, sign, accumulate, indices);
                return true;
              }
              return false;
            }))
          return;
        if constexpr (p1)
        {
          if (vertices == 4 && quadrature.getMinimumVertexCount() == 4 &&
            vdim == VectorCapacity)
          {
            integratePolytopes<4, true>(fes.getMesh(), m_polytopes, test, dofs, vdim,
              values, offset, sign, accumulate, indices);
            return;
          }
        }
        if (vertices <= 4)
        {
          integratePolytopes<4>(fes.getMesh(), m_polytopes, test, dofs, vdim, values,
            offset, sign, accumulate, indices);
          return;
        }
        integratePolytopes<IntegrationPointView::Capacity>(fes.getMesh(), m_polytopes,
          test, dofs, vdim, values, offset, sign, accumulate, indices);
      }

      /**
       * @brief Integrates the polytopes of at most Vertices vertices.
       *
       * When Exact, the space is P1, every polytope has Vertices vertices and
       * the vector dimension is maximal, so the local loop has a static bound.
       * Shape is the geometry tag of the integration points. Public as the
       * enclosing function of a device lambda.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Test>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const Test& test, const DOFsType& dofs, Index vdim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& indices)
      {
        constexpr Index Capacity = DOFsType::getCapacity(Vertices) * VectorCapacity;
        const auto& execution = mesh.getContext().getExecutionSpace();
        const auto view = mesh.getDeviceView();
        const auto quadrature = polytopes.getQuadrature();
        const auto attributes = polytopes.getAttributes();
        const auto offsets = polytopes.getOffsets();
        const size_t spacing = polytopes.getSpacing();
        const auto region = polytopes.getRegion();
        const Index dimension = polytopes.getDimension();
        const bool index = indices.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin linear form integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()),
          KOKKOS_LAMBDA(Index p) {
            const auto polytope = view.getPolytope(dimension, p);
            if (!Rodin::Kokkos::Variational::isIntegrated(
                  view, region, polytope, attributes))
              return;
            KOKKOS_ASSERT(Exact ? polytope.getVertices().size() == Vertices
                                : polytope.getVertices().size() <= Vertices);
            const auto testIndices = dofs.bind(polytope);
            const Index nte = Exact ? Capacity : testIndices.size() * vdim;
            const Index stride = Exact ? VectorCapacity : vdim;
            const Index g = static_cast<Index>(polytope.getGeometry());
            ::Kokkos::Array<Real, Capacity> local{};
            const Index begin = quadrature.begin(g);
            const Index end = quadrature.end(g);
            if (begin < end)
            {
              IntegrationPointView ip(
                polytope, quadrature.getReference(polytope, begin), Shape{});
              for (Index q = begin; q < end; ++q)
              {
                if (q > begin)
                  ip.setPoint(quadrature.getReference(polytope, q), Shape{});
                const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
                const auto basis = Rodin::Kokkos::Variational::bindBasis(test, ip);
                RODIN_KOKKOS_UNROLL
                for (Index te = 0; te < nte; ++te)
                  local[te] += wdet * basis.getBasis(te / stride, te % stride);
              }
            }
            // Local index a * vdim + r is component r of the scalar degree of
            // freedom a, numbered globally as dof + r * stride.
            const size_t base = offset + offsets(p);
            RODIN_KOKKOS_UNROLL
            for (Index te = 0; te < nte; ++te)
            {
              const size_t e = base + te * spacing;
              const Real v = sign * local[te];
              if (accumulate)
                values(e) += v;
              else
                values(e) = v;
              if (index)
                indices(e) = testIndices[te / stride] + (te % stride) * dofs.getStride();
            }
          });
      }

      DeviceLinearIntegral* copy() const noexcept override = 0;

    private:
      /// @brief Maximum number of components of a value of the range.
      static constexpr Index VectorCapacity =
        std::is_same_v<Range, Real> ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION;

      /**
       * @brief Selects the quadrature order on a polytope.
       * @param integrand Integrated test operator.
       * @param polytope Host polytope.
       * @return Order chosen as by the generic host QuadratureRule: the
       * integrand order, else the test function order, else the element
       * degree.
       */
      size_t getOrder(
        const IntegrandType& integrand, const ::Rodin::Geometry::Polytope& polytope) const
      {
        auto order = this->Parent::getOrder(polytope);
        if (!order)
          order = integrand.getOrder(polytope);
        if (!order)
          order = integrand.getLeaf().getOrder(polytope);
        if (!order)
        {
          order = integrand.getFiniteElementSpace()
                    .getFiniteElement(polytope.getDimension(), polytope.getIndex())
                    .getDegree();
        }
        if (!order)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic quadrature order is established. "
            << "Specify an integration order with setOrder()." << Alert::Raise;
        }
        return *order;
      }

      std::unique_ptr<IntegrandType> m_integrand;
      mutable PolytopesType m_polytopes;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the cells of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(u) \cdot B(v) \ dx
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::DeviceBilinearIntegral.
   */
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
    requires Rodin::Kokkos::Variational::DeviceSpace<TrialFES> &&
    Rodin::Kokkos::Variational::DeviceSpace<TestFES>
  class Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
              RHSDerived, TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
        RHSDerived, TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      Integral(const LHSType& lhs, const RHSType& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      Integral(const Integral&) = default;

      /// @brief Move constructor.
      Integral(Integral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the cells of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(v) \ dx
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::DeviceLinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::DeviceSpace<FES>
  class Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      Integral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      Integral(const Integral&) = default;

      /// @brief Move constructor.
      Integral(Integral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the boundary faces of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{B}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::DeviceBilinearIntegral.
   */
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
    requires Rodin::Kokkos::Variational::DeviceSpace<TrialFES> &&
    Rodin::Kokkos::Variational::DeviceSpace<TestFES>
  class BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
              RHSDerived, TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
        RHSDerived, TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      BoundaryIntegral(const LHSType& lhs, const RHSType& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      BoundaryIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      BoundaryIntegral(const BoundaryIntegral&) = default;

      /// @brief Move constructor.
      BoundaryIntegral(BoundaryIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Boundary, the integral running over boundary faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the boundary faces of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{B}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::DeviceLinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::DeviceSpace<FES>
  class BoundaryIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      BoundaryIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : BoundaryIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      BoundaryIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      BoundaryIntegral(const BoundaryIntegral&) = default;

      /// @brief Move constructor.
      BoundaryIntegral(BoundaryIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Boundary, the integral running over boundary faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Boundary;
      }

      BoundaryIntegral* copy() const noexcept override
      {
        return new BoundaryIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over every face of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{F}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::DeviceBilinearIntegral.
   */
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
    requires Rodin::Kokkos::Variational::DeviceSpace<TrialFES> &&
    Rodin::Kokkos::Variational::DeviceSpace<TestFES>
  class FaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
              RHSDerived, TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
        RHSDerived, TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      FaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      FaceIntegral(const FaceIntegral&) = default;

      /// @brief Move constructor.
      FaceIntegral(FaceIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Faces, the integral running over every face.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over every face of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{F}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::DeviceLinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::DeviceSpace<FES>
  class FaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      FaceIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : FaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      FaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      FaceIntegral(const FaceIntegral&) = default;

      /// @brief Move constructor.
      FaceIntegral(FaceIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Faces, the integral running over every face.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Faces;
      }

      FaceIntegral* copy() const noexcept override
      {
        return new FaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product over the interface faces of a
   * Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{I}_h} A(u) \cdot B(v) \ d\sigma(x)
   * @f]
   * for device trial and test spaces on a resident Kokkos mesh, as described in
   * Rodin::Kokkos::Variational::DeviceBilinearIntegral.
   */
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
    requires Rodin::Kokkos::Variational::DeviceSpace<TrialFES> &&
    Rodin::Kokkos::Variational::DeviceSpace<TestFES>
  class InterfaceIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
              RHSDerived, TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceBilinearIntegral<LHSDerived,
        RHSDerived, TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of two operators.
       * @param[in] lhs Trial operator @f$ A(u) @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      InterfaceIntegral(const LHSType& lhs, const RHSType& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      InterfaceIntegral(const InterfaceIntegral&) = default;

      /// @brief Move constructor.
      InterfaceIntegral(InterfaceIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Interface, the integral running over interface faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator over the interface faces of a Kokkos
   * mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{I}_h} A(v) \ d\sigma(x)
   * @f]
   * for a scalar valued operator of a device test function on a resident Kokkos
   * mesh, as described in Rodin::Kokkos::Variational::DeviceLinearIntegral.
   */
  template <class NestedDerived, class FES>
    requires ScalarValuedIntegrand<ShapeFunctionBase<NestedDerived, FES, TestSpace>> &&
    Rodin::Kokkos::Variational::DeviceSpace<FES>
  class InterfaceIntegral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DeviceLinearIntegral<NestedDerived, FES>;

      /// @brief Finite element space of the test functions.
      using FESType = typename Parent::FESType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of the dot product of a function and a
       * test operator.
       * @param[in] lhs Function @f$ f @f$.
       * @param[in] rhs Test operator @f$ B(v) @f$.
       */
      template <class LHSDerived, class RHSDerived>
      InterfaceIntegral(const FunctionBase<LHSDerived>& lhs,
        const ShapeFunctionBase<RHSDerived, FESType, TestSpace>& rhs)
        : InterfaceIntegral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      InterfaceIntegral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      InterfaceIntegral(const InterfaceIntegral&) = default;

      /// @brief Move constructor.
      InterfaceIntegral(InterfaceIntegral&&) = default;

      /**
       * @brief Gets the integration region.
       * @return Region::Interface, the integral running over interface faces.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Interface;
      }

      InterfaceIntegral* copy() const noexcept override
      {
        return new InterfaceIntegral(*this);
      }
  };
}

#endif
