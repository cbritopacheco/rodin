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
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/LinearFormIntegrator.h"

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

    /// @brief Selected attributes; empty selects every cell.
    std::vector<::Rodin::Geometry::Attribute> attributes;

    bool operator==(const DeviceLayout&) const = default;
  };

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
       * @brief Integrates the form over its region.
       * @param values Entry values, entry e being written at offset + e.
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
       * indexing trial degrees of freedom.
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
       * @brief Integrates the form over its region.
       * @param values Entry values, entry e being written at offset + e.
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
       * @return Local contributions, indexed by test degrees of freedom.
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
}

namespace Rodin::Kokkos::Variational::Internal
{
  /**
   * @brief Whether the host integrates a trial/test pair with a P1 derivative
   * handler.
   *
   * The host P1 handlers of forms with derivatives on both sides, such as
   * @f$ \nabla u \cdot \nabla v @f$, @f$ \nabla u : \nabla v @f$ or
   * @f$ (\nabla \cdot u) q @f$, fall back to the quadrature order
   * @f$ k_u + k_v - 2 @f$, which the device reproduces.
   */
  template <class LHSDerived, class RHSDerived>
  struct IsP1DerivativeForm
  {
    static constexpr bool Value = false;
  };

  /// @brief Gradient trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Scaled gradient trial operator and gradient test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Divergence trial operator and scalar test function.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    static constexpr bool Value = std::is_same_v<typename TestFES::RangeType, Real>;
  };

  /// @brief Scalar trial function and divergence test operator.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct IsP1DerivativeForm<::Rodin::Variational::ShapeFunction<LHSDerived, TrialFES,
                              ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
    ::Rodin::Variational::Div<::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = std::is_same_v<typename TrialFES::RangeType, Real>;
  };

  /// @brief Jacobian trial and test operators.
  template <class LHSDerived, class RHSDerived, class TrialFES, class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
      TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Scaled Jacobian trial operator and Jacobian test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<RHSDerived,
      TestFES, ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Jacobian trial operator times a function, against a test function.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class TrialFES,
    class TestFES>
  struct IsP1DerivativeForm<
    ::Rodin::Variational::Mult<
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Jacobian<::Rodin::Variational::ShapeFunction<LHSDerived,
          TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        TrialFES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>,
      ::Rodin::Variational::FunctionBase<CoefficientDerived>>,
    ::Rodin::Variational::ShapeFunction<RHSDerived, TestFES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>
  {
    static constexpr bool Value = true;
  };

  /**
   * @brief Quadrature tables of a resident mesh, one formula per reference
   * geometry.
   *
   * The quadrature order is selected on the host, once per reference
   * geometry, using the first cell of each geometry. It must therefore
   * depend only on the geometry, which holds for meshes whose cells of a
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
       * @param order Callable giving the quadrature order on a host cell.
       */
      template <class Order>
      DeviceQuadrature(const Mesh& mesh, const Order& order)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t dimension = mesh.getDimension();
        m_representatives = getRepresentatives(mesh);
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
       * @param order Callable giving the quadrature order on a host cell.
       * @return Whether every geometry present has the same order.
       */
      template <class Order>
      bool hasOrders(const Mesh& mesh, const Order& order) const
      {
        const size_t dimension = mesh.getDimension();
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

      /// @brief Gets the largest vertex count of the cells. @return Vertex count.
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
       * @brief Maps a quadrature point onto a cell.
       * @param cell Device cell.
       * @param q Quadrature point of the cell geometry.
       * @return Point on the cell.
       */
      KOKKOS_INLINE_FUNCTION
      Rodin::Kokkos::Geometry::PointView getPoint(
        const Rodin::Kokkos::Geometry::PolytopeView& cell, Index q) const
      {
        const Index rdim = cell.getDimension();
        Rodin::Kokkos::Geometry::Vector reference(rdim);
        for (Index k = 0; k < rdim; ++k)
          reference(k) = m_points(q * ReferenceStride + k);
        return cell.getPoint(reference);
      }

      /**
       * @brief Maps a quadrature point onto a cell.
       * @param cell Device cell.
       * @param q Quadrature point of the cell geometry.
       * @return Integration point on the cell.
       */
      KOKKOS_INLINE_FUNCTION
      IntegrationPointView getIntegrationPoint(
        const Rodin::Kokkos::Geometry::PolytopeView& cell, Index q) const
      {
        return IntegrationPointView(getPoint(cell, q));
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
       * @brief Finds the first cell of every reference geometry.
       * @param mesh Resident mesh.
       * @return Index of the first cell per geometry, empty when absent.
       */
      static std::array<Optional<Index>, GeometryCount> getRepresentatives(
        const Mesh& mesh)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        constexpr Index None = std::numeric_limits<Index>::max();
        ::Kokkos::View<Index[GeometryCount], MemorySpace> first(
          ::Kokkos::view_alloc(execution, std::string("Rodin geometry representatives")));
        ::Kokkos::deep_copy(execution, first, None);
        const auto view = mesh.getDeviceView();
        findRepresentatives(execution, view, first, mesh.getCellCount());
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
       * @brief Records the minimum cell index of every geometry.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param first Per-geometry minimum, initialized to the maximum index.
       * @param count Cell count.
       *
       * Public as the enclosing function of a device lambda.
       */
      static void findRepresentatives(const ExecutionSpace& execution,
        const Rodin::Kokkos::Geometry::MeshView& view,
        const ::Kokkos::View<Index[GeometryCount], MemorySpace>& first, Index count)
      {
        ::Kokkos::parallel_for(
          "Rodin geometry representatives", policy(execution, count),
          KOKKOS_LAMBDA(Index c) {
            const Index g = static_cast<Index>(view.getCell(c).getGeometry());
            ::Kokkos::atomic_min(&first(g), c);
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
   * @brief Cells integrated by a device integrator and the position of their
   * entries.
   *
   * The quadrature tables and the attribute selection are uploaded on the
   * first integration and again only when the mesh, the quadrature orders
   * or the selection change, so that repeated integrations do not allocate
   * device memory.
   *
   * @tparam Mesh Resident Kokkos mesh type.
   */
  template <class Mesh>
  class DeviceCells
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
       * @brief Locates the entries of every selected cell.
       * @param mesh Resident mesh.
       * @param attributes Selected attributes; empty selects every cell.
       * @param stride Entry count of a cell is stride times the square of its
       * vertex count for bilinear forms, or stride times its vertex count.
       * @param square Whether the entry count grows with the square of the
       * vertex count.
       * @param order Callable giving the quadrature order on a host cell.
       * @return Number of entries.
       */
      template <class Order>
      size_t prepare(const Mesh& mesh,
        const FlatSet<::Rodin::Geometry::Attribute>& attributes, Index stride,
        bool square, const Order& order)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const Index cellCount = mesh.getCellCount();
        if (m_mesh != &mesh || m_cellCount != cellCount || !m_quadrature ||
          !m_quadrature->hasOrders(mesh, order))
        {
          m_quadrature.emplace(mesh, order);
          m_offsets = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin entry offsets")),
            cellCount);
          m_mesh = &mesh;
          m_cellCount = cellCount;
        }
        if (!m_selection || *m_selection != attributes)
        {
          m_attributes = uploadAttributes(execution, attributes);
          m_selection = attributes;
        }
        m_count = locate(execution, mesh.getDeviceView(), m_attributes, m_offsets,
          cellCount, stride, square);
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

      /// @brief Gets the position of the first entry of every cell. @return Offsets.
      const ::Kokkos::View<Index*, MemorySpace>& getOffsets() const
      {
        return m_offsets;
      }

      /// @brief Gets the entry count of the last integration. @return Count.
      size_t getCount() const
      {
        return m_count;
      }

      /**
       * @brief Computes the entry offsets of the cells.
       *
       * Public as the enclosing function of a device lambda.
       */
      static size_t locate(const ExecutionSpace& execution,
        const Rodin::Kokkos::Geometry::MeshView& view, const AttributesType& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& offsets, Index cellCount, Index stride,
        bool square)
      {
        Index total = 0;
        ::Kokkos::parallel_scan(
          "Rodin entry offsets", QuadratureType::policy(execution, cellCount),
          KOKKOS_LAMBDA(Index c, Index & update, bool final) {
            if (final)
              offsets(c) = update;
            const auto cell = view.getCell(c);
            if (Rodin::Kokkos::Variational::isSelected(cell, attributes))
            {
              const Index n = cell.getVertices().size();
              update += (square ? n * n : n) * stride;
            }
          },
          total);
        return total;
      }

    private:
      const Mesh* m_mesh = nullptr;
      Index m_cellCount = 0;
      Optional<QuadratureType> m_quadrature;
      Optional<FlatSet<::Rodin::Geometry::Attribute>> m_selection;
      AttributesType m_attributes;
      ::Kokkos::View<Index*, MemorySpace> m_offsets;
      size_t m_count = 0;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a trial/test dot product on a Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(u) \cdot B(v) \ dx
   * @f]
   * for scalar or vector P1 trial and test spaces on a resident Kokkos mesh,
   * possibly different, such as the velocity and pressure spaces of a
   * Stokes problem. The
   * operators @f$ A @f$ and @f$ B @f$ are lowered once to
   * Rodin::Kokkos::Variational::DeviceView, then a kernel integrates every
   * cell with one thread and writes its local matrix as device triplets.
   *
   * The quadrature order is selected on the host as by the QuadratureRule
   * the host would use, once per reference geometry, as described in
   * Rodin::Kokkos::Variational::Internal::DeviceQuadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   */
  template <class LHSDerived, class RHSDerived, class TrialRange, class TestRange,
    class Backend>
  class Integral<
    Dot<ShapeFunctionBase<LHSDerived,
          Rodin::Kokkos::Variational::DeviceP1Space<TrialRange, Backend>, TrialSpace>,
      ShapeFunctionBase<RHSDerived,
        Rodin::Kokkos::Variational::DeviceP1Space<TestRange, Backend>, TestSpace>>>
    final : public LocalBilinearFormIntegratorBase<Real>,
            public Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>
  {
    public:
      /// @brief Finite element space of the trial functions.
      using TrialFESType = Rodin::Kokkos::Variational::DeviceP1Space<TrialRange, Backend>;

      /// @brief Finite element space of the test functions.
      using TestFESType = Rodin::Kokkos::Variational::DeviceP1Space<TestRange, Backend>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFESType, TrialSpace>;

      /// @brief Type of the right operand of the dot product.
      using RHSType = ShapeFunctionBase<RHSDerived, TestFESType, TestSpace>;

      /// @brief Type of the integrand.
      using IntegrandType = Dot<LHSType, RHSType>;

      /// @brief Host integrator interface.
      using Parent = LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent =
        Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>;

      /// @brief Device triplet type.
      using TripletsType = typename DeviceParent::TripletsType;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos execution space.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType =
        Rodin::Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated cells.
      using CellsType = Rodin::Kokkos::Variational::Internal::DeviceCells<MeshType>;

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
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      Integral(const Integral& other)
        : Parent(other),
          DeviceParent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      Integral(Integral&& other)
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
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      Integral& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope) override
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

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Rodin::Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&trialfes.getMesh(), static_cast<Index>(trialfes.getVectorDimension()),
          static_cast<Index>(testfes.getVectorDimension()),
          {attributes.begin(), attributes.end()}};
      }

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries, each selected cell with n vertices owning
       * n * trialVDim * n * testVDim entries.
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
        return m_cells.prepare(mesh, this->getAttributes(),
          trialfes.getVectorDimension() * testfes.getVectorDimension(), true,
          [&](const Geometry::Polytope& polytope) {
            return getOrder(integrand, polytope);
          });
      }

      /**
       * @brief Integrates every selected cell on the device.
       *
       * Each selected cell writes its local matrix, the entry of local test
       * index te and local trial index tr at te * ntr + tr past the offset of
       * the cell.
       */
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        using Rodin::Kokkos::Variational::DeviceView;
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const MeshType& mesh = trialfes.getMesh();
        const DeviceView<LHSType> trial(integrand.getLHS());
        const DeviceView<RHSType> test(integrand.getRHS());
        const auto& quadrature = m_cells.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        if (vertices == 4 && quadrature.getMinimumVertexCount() == 4 &&
          trialfes.getVectorDimension() == VectorCapacity<TrialRange> &&
          testfes.getVectorDimension() == VectorCapacity<TestRange>)
        {
          integrateCells<4, true>(mesh, m_cells, trial, test,
            trialfes.getVectorDimension(), testfes.getVectorDimension(), values, offset,
            sign, accumulate, rows, cols);
        }
        else if (vertices <= 4)
        {
          integrateCells<4>(mesh, m_cells, trial, test, trialfes.getVectorDimension(),
            testfes.getVectorDimension(), values, offset, sign, accumulate, rows, cols);
        }
        else
        {
          integrateCells<Rodin::Kokkos::Variational::IntegrationPointView::Capacity>(mesh,
            m_cells, trial, test, trialfes.getVectorDimension(),
            testfes.getVectorDimension(), values, offset, sign, accumulate, rows, cols);
        }
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * When Exact, every cell has Vertices vertices and both vector
       * dimensions are maximal, so the local loops have static bounds.
       * Public as the enclosing function of a device lambda.
       */
      template <Index Vertices, bool Exact = false, class Trial, class Test>
      static void integrateCells(const MeshType& mesh, const CellsType& cells,
        const Trial& trial, const Test& test, Index trialVDim, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        using Rodin::Kokkos::Variational::IntegrationPointView;
        using TestValue = std::decay_t<decltype(test.getBasis(
          std::declval<const IntegrationPointView&>(), Index(0)))>;
        constexpr Index TrialCapacity = Vertices * VectorCapacity<TrialRange>;
        constexpr Index TestCapacity = Vertices * VectorCapacity<TestRange>;
        const Index vertexCount = mesh.getVertexCount();
        const auto& execution = mesh.getContext().getExecutionSpace();
        const auto view = mesh.getDeviceView();
        const auto quadrature = cells.getQuadrature();
        const auto attributes = cells.getAttributes();
        const auto offsets = cells.getOffsets();
        const bool index = rows.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin bilinear form integration",
          QuadratureType::policy(execution, mesh.getCellCount()), KOKKOS_LAMBDA(Index c) {
            const auto cell = view.getCell(c);
            if (!Rodin::Kokkos::Variational::isSelected(cell, attributes))
              return;
            const auto vertices = cell.getVertices();
            KOKKOS_ASSERT(
              Exact ? vertices.size() == Vertices : vertices.size() <= Vertices);
            const Index ntr = Exact ? TrialCapacity : vertices.size() * trialVDim;
            const Index nte = Exact ? TestCapacity : vertices.size() * testVDim;
            const Index trialStride = Exact ? VectorCapacity<TrialRange> : trialVDim;
            const Index testStride = Exact ? VectorCapacity<TestRange> : testVDim;
            const Index g = static_cast<Index>(cell.getGeometry());
            ::Kokkos::Array<Real, TrialCapacity * TestCapacity> local{};
            ::Kokkos::Array<TestValue, TestCapacity> phi;
            const Index begin = quadrature.begin(g);
            const Index end = quadrature.end(g);
            if (begin < end)
            {
              IntegrationPointView ip = quadrature.getIntegrationPoint(cell, begin);
              for (Index q = begin; q < end; ++q)
              {
                if (q > begin)
                  ip.setPoint(quadrature.getPoint(cell, q));
                const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
                const auto testBasis = Rodin::Kokkos::Variational::bindBasis(test, ip);
                const auto trialBasis = Rodin::Kokkos::Variational::bindBasis(trial, ip);
                for (Index te = 0; te < nte; ++te)
                  phi[te] = testBasis.getBasis(te);
                for (Index tr = 0; tr < ntr; ++tr)
                {
                  const auto psi = trialBasis.getBasis(tr);
                  for (Index te = 0; te < nte; ++te)
                  {
                    if (Rodin::Kokkos::Variational::overlaps(psi, phi[te]))
                    {
                      local[te * ntr + tr] +=
                        wdet * Rodin::Kokkos::Variational::dot(psi, phi[te]);
                    }
                  }
                }
              }
            }
            // Local index a * vdim + r is component r of local vertex a,
            // numbered globally as vertex + r * vertexCount.
            const size_t base = offset + offsets(c);
            for (Index te = 0; te < nte; ++te)
            {
              const Index row =
                vertices[te / testStride] + (te % testStride) * vertexCount;
              for (Index tr = 0; tr < ntr; ++tr)
              {
                const size_t e = base + te * ntr + tr;
                const Real v = sign * local[te * ntr + tr];
                if (accumulate)
                  values(e) += v;
                else
                  values(e) = v;
                if (index)
                {
                  rows(e) = row;
                  cols(e) = vertices[tr / trialStride] + (tr % trialStride) * vertexCount;
                }
              }
            }
          });
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      /// @brief Maximum number of components of a value of the given range.
      template <class Range>
      static constexpr Index VectorCapacity =
        std::is_same_v<Range, Real> ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION;

      /**
       * @brief Selects the quadrature order on a polytope.
       * @param integrand Integrated dot product.
       * @param polytope Host polytope.
       * @return Order chosen as by the host QuadratureRule handling the same
       * integrand on a local mesh.
       *
       * The P1 derivative handlers fall back to @f$ k_u + k_v - 2 @f$, the
       * generic handler to the sum of the leaf orders.
       */
      size_t getOrder(
        const IntegrandType& integrand, const Geometry::Polytope& polytope) const
      {
        auto order = this->Parent::getOrder(polytope);
        if (!order)
          order = integrand.getOrder(polytope);
        if (!order &&
          Rodin::Kokkos::Variational::Internal::IsP1DerivativeForm<LHSDerived,
            RHSDerived>::Value)
        {
          const size_t d = polytope.getDimension();
          const Index i = polytope.getIndex();
          const auto& trialfe =
            integrand.getLHS().getFiniteElementSpace().getFiniteElement(d, i);
          const auto& testfe =
            integrand.getRHS().getFiniteElementSpace().getFiniteElement(d, i);
          const size_t trialDegree = trialfe.getDegree().value_or(2 * trialfe.getOrder());
          const size_t testDegree = testfe.getDegree().value_or(2 * testfe.getOrder());
          order =
            (trialDegree == 0 || testDegree == 0) ? 0 : trialDegree + testDegree - 2;
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
      mutable CellsType m_cells;
  };

  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a test operator on a Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} A(v) \ dx
   * @f]
   * for a scalar valued operator @f$ A @f$ of a scalar or vector P1 test
   * function on a resident Kokkos mesh, such as @f$ f \cdot v @f$ or
   * @f$ g \nabla \cdot v @f$. The operator is lowered once to
   * Rodin::Kokkos::Variational::DeviceView, then a kernel integrates every
   * cell with one thread and writes its local vector as device entries.
   *
   * The quadrature order is selected on the host as by the generic
   * QuadratureRule, once per reference geometry, as described in
   * Rodin::Kokkos::Variational::Internal::DeviceQuadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   */
  template <class NestedDerived, class Range, class Backend>
    requires(std::is_same_v<typename FormLanguage::Traits<ShapeFunctionBase<NestedDerived,
                              Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>,
                              TestSpace>>::RangeType,
              Real>)
  class Integral<ShapeFunctionBase<NestedDerived,
    Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>, TestSpace>>
    final : public LinearFormIntegratorBase<Real>,
            public Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>
  {
    public:
      /// @brief Finite element space of the test functions.
      using FESType = Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>;

      /// @brief Type of the integrand.
      using IntegrandType = ShapeFunctionBase<NestedDerived, FESType, TestSpace>;

      /// @brief Host integrator interface.
      using Parent = LinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent =
        Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>;

      /// @brief Device entry type.
      using EntriesType = typename DeviceParent::EntriesType;

      /// @brief Mesh type.
      using MeshType = typename FESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType =
        Rodin::Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated cells.
      using CellsType = Rodin::Kokkos::Variational::Internal::DeviceCells<MeshType>;

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
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      Integral(const Integral& other)
        : Parent(other),
          DeviceParent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      Integral(Integral&& other)
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
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integral has no host evaluation. @throws Alert::Exception Always.
      Integral& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope) override
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

      /**
       * @brief Gets the integration region.
       * @return Region::Cells, the integral running over mesh cells.
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Rodin::Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& fes = getIntegrand().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&fes.getMesh(), 0, static_cast<Index>(fes.getVectorDimension()),
          {attributes.begin(), attributes.end()}};
      }

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries, each selected cell with n vertices owning
       * n * vdim entries.
       */
      size_t prepare() const override
      {
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        return m_cells.prepare(fes.getMesh(), this->getAttributes(),
          fes.getVectorDimension(), false, [&](const Geometry::Polytope& polytope) {
            return getOrder(integrand, polytope);
          });
      }

      /// @brief Integrates every selected cell on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const override
      {
        using Rodin::Kokkos::Variational::DeviceView;
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const DeviceView<IntegrandType> test(integrand);
        const auto& quadrature = m_cells.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        if (vertices == 4 && quadrature.getMinimumVertexCount() == 4 &&
          fes.getVectorDimension() ==
            (std::is_same_v<Range, Real> ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION))
        {
          integrateCells<4, true>(fes.getMesh(), m_cells, test, fes.getVectorDimension(),
            values, offset, sign, accumulate, indices);
        }
        else if (vertices <= 4)
        {
          integrateCells<4>(fes.getMesh(), m_cells, test, fes.getVectorDimension(),
            values, offset, sign, accumulate, indices);
        }
        else
        {
          integrateCells<Rodin::Kokkos::Variational::IntegrationPointView::Capacity>(
            fes.getMesh(), m_cells, test, fes.getVectorDimension(), values, offset, sign,
            accumulate, indices);
        }
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * When Exact, every cell has Vertices vertices and the vector dimension
       * is maximal, so the local loop has a static bound. Public as the
       * enclosing function of a device lambda.
       */
      template <Index Vertices, bool Exact = false, class Test>
      static void integrateCells(const MeshType& mesh, const CellsType& cells,
        const Test& test, Index vdim, const ::Kokkos::View<Real*, MemorySpace>& values,
        size_t offset, Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices)
      {
        using Rodin::Kokkos::Variational::IntegrationPointView;
        constexpr Index Capacity =
          Vertices * (std::is_same_v<Range, Real> ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION);
        const Index vertexCount = mesh.getVertexCount();
        const auto& execution = mesh.getContext().getExecutionSpace();
        const auto view = mesh.getDeviceView();
        const auto quadrature = cells.getQuadrature();
        const auto attributes = cells.getAttributes();
        const auto offsets = cells.getOffsets();
        const bool index = indices.extent(0) > 0;
        ::Kokkos::parallel_for(
          "Rodin linear form integration",
          QuadratureType::policy(execution, mesh.getCellCount()), KOKKOS_LAMBDA(Index c) {
            const auto cell = view.getCell(c);
            if (!Rodin::Kokkos::Variational::isSelected(cell, attributes))
              return;
            const auto vertices = cell.getVertices();
            KOKKOS_ASSERT(
              Exact ? vertices.size() == Vertices : vertices.size() <= Vertices);
            const Index nte = Exact ? Capacity : vertices.size() * vdim;
            const Index stride = Exact ? Capacity / Vertices : vdim;
            const Index g = static_cast<Index>(cell.getGeometry());
            ::Kokkos::Array<Real, Capacity> local{};
            const Index begin = quadrature.begin(g);
            const Index end = quadrature.end(g);
            if (begin < end)
            {
              IntegrationPointView ip = quadrature.getIntegrationPoint(cell, begin);
              for (Index q = begin; q < end; ++q)
              {
                if (q > begin)
                  ip.setPoint(quadrature.getPoint(cell, q));
                const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
                const auto basis = Rodin::Kokkos::Variational::bindBasis(test, ip);
                for (Index te = 0; te < nte; ++te)
                  local[te] += wdet * basis.getBasis(te);
              }
            }
            // Local index a * vdim + r is component r of local vertex a,
            // numbered globally as vertex + r * vertexCount.
            const size_t base = offset + offsets(c);
            for (Index te = 0; te < nte; ++te)
            {
              const Real v = sign * local[te];
              if (accumulate)
                values(base + te) += v;
              else
                values(base + te) = v;
              if (index)
                indices(base + te) = vertices[te / stride] + (te % stride) * vertexCount;
            }
          });
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      /**
       * @brief Selects the quadrature order on a polytope.
       * @param integrand Integrated test operator.
       * @param polytope Host polytope.
       * @return Order chosen as by the generic host QuadratureRule: the
       * integrand order, else the test function order, else the element
       * degree.
       */
      size_t getOrder(
        const IntegrandType& integrand, const Geometry::Polytope& polytope) const
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
      mutable CellsType m_cells;
  };
}

#endif
