/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_GRIDFUNCTION_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_GRIDFUNCTION_H

/**
 * @file
 * @brief Device evaluation of grid-function values and derivatives.
 */

#include "Common.h"
#include "FiniteElementSpace.h"
#include "Function.h"
#include "Rodin/Kokkos/Variational/GridFunctionData.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Trace domain of a function, in device memory.
   *
   * A derivative of a grid function evaluated on a face is that of an
   * incident cell, as on the host. A boundary face has a single incident
   * cell. On an interface, the cell is the one whose attribute is in the
   * trace domain, and a kernel evaluating it with no such cell aborts.
   */
  class TraceDomainView
  {
    public:
      /// @brief Maximum number of attributes.
      static constexpr Index Capacity = 16;

      /// @brief Constructs the empty trace domain.
      TraceDomainView() = default;

      /**
       * @brief Copies a trace domain.
       * @param domain Attributes of the cells giving the trace.
       */
      explicit TraceDomainView(const FlatSet<::Rodin::Geometry::Attribute>& domain)
        : m_count(domain.size())
      {
        if (m_count > Capacity)
        {
          Alert::Exception() << "Device trace domains hold at most " << Capacity
                             << " attributes." << Alert::Raise;
        }
        Index i = 0;
        for (const auto& attribute : domain)
          m_attributes[i++] = attribute;
      }

      /**
       * @brief Tests whether a point lies on a face.
       * @param p Device point. @return Whether the point is bound to a face.
       */
      KOKKOS_INLINE_FUNCTION static bool isFace(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p)
      {
        const auto& polytope = p.getPolytope();
        return polytope.getDimension() + 1 == polytope.mesh.getDimension();
      }

      /**
       * @brief Carries a point of a face to the cell giving the trace.
       * @param p Point on a face; face-to-cell incidence must be computed.
       * @return Same physical point, bound to the cell.
       */
      KOKKOS_INLINE_FUNCTION ::Rodin::Kokkos::View<::Rodin::Geometry::Point> getCellPoint(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
      {
        const auto& face = p.getPolytope();
        const Index d = face.mesh.getDimension();
        const auto connectivity = face.mesh.getConnectivity();
        KOKKOS_ASSERT(connectivity.hasIncidence(d - 1, d));
        const auto cells = connectivity.getIncidence(d - 1, d, face.getIndex());
        if (cells.size() == 1)
          return p.getCellPoint(cells[0]);
        for (Index c = 0; c < cells.size(); ++c)
        {
          const auto cell = face.mesh.getCell(cells[c]);
          if (!cell.hasAttribute())
            continue;
          for (Index i = 0; i < m_count; ++i)
          {
            if (m_attributes[i] == cell.getAttribute())
              return p.getCellPoint(cells[c]);
          }
        }
        ::Kokkos::abort(
          "Rodin: derivative at an interface with no cell in the trace domain.");
        return p;
      }

    private:
      Index m_count = 0;
      ::Kokkos::Array<::Rodin::Geometry::Attribute, Capacity> m_attributes = {};
  };

  /**
   * @brief Coefficients of a grid function of a device space in device
   * memory.
   *
   * Backend coefficient vectors are borrowed directly. Explicit host coefficient
   * vectors are staged when the view is built, and PETSc vectors retain their
   * prepared coefficients through their state-aware adapter. Values and
   * derivatives are interpolated from the degrees of freedom of the polytope
   * holding the point, as on the host.
   *
   * @tparam FES Device space of the grid function.
   */
  template <Space FES>
  class GridFunctionView
  {
    public:
      /// @brief Finite element space.
      using FESType = FES;

      /// @brief Device degrees of freedom.
      using DOFsType = typename FESType::DOFsType;

      /// @brief Physical element evaluation and numbering.
      using ElementType = ::Rodin::Kokkos::View<FESType>;

      /// @brief Memory space of the coefficients.
      using MemorySpace = typename FESType::MeshType::MemorySpace;

      /**
       * @brief Prepares the coefficients of a grid function for device evaluation.
       * @param fes Space of the grid function.
       * @param data Coefficients, numbered as by the space.
       */
      template <class Data>
      GridFunctionView(const FESType& fes, const Data& data)
        : m_dofs(fes),
          m_data(prepareData(fes, data))
      {}

      /// @brief Gets the device degrees of freedom. @return Degrees of freedom.
      KOKKOS_INLINE_FUNCTION const DOFsType& getDOFs() const
      {
        return m_dofs.getDOFs();
      }

      /// @brief Gets the coefficients. @return Device coefficients.
      KOKKOS_INLINE_FUNCTION const ::Kokkos::View<const Real*, MemorySpace>&
      getData() const
      {
        return m_data;
      }

      /// @brief Returns the vector dimension. @return Component count.
      KOKKOS_INLINE_FUNCTION Index getVectorDimension() const
      {
        return m_dofs.getVectorDimension();
      }

      /**
       * @brief Interpolates one component.
       * @param p Device point.
       * @param r Component index.
       * @return Component @f$ u_r(p) @f$.
       *
       * The P1 basis is evaluated directly at the point, the other bases
       * through an IntegrationPointView of the point.
       */
      KOKKOS_INLINE_FUNCTION Real interpolate(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p, Index r) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          const auto vertices = p.getPolytope().getVertices();
          Geometry::Vector reference(p.getReferenceCoordinates().size());
          Real result = 0;
          for (Index a = 0; a < vertices.size(); ++a)
          {
            Real phi;
            p.evaluateBasis(a, phi, reference);
            result += phi * m_data(vertices[a] + r * m_dofs.getStride());
          }
          return result;
        }
        else
        {
          return interpolate(IntegrationPointView(p), r);
        }
      }

      /**
       * @brief Interpolates the physical Jacobian.
       * @param p Device point.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       *
       * On a face, the Jacobian is that of an incident cell.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix interpolateJacobian(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p,
        const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(p) ? getCellJacobian(trace.getCellPoint(p))
                                          : getCellJacobian(p);
      }

      /**
       * @brief Interpolates the divergence.
       * @param p Device point.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Trace of interpolateJacobian(p), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION Real interpolateDivergence(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p,
        const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(p) ? getCellDivergence(trace.getCellPoint(p))
                                          : getCellDivergence(p);
      }

      /**
       * @brief Interpolates one component with the basis of an integration
       * point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param r Component index.
       * @return Component @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Real interpolate(
        const IntegrationPointView& ip, Index r) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        Real result = 0;
        forEachDOF(dofs.size(), [&](Index q) {
          if constexpr (ElementType::Componentwise)
          {
            result +=
              m_dofs.getGenerator(ip, q) * m_data(dofs[q] + r * m_dofs.getStride());
          }
          else
            result += m_dofs.getGenerator(ip, q)(r) * m_data(dofs[q]);
        });
        return result;
      }

      /**
       * @brief Interpolates every component with the basis of an integration
       * point.
       * @param ip Integration point on a polytope of the space mesh.
       * @return Vector of the components @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Vector interpolateVector(
        const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        Geometry::Vector result(m_dofs.getVectorDimension());
        forEachDOF(dofs.size(), [&](Index q) {
          const auto basis = m_dofs.getGenerator(ip, q);
          if constexpr (ElementType::Componentwise)
          {
            forEachComponent(m_dofs.getVectorDimension(), [&](Index r) {
              result(r) += basis * m_data(dofs[q] + r * m_dofs.getStride());
            });
          }
          else
          {
            const Real coefficient = m_data(dofs[q]);
            forEachComponent(m_dofs.getVectorDimension(),
              [&](Index r) { result(r) += basis(r) * coefficient; });
          }
        });
        return result;
      }

      /**
       * @brief Interpolates every component.
       * @param p Device point.
       * @return Vector of the components @f$ u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Vector interpolateVector(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
      {
        Geometry::Vector result(m_dofs.getVectorDimension());
        for (Index r = 0; r < m_dofs.getVectorDimension(); ++r)
          result(r) = interpolate(p, r);
        return result;
      }

      /**
       * @brief Interpolates the physical Jacobian with the basis gradients of
       * an integration point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix interpolateJacobian(
        const IntegrationPointView& ip, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(ip.getPoint())
          ? getCellJacobian(trace.getCellPoint(ip.getPoint()))
          : getCellJacobian(ip);
      }

      /**
       * @brief Interpolates the divergence with the basis gradients of an
       * integration point.
       * @param ip Integration point on a polytope of the space mesh.
       * @param trace Trace domain, used when the point lies on an interface.
       * @return Trace of interpolateJacobian(ip), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION Real interpolateDivergence(
        const IntegrationPointView& ip, const TraceDomainView& trace = {}) const
      {
        return TraceDomainView::isFace(ip.getPoint())
          ? getCellDivergence(trace.getCellPoint(ip.getPoint()))
          : getCellDivergence(ip);
      }

      /**
       * @brief Interpolates the physical Jacobian at an integration point of
       * a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @param trace Trace domain, unused since a cell is not an interface.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix interpolateJacobian(
        const CellIntegrationPointView& ip,
        [[maybe_unused]] const TraceDomainView& trace = {}) const
      {
        return getCellJacobian(ip);
      }

      /**
       * @brief Interpolates the divergence at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @param trace Trace domain, unused since a cell is not an interface.
       * @return Trace of interpolateJacobian(ip), summed in the same order.
       */
      KOKKOS_INLINE_FUNCTION Real interpolateDivergence(
        const CellIntegrationPointView& ip,
        [[maybe_unused]] const TraceDomainView& trace = {}) const
      {
        return getCellDivergence(ip);
      }

    private:
      /**
       * @brief Calls a function on each local degree of freedom.
       * @param count Local degree of freedom count.
       * @param f Function of the local index.
       *
       * The P1 loop has a static bound, so that the basis cached by the
       * integration point is indexed statically and stays in registers.
       */
      template <class F>
      KOKKOS_INLINE_FUNCTION static void forEachDOF(Index count, const F& f)
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          RODIN_KOKKOS_UNROLL
          for (Index q = 0; q < IntegrationPointView::Capacity; ++q)
          {
            if (q < count)
              f(q);
          }
        }
        else
        {
          for (Index q = 0; q < count; ++q)
            f(q);
        }
      }

      /**
       * @brief Calls a function on each component, or each coordinate.
       * @param count Component count, at most three.
       * @param f Function of the index.
       *
       * The static bound keeps the device vectors and matrices in registers.
       */
      template <class F>
      KOKKOS_INLINE_FUNCTION static void forEachComponent(Index count, const F& f)
      {
        KOKKOS_ASSERT(count <= 3);
        RODIN_KOKKOS_UNROLL
        for (Index r = 0; r < 3; ++r)
        {
          if (r < count)
            f(r);
        }
      }

      /**
       * @brief Interpolates the physical Jacobian at a point of a cell.
       *
       * The P1 basis is evaluated directly at the point, the other bases
       * through an IntegrationPointView of the point.
       *
       * @param p Point on a cell of the space mesh.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix getCellJacobian(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          const auto vertices = p.getPolytope().getVertices();
          const Geometry::Matrix inverse = p.getJacobianInverse();
          const Index rdim = inverse.rows();
          const Index sdim = inverse.cols();
          Geometry::Matrix result(m_dofs.getVectorDimension(), sdim);
          for (Index a = 0; a < vertices.size(); ++a)
          {
            Real phi;
            Geometry::Vector reference(rdim);
            p.evaluateBasis(a, phi, reference);
            for (Index c = 0; c < sdim; ++c)
            {
              Real gradient = 0;
              for (Index k = 0; k < rdim; ++k)
                gradient += inverse(k, c) * reference(k);
              for (Index r = 0; r < m_dofs.getVectorDimension(); ++r)
                result(r, c) += gradient * m_data(vertices[a] + r * m_dofs.getStride());
            }
          }
          return result;
        }
        else
        {
          return getCellJacobian(IntegrationPointView(p));
        }
      }

      /**
       * @brief Interpolates the divergence at a point of a cell.
       * @param p Point on a cell of the space mesh.
       * @return Trace of getCellJacobian(p).
       */
      KOKKOS_INLINE_FUNCTION Real getCellDivergence(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p) const
      {
        if constexpr (std::is_same_v<DOFsType, VertexDOFs>)
        {
          const auto vertices = p.getPolytope().getVertices();
          const Geometry::Matrix inverse = p.getJacobianInverse();
          const Index rdim = inverse.rows();
          KOKKOS_ASSERT(m_dofs.getVectorDimension() == inverse.cols());
          Real result = 0;
          for (Index r = 0; r < m_dofs.getVectorDimension(); ++r)
          {
            Real diagonal = 0;
            for (Index a = 0; a < vertices.size(); ++a)
            {
              Real phi;
              Geometry::Vector reference(rdim);
              p.evaluateBasis(a, phi, reference);
              Real gradient = 0;
              for (Index k = 0; k < rdim; ++k)
                gradient += inverse(k, r) * reference(k);
              diagonal += gradient * m_data(vertices[a] + r * m_dofs.getStride());
            }
            result += diagonal;
          }
          return result;
        }
        else
        {
          return getCellDivergence(IntegrationPointView(p));
        }
      }

      /**
       * @brief Interpolates the physical Jacobian at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @return Matrix with one row per component, row @f$ r @f$ holding
       * @f$ \nabla u_r(p) @f$.
       */
      KOKKOS_INLINE_FUNCTION Geometry::Matrix getCellJacobian(
        const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        const Index sdim = ip.getGradient(0).size();
        Geometry::Matrix result(m_dofs.getVectorDimension(), sdim);
        forEachDOF(dofs.size(), [&](Index q) {
          if constexpr (ElementType::Componentwise)
          {
            const auto& gradient = m_dofs.getGradient(ip, q);
            forEachComponent(m_dofs.getVectorDimension(), [&](Index r) {
              const Real value = m_data(dofs[q] + r * m_dofs.getStride());
              forEachComponent(
                sdim, [&](Index c) { result(r, c) += gradient(c) * value; });
            });
          }
          else
          {
            const auto jacobian = m_dofs.getJacobian(ip, q, 0);
            const Real value = m_data(dofs[q]);
            forEachComponent(m_dofs.getVectorDimension(), [&](Index r) {
              forEachComponent(
                sdim, [&](Index c) { result(r, c) += jacobian(r, c) * value; });
            });
          }
        });
        return result;
      }

      /**
       * @brief Interpolates the divergence at an integration point of a cell.
       * @param ip Integration point on a cell of the space mesh.
       * @return Trace of getCellJacobian(ip).
       */
      KOKKOS_INLINE_FUNCTION Real getCellDivergence(const IntegrationPointView& ip) const
      {
        const auto dofs = m_dofs.bind(ip.getPoint().getPolytope());
        KOKKOS_ASSERT(m_dofs.getVectorDimension() == ip.getGradient(0).size());
        Real result = 0;
        if constexpr (ElementType::Componentwise)
        {
          forEachComponent(m_dofs.getVectorDimension(), [&](Index r) {
            Real diagonal = 0;
            forEachDOF(dofs.size(), [&](Index q) {
              diagonal +=
                m_dofs.getGradient(ip, q)(r) * m_data(dofs[q] + r * m_dofs.getStride());
            });
            result += diagonal;
          });
        }
        else
        {
          forEachDOF(dofs.size(),
            [&](Index q) { result += m_dofs.getDivergence(ip, q, 0) * m_data(dofs[q]); });
        }
        return result;
      }

      /**
       * @tparam Data Coefficient storage type.
       * @param fes Space of the grid function.
       * @param data Coefficients in the space's numbering.
       * @returns Read-only coefficients prepared for execution.
       */
      template <class Data>
      static auto prepareData(const FESType& fes, const Data& data)
      {
        using Access = GridFunctionData<std::decay_t<Data>>;
        if constexpr (requires { Access::prepareView(fes, data); })
          return Access::prepareView(fes, data);
        else
          return Access::prepare(fes, data);
      }

      ElementType m_dofs;
      ::Kokkos::View<const Real*, MemorySpace> m_data;
  };

  /// @brief Grid function of a device space.
  template <Space FES, class Data>
  struct View<::Rodin::Variational::GridFunction<FES, Data>>
  {
    /// @param node Grid function.
    View(const ::Rodin::Variational::GridFunction<FES, Data>& node)
      : u(node.getFiniteElementSpace(), node.getData())
    {}

    /// @param fes Space of the grid function. @param data Coefficients.
    View(const FES& fes, const Data& data)
      : u(fes, data)
    {}

    /// @param p Device or integration point. @return Function value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getValue(const Point& p) const
    {
      if constexpr (std::is_same_v<typename FES::RangeType, Real>)
      {
        return u.interpolate(p, 0);
      }
      else
      {
        return u.interpolateVector(p);
      }
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;
  };

  /// @brief Grid function referenced by an expression.
  template <Space FES, class Data>
  struct View<::Rodin::Variational::GridFunctionBaseReference<
    ::Rodin::Variational::GridFunction<FES, Data>>>
    : View<::Rodin::Variational::GridFunction<FES, Data>>
  {
    /// @param node Reference to a grid function.
    View(const ::Rodin::Variational::GridFunctionBaseReference<
      ::Rodin::Variational::GridFunction<FES, Data>>& node)
      : View<::Rodin::Variational::GridFunction<FES, Data>>(
          node.getFiniteElementSpace(), node.getData())
    {}
  };

  /// @brief Physical gradient of a scalar grid function.
  template <ScalarSpace FES, class Data>
  struct View<::Rodin::Variational::Grad<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Gradient expression.
    View(
      const ::Rodin::Variational::Grad<::Rodin::Variational::GridFunction<FES, Data>>&
        node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Gradient @f$ \nabla u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Vector getValue(const Point& p) const
    {
      const Geometry::Matrix jacobian = u.interpolateJacobian(p, trace);
      Geometry::Vector result(jacobian.cols());
      for (Index c = 0; c < jacobian.cols(); ++c)
        result(c) = jacobian(0, c);
      return result;
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
  };

  /// @brief Physical Jacobian of a vector grid function.
  template <VectorSpace FES, class Data>
  struct View<
    ::Rodin::Variational::Jacobian<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Jacobian expression.
    View(const ::Rodin::Variational::Jacobian<
      ::Rodin::Variational::GridFunction<FES, Data>>& node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Jacobian @f$ \nabla u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Geometry::Matrix getValue(const Point& p) const
    {
      return u.interpolateJacobian(p, trace);
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
  };

  /// @brief Divergence of a vector grid function.
  template <VectorSpace FES, class Data>
  struct View<::Rodin::Variational::Div<::Rodin::Variational::GridFunction<FES, Data>>>
  {
    /// @param node Divergence expression.
    View(const ::Rodin::Variational::Div<::Rodin::Variational::GridFunction<FES, Data>>&
        node)
      : u(node.getOperand().getFiniteElementSpace(), node.getOperand().getData()),
        trace(node.getTraceDomain())
    {}

    /// @param p Device or integration point. @return Divergence @f$ \nabla \cdot u(p) @f$.
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getValue(const Point& p) const
    {
      return u.interpolateDivergence(p, trace);
    }

    /// @brief Device coefficients.
    GridFunctionView<FES> u;

    /// @brief Trace domain of the derivative on interfaces.
    TraceDomainView trace;
  };

}

#endif
