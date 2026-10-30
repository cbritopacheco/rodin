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
 * @brief Device integration of bilinear forms on resident Kokkos meshes.
 */

#include <array>
#include <cstdlib>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/Dot.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

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
   * @brief Interface of bilinear form integrators evaluated on the device.
   *
   * Assembly::Kokkos discovers this interface on the local integrators of a
   * bilinear form and asks each one for its contributions as device triplets.
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

      /**
       * @brief Integrates the form over its region.
       * @return Local contributions, with rows indexing test and columns
       * indexing trial degrees of freedom.
       */
      virtual TripletsType assemble() const = 0;
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
   * @brief Whether the host integrates a trial/test pair with a P1 gradient
   * handler.
   *
   * The host P1 handlers of @f$ \nabla u \cdot \nabla v @f$ and
   * @f$ f \nabla u \cdot \nabla v @f$ have their own quadrature order
   * fallback, which the device reproduces.
   */
  template <class LHSDerived, class RHSDerived>
  struct IsP1GradientForm
  {
    static constexpr bool Value = false;
  };

  /// @brief Gradient trial and test operators.
  template <class LHSDerived, class RHSDerived, class FES>
  struct IsP1GradientForm<
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived, FES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, FES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
  };

  /// @brief Scaled gradient trial operator and gradient test operator.
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class FES>
  struct IsP1GradientForm<
    ::Rodin::Variational::Mult<::Rodin::Variational::FunctionBase<CoefficientDerived>,
      ::Rodin::Variational::ShapeFunctionBase<
        ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<LHSDerived, FES,
          ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
        FES, ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>,
    ::Rodin::Variational::Grad<::Rodin::Variational::ShapeFunction<RHSDerived, FES,
      ::Rodin::Variational::ShapeFunctionSpaceType::Test>>>
  {
    static constexpr bool Value = true;
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
   * for scalar P1 trial and test spaces on a resident Kokkos mesh. The
   * operators @f$ A @f$ and @f$ B @f$ are lowered once to
   * Rodin::Kokkos::Variational::DeviceView, then a kernel integrates every
   * cell with one thread and writes its local matrix as device triplets.
   *
   * The quadrature order is selected on the host as by the QuadratureRule
   * the host would use, once per reference geometry, using the first cell of
   * each geometry. It must therefore depend only on the geometry, which
   * holds for meshes whose cells of a given geometry share one
   * transformation degree.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   */
  template <class LHSDerived, class RHSDerived, class Backend>
  class Integral<Dot<ShapeFunctionBase<LHSDerived,
                       Rodin::Kokkos::Variational::DeviceP1<Backend>, TrialSpace>,
    ShapeFunctionBase<RHSDerived, Rodin::Kokkos::Variational::DeviceP1<Backend>,
      TestSpace>>>
    final : public LocalBilinearFormIntegratorBase<Real>,
            public Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>
  {
    public:
      /// @brief Finite element space of the trial and test functions.
      using FESType = Rodin::Kokkos::Variational::DeviceP1<Backend>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = ShapeFunctionBase<LHSDerived, FESType, TrialSpace>;

      /// @brief Type of the right operand of the dot product.
      using RHSType = ShapeFunctionBase<RHSDerived, FESType, TestSpace>;

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
      using MeshType = typename FESType::MeshType;

      /// @brief Kokkos execution space.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

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

      /**
       * @brief Integrates every selected cell on the device.
       * @return Local matrices as triplets, with one entry per pair of local
       * trial and test degrees of freedom of each cell.
       */
      TripletsType assemble() const override
      {
        using Rodin::Kokkos::Variational::DeviceView;
        using Rodin::Kokkos::Variational::IntegrationPointView;

        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getLHS().getFiniteElementSpace();
        const MeshType& mesh = fes.getMesh();
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t dimension = mesh.getDimension();
        const Index cellCount = mesh.getCellCount();
        const auto view = mesh.getDeviceView();

        // Quadrature tables per reference geometry, flattened with the
        // reference coordinates padded to three components.
        const auto representatives = getRepresentatives(mesh);
        std::vector<Index> quadratureOffsets(GeometryCount + 1, 0);
        std::vector<Real> weights;
        std::vector<Real> points;
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          quadratureOffsets[g] = weights.size();
          if (!representatives[g])
            continue;
          const auto it = mesh.getPolytope(dimension, *representatives[g]);
          const Geometry::Polytope& polytope = *it;
          const auto& qf = QF::QuadratureFormula::get(
            getOrder(integrand, polytope), polytope.getGeometry());
          for (size_t qp = 0; qp < qf.getSize(); ++qp)
          {
            weights.push_back(qf.getWeight(qp));
            const auto& point = qf.getPoint(qp);
            for (size_t k = 0; k < ReferenceStride; ++k)
              points.push_back(
                k < static_cast<size_t>(point.size()) ? point(k) : Real(0));
          }
        }
        quadratureOffsets[GeometryCount] = weights.size();

        const auto& attributeSet = this->getAttributes();
        const std::vector<Geometry::Attribute> attributeList(
          attributeSet.begin(), attributeSet.end());

        const auto qoffsets = Rodin::Kokkos::Internal::upload(
          execution, quadratureOffsets, "Rodin quadrature offsets");
        const auto qweights =
          Rodin::Kokkos::Internal::upload(execution, weights, "Rodin quadrature weights");
        const auto qpoints =
          Rodin::Kokkos::Internal::upload(execution, points, "Rodin quadrature points");
        const auto attributes = Rodin::Kokkos::Internal::upload(
          execution, attributeList, "Rodin integration attributes");

        // Entry offsets of the cells, each selected cell owning n * n entries.
        ::Kokkos::View<Index*, MemorySpace> offsets(
          ::Kokkos::view_alloc(execution, std::string("Rodin triplet offsets")),
          cellCount);
        Index total = 0;
        ::Kokkos::parallel_scan(
          "Rodin triplet offsets", policy(execution, cellCount),
          KOKKOS_LAMBDA(Index c, Index & update, bool final) {
            if (final)
              offsets(c) = update;
            const auto cell = view.getCell(c);
            if (Rodin::Kokkos::Variational::isSelected(cell, attributes))
            {
              const Index n = cell.getVertices().size();
              update += n * n;
            }
          },
          total);

        TripletsType result;
        result.rows = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin triplet rows")), total);
        result.cols = ::Kokkos::View<Index*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin triplet columns")), total);
        result.values = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin triplet values")), total);
        const auto rows = result.rows;
        const auto cols = result.cols;
        const auto values = result.values;

        const DeviceView<LHSType> trial(integrand.getLHS());
        const DeviceView<RHSType> test(integrand.getRHS());
        using TestValue = std::decay_t<decltype(test.getBasis(
          std::declval<const IntegrationPointView&>(), Index(0)))>;
        constexpr Index Capacity = IntegrationPointView::Capacity;

        ::Kokkos::parallel_for(
          "Rodin bilinear form integration", policy(execution, cellCount),
          KOKKOS_LAMBDA(Index c) {
            const auto cell = view.getCell(c);
            if (!Rodin::Kokkos::Variational::isSelected(cell, attributes))
              return;
            const auto vertices = cell.getVertices();
            const Index n = vertices.size();
            const Index rdim = cell.getDimension();
            const Index g = static_cast<Index>(cell.getGeometry());
            ::Kokkos::Array<Real, Capacity * Capacity> local{};
            ::Kokkos::Array<TestValue, Capacity> phi;
            for (Index q = qoffsets(g); q < qoffsets(g + 1); ++q)
            {
              Rodin::Kokkos::Geometry::Vector reference(rdim);
              for (Index k = 0; k < rdim; ++k)
                reference(k) = qpoints(q * ReferenceStride + k);
              const IntegrationPointView ip(cell.getPoint(reference));
              const Real wdet = qweights(q) * ip.getPoint().getDistortion();
              for (Index te = 0; te < n; ++te)
                phi[te] = test.getBasis(ip, te);
              for (Index tr = 0; tr < n; ++tr)
              {
                const auto psi = trial.getBasis(ip, tr);
                for (Index te = 0; te < n; ++te)
                  local[te * n + tr] +=
                    wdet * Rodin::Kokkos::Variational::dot(psi, phi[te]);
              }
            }
            const Index offset = offsets(c);
            for (Index te = 0; te < n; ++te)
            {
              for (Index tr = 0; tr < n; ++tr)
              {
                const Index e = offset + te * n + tr;
                rows(e) = vertices[te];
                cols(e) = vertices[tr];
                values(e) = local[te * n + tr];
              }
            }
          });
        execution.fence("Rodin bilinear form integration");
        return result;
      }

      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      /// @brief Number of reference geometries.
      static constexpr size_t GeometryCount = Geometry::Polytope::Types.size();

      /// @brief Stride of the flattened reference coordinates.
      static constexpr Index ReferenceStride = 3;

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
        const MeshType& mesh)
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

    public:
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
      /**
       * @brief Selects the quadrature order on a polytope.
       * @param integrand Integrated dot product.
       * @param polytope Host polytope.
       * @return Order chosen as by the host QuadratureRule handling the same
       * integrand on a local mesh.
       *
       * The P1 gradient handlers fall back to @f$ k_u + k_v - 2 @f$, the
       * generic handler to the sum of the leaf orders.
       */
      size_t getOrder(
        const IntegrandType& integrand, const Geometry::Polytope& polytope) const
      {
        auto order = this->Parent::getOrder(polytope);
        if (!order)
          order = integrand.getOrder(polytope);
        if (!order &&
          Rodin::Kokkos::Variational::Internal::IsP1GradientForm<LHSDerived,
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
  };
}

#endif
