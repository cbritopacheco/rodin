/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_LINEARINTEGRAL_H
#define RODIN_KOKKOS_VARIATIONAL_LINEARINTEGRAL_H

/**
 * @file
 * @brief Device integration of linear finite-element forms.
 */

#include <algorithm>
#include <array>
#include <cstdint>
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

#include "Rodin/Kokkos/Variational/View.h"
#include "Rodin/Kokkos/Variational/BilinearFormIntegrator.h"
#include "Rodin/Kokkos/Variational/LinearFormIntegrator.h"

#include "Rodin/Kokkos/Geometry/Region.h"
#include "Rodin/Kokkos/Assembly/Polytopes.h"
#include "Integral/Traits.h"

namespace Rodin::Kokkos::Variational
{
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
   * getRegion(). The operator is lowered once to View, then a kernel
   * integrates every polytope of the region with one thread and writes its
   * local vector as device entries.
   *
   * A face is integrated with the basis of the face, whose degrees of
   * freedom are the ones of the face, as on the host.
   *
   * The quadrature order is selected on the host as by the generic
   * QuadratureRule, once per reference geometry, as described in
   * Quadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   *
   * Rodin::Variational::Integral, BoundaryIntegral, FaceIntegral and
   * InterfaceIntegral specialize this class for device spaces.
   */
  template <class NestedDerived, Space FES>
  class LinearIntegral : public LinearFormIntegrator<Internal::BackendOf<FES>>
  {
    public:
      /// @brief Finite element space of the test functions.
      using FESType = FES;

      /// @brief Range type of the test space.
      using Range = typename FESType::RangeType;

      /// @brief Device degrees of freedom of the test space.
      using DOFsType = ::Rodin::Kokkos::View<FESType>;

      /// @brief Type of the integrand.
      using IntegrandType = ::Rodin::Variational::ShapeFunctionBase<NestedDerived,
        FESType, ::Rodin::Variational::ShapeFunctionSpaceType::Test>;

      /// @brief Kokkos integrator base.
      using Parent = LinearFormIntegrator<Internal::BackendOf<FES>>;

      /// @brief Bulk integration interface.
      using ViewParent = LinearFormIntegrator<Internal::BackendOf<FES>>;

      /// @brief Device entry type.
      using EntriesType = typename ViewParent::EntriesType;

      /// @brief Mesh type.
      using MeshType = typename FESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = Quadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = ::Rodin::Assembly::Polytopes<MeshType>;

      /**
       * @brief Constructs the integral of a test operator.
       * @param[in] integrand Scalar valued test operator.
       */
      LinearIntegral(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      LinearIntegral(const LinearIntegral& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      LinearIntegral(LinearIntegral&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      /// @brief Gets the integrand. @return Test operator being integrated.
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @return Never returns.
       * @throws Alert::Exception Always.
       */
      const ::Rodin::Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @param polytope Ignored host polytope.
       * @return Never returns.
       * @throws Alert::Exception Always.
       */
      LinearIntegral& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @param local Ignored local test index.
       * @return Never returns.
       * @throws Alert::Exception Always.
       */
      Real integrate([[maybe_unused]] size_t local) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Layout getLayout() const override
      {
        const auto& fes = getIntegrand().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        Layout layout{&fes.getMesh(), 0,
          static_cast<Index>(DOFsType(fes).getComponentCount()),
          {attributes.begin(), attributes.end()}, this->getRegion(), 0,
          DOFsType(fes).getStride(), 0,
          DOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
        if constexpr (!DOFsType::Componentwise)
          layout.testSpace = &fes;
        return layout;
      }

      /**
       * @brief Locates the entries of every selected polytope.
       * @return Number of entries, each selected polytope with n scalar
       * degrees of freedom owning n * vdim entries.
       */
      size_t prepare() override
      {
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        return m_polytopes.prepare(fes.getMesh(), this->getRegion(),
          this->getAttributes(), ConstantDOFs{}, 1, DOFsType(fes),
          DOFsType(fes).getComponentCount(),
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
       * @param values Entry values, the entry of rank r being written at
       * offset + getEntryPosition(r, count, getSpacing()).
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values
       * instead of replacing them.
       * @param indices Entry indices, test degrees of freedom, written as the
       * values unless empty.
       */
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const override
      {
        const auto& integrand = getIntegrand();
        const auto& fes = integrand.getFiniteElementSpace();
        const View<IntegrandType> test(integrand);
        const DOFsType dofs = DOFsType(fes);
        const Index vdim = DOFsType(fes).getComponentCount();
        const auto& quadrature = m_polytopes.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        constexpr bool p1 = std::is_same_v<typename FESType::DOFsType, VertexDOFs>;
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
       * Shape is the geometry tag of the integration points.
       *
       * @param mesh Resident mesh.
       * @param polytopes Selected polytopes, located by prepare().
       * @param test Device view of the test operator.
       * @param dofs Device degrees of freedom of the test space.
       * @param vdim Vector dimension of the test space.
       * @param values Entry values.
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values.
       * @param indices Entry indices, written unless empty.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Test>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const Test& test, const DOFsType& dofs, Index vdim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& indices)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const PolytopeKernel<Vertices, Exact, Shape, Test> kernel{mesh.getView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(), test,
          dofs, vdim, values, offset, sign, accumulate, indices.extent(0) > 0, indices};
        ::Kokkos::parallel_for("Rodin linear form integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      LinearIntegral* copy() const noexcept override = 0;

    private:
      /// @brief Maximum number of independently numbered component copies.
      static constexpr Index VectorCapacity = DOFsType::ComponentCapacity;

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
      PolytopesType m_polytopes;

    public:
      /// @brief Kernel of integratePolytopes(), one thread per polytope.
      template <Index Vertices, bool Exact, class Shape, class Test>
      struct PolytopeKernel
      {
        /// @brief Maximum number of local test indices of a polytope.
        static constexpr Index Capacity =
          DOFsType::getCapacity(Vertices) * VectorCapacity;

        /// @brief Device mesh.
        ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase> view;

        /// @brief Quadrature tables of the polytopes.
        QuadratureType quadrature;

        /// @brief Selected attributes; empty selects every polytope of the region.
        typename PolytopesType::AttributesType attributes;

        /// @brief Entry offset of every polytope, relative to offset.
        ::Kokkos::View<Index*, MemorySpace> offsets;

        /// @brief Distance between consecutive entries of a polytope.
        size_t spacing;

        /// @brief Integration region.
        ::Rodin::Geometry::Region region;

        /// @brief Dimension of the integrated polytopes.
        Index dimension;

        /// @brief Device view of the test operator.
        Test test;

        /// @brief Device degrees of freedom of the test space.
        DOFsType dofs;

        /// @brief Vector dimension of the test space.
        Index vdim;

        /// @brief Entry values.
        ::Kokkos::View<Real*, MemorySpace> values;

        /// @brief Position of the first entry of the integral.
        size_t offset;

        /// @brief Factor of the contributions.
        Real sign;

        /// @brief Whether the contributions are added to the values.
        bool accumulate;

        /// @brief Whether the indices of the entries are written.
        bool index;

        /// @brief Entry indices, test degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> indices;

        /**
         * @brief Integrates a polytope if it is selected.
         * @param p Polytope index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Geometry::isIntegrated(
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
            Rodin::Kokkos::Variational::KernelPointView<Shape> ip(
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
        }
      };
  };
}

#endif
