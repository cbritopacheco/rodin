/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_POTENTIAL_H
#define RODIN_KOKKOS_VARIATIONAL_POTENTIAL_H

/**
 * @file
 * @brief Device integration of potentials tested against shape functions on
 * resident Kokkos meshes.
 */

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/Variational/Potential.h"

#include "Integral.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Device integration of a potential tested against a shape function.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} \int_{\mathcal{T}_h} k(x, y) \, u(x) \cdot v(y)
   *   \ dx \ dy
   * @f]
   * for P1 trial and test spaces on a resident Kokkos mesh, with the
   * quadrature of the host QuadratureRule. A thread integrates one pair of a
   * selected test cell and a selected trial cell. A pair of coincident
   * triangles is integrated by the Sauter-Schwab transformation with a tensor
   * Gauss-Legendre rule, any other pair by the centroid rules of both cells.
   *
   * The kernel is a KOKKOS_LAMBDA. For scalar spaces it takes the trial and
   * the test points, as const Rodin::Kokkos::Geometry::PointView&, and
   * returns a Real. For vector spaces it takes first a
   * Rodin::Kokkos::Geometry::Matrix&, sized to the vector dimension and set
   * to zero, which it fills.
   *
   * The entries of a pair follow the entries of the previous pairs, the
   * pairs being ordered by test cell, then by trial cell, as the host
   * assemblies visit them.
   *
   * The order of the Gauss-Legendre rule is read once, on the first
   * triangle of the mesh, so it must depend only on the geometry.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   */
  template <class Kernel, class LHSDerived, class RHSDerived, DeviceSpace TrialFES,
    DeviceSpace TestFES>
  class DevicePotentialIntegral
    : public ::Rodin::Variational::GlobalBilinearFormIntegratorBase<Real>,
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

      /// @brief Kernel type.
      using KernelType = Kernel;

      /// @brief Type of the left operand of the dot product.
      using LHSType = ::Rodin::Variational::Potential<KernelType,
        ::Rodin::Variational::ShapeFunctionBase<LHSDerived, TrialFESType,
          ::Rodin::Variational::ShapeFunctionSpaceType::Trial>>;

      /// @brief Type of the right operand of the dot product.
      using RHSType = ::Rodin::Variational::ShapeFunctionBase<RHSDerived, TestFESType,
        ::Rodin::Variational::ShapeFunctionSpaceType::Test>;

      /// @brief Type of the integrand.
      using IntegrandType = ::Rodin::Variational::Dot<LHSType, RHSType>;

      /// @brief Host integrator interface.
      using Parent = ::Rodin::Variational::GlobalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = DeviceBilinearFormIntegrator<Internal::BackendOf<TrialFES>>;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos execution space.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device attribute selection.
      using AttributesType = ::Kokkos::View<::Rodin::Geometry::Attribute*, MemorySpace>;

      /// @brief Device quadrature table.
      using TableType = ::Kokkos::View<Real*, MemorySpace>;

      /// @brief Whether the spaces are scalar.
      static constexpr bool Scalar = std::is_same_v<TrialRange, Real>;

      static_assert(std::is_same_v<TrialRange, TestRange>,
        "Kokkos potentials need trial and test spaces of the same range.");

      static_assert(std::is_same_v<TrialDOFsType, VertexDOFs> &&
          std::is_same_v<TestDOFsType, VertexDOFs>,
        "Kokkos potentials are implemented for P1 spaces only.");

      static_assert(Scalar ? std::is_invocable_r_v<Real, const KernelType&,
                               const Geometry::PointView&, const Geometry::PointView&>
                           : std::is_invocable_v<const KernelType&, Geometry::Matrix&,
                               const Geometry::PointView&, const Geometry::PointView&>,
        "The kernel of a Kokkos potential must be a KOKKOS_LAMBDA taking two const "
        "Rodin::Kokkos::Geometry::PointView&, after a Rodin::Kokkos::Geometry::Matrix& "
        "for vector spaces.");

      /// @brief Maximum number of components of the spaces.
      static constexpr Index VectorCapacity = Scalar ? 1 : RODIN_MAXIMAL_SPACE_DIMENSION;

      /// @brief Maximum local count of a triangle, integrated in thread memory.
      static constexpr Index TriangleCapacity = 3 * VectorCapacity;

      /// @brief Number of reference geometries.
      static constexpr size_t GeometryCount = ::Rodin::Geometry::Polytope::Types.size();

      /// @brief Stride of the flattened reference coordinates of the centroids.
      static constexpr Index ReferenceStride = 3;

      /// @brief Selected cells of one side of the pairs.
      struct Selection
      {
        /// @brief Selected cells, in increasing order.
        ::Kokkos::View<Index*, MemorySpace> cells;

        /// @brief Sum of the local counts of the previous selected cells.
        ::Kokkos::View<Index*, MemorySpace> offsets;

        /// @brief Number of selected cells.
        Index count = 0;

        /// @brief Sum of the local counts of the selected cells.
        Index total = 0;

        /// @brief Smallest local count.
        Index smallest = 0;

        /// @brief Largest local count.
        Index largest = 0;
      };

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a potential and a test operator.
       */
      DevicePotentialIntegral(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getOperand().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      DevicePotentialIntegral(const DevicePotentialIntegral& other)
        : Parent(other),
          DeviceParent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      DevicePotentialIntegral(DevicePotentialIntegral&& other)
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

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @param tau Trial polytope.
       * @param t Test polytope.
       * @returns Never returns.
       * @throws Alert::Exception Always.
       */
      DevicePotentialIntegral& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& tau,
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& t) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @param tr Local trial degree of freedom.
       * @param te Local test degree of freedom.
       * @returns Never returns.
       * @throws Alert::Exception Always.
       */
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the trial region. @return Region::Cells.
      ::Rodin::Geometry::Region getTrialRegion() const override
      {
        return getIntegrand().getLHS().getRegion();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      DeviceLayout getLayout() const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getOperand().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const auto& testAttributes = this->getTestAttributes();
        const auto& trialAttributes = this->getTrialAttributes();
        DeviceLayout layout{&trialfes.getMesh(),
          static_cast<Index>(trialfes.getVectorDimension()),
          static_cast<Index>(testfes.getVectorDimension()),
          {testAttributes.begin(), testAttributes.end()}, this->getTestRegion(),
          trialfes.getDeviceDOFs().getStride(), testfes.getDeviceDOFs().getStride(),
          TrialDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES),
          TestDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
        layout.global = true;
        layout.trialAttributes.assign(trialAttributes.begin(), trialAttributes.end());
        return layout;
      }

      /**
       * @brief Locates the entries of every pair of selected cells.
       * @return Number of entries, the pair of a test cell with n test and a
       * trial cell with m trial local degrees of freedom owning n * m
       * entries.
       */
      size_t prepare(bool) const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getOperand().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const MeshType& mesh = trialfes.getMesh();
        if (&testfes.getMesh() != &mesh)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos integrals require trial and test spaces on the same mesh."
            << Alert::Raise;
        }
        const auto& execution = mesh.getContext().getExecutionSpace();
        const Index dimension = mesh.getDimension();
        const Index count = mesh.getPolytopeCount(dimension);
        if (m_mesh != &mesh || m_cellCount != count)
        {
          upload(mesh, dimension);
          m_mesh = &mesh;
          m_cellCount = count;
        }
        if (!m_testSelection || *m_testSelection != this->getTestAttributes())
        {
          m_testAttributes =
            Internal::uploadAttributes(execution, this->getTestAttributes());
          m_testSelection = this->getTestAttributes();
        }
        if (!m_trialSelection || *m_trialSelection != this->getTrialAttributes())
        {
          m_trialAttributes =
            Internal::uploadAttributes(execution, this->getTrialAttributes());
          m_trialSelection = this->getTrialAttributes();
        }
        const auto view = mesh.getDeviceView();
        select(execution, view, dimension, count, m_testAttributes,
          testfes.getDeviceDOFs(), testfes.getVectorDimension(), m_test);
        select(execution, view, dimension, count, m_trialAttributes,
          trialfes.getDeviceDOFs(), trialfes.getVectorDimension(), m_trial);
        const size_t pairs = size_t(m_test.count) * m_trial.count;
        m_uniform = pairs > 1 && m_test.smallest == m_test.largest &&
          m_trial.smallest == m_trial.largest;
        m_spacing = m_uniform ? pairs : 1;
        m_count = size_t(m_test.total) * m_trial.total;
        return m_count;
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_spacing;
      }

      /**
       * @brief Integrates every pair of selected cells on the device.
       *
       * The entry of local test index te and local trial index tr of a pair
       * has rank te * ntr + tr among the entries of the pair.
       *
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
       */
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getOperand().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        integratePairs(trialfes.getMesh(), integrand.getLHS().getKernel(), m_test,
          m_trial, m_uniform, m_spacing, m_nodes, m_weights, m_centroids,
          m_centroidWeights, trialfes.getDeviceDOFs(), trialfes.getVectorDimension(),
          testfes.getDeviceDOFs(), testfes.getVectorDimension(), values, offset, sign,
          accumulate, rows, cols);
      }

      /// @brief Counts the selected cells and bounds their local counts.
      template <class DOFs>
      struct CellRange
      {
        /// @brief Device mesh.
        Geometry::MeshView view;
        /// @brief Integration region of the cells.
        ::Rodin::Geometry::Region region;
        /// @brief Dimension of the cells.
        Index dimension;
        /// @brief Selected attributes; empty selects every cell.
        AttributesType attributes;
        /// @brief Device degrees of freedom of the side.
        DOFs dofs;
        /// @brief Vector dimension of the side.
        Index vdim;

        /**
         * @brief Counts one cell if it is selected and bounds its local count.
         * @param p Cell index.
         * @param m Number of selected cells.
         * @param lo Smallest local count.
         * @param hi Largest local count.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          Index p, Index& m, Index& lo, Index& hi) const
        {
          const auto cell = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(view, region, cell, attributes))
            return;
          const Index n = dofs.bind(cell).size() * vdim;
          m += 1;
          lo = n < lo ? n : lo;
          hi = n > hi ? n : hi;
        }
      };

      /// @brief Writes the selected cells, in increasing order.
      struct SelectedCells
      {
        /// @brief Device mesh.
        Geometry::MeshView view;
        /// @brief Integration region of the cells.
        ::Rodin::Geometry::Region region;
        /// @brief Dimension of the cells.
        Index dimension;
        /// @brief Selected attributes; empty selects every cell.
        AttributesType attributes;
        /// @brief Selected cells, in increasing order.
        ::Kokkos::View<Index*, MemorySpace> cells;

        /**
         * @brief Writes one cell at its rank among the selected cells if it is
         * selected.
         * @param p Cell index.
         * @param update Number of selected cells before the cell.
         * @param final Whether the scan is in its final pass.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          Index p, Index& update, const bool final) const
        {
          const auto cell = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(view, region, cell, attributes))
            return;
          if (final)
            cells(update) = p;
          update += 1;
        }
      };

      /// @brief Writes the sum of the local counts of the previous selected cells.
      template <class DOFs>
      struct CellOffsets
      {
        /// @brief Device mesh.
        Geometry::MeshView view;
        /// @brief Dimension of the cells.
        Index dimension;
        /// @brief Device degrees of freedom of the side.
        DOFs dofs;
        /// @brief Vector dimension of the side.
        Index vdim;
        /// @brief Selected cells, in increasing order.
        ::Kokkos::View<Index*, MemorySpace> cells;
        /// @brief Sum of the local counts of the previous selected cells.
        ::Kokkos::View<Index*, MemorySpace> offsets;

        /**
         * @brief Writes the offset of one selected cell.
         * @param a Rank of the cell among the selected cells.
         * @param update Sum of the local counts of the previous selected cells.
         * @param final Whether the scan is in its final pass.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          Index a, Index& update, const bool final) const
        {
          if (final)
            offsets(a) = update;
          update += dofs.bind(view.getPolytope(dimension, cells(a))).size() * vdim;
        }
      };

      /**
       * @brief Selects the cells of one side of the pairs.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param dimension Dimension of the cells.
       * @param count Number of cells.
       * @param attributes Selected attributes; empty selects every cell.
       * @param dofs Device degrees of freedom of the side.
       * @param vdim Vector dimension of the side.
       * @param selection Selected cells and their offsets.
       */
      template <class DOFs>
      static void select(const ExecutionSpace& execution, const Geometry::MeshView& view,
        Index dimension, Index count, const AttributesType& attributes, const DOFs& dofs,
        Index vdim, Selection& selection)
      {
        const auto policy = Internal::DeviceQuadrature<MeshType>::policy;
        const ::Rodin::Geometry::Region region = ::Rodin::Geometry::Region::Cells;
        Index selected = 0;
        Index smallest = 0;
        Index largest = 0;
        ::Kokkos::parallel_reduce(
          "Rodin potential cell counts", policy(execution, count),
          CellRange<DOFs>{view, region, dimension, attributes, dofs, vdim},
          ::Kokkos::Sum<Index>(selected), ::Kokkos::Min<Index>(smallest),
          ::Kokkos::Max<Index>(largest));
        if (selection.cells.extent(0) != selected)
        {
          selection.cells = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin potential cells")),
            selected);
          selection.offsets = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin potential offsets")),
            selected);
        }
        const auto cells = selection.cells;
        const auto offsets = selection.offsets;
        ::Kokkos::parallel_scan(
          "Rodin potential cells", policy(execution, count),
          SelectedCells{view, region, dimension, attributes, cells});
        Index total = 0;
        ::Kokkos::parallel_scan(
          "Rodin potential offsets", policy(execution, selected),
          CellOffsets<DOFs>{view, dimension, dofs, vdim, cells, offsets}, total);
        selection.count = selected;
        selection.total = total;
        selection.smallest = selected > 0 ? smallest : 0;
        selection.largest = largest;
      }

      /**
       * @brief Integrates the pairs of selected cells.
       * @param mesh Resident mesh.
       * @param kernel Device kernel @f$ k(x, y) @f$.
       * @param test Selected test cells.
       * @param trial Selected trial cells.
       * @param uniform Whether all the pairs have the same local counts, the
       * entries then being interleaved with the given spacing.
       * @param spacing Spacing of the entries.
       * @param nodes Gauss-Legendre nodes on the unit interval.
       * @param weights Gauss-Legendre weights on the unit interval.
       * @param centroids Reference coordinates of the centroid of each
       * reference geometry, with stride ReferenceStride.
       * @param centroidWeights Centroid rule weight of each reference geometry.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param values Entry values.
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values
       * instead of replacing them.
       * @param rows Row indices, written as the values unless empty.
       * @param cols Column indices, written as the values unless empty.
       */
      static void integratePairs(const MeshType& mesh, const KernelType& kernel,
        const Selection& test, const Selection& trial, bool uniform, size_t spacing,
        const TableType& nodes, const TableType& weights, const TableType& centroids,
        const TableType& centroidWeights, const TrialDOFsType& trialDOFs, Index trialVDim,
        const TestDOFsType& testDOFs, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        const size_t pairs = size_t(test.count) * trial.count;
        if (pairs == 0)
          return;
        const auto& execution = mesh.getContext().getExecutionSpace();
        ::Kokkos::parallel_for(
          "Rodin potential integration",
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<size_t>>(
            execution, 0, pairs),
          PairKernel{kernel, mesh.getDeviceView(), mesh.getDimension(), test.cells,
            test.offsets, trial.cells, trial.offsets, trial.count, trial.total, uniform,
            spacing, nodes, weights, centroids, centroidWeights, nodes.extent(0),
            trialDOFs, trialVDim, testDOFs, testVDim, values, offset, sign, accumulate,
            rows.extent(0) > 0, rows, cols});
      }

      /**
       * @brief Creates a copy of the integral.
       * @returns Pointer to a new instance of the derived integral.
       */
      DevicePotentialIntegral* copy() const noexcept override = 0;

    private:
      /**
       * @brief Uploads the quadrature tables.
       * @param mesh Resident mesh.
       * @param dimension Dimension of the cells.
       */
      void upload(const MeshType& mesh, Index dimension) const
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const auto representatives =
          Internal::DeviceQuadrature<MeshType>::getRepresentatives(mesh, dimension);
        size_t order = 7;
        const auto& triangle = representatives[static_cast<size_t>(
          ::Rodin::Geometry::Polytope::Type::Triangle)];
        if (triangle)
        {
          const auto it = mesh.getPolytope(dimension, *triangle);
          order = this->getOrder(*it).value_or(7);
        }
        std::vector<Real> nodes;
        std::vector<Real> weights;
        QF::GaussLegendre::gl1dUnit(std::max<size_t>(1, (order + 2) / 2), nodes, weights);
        m_nodes =
          Rodin::Kokkos::Internal::upload(execution, nodes, "Rodin potential nodes");
        m_weights =
          Rodin::Kokkos::Internal::upload(execution, weights, "Rodin potential weights");
        std::vector<Real> centroids(GeometryCount * ReferenceStride, Real(0));
        std::vector<Real> centroidWeights(GeometryCount, Real(0));
        for (size_t g = 0; g < GeometryCount; ++g)
        {
          if (!representatives[g])
            continue;
          const QF::Centroid qf(::Rodin::Geometry::Polytope::Types[g]);
          const auto& point = qf.getPoint(0);
          for (Index d = 0; d < static_cast<Index>(point.size()); ++d)
            centroids[g * ReferenceStride + d] = point(d);
          centroidWeights[g] = qf.getWeight(0);
        }
        m_centroids = Rodin::Kokkos::Internal::upload(
          execution, centroids, "Rodin potential centroids");
        m_centroidWeights = Rodin::Kokkos::Internal::upload(
          execution, centroidWeights, "Rodin potential centroid weights");
      }

      std::unique_ptr<IntegrandType> m_integrand;
      mutable const MeshType* m_mesh = nullptr;
      mutable Index m_cellCount = 0;
      mutable Optional<FlatSet<::Rodin::Geometry::Attribute>> m_testSelection;
      mutable Optional<FlatSet<::Rodin::Geometry::Attribute>> m_trialSelection;
      mutable AttributesType m_testAttributes;
      mutable AttributesType m_trialAttributes;
      mutable TableType m_nodes;
      mutable TableType m_weights;
      mutable TableType m_centroids;
      mutable TableType m_centroidWeights;
      mutable Selection m_test;
      mutable Selection m_trial;
      mutable bool m_uniform = false;
      mutable size_t m_spacing = 1;
      mutable size_t m_count = 0;

    public:
      /// @brief Kernel of integratePairs(), one pair of selected cells per thread.
      struct PairKernel
      {
        /// @brief Device kernel @f$ k(x, y) @f$.
        KernelType kernel;
        /// @brief Device mesh.
        Geometry::MeshView view;
        /// @brief Dimension of the cells.
        Index dimension;
        /// @brief Selected test cells.
        ::Kokkos::View<Index*, MemorySpace> testCells;
        /// @brief Sum of the local counts of the previous selected test cells.
        ::Kokkos::View<Index*, MemorySpace> testOffsets;
        /// @brief Selected trial cells.
        ::Kokkos::View<Index*, MemorySpace> trialCells;
        /// @brief Sum of the local counts of the previous selected trial cells.
        ::Kokkos::View<Index*, MemorySpace> trialOffsets;
        /// @brief Number of selected trial cells.
        Index trialCount;
        /// @brief Sum of the local counts of the selected trial cells.
        Index trialTotal;
        /// @brief Whether all the pairs have the same local counts.
        bool uniform;
        /// @brief Spacing of the entries.
        size_t spacing;
        /// @brief Gauss-Legendre nodes on the unit interval.
        TableType nodes;
        /// @brief Gauss-Legendre weights on the unit interval.
        TableType weights;
        /// @brief Reference coordinates of the centroid of each reference geometry.
        TableType centroids;
        /// @brief Centroid rule weight of each reference geometry.
        TableType centroidWeights;
        /// @brief Number of Gauss-Legendre nodes.
        Index n;
        /// @brief Device degrees of freedom of the trial space.
        TrialDOFsType trialDOFs;
        /// @brief Vector dimension of the trial space.
        Index trialVDim;
        /// @brief Device degrees of freedom of the test space.
        TestDOFsType testDOFs;
        /// @brief Vector dimension of the test space.
        Index testVDim;
        /// @brief Entry values.
        ::Kokkos::View<Real*, MemorySpace> values;
        /// @brief Position of the first entry.
        size_t offset;
        /// @brief Factor of the contributions.
        Real sign;
        /// @brief Whether the contributions are added to the values.
        bool accumulate;
        /// @brief Whether the row and column indices are written.
        bool index;
        /// @brief Row indices, test degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> rows;
        /// @brief Column indices, trial degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> cols;

        /**
         * @brief Integrates one pair of selected cells.
         * @param k Rank of the pair, ordered by test cell, then by trial cell.
         */
        KOKKOS_INLINE_FUNCTION void operator()(size_t k) const
        {
          const Index a = k / trialCount;
          const Index b = k % trialCount;
          const Index te = testCells(a);
          const Index tr = trialCells(b);
          const auto testCell = view.getPolytope(dimension, te);
          const auto trialCell = view.getPolytope(dimension, tr);
          const auto testIndices = testDOFs.bind(testCell);
          const auto trialIndices = trialDOFs.bind(trialCell);
          const Index nte = testIndices.size() * testVDim;
          const Index ntr = trialIndices.size() * trialVDim;
          const size_t base = offset +
            (uniform
                ? k
                : size_t(testOffsets(a)) * trialTotal + size_t(nte) * trialOffsets(b));
          // Local index q * vdim + r is component r of the vertex q,
          // numbered globally as vertex + r * stride.
          const auto write = [&](Index l, Index m, Real v) {
            const size_t e = base + size_t(l * ntr + m) * spacing;
            if (accumulate)
              values(e) += sign * v;
            else
              values(e) = sign * v;
            if (index)
            {
              rows(e) =
                testIndices[l / testVDim] + (l % testVDim) * testDOFs.getStride();
              cols(e) =
                trialIndices[m / trialVDim] + (m % trialVDim) * trialDOFs.getStride();
            }
          };
          if (te == tr &&
            testCell.getGeometry() == ::Rodin::Geometry::Polytope::Type::Triangle)
          {
            ::Kokkos::Array<Real, TriangleCapacity * TriangleCapacity> local{};
            Geometry::Vector gradient(2);
            for (Index i3 = 0; i3 < n; ++i3)
            {
              const Real eta3 = nodes(i3);
              const Real w3 = weights(i3);
              for (Index i2 = 0; i2 < n; ++i2)
              {
                const Real eta2 = nodes(i2);
                const Real w2 = weights(i2);
                for (Index i1 = 0; i1 < n; ++i1)
                {
                  const Real eta1 = nodes(i1);
                  const Real w1 = weights(i1);
                  for (Index i0 = 0; i0 < n; ++i0)
                  {
                    const Real xi = nodes(i0);
                    const Real w0 = weights(i0);
                    const Real jacobian = xi * xi * xi * eta1 * eta1 * eta2;
                    const Real factor = w0 * w1 * w2 * w3 * jacobian;
                    // The six sub-domains of the identical-panel case of
                    // Sauter and Schwab, as on the host.
                    const Real ca = xi - xi * eta1;
                    const Real cb = xi - xi * eta1 * eta2 * eta3;
                    const Real cc = xi - xi * eta1 * eta2;
                    const Real cd = xi * eta1 * (1 - eta2);
                    const Real ce = xi * eta1 * (1 - eta2 * eta3);
                    const Real cf = xi * eta1 * (1 - eta2 + eta2 * eta3);
                    const Real cg = xi - xi * eta1 + xi * eta1 * eta2;
                    const ::Kokkos::Array<Real, 24> r{1 - xi, cg, 1 - cb, ca, 1 - cb,
                      ca, 1 - xi, cg, 1 - xi, cf, 1 - cc, cd, 1 - cc, cd, 1 - xi, cf,
                      1 - cb, ce, 1 - xi, cd, 1 - xi, cd, 1 - cb, ce};
                    ::Kokkos::Array<Real, 6> s;
                    ::Kokkos::Array<Real, 18> px;
                    ::Kokkos::Array<Real, 18> pz;
                    ::Kokkos::Array<std::conditional_t<Scalar, Real, Geometry::Matrix>,
                      6>
                      kv;
                    for (Index i = 0; i < 6; ++i)
                    {
                      Geometry::Vector rx(2);
                      rx(0) = r[4 * i];
                      rx(1) = r[4 * i + 1];
                      Geometry::Vector rz(2);
                      rz(0) = r[4 * i + 2];
                      rz(1) = r[4 * i + 3];
                      const auto x = trialCell.getPoint(rx);
                      const auto z = testCell.getPoint(rz);
                      s[i] = x.getDistortion() * z.getDistortion();
                      if constexpr (Scalar)
                      {
                        kv[i] = kernel(x, z);
                      }
                      else
                      {
                        kv[i] = Geometry::Matrix(testVDim, trialVDim);
                        kernel(kv[i], x, z);
                      }
                      for (Index q = 0; q < 3; ++q)
                      {
                        x.evaluateBasis(q, px[3 * i + q], gradient);
                        z.evaluateBasis(q, pz[3 * i + q], gradient);
                      }
                    }
                    for (Index l = 0; l < nte; ++l)
                    {
                      for (Index m = 0; m < ntr; ++m)
                      {
                        for (Index i = 0; i < 6; ++i)
                        {
                          const Real phi = px[3 * i + m / trialVDim];
                          const Real psi = pz[3 * i + l / testVDim];
                          if constexpr (Scalar)
                            local[l * ntr + m] += factor * s[i] * kv[i] * (phi * psi);
                          else
                            local[l * ntr + m] += factor * s[i] *
                              (kv[i](l % testVDim, m % trialVDim) * phi * psi);
                        }
                      }
                    }
                  }
                }
              }
            }
            for (Index l = 0; l < nte; ++l)
            {
              for (Index m = 0; m < ntr; ++m)
                write(l, m, local[l * ntr + m]);
            }
            return;
          }
          const auto centroid = [&](const Geometry::PolytopeView& cell) {
            const Index g = static_cast<Index>(cell.getGeometry());
            Geometry::Vector reference(cell.getDimension());
            for (Index d = 0; d < cell.getDimension(); ++d)
              reference(d) = centroids(g * ReferenceStride + d);
            return cell.getPoint(reference);
          };
          const auto x = centroid(trialCell);
          const auto y = centroid(testCell);
          const Real scale = (x.getDistortion() * y.getDistortion()) *
            (centroidWeights(static_cast<Index>(trialCell.getGeometry())) *
              centroidWeights(static_cast<Index>(testCell.getGeometry())));
          Geometry::Vector gradient(dimension);
          ::Kokkos::Array<Real, RODIN_MAXIMUM_POLYTOPE_VERTICES> px;
          ::Kokkos::Array<Real, RODIN_MAXIMUM_POLYTOPE_VERTICES> py;
          for (Index q = 0; q < trialIndices.size(); ++q)
            x.evaluateBasis(q, px[q], gradient);
          for (Index q = 0; q < testIndices.size(); ++q)
            y.evaluateBasis(q, py[q], gradient);
          if constexpr (Scalar)
          {
            const Real sk = kernel(x, y);
            for (Index l = 0; l < nte; ++l)
            {
              for (Index m = 0; m < ntr; ++m)
                write(l, m, scale * (sk * (px[m] * py[l])));
            }
          }
          else
          {
            Geometry::Matrix mk(testVDim, trialVDim);
            kernel(mk, x, y);
            for (Index l = 0; l < nte; ++l)
            {
              for (Index m = 0; m < ntr; ++m)
              {
                write(l, m,
                  scale *
                    (mk(l % testVDim, m % trialVDim) * px[m / trialVDim] *
                      py[l / testVDim]));
              }
            }
          }
        }
      };
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup IntegralSpecializations
   * @brief Device integration of a potential tested against a shape function
   * on a Kokkos mesh.
   *
   * Represents
   * @f[
   *   \int_{\mathcal{T}_h} (K u)(y) \cdot v(y) \ dy
   * @f]
   * for P1 device trial and test spaces on a resident Kokkos mesh, as
   * described in Rodin::Kokkos::Variational::DevicePotentialIntegral.
   */
  template <class Kernel, class LHSDerived, class TrialFES, class RHSDerived,
    class TestFES>
    requires Rodin::Kokkos::Variational::DeviceSpace<TrialFES> &&
    Rodin::Kokkos::Variational::DeviceSpace<TestFES>
  class Integral<
    Dot<Potential<Kernel, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    final : public Rodin::Kokkos::Variational::DevicePotentialIntegral<Kernel, LHSDerived,
              RHSDerived, TrialFES, TestFES>
  {
    public:
      /// @brief Device integral type.
      using Parent = Rodin::Kokkos::Variational::DevicePotentialIntegral<Kernel,
        LHSDerived, RHSDerived, TrialFES, TestFES>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = typename Parent::LHSType;

      /// @brief Type of the right operand of the dot product.
      using RHSType = typename Parent::RHSType;

      /// @brief Type of the integrand.
      using IntegrandType = typename Parent::IntegrandType;

      /**
       * @brief Constructs the integral of a potential against a test operator.
       * @param[in] lhs Potential @f$ K u @f$.
       * @param[in] rhs Test operator @f$ v @f$.
       */
      Integral(const LHSType& lhs, const RHSType& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a potential and a test operator.
       */
      Integral(const IntegrandType& integrand)
        : Parent(integrand)
      {}

      /// @brief Copy constructor.
      Integral(const Integral&) = default;

      /// @brief Move constructor.
      Integral(Integral&&) = default;

      /// @brief Gets the test region. @return Region::Cells.
      Geometry::Region getTestRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Creates a copy of the integral.
       * @returns Pointer to a new Integral instance.
       */
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }
  };
}

#endif
