/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_BILINEARINTEGRAL_H
#define RODIN_KOKKOS_VARIATIONAL_BILINEARINTEGRAL_H

/**
 * @file
 * @brief Device integration of bilinear finite-element forms.
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
   * @f$ A @f$ and @f$ B @f$ are lowered once to View, then a kernel
   * integrates every polytope of the region with one thread and writes its
   * local matrix as device triplets.
   *
   * A face is integrated with the basis of the face, whose degrees of
   * freedom are the ones of the face, as on the host.
   *
   * The quadrature order is selected on the host as by the QuadratureRule
   * the host would use, once per reference geometry, as described in
   * Quadrature.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   *
   * Rodin::Variational::Integral, BoundaryIntegral, FaceIntegral and
   * InterfaceIntegral specialize this class for device spaces.
   */
  template <class LHSDerived, class RHSDerived, Space TrialFES, Space TestFES>
  class BilinearIntegral
    : public LocalBilinearFormIntegrator<Internal::BackendOf<TrialFES>>
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
      using TrialDOFsType = ::Rodin::Kokkos::View<TrialFESType>;

      /// @brief Device degrees of freedom of the test space.
      using TestDOFsType = ::Rodin::Kokkos::View<TestFESType>;

      /// @brief Type of the left operand of the dot product.
      using LHSType = ::Rodin::Variational::ShapeFunctionBase<LHSDerived, TrialFESType,
        ::Rodin::Variational::ShapeFunctionSpaceType::Trial>;

      /// @brief Type of the right operand of the dot product.
      using RHSType = ::Rodin::Variational::ShapeFunctionBase<RHSDerived, TestFESType,
        ::Rodin::Variational::ShapeFunctionSpaceType::Test>;

      /// @brief Type of the integrand.
      using IntegrandType = ::Rodin::Variational::Dot<LHSType, RHSType>;

      /// @brief Kokkos integrator base.
      using Parent = LocalBilinearFormIntegrator<Internal::BackendOf<TrialFES>>;

      /// @brief Bulk integration interface.
      using ViewParent = BilinearFormIntegrator<Internal::BackendOf<TrialFES>>;

      /// @brief Device triplet type.
      using TripletsType = typename ViewParent::TripletsType;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos execution space.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = Quadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = ::Rodin::Assembly::Polytopes<MeshType>;

      /**
       * @brief Constructs the integral of a dot product.
       * @param[in] integrand Dot product of a trial and a test operator.
       */
      BilinearIntegral(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integral to copy.
       */
      BilinearIntegral(const BilinearIntegral& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integral to move from.
       */
      BilinearIntegral(BilinearIntegral&& other)
        : Parent(std::move(other)),
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
      BilinearIntegral& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integral has no host evaluation.
       * @param tr Ignored local trial index.
       * @param te Ignored local test index.
       * @return Never returns.
       * @throws Alert::Exception Always.
       */
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Layout getLayout() const override
      {
        const auto& integrand = getIntegrand();
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        Layout layout{&trialfes.getMesh(),
          static_cast<Index>(TrialDOFsType(trialfes).getComponentCount()),
          static_cast<Index>(TestDOFsType(testfes).getComponentCount()),
          {attributes.begin(), attributes.end()}, this->getRegion(),
          TrialDOFsType(trialfes).getStride(), TestDOFsType(testfes).getStride(),
          TrialDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES),
          TestDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
        if constexpr (!TrialDOFsType::Componentwise)
          layout.trialSpace = &trialfes;
        if constexpr (!TestDOFsType::Componentwise)
          layout.testSpace = &testfes;
        return layout;
      }

      /**
       * @brief Tells whether the entries vanish off the diagonal of the
       * component blocks.
       * @return Whether the integrand is component diagonal
       * (isComponentDiagonal()).
       */
      bool isBlockDiagonal() const override
      {
        return isComponentDiagonal<View<LHSType>, View<RHSType>>();
      }

      /**
       * @brief Locates the entries of every selected polytope.
       * @param diagonal Whether to list the diagonals of the component
       * blocks only, honored when isBlockDiagonal().
       * @return Number of entries, each selected polytope with m trial and
       * n test scalar degrees of freedom owning m * trialVDim * n * testVDim
       * entries, or m * n * testVDim for the diagonals only.
       */
      size_t prepare(bool diagonal) override
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
        m_compact = diagonal && isBlockDiagonal();
        // The diagonals hold one entry per component of each pair of scalar
        // degrees of freedom.
        return m_polytopes.prepare(mesh, this->getRegion(), this->getAttributes(),
          TrialDOFsType(trialfes),
          m_compact ? 1 : TrialDOFsType(trialfes).getComponentCount(),
          TestDOFsType(testfes), TestDOFsType(testfes).getComponentCount(),
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
       * among the entries of the polytope, or the diagonals of its
       * component blocks as described by writeEntries() when
       * prepared so.
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
        const auto& trialfes = integrand.getLHS().getFiniteElementSpace();
        const auto& testfes = integrand.getRHS().getFiniteElementSpace();
        const MeshType& mesh = trialfes.getMesh();
        const View<LHSType> trial(integrand.getLHS());
        const View<RHSType> test(integrand.getRHS());
        const TrialDOFsType trialDOFs = TrialDOFsType(trialfes);
        const TestDOFsType testDOFs = TestDOFsType(testfes);
        const Index trialVDim = TrialDOFsType(trialfes).getComponentCount();
        const Index testVDim = TestDOFsType(testfes).getComponentCount();
        const auto& quadrature = m_polytopes.getQuadrature();
        const Index vertices = quadrature.getMaximumVertexCount();
        const bool compact = m_compact;
        constexpr bool p1 = std::is_same_v<typename TrialFESType::DOFsType, VertexDOFs> &&
          std::is_same_v<typename TestFESType::DOFsType, VertexDOFs>;
        const bool maximal = trialVDim == TrialDOFsType::ComponentCapacity &&
          testVDim == TestDOFsType::ComponentCapacity;
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
                trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, compact,
                rows, cols);
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
              trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, compact,
              rows, cols);
            return;
          }
        }
        constexpr Index Capacity = IntegrationPointView::Capacity;
        if (vertices <= 4)
        {
          integratePolytopes<4>(mesh, m_polytopes, trial, test, trialDOFs, trialVDim,
            testDOFs, testVDim, values, offset, sign, accumulate, compact, rows, cols);
        }
        else
        {
          integratePolytopes<Capacity>(mesh, m_polytopes, trial, test, trialDOFs,
            trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, compact,
            rows, cols);
        }
      }

      /**
       * @brief Gets the number of entries of the local matrix of a polytope.
       * @tparam Vertices Maximum vertex count of the polytopes.
       * @return Entry count of the local matrix in the kernel.
       */
      template <Index Vertices, class Trial, class Test>
      static constexpr Index getLocalSize()
      {
        return TrialDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TrialFESType, Trial, Test>() *
          TestDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TestFESType, Trial, Test>();
      }

      /**
       * @brief Tells whether the kernels integrate one scalar matrix
       * replicated on the diagonal of the component blocks.
       * @return Whether the test and trial values are sparse by component
       * (Rodin::Kokkos::Variational::IsComponentSparse).
       */
      template <class Trial, class Test>
      static constexpr bool isComponentDiagonal()
      {
        using TestValue = std::decay_t<decltype(std::declval<const Test&>().getBasis(
          std::declval<const IntegrationPointView&>(), Index(0), Index(0)))>;
        using TrialValue = std::decay_t<decltype(std::declval<const Trial&>().getBasis(
          std::declval<const IntegrationPointView&>(), Index(0), Index(0)))>;
        return TrialDOFsType::Componentwise && TestDOFsType::Componentwise &&
          Rodin::Kokkos::Variational::IsComponentSparse<TestValue> &&
          Rodin::Kokkos::Variational::IsComponentSparse<TrialValue>;
      }

      /**
       * @brief Gets the maximum number of components integrated per scalar
       * degree of freedom.
       * @return One for component diagonal integrands, the vector capacity
       * of the range otherwise.
       */
      template <class FES, class Trial, class Test>
      static constexpr Index getComponentCapacity()
      {
        return isComponentDiagonal<Trial, Test>()
          ? 1
          : ::Rodin::Kokkos::View<FES>::ComponentCapacity;
      }

      /**
       * @brief Value of a basis function stored by the kernels.
       *
       * The payload of a component sparse value for component diagonal
       * integrands (Rodin::Kokkos::Variational::getStored()), the basis value
       * otherwise.
       *
       * @tparam Space Test or trial function.
       */
      template <class Space, class Trial, class Test>
      using StoredValue = std::decay_t<decltype(Rodin::Kokkos::Variational::getStored<
        isComponentDiagonal<Trial, Test>()>(std::declval<const Space&>().getBasis(
        std::declval<const IntegrationPointView&>(), Index(0), Index(0))))>;

      /**
       * @brief Gets the team scratch bytes of one quadrature point.
       * @return Bytes of an integration point, its weight and the test and
       * trial basis values at this point.
       */
      template <Index Vertices, class Trial, class Test>
      static constexpr size_t getTeamPointBytes()
      {
        using TestValue = StoredValue<Test, Trial, Test>;
        using TrialValue = StoredValue<Trial, Trial, Test>;
        return sizeof(IntegrationPointView) + sizeof(Real) +
          TestDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TestFESType, Trial, Test>() * sizeof(TestValue) +
          TrialDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TrialFESType, Trial, Test>() * sizeof(TrialValue);
      }

      /**
       * @brief Gets the number of quadrature points of a group of the team
       * kernel that fit in TeamScratch.
       * @return TeamScratch divided by getTeamPointBytes().
       */
      template <Index Vertices, class Trial, class Test>
      static constexpr Index getTeamPoints()
      {
        return TeamScratch / getTeamPointBytes<Vertices, Trial, Test>();
      }

      /**
       * @brief Integrates the polytopes of at most Vertices vertices.
       *
       * When Exact, both spaces are P1, every polytope has Vertices vertices
       * and both vector dimensions are maximal, so the local loops have
       * static bounds. Shape is the geometry tag of the integration points.
       *
       * @param mesh Resident mesh.
       * @param polytopes Selected polytopes, located by prepare().
       * @param trial Device view of the trial operator.
       * @param test Device view of the test operator.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param values Entry values.
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values.
       * @param compact Whether to list the diagonals of the component blocks
       * only.
       * @param rows Row indices, written unless empty.
       * @param cols Column indices, written unless empty.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Trial, class Test>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const Trial& trial, const Test& test, const TrialDOFsType& trialDOFs,
        Index trialVDim, const TestDOFsType& testDOFs, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, bool compact, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        constexpr Index size = getLocalSize<Vertices, Trial, Test>();
        // One thread per polytope bounds its local matrix by LocalCapacity. On
        // the device, teams are faster beyond TeamThreshold when the basis
        // values of a group of points fit TeamScratch.
        if constexpr (size > LocalCapacity ||
          (!OnHost && size > TeamThreshold &&
            getTeamPoints<Vertices, Trial, Test>() >= TeamMinimumPoints))
        {
          integrateTeams<Vertices, Exact, Shape>(mesh, polytopes, trial, test, trialDOFs,
            trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, compact,
            rows, cols);
        }
        else
        {
          integrateThreads<Vertices, Exact, Shape>(mesh, polytopes, trial, test,
            trialDOFs, trialVDim, testDOFs, testVDim, values, offset, sign, accumulate,
            compact, rows, cols);
        }
      }

      /**
       * @brief Integrates the polytopes with one team per polytope.
       *
       * The quadrature points are taken by groups of P points. The threads
       * of the team build the integration points of a group in team scratch
       * memory, one point per thread, then evaluate the test and trial basis
       * functions at these points, one pair of a point and a basis function
       * per thread. Each thread then accumulates in registers the entries of
       * rank @f$ k @f$ such that @f$ k \equiv r \pmod T @f$, for its rank r
       * in a team of T threads. The local matrix is covered by blocks of
       * TeamEntries entries per thread, or TeamWideEntries beyond
       * LocalCapacity, so its size is not bounded. A component diagonal
       * integrand (isComponentDiagonal()) integrates the pairs of scalar
       * degrees of freedom only, with the payloads of the values of component 0
       * (StoredValue), and writes
       * each entry on the diagonal of its component block
       * (writeEntries()).
       *
       * @param mesh Resident mesh.
       * @param polytopes Selected polytopes, located by prepare().
       * @param trial Device view of the trial operator.
       * @param test Device view of the test operator.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param values Entry values.
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values.
       * @param compact Whether to list the diagonals of the component blocks
       * only.
       * @param rows Row indices, written unless empty.
       * @param cols Column indices, written unless empty.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Trial, class Test>
      static void integrateTeams(const MeshType& mesh, const PolytopesType& polytopes,
        const Trial& trial, const Test& test, const TrialDOFsType& trialDOFs,
        Index trialVDim, const TestDOFsType& testDOFs, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, bool compact, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        using Kernel = TeamKernel<Vertices, Exact, Shape, Trial, Test>;
        using Policy = typename Kernel::Policy;
        static_assert(std::is_trivially_copyable_v<typename Kernel::TestValue> &&
          std::is_trivially_copyable_v<typename Kernel::TrialValue> &&
          std::is_trivially_destructible_v<IntegrationPointView>);
        constexpr Index TrialCapacity = Kernel::TrialCapacity;
        constexpr Index TestCapacity = Kernel::TestCapacity;
        constexpr Index E = Kernel::E;
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t scratch = Kernel::P * Kernel::PointBytes + 4 * Kernel::Alignment;
        const int level = scratch <= size_t(Policy::scratch_size_max(0)) ? 0 : 1;
        const Kernel kernel{mesh.getView(), polytopes.getQuadrature(),
          polytopes.getAttributes(), polytopes.getOffsets(), polytopes.getSpacing(),
          polytopes.getRegion(), polytopes.getDimension(), trial, test, trialDOFs,
          trialVDim, testDOFs, testVDim, values, offset, sign, accumulate, compact,
          rows.extent(0) > 0, rows, cols, level};
        Policy policy(execution, polytopes.getPolytopeCount(), ::Kokkos::AUTO);
        policy.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        // About TeamEntries / 2 entries per thread balance the occupancy against
        // the passes over the quadrature.
        const Index wanted =
          (2 * TrialCapacity * TestCapacity + 32 * E - 1) / (32 * E) * 32;
        const int maximum = policy.team_size_max(kernel, ::Kokkos::ParallelForTag());
        int threads = std::max(1, std::min(maximum, int(std::min<Index>(wanted, 1024))));
        Policy launch(execution, polytopes.getPolytopeCount(), threads);
        launch.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        ::Kokkos::parallel_for("Rodin bilinear form team integration", launch, kernel);
      }

      /**
       * @brief Integrates the polytopes with one thread per polytope.
       *
       * Each thread holds the local matrix of its polytope, reduced to the
       * pairs of scalar degrees of freedom for a component diagonal
       * integrand as in integrateTeams().
       *
       * @param mesh Resident mesh.
       * @param polytopes Selected polytopes, located by prepare().
       * @param trial Device view of the trial operator.
       * @param test Device view of the test operator.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param values Entry values.
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions are added to the values.
       * @param compact Whether to list the diagonals of the component blocks
       * only.
       * @param rows Row indices, written unless empty.
       * @param cols Column indices, written unless empty.
       */
      template <Index Vertices, bool Exact = false, class Shape = DynamicGeometry,
        class Trial, class Test>
      static void integrateThreads(const MeshType& mesh, const PolytopesType& polytopes,
        const Trial& trial, const Test& test, const TrialDOFsType& trialDOFs,
        Index trialVDim, const TestDOFsType& testDOFs, Index testVDim,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, bool compact, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const ThreadKernel<Vertices, Exact, Shape, Trial, Test> kernel{mesh.getView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(), trial,
          test, trialDOFs, trialVDim, testDOFs, testVDim, values, offset, sign,
          accumulate, compact, rows.extent(0) > 0, rows, cols};
        ::Kokkos::parallel_for("Rodin bilinear form integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

      /**
       * @brief Creates a polymorphic copy of this integral.
       * @return Pointer to a new copy.
       */
      BilinearIntegral* copy() const noexcept override = 0;

    private:
      /**
       * @brief Writes an integrated entry of a local matrix and its copies.
       *
       * The integrated entry (te, tr) stands for copies by copies entries of
       * the local matrix: test index te * copies + t and trial index
       * tr * copies + s hold the value when t equals s, zero otherwise. Local
       * index a * vdim + r is component r of the scalar degree of freedom a,
       * numbered globally as dof + r * stride, with vdim the product of the
       * integrated components and the copies. The entry of local test index i
       * and local trial index j has rank i * ntr * copies + j. When @p compact,
       * only the entries on the diagonal of the blocks are listed, the copy t
       * having rank (te * ntr + tr) * copies + t.
       *
       * @param values Values of the entries.
       * @param rows Rows of the entries, written when @p index.
       * @param cols Columns of the entries, written when @p index.
       * @param base Position of the first entry of the polytope.
       * @param spacing Distance between consecutive entries of the polytope.
       * @param accumulate Whether to add to the values.
       * @param index Whether to write the rows and columns.
       * @param testIndices Global scalar test degrees of freedom of the polytope.
       * @param testStride Global distance between the components of a test
       * degree of freedom.
       * @param testComponents Integrated test components per scalar degree of
       * freedom.
       * @param trialIndices Global scalar trial degrees of freedom of the
       * polytope.
       * @param trialStride Global distance between the components of a trial
       * degree of freedom.
       * @param trialComponents Integrated trial components per scalar degree
       * of freedom.
       * @param copies Copies of each entry on the diagonal of its block.
       * @param compact Whether to list the diagonal entries of the blocks only.
       * @param ntr Number of integrated trial indices.
       * @param te Integrated test index.
       * @param tr Integrated trial index.
       * @param v Integrated value.
       */
      template <class Values, class Indices, class TestIndices, class TrialIndices>
      KOKKOS_INLINE_FUNCTION static void writeEntries(const Values& values,
        const Indices& rows, const Indices& cols, size_t base, size_t spacing,
        bool accumulate, bool index, const TestIndices& testIndices, Index testStride,
        Index testComponents, const TrialIndices& trialIndices, Index trialStride,
        Index trialComponents, Index copies, bool compact, Index ntr, Index te, Index tr,
        Real v)
      {
        const Index testVDim = testComponents * copies;
        const Index trialVDim = trialComponents * copies;
        for (Index t = 0; t < copies; ++t)
        {
          const Index i = te * copies + t;
          // The compact listing visits the diagonal copy only.
          for (Index s = compact ? t : 0; s < (compact ? t + 1 : copies); ++s)
          {
            const Index j = tr * copies + s;
            const Index rank =
              compact ? (te * ntr + tr) * copies + t : i * ntr * copies + j;
            const size_t e = base + size_t(rank) * spacing;
            // The entries off the diagonal of a block vanish.
            if (t == s)
            {
              if (accumulate)
                values(e) += v;
              else
                values(e) = v;
            }
            else if (!accumulate)
            {
              values(e) = 0;
            }
            if (index)
            {
              rows(e) = testIndices[i / testVDim] + (i % testVDim) * testStride;
              cols(e) = trialIndices[j / trialVDim] + (j % trialVDim) * trialStride;
            }
          }
        }
      }

      /**
       * @brief Maximum number of entries of the local matrix of a kernel.
       *
       * Each thread holds its local matrix, which spills beyond the
       * registers, so larger local matrices take one team per polytope.
       */
      static constexpr Index LocalCapacity = 8192;

      /// @brief Local matrix size above which device kernels use one team per polytope.
      static constexpr Index TeamThreshold = 256;

      /// @brief Number of entries of the local matrix accumulated by each team thread.
      static constexpr Index TeamEntries = 8;

      /// @brief Maximum number of quadrature points of a group of the team kernel.
      static constexpr Index TeamPoints = 32;

      /// @brief Team scratch bytes targeted by the groups of quadrature points.
      static constexpr size_t TeamScratch = 24576;

      /**
       * @brief Minimum number of quadrature points of a group.
       *
       * Smaller groups leave the threads waiting on the integration points,
       * so below LocalCapacity the forms whose groups would be smaller keep
       * one thread per polytope.
       */
      static constexpr Index TeamMinimumPoints = 4;

      /// @brief Number of points of a group whose values exceed TeamScratch, in the slower scratch level.
      static constexpr Index TeamSlowPoints = 8;

      /// @brief Number of entries per team thread beyond LocalCapacity.
      static constexpr Index TeamWideEntries = 32;

      /// @brief Whether the kernels run on the host.
      static constexpr bool OnHost =
        ::Kokkos::SpaceAccessibility<ExecutionSpace, ::Kokkos::HostSpace>::accessible;

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
        if constexpr (kind != Internal::DerivativeKind::None && (p1 || h1))
        {
          if (!order)
          {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto& trialfe =
              integrand.getLHS().getFiniteElementSpace().getFiniteElement(d, i);
            const auto& testfe =
              integrand.getRHS().getFiniteElementSpace().getFiniteElement(d, i);
            const size_t trialDegree =
              trialfe.getDegree().value_or(2 * trialfe.getOrder());
            const size_t testDegree = testfe.getDegree().value_or(2 * testfe.getOrder());
            if (p1 || kind == Internal::DerivativeKind::Gradients)
            {
              order =
                (trialDegree == 0 || testDegree == 0) ? 0 : trialDegree + testDegree - 2;
            }
            else if (kind == Internal::DerivativeKind::TrialDivergence)
              order = trialDegree == 0 ? 0 : trialDegree + testDegree - 1;
            else if (kind == Internal::DerivativeKind::TestDivergence)
              order = testDegree == 0 ? 0 : trialDegree + testDegree - 1;
            else
              order = trialDegree == 0 ? testDegree : trialDegree + testDegree - 1;
          }
        }
        if (!order)
        {
          const size_t d = polytope.getDimension();
          const Index i = polytope.getIndex();
          auto trialOrder = integrand.getLHS().getLeaf().getOrder(polytope);
          auto testOrder = integrand.getRHS().getLeaf().getOrder(polytope);
          if (!trialOrder)
          {
            trialOrder = integrand.getLHS()
                           .getFiniteElementSpace()
                           .getFiniteElement(d, i)
                           .getDegree();
          }
          if (!testOrder)
          {
            testOrder = integrand.getRHS()
                          .getFiniteElementSpace()
                          .getFiniteElement(d, i)
                          .getDegree();
          }
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
      PolytopesType m_polytopes;
      bool m_compact = false;

    public:
      /// @brief Kernel of integrateTeams(), one team per polytope.
      template <Index Vertices, bool Exact, class Shape, class Trial, class Test>
      struct TeamKernel
      {
        /// @brief Stored value of a test basis function at a point (StoredValue).
        using TestValue = StoredValue<Test, Trial, Test>;

        /// @brief Stored value of a trial basis function at a point (StoredValue).
        using TrialValue = StoredValue<Trial, Trial, Test>;

        /// @brief Integration point type of the geometry tag.
        using PointType = Rodin::Kokkos::Variational::KernelPointView<Shape>;

        /// @brief Team policy of the kernel launch.
        using Policy = ::Kokkos::TeamPolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>;

        /// @brief Team handle given to the kernel.
        using Member = typename Policy::member_type;

        /**
         * @brief Whether the integrand is component diagonal, so that one entry
         * is integrated per pair of scalar degrees of freedom.
         */
        static constexpr bool Diagonal = isComponentDiagonal<Trial, Test>();

        /// @brief Maximum number of integrated trial indices of a polytope.
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TrialFESType, Trial, Test>();

        /// @brief Maximum number of integrated test indices of a polytope.
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TestFESType, Trial, Test>();

        /// @brief Team scratch bytes of one quadrature point of a group.
        static constexpr size_t PointBytes = getTeamPointBytes<Vertices, Trial, Test>();
        // Groups that do not fit TeamScratch go to the slower scratch level.
        /// @brief Number of quadrature points of a group that fit TeamScratch.
        static constexpr Index Fit = getTeamPoints<Vertices, Trial, Test>();

        /**
         * @brief Number of quadrature points of a group, at most TeamPoints, or
         * TeamSlowPoints when fewer than TeamMinimumPoints fit.
         */
        static constexpr Index P = Fit < TeamMinimumPoints ? TeamSlowPoints
          : Fit > TeamPoints                               ? TeamPoints
                                                           : Fit;
        // Wider blocks reduce the passes over the quadrature of large matrices.
        /**
         * @brief Number of local matrix entries accumulated by each thread per
         * pass over the quadrature.
         */
        static constexpr Index E =
          TrialCapacity * TestCapacity > LocalCapacity ? TeamWideEntries : TeamEntries;

        /// @brief Alignment of the arrays in team scratch memory.
        static constexpr size_t Alignment = std::max({alignof(TestValue),
          alignof(TrialValue), alignof(IntegrationPointView), alignof(Real)});

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

        /// @brief Device view of the trial operator.
        Trial trial;

        /// @brief Device view of the test operator.
        Test test;

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

        /// @brief Position of the first entry of the integral.
        size_t offset;

        /// @brief Factor of the contributions.
        Real sign;

        /// @brief Whether the contributions are added to the values.
        bool accumulate;

        /// @brief Whether only the diagonals of the component blocks are listed.
        bool compact;

        /// @brief Whether the indices of the entries are written.
        bool index;

        /// @brief Row indices, test degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> rows;

        /// @brief Column indices, trial degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> cols;

        /// @brief Team scratch level holding the groups of quadrature points.
        int level;

        /**
         * @brief Integrates the polytope of a team if it is selected.
         * @param team Team handle, whose league rank is the polytope index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(const Member& team) const
        {
          const Index p = team.league_rank();
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Geometry::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(Exact ? polytope.getVertices().size() == Vertices
                              : polytope.getVertices().size() <= Vertices);
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const Index trialStride =
            Exact ? TrialDOFsType::ComponentCapacity : trialVDim;
          const Index testStride = Exact ? TestDOFsType::ComponentCapacity : testVDim;
          KOKKOS_ASSERT(!Diagonal || trialStride == testStride);
          // Integrated components per scalar degree of freedom, and copies of
          // each integrated entry on the diagonal of the component blocks.
          const Index trialComponents = Diagonal ? 1 : trialStride;
          const Index testComponents = Diagonal ? 1 : testStride;
          const Index copies = Diagonal ? trialStride : 1;
          const Index ntr =
            Exact ? TrialCapacity : trialIndices.size() * trialComponents;
          const Index nte = Exact ? TestCapacity : testIndices.size() * testComponents;
          const Index g = static_cast<Index>(polytope.getGeometry());
          const Index rank = team.team_rank();
          const Index threads = team.team_size();
          auto* const points =
            static_cast<PointType*>(team.team_scratch(level).get_shmem_aligned(
              P * sizeof(PointType), Alignment));
          auto* const weights = static_cast<Real*>(
            team.team_scratch(level).get_shmem_aligned(P * sizeof(Real), Alignment));
          auto* const phi =
            static_cast<TestValue*>(team.team_scratch(level).get_shmem_aligned(
              P * TestCapacity * sizeof(TestValue), Alignment));
          auto* const psi =
            static_cast<TrialValue*>(team.team_scratch(level).get_shmem_aligned(
              P * TrialCapacity * sizeof(TrialValue), Alignment));
          const Index entries = nte * ntr;
          const Index begin = quadrature.begin(g);
          const Index end = quadrature.end(g);
          const size_t base = offset + offsets(p);
          // Groups of balanced sizes, at most P points each.
          const Index groups = (end - begin + P - 1) / P;
          const Index size = groups > 0 ? (end - begin + groups - 1) / groups : 0;
          const auto write = [&](Index te, Index tr, Real v) {
            BilinearIntegral::writeEntries(values, rows, cols, base, spacing,
              accumulate, index, testIndices, testDOFs.getStride(), testComponents,
              trialIndices, trialDOFs.getStride(), trialComponents, copies, compact,
              ntr, te, tr, v);
          };
          for (Index first = 0; first < entries; first += threads * E)
          {
            Real local[E] = {};
            for (Index q0 = begin; q0 < end; q0 += size)
            {
              const Index count = end - q0 < size ? end - q0 : size;
              // The points and the values of the previous group are no longer read.
              team.team_barrier();
              for (Index c = rank; c < count; c += threads)
              {
                new (points + c) PointType(
                  polytope, quadrature.getReference(polytope, q0 + c), Shape{});
                weights[c] = quadrature.getWeight(q0 + c) * points[c].getDistortion();
              }
              team.team_barrier();
              // Separate loops keep the lanes of a warp on the same basis.
              for (Index i = rank; i < count * nte; i += threads)
              {
                const Index c = i / nte;
                const Index te = i - c * nte;
                const auto basis =
                  Rodin::Kokkos::Variational::bindBasis(test, points[c]);
                phi[c * TestCapacity + te] =
                  Rodin::Kokkos::Variational::getStored<Diagonal>(
                    basis.getBasis(te / testComponents, te % testComponents));
              }
              for (Index i = rank; i < count * ntr; i += threads)
              {
                const Index c = i / ntr;
                const Index tr = i - c * ntr;
                const auto basis =
                  Rodin::Kokkos::Variational::bindBasis(trial, points[c]);
                psi[c * TrialCapacity + tr] =
                  Rodin::Kokkos::Variational::getStored<Diagonal>(
                    basis.getBasis(tr / trialComponents, tr % trialComponents));
              }
              team.team_barrier();
              RODIN_KOKKOS_UNROLL
              for (Index j = 0; j < E; ++j)
              {
                const Index k = first + j * threads + rank;
                if (k < entries)
                {
                  const Index te = k / ntr;
                  const Index tr = k - te * ntr;
                  for (Index c = 0; c < count; ++c)
                  {
                    const auto& u = psi[c * TrialCapacity + tr];
                    const auto& v = phi[c * TestCapacity + te];
                    if (Rodin::Kokkos::Variational::overlaps(u, v))
                      local[j] += weights[c] * Rodin::Kokkos::Variational::dot(u, v);
                  }
                }
              }
            }
            RODIN_KOKKOS_UNROLL
            for (Index j = 0; j < E; ++j)
            {
              const Index k = first + j * threads + rank;
              if (k < entries)
              {
                const Index te = k / ntr;
                const Index tr = k - te * ntr;
                write(te, tr, sign * local[j]);
              }
            }
          }
        }
      };

      /// @brief Kernel of integrateThreads(), one thread per polytope.
      template <Index Vertices, bool Exact, class Shape, class Trial, class Test>
      struct ThreadKernel
      {
        /// @brief Stored value of a test basis function at a point (StoredValue).
        using TestValue = StoredValue<Test, Trial, Test>;

        /**
         * @brief Whether the integrand is component diagonal, so that one entry
         * is integrated per pair of scalar degrees of freedom.
         */
        static constexpr bool Diagonal = isComponentDiagonal<Trial, Test>();

        /// @brief Maximum number of integrated trial indices of a polytope.
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TrialFESType, Trial, Test>();

        /// @brief Maximum number of integrated test indices of a polytope.
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices) *
          getComponentCapacity<TestFESType, Trial, Test>();

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

        /// @brief Device view of the trial operator.
        Trial trial;

        /// @brief Device view of the test operator.
        Test test;

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

        /// @brief Position of the first entry of the integral.
        size_t offset;

        /// @brief Factor of the contributions.
        Real sign;

        /// @brief Whether the contributions are added to the values.
        bool accumulate;

        /// @brief Whether only the diagonals of the component blocks are listed.
        bool compact;

        /// @brief Whether the indices of the entries are written.
        bool index;

        /// @brief Row indices, test degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> rows;

        /// @brief Column indices, trial degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> cols;

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
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const Index trialStride =
            Exact ? TrialDOFsType::ComponentCapacity : trialVDim;
          const Index testStride = Exact ? TestDOFsType::ComponentCapacity : testVDim;
          KOKKOS_ASSERT(!Diagonal || trialStride == testStride);
          // Integrated components per scalar degree of freedom, and copies of
          // each integrated entry on the diagonal of the component blocks.
          const Index trialComponents = Diagonal ? 1 : trialStride;
          const Index testComponents = Diagonal ? 1 : testStride;
          const Index copies = Diagonal ? trialStride : 1;
          const Index ntr =
            Exact ? TrialCapacity : trialIndices.size() * trialComponents;
          const Index nte = Exact ? TestCapacity : testIndices.size() * testComponents;
          const Index g = static_cast<Index>(polytope.getGeometry());
          ::Kokkos::Array<Real, TrialCapacity * TestCapacity> local{};
          ::Kokkos::Array<TestValue, TestCapacity> phi;
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
              const auto testBasis = Rodin::Kokkos::Variational::bindBasis(test, ip);
              const auto trialBasis = Rodin::Kokkos::Variational::bindBasis(trial, ip);
              RODIN_KOKKOS_UNROLL
              for (Index te = 0; te < nte; ++te)
              {
                phi[te] = Rodin::Kokkos::Variational::getStored<Diagonal>(
                  testBasis.getBasis(te / testComponents, te % testComponents));
              }
              RODIN_KOKKOS_UNROLL
              for (Index tr = 0; tr < ntr; ++tr)
              {
                const auto psi = Rodin::Kokkos::Variational::getStored<Diagonal>(
                  trialBasis.getBasis(tr / trialComponents, tr % trialComponents));
                RODIN_KOKKOS_UNROLL
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
          const size_t base = offset + offsets(p);
          RODIN_KOKKOS_UNROLL
          for (Index te = 0; te < nte; ++te)
          {
            RODIN_KOKKOS_UNROLL
            for (Index tr = 0; tr < ntr; ++tr)
            {
              BilinearIntegral::writeEntries(values, rows, cols, base, spacing,
                accumulate, index, testIndices, testDOFs.getStride(), testComponents,
                trialIndices, trialDOFs.getStride(), trialComponents, copies, compact,
                ntr, te, tr, sign * local[te * ntr + tr]);
            }
          }
        }
      };
  };
}

#endif
