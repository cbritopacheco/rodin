/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_INTERNALVIRTUALWORK_H
#define RODIN_KOKKOS_VARIATIONAL_INTERNALVIRTUALWORK_H

/**
 * @file InternalVirtualWork.h
 * @brief Device internal virtual work of hyperelastic laws.
 *
 * Specializes Solid::InternalVirtualWorkTangent and
 * Solid::InternalVirtualWorkResidual on vector device spaces for the
 * Neo-Hookean law, with closed form entries, and for the laws of
 * Kokkos::Variational::DeviceLaw, through their elasticity tensors, so that
 * Solid::InternalVirtualWork assembles Newton problems of resident Kokkos
 * meshes on the device.
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>

#include "Rodin/Solid/Constitutive/NeoHookean.h"
#include "Rodin/Solid/Integrators/InternalVirtualWork.h"

#include "HyperElasticLaw.h"
#include "Integral.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Neo-Hookean kinematics at a displacement gradient.
   *
   * Holds @f$ F = I + H @f$, @f$ F^{-1} @f$ and @f$ \ln J @f$, as
   * Solid::KinematicState does on the host.
   */
  struct NeoHookeanPoint
  {
    /**
     * @param h Displacement gradient @f$ H @f$, a square matrix.
     */
    KOKKOS_INLINE_FUNCTION
    explicit NeoHookeanPoint(const Geometry::Matrix& h)
      : F(h),
        Finv(h.rows(), h.cols())
    {
      KOKKOS_ASSERT(h.rows() == h.cols());
      for (Index i = 0; i < F.rows(); ++i)
        F(i, i) += 1;
      const Real J = F.determinant();
      KOKKOS_ASSERT(J > 0);
      Finv = F.inverse();
      logJ = ::Kokkos::log(J);
    }

    /**
     * @param g Device vector of the space dimension.
     * @return @f$ F^{-T} g @f$.
     */
    KOKKOS_INLINE_FUNCTION
    Geometry::Vector pullBack(const Geometry::Vector& g) const
    {
      Geometry::Vector result(g.size());
      for (Index s = 0; s < g.size(); ++s)
      {
        Real v = 0;
        for (Index k = 0; k < g.size(); ++k)
          v += Finv(k, s) * g(k);
        result(s) = v;
      }
      return result;
    }

    /// @brief Deformation gradient @f$ F @f$.
    Geometry::Matrix F;

    /// @brief Inverse @f$ F^{-1} @f$.
    Geometry::Matrix Finv;

    /// @brief Logarithm of @f$ J = \det F @f$.
    Real logJ = 0;
  };

  /**
   * @brief Selects the kernel size from the vertex count of the polytopes.
   * @param vertices Maximum vertex count of the integrated polytopes.
   * @param f Callable taking a std::integral_constant of the vertex
   * capacity, returning false when the capacity is unsupported.
   * @return Whether a kernel ran.
   */
  template <class F>
  bool visitVertexCapacity(Index vertices, const F& f)
  {
    if (vertices <= 4)
      return f(std::integral_constant<Index, 4>{});
    return f(std::integral_constant<Index, IntegrationPointView::Capacity>{});
  }
}

namespace Rodin::Solid
{
  /**
   * @brief Device tangent of the internal virtual work of the Neo-Hookean
   * law.
   *
   * Represents
   * @f[
   *   \int_\Omega \frac{\partial P}{\partial F}(F(u)) : \nabla \delta u :
   *   \nabla v \ dx
   * @f]
   * for vector device spaces, with the Neo-Hookean material tangent of
   * NeoHookean::getMaterialTangent(). For the vector Lagrange bases
   * @f$ \phi_a e_s @f$ and @f$ \phi_b e_c @f$, with
   * @f$ \alpha_a = F^{-T} \nabla \phi_a @f$, the entry is
   * @f[
   *   \mu \delta_{cs} \nabla \phi_a \cdot \nabla \phi_b
   *   + \lambda (\alpha_a)_s (\alpha_b)_c
   *   + (\mu - \lambda \ln J) (\alpha_a)_c (\alpha_b)_s,
   * @f]
   * the contraction of the host dense tangent.
   *
   * The quadrature order is that of the host,
   * @f$ 2 \max(k_{\delta u}, k_v, k_u) @f$ unless set by
   * setQuadratureOrder(). The displacement is copied to the device at
   * every assembly. Constitutive inputs are not available on the device.
   */
  template <class Solution, Kokkos::Variational::VectorDeviceSpace TrialFES,
    Kokkos::Variational::VectorDeviceSpace TestFES,
    Kokkos::Variational::VectorDeviceSpace StateFES, class Data>
  class InternalVirtualWorkTangent<NeoHookean,
    Variational::TrialFunction<Solution, TrialFES>, Variational::TestFunction<TestFES>,
    Variational::GridFunction<StateFES, Data>>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>,
            public Kokkos::Variational::DeviceBilinearFormIntegrator<
              Kokkos::Variational::Internal::BackendOf<TrialFES>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Constitutive law type.
      using LawType = NeoHookean;

      /// @brief Trial function type.
      using TrialType = Variational::TrialFunction<Solution, TrialFES>;

      /// @brief Test function type.
      using TestType = Variational::TestFunction<TestFES>;

      /// @brief Current displacement state type.
      using StateType = Variational::GridFunction<StateFES, Data>;

      /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Displacement state finite element space type.
      using StateFESType = StateFES;

      /// @brief Host integrator interface.
      using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = Kokkos::Variational::DeviceBilinearFormIntegrator<
        Kokkos::Variational::Internal::BackendOf<TrialFES>>;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device degrees of freedom of the trial space.
      using TrialDOFsType = typename TrialFESType::DOFsType;

      /// @brief Device degrees of freedom of the test space.
      using TestDOFsType = typename TestFESType::DOFsType;

      /// @brief Device quadrature tables.
      using QuadratureType = Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Kokkos::Variational::Internal::DevicePolytopes<MeshType>;

      /// @brief Device coefficients of the displacement.
      using StateViewType = Kokkos::Variational::GridFunctionView<StateFESType>;

      /**
       * @brief Constructs the tangent at a displacement.
       * @param law Neo-Hookean law.
       * @param u Trial displacement increment.
       * @param v Test displacement.
       * @param displacement Current displacement.
       */
      InternalVirtualWorkTangent(const LawType& law, const TrialType& u,
        const TestType& v, const StateType& displacement)
        : Parent(u, v),
          m_law(law),
          m_trial(u),
          m_test(v),
          m_displacement(displacement),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      /// @brief Copy constructor. @param other Tangent to copy.
      InternalVirtualWorkTangent(const InternalVirtualWorkTangent& other)
        : Parent(other),
          DeviceParent(other),
          m_law(other.m_law),
          m_trial(other.m_trial),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_quadOrder(other.m_quadOrder)
      {}

      /// @brief Sets the displacement. @param gf Current displacement. @return Reference to this.
      InternalVirtualWorkTangent& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        return *this;
      }

      /// @brief Sets the quadrature order. @param order Order. @return Reference to this.
      InternalVirtualWorkTangent& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /// @brief Unavailable on the device. @throws Alert::Exception Always.
      InternalVirtualWorkTangent& setInput([[maybe_unused]] InputFunction input)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Constitutive inputs are not available on Kokkos devices." << Alert::Raise;
        return *this;
      }

      /// @brief Gets the law. @return Neo-Hookean law.
      const LawType& getLaw() const
      {
        return m_law;
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      InternalVirtualWorkTangent& setPolytope(
        [[maybe_unused]] const Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&trialfes.getMesh(), static_cast<Index>(trialfes.getVectorDimension()),
          static_cast<Index>(testfes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(),
          trialfes.getDeviceDOFs().getStride(), testfes.getDeviceDOFs().getStride(),
          TrialDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES),
          TestDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /**
       * @brief Locates the entries of every cell.
       * @return Number of entries, laid out as by the device integrals.
       */
      size_t prepare(bool) const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        return m_polytopes.prepare(trialfes.getMesh(), this->getRegion(),
          this->getAttributes(), trialfes.getDeviceDOFs(), trialfes.getVectorDimension(),
          testfes.getDeviceDOFs(), testfes.getVectorDimension(),
          [&](const Geometry::Polytope& polytope) { return getOrder(polytope); });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /// @brief Integrates every cell on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& displacement = m_displacement.get();
        const StateViewType state(
          displacement.getFiniteElementSpace(), displacement.getData());
        const Index vertices = m_polytopes.getQuadrature().getMaximumVertexCount();
        Kokkos::Variational::visitVertexCapacity(vertices, [&](auto n) {
          constexpr Index Vertices = decltype(n)::value;
          // Larger local matrices take one team per cell.
          if constexpr (getLocalSize<Vertices>() <= LocalCapacity)
            integratePolytopes<Vertices>(trialfes.getMesh(), m_polytopes,
              trialfes.getDeviceDOFs(), testfes.getDeviceDOFs(),
              testfes.getVectorDimension(), state, m_law.getLameFirstParameter(),
              m_law.getShearModulus(), values, offset, sign, accumulate, rows, cols);
          else
            integrateTeams<Vertices>(trialfes.getMesh(), m_polytopes,
              trialfes.getDeviceDOFs(), testfes.getDeviceDOFs(),
              testfes.getVectorDimension(), state, m_law.getLameFirstParameter(),
              m_law.getShearModulus(), values, offset, sign, accumulate, rows, cols);
          return true;
        });
      }

      /// @brief Gets the integration region. @return Cells.
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      InternalVirtualWorkTangent* copy() const noexcept override
      {
        return new InternalVirtualWorkTangent(*this);
      }

      /**
       * @brief Gets the number of entries of the local matrix of a cell.
       * @tparam Vertices Maximum vertex count of the cells.
       * @return Entry count of the local matrix in the kernel.
       */
      template <Index Vertices>
      static constexpr Index getLocalSize()
      {
        return TrialDOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION *
          TestDOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION;
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * The entry of local test index te and local trial index tr has rank
       * te * ntr + tr among the entries of the cell, the local index
       * a * vdim + r being component r of the scalar degree of freedom a.
       */
      template <Index Vertices>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const TrialDOFsType& trialDOFs, const TestDOFsType& testDOFs, Index vdim,
        const StateViewType& state, Real lambda, Real mu,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const PolytopeKernel<Vertices> kernel{mesh.getDeviceView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(),
          trialDOFs, testDOFs, vdim, state, lambda, mu, values, offset, sign, accumulate,
          rows.extent(0) > 0, rows, cols};
        ::Kokkos::parallel_for(
          "Rodin hyperelastic tangent integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices with one
       * team per cell.
       *
       * The quadrature points are taken by groups of P points. The threads
       * of the team build the integration points of a group in team scratch
       * memory, one point per thread, then the gradients
       * @f$ \nabla \phi_a @f$, one pair of a point and a scalar basis
       * function per thread, which the test and state spaces share when
       * they have the type of the trial space. They then sum the
       * displacement gradients, one pair of a point and a basis function
       * of the state per thread, build the kinematics, one point per
       * thread, and the pullbacks @f$ \alpha_a @f$, one pair of a point
       * and a scalar basis function per thread. Each thread then takes a
       * slot, a trial function and a block of TeamPairs test functions,
       * and accumulates in registers the vdim by vdim entries of each of
       * its pairs over the group, reusing the products of the trial
       * function across the block. The entries are laid out as by
       * integratePolytopes().
       */
      template <Index Vertices>
      static void integrateTeams(const MeshType& mesh, const PolytopesType& polytopes,
        const TrialDOFsType& trialDOFs, const TestDOFsType& testDOFs, Index vdim,
        const StateViewType& state, Real lambda, Real mu,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        using Kernel = TeamKernel<Vertices>;
        using Policy = typename Kernel::Policy;
        static_assert(std::is_trivially_destructible_v<typename Kernel::Point> &&
          std::is_trivially_destructible_v<typename Kernel::Kinematics> &&
          std::is_trivially_copyable_v<typename Kernel::Vector>);
        constexpr Index TrialCapacity = Kernel::TrialCapacity;
        constexpr Index TestCapacity = Kernel::TestCapacity;
        constexpr Index E = Kernel::E;
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t scratch = Kernel::P * Kernel::PointBytes + 8 * Kernel::Alignment;
        const int level = scratch <= size_t(Policy::scratch_size_max(0)) ? 0 : 1;
        const Kernel kernel{mesh.getDeviceView(), polytopes.getQuadrature(),
          polytopes.getAttributes(), polytopes.getOffsets(), polytopes.getSpacing(),
          polytopes.getRegion(), polytopes.getDimension(), trialDOFs, testDOFs, vdim,
          state, lambda, mu, values, offset, sign, accumulate, rows.extent(0) > 0, rows,
          cols, level};
        Policy policy(execution, polytopes.getPolytopeCount(), ::Kokkos::AUTO);
        policy.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        // One slot per thread, the fewest passes over the quadrature.
        constexpr Index Slots = TrialCapacity * ((TestCapacity + E - 1) / E);
        const Index wanted = (Slots + 31) / 32 * 32;
        const int maximum = policy.team_size_max(kernel, ::Kokkos::ParallelForTag());
        const int threads =
          std::max(1, std::min(maximum, int(std::min<Index>(wanted, TeamThreads))));
        Policy launch(execution, polytopes.getPolytopeCount(), threads);
        launch.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        ::Kokkos::parallel_for(
          "Rodin hyperelastic tangent team integration", launch, kernel);
      }

    private:
      /// @brief Maximum number of entries of the local matrix of a thread,
      /// so larger local matrices take one team per cell.
      static constexpr Index LocalCapacity = 1024;

      /// @brief Largest group of quadrature points of a team.
      static constexpr Index TeamPoints = 32;

      /// @brief Bytes of team scratch memory of a group in the fast level.
      static constexpr size_t TeamScratch = 47104;

      /// @brief Fewest points of a group in the fast scratch level.
      static constexpr Index TeamMinimumPoints = 4;

      /// @brief Number of points of a group whose values exceed TeamScratch,
      /// in the slower scratch level.
      static constexpr Index TeamSlowPoints = 8;

      /// @brief Largest team, whose register bound the kernel is compiled for.
      static constexpr unsigned int TeamThreads = 512;

      /// @brief Test functions of the block of a thread, whose entries
      /// with one trial function it accumulates per pass over the quadrature.
      static constexpr Index TeamPairs = 4;

      /// @brief Selects the quadrature order of a cell, as on the host.
      size_t getOrder(const Geometry::Polytope& polytope) const
      {
        if (m_quadOrder > 0)
          return m_quadOrder;
        const auto trialOrder = m_trial.get().getOrder(polytope);
        const auto testOrder = m_test.get().getOrder(polytope);
        const auto stateOrder = m_displacement.get().getOrder(polytope);
        if (!trialOrder || !testOrder || !stateOrder)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic constitutive quadrature order is established. "
            << "Specify it with setQuadratureOrder()." << Alert::Raise;
        }
        return 2 * std::max({*trialOrder, *testOrder, *stateOrder});
      }

      void checkCompatibility(const StateType& displacement) const
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& statefes = displacement.getFiniteElementSpace();
        const size_t d = trialfes.getMesh().getSpaceDimension();
        if (&trialfes.getMesh() != &testfes.getMesh() ||
          &statefes.getMesh() != &testfes.getMesh() ||
          trialfes.getVectorDimension() != d || testfes.getVectorDimension() != d ||
          statefes.getVectorDimension() != d || trialfes.getMesh().getDimension() != d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos hyperelastic tangents require displacement spaces of the "
            << "space dimension on one mesh of full dimension." << Alert::Raise;
        }
      }

      LawType m_law;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      size_t m_quadOrder;
      mutable PolytopesType m_polytopes;

    public:
      /**
       * @brief Kernel of integratePolytopes(), the local matrix of one cell
       * per thread.
       */
      template <Index Vertices>
      struct PolytopeKernel
      {
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices);
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices);

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        TrialDOFsType trialDOFs;
        TestDOFsType testDOFs;
        Index vdim;
        StateViewType state;
        Real lambda;
        Real mu;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> rows;
        ::Kokkos::View<Index*, MemorySpace> cols;

        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const Index na = trialIndices.size();
          const Index nb = testIndices.size();
          const Index ntr = na * vdim;
          const Index nte = nb * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          ::Kokkos::Array<Real,
            TrialCapacity * RODIN_MAXIMAL_SPACE_DIMENSION * TestCapacity *
              RODIN_MAXIMAL_SPACE_DIMENSION>
            local{};
          ::Kokkos::Array<Rodin::Kokkos::Geometry::Vector, TrialCapacity> ga;
          ::Kokkos::Array<Rodin::Kokkos::Geometry::Vector, TrialCapacity> aa;
          ::Kokkos::Array<Rodin::Kokkos::Geometry::Vector, TestCapacity> gb;
          ::Kokkos::Array<Rodin::Kokkos::Geometry::Vector, TestCapacity> ab;
          for (Index q = quadrature.begin(g); q < quadrature.end(g); ++q)
          {
            const Rodin::Kokkos::Variational::CellIntegrationPointView ip(polytope,
              quadrature.getReference(polytope, q),
              Rodin::Kokkos::Variational::DynamicGeometry{});
            const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
            const Rodin::Kokkos::Variational::NeoHookeanPoint point(
              state.interpolateJacobian(ip));
            for (Index a = 0; a < na; ++a)
            {
              ga[a] = trialDOFs.getGradient(ip, a);
              aa[a] = point.pullBack(ga[a]);
            }
            for (Index b = 0; b < nb; ++b)
            {
              gb[b] = testDOFs.getGradient(ip, b);
              ab[b] = point.pullBack(gb[b]);
            }
            const Real c1 = wdet * mu;
            const Real c2 = wdet * lambda;
            const Real c3 = wdet * (mu - lambda * point.logJ);
            for (Index b = 0; b < nb; ++b)
            {
              for (Index a = 0; a < na; ++a)
              {
                const Real shear = c1 * Rodin::Kokkos::Variational::dot(ga[a], gb[b]);
                for (Index c = 0; c < vdim; ++c)
                {
                  for (Index s = 0; s < vdim; ++s)
                  {
                    Real v = c2 * aa[a](s) * ab[b](c) + c3 * aa[a](c) * ab[b](s);
                    if (c == s)
                      v += shear;
                    local[(b * vdim + c) * ntr + a * vdim + s] += v;
                  }
                }
              }
            }
          }
          const size_t base = offset + offsets(p);
          for (Index te = 0; te < nte; ++te)
          {
            const Index row =
              testIndices[te / vdim] + (te % vdim) * testDOFs.getStride();
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
                cols(e) = trialIndices[tr / vdim] + (tr % vdim) * trialDOFs.getStride();
              }
            }
          }
        }
      };

      /**
       * @brief Kernel of integrateTeams(), the local matrix of one cell per
       * team.
       */
      template <Index Vertices>
      struct TeamKernel
      {
        using Point = Kokkos::Variational::CellIntegrationPointView;
        using Kinematics = Kokkos::Variational::NeoHookeanPoint;
        using Vector = Kokkos::Geometry::Vector;
        using ExecutionSpace = typename MeshType::ExecutionSpace;
        // The register bound keeps enough threads for one pass over the
        // quadrature, at the price of a few spills.
        using Policy = ::Kokkos::TeamPolicy<ExecutionSpace, ::Kokkos::IndexType<Index>,
          ::Kokkos::LaunchBounds<TeamThreads, 1>>;
        using Member = typename Policy::member_type;
        using StateDOFsType =
          std::decay_t<decltype(std::declval<const StateViewType&>().getDOFs())>;
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices);
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices);
        static constexpr Index D = RODIN_MAXIMAL_SPACE_DIMENSION;
        static constexpr Index Square = D * D;
        // Spaces of one type share their basis gradients.
        static constexpr bool SharedTest = std::is_same_v<TrialDOFsType, TestDOFsType>;
        static constexpr bool SharedState = std::is_same_v<TrialDOFsType, StateDOFsType>;
        static constexpr Index TestStorage = SharedTest ? 0 : TestCapacity;
        static constexpr size_t PointBytes = sizeof(Point) + sizeof(Kinematics) +
          (3 + Square) * sizeof(Real) +
          2 * (TrialCapacity + TestStorage) * sizeof(Vector);
        // Groups that do not fit TeamScratch go to the slower scratch level.
        static constexpr Index Fit = TeamScratch / PointBytes;
        static constexpr Index P = Fit < TeamMinimumPoints ? TeamSlowPoints
          : Fit > TeamPoints                               ? TeamPoints
                                                           : Fit;
        static constexpr Index E = TeamPairs;
        static constexpr size_t Alignment =
          std::max({alignof(Point), alignof(Kinematics), alignof(Vector), alignof(Real)});

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        TrialDOFsType trialDOFs;
        TestDOFsType testDOFs;
        Index vdim;
        StateViewType state;
        Real lambda;
        Real mu;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> rows;
        ::Kokkos::View<Index*, MemorySpace> cols;
        int level;

        KOKKOS_INLINE_FUNCTION void operator()(const Member& team) const
        {
          const Index p = team.league_rank();
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const auto& stateDOFs = state.getDOFs();
          const auto& data = state.getData();
          const auto stateIndices = stateDOFs.bind(polytope);
          const Index ns = stateIndices.size();
          const Index na = trialIndices.size();
          const Index nb = testIndices.size();
          const Index ntr = na * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          const Index rank = team.team_rank();
          const Index threads = team.team_size();
          auto& memory = team.team_scratch(level);
          auto* const points =
            static_cast<Point*>(memory.get_shmem_aligned(P * sizeof(Point), Alignment));
          auto* const kinematics = static_cast<Kinematics*>(
            memory.get_shmem_aligned(P * sizeof(Kinematics), Alignment));
          auto* const coefficients =
            static_cast<Real*>(memory.get_shmem_aligned(3 * P * sizeof(Real), Alignment));
          auto* const jacobians = static_cast<Real*>(
            memory.get_shmem_aligned(Square * P * sizeof(Real), Alignment));
          auto* const ga = static_cast<Vector*>(
            memory.get_shmem_aligned(P * TrialCapacity * sizeof(Vector), Alignment));
          auto* const aa = static_cast<Vector*>(
            memory.get_shmem_aligned(P * TrialCapacity * sizeof(Vector), Alignment));
          Vector* gb = ga;
          Vector* ab = aa;
          if (!SharedTest)
          {
            gb = static_cast<Vector*>(
              memory.get_shmem_aligned(P * TestCapacity * sizeof(Vector), Alignment));
            ab = static_cast<Vector*>(
              memory.get_shmem_aligned(P * TestCapacity * sizeof(Vector), Alignment));
          }
          // A slot is a trial function and a block of E test functions.
          const Index slots = na * ((nb + E - 1) / E);
          const Index begin = quadrature.begin(g);
          const Index end = quadrature.end(g);
          const size_t base = offset + offsets(p);
          // Groups of balanced sizes, at most P points each.
          const Index groups = (end - begin + P - 1) / P;
          const Index size = groups > 0 ? (end - begin + groups - 1) / groups : 0;
          for (Index first = 0; first < slots; first += threads)
          {
            const Index slot = first + rank;
            const Index a = slot % na;
            const Index b0 = (slot / na) * E;
            Real local[E][D][D] = {};
            for (Index q0 = begin; q0 < end; q0 += size)
            {
              const Index count = end - q0 < size ? end - q0 : size;
              // The values of the previous group are no longer read.
              team.team_barrier();
              for (Index c = rank; c < count; c += threads)
              {
                new (points + c)
                  Point(polytope, quadrature.getReference(polytope, q0 + c),
                    Rodin::Kokkos::Variational::DynamicGeometry{});
                for (Index k = 0; k < Square; ++k)
                  jacobians[c * Square + k] = 0;
              }
              team.team_barrier();
              for (Index i = rank; i < count * na; i += threads)
              {
                const Index c = i / na;
                const Index a = i - c * na;
                ga[c * TrialCapacity + a] = trialDOFs.getGradient(points[c], a);
              }
              if (!SharedTest)
              {
                for (Index i = rank; i < count * nb; i += threads)
                {
                  const Index c = i / nb;
                  const Index b = i - c * nb;
                  gb[c * TestCapacity + b] = testDOFs.getGradient(points[c], b);
                }
              }
              // The displacement gradient sums one term per pair of a point
              // and a basis function of the state, after the trial gradients
              // when the state shares them.
              const auto sum = [&](Index i) {
                const Index c = i / ns;
                const Index q = i - c * ns;
                Vector gradient;
                if (SharedState)
                  gradient = ga[c * TrialCapacity + q];
                else
                  gradient = stateDOFs.getGradient(points[c], q);
                for (Index r = 0; r < vdim; ++r)
                {
                  const Real value = data(stateIndices[q] + r * stateDOFs.getStride());
                  for (Index k = 0; k < vdim; ++k)
                    ::Kokkos::atomic_add(
                      &jacobians[c * Square + r * vdim + k], gradient(k) * value);
                }
              };
              if (!SharedState)
              {
                for (Index i = rank; i < count * ns; i += threads)
                  sum(i);
              }
              team.team_barrier();
              if (SharedState)
              {
                for (Index i = rank; i < count * ns; i += threads)
                  sum(i);
                team.team_barrier();
              }
              for (Index c = rank; c < count; c += threads)
              {
                Kokkos::Geometry::Matrix h(vdim, vdim);
                for (Index r = 0; r < vdim; ++r)
                  for (Index k = 0; k < vdim; ++k)
                    h(r, k) = jacobians[c * Square + r * vdim + k];
                new (kinematics + c) Kinematics(h);
                const Real wdet =
                  quadrature.getWeight(q0 + c) * points[c].getDistortion();
                coefficients[3 * c] = wdet * mu;
                coefficients[3 * c + 1] = wdet * lambda;
                coefficients[3 * c + 2] = wdet * (mu - lambda * kinematics[c].logJ);
              }
              team.team_barrier();
              for (Index i = rank; i < count * na; i += threads)
              {
                const Index c = i / na;
                const Index a = i - c * na;
                aa[c * TrialCapacity + a] =
                  kinematics[c].pullBack(ga[c * TrialCapacity + a]);
              }
              if (!SharedTest)
              {
                for (Index i = rank; i < count * nb; i += threads)
                {
                  const Index c = i / nb;
                  const Index b = i - c * nb;
                  ab[c * TestCapacity + b] =
                    kinematics[c].pullBack(gb[c * TestCapacity + b]);
                }
              }
              team.team_barrier();
              if (slot < slots)
              {
                for (Index c = 0; c < count; ++c)
                {
                  // The products of the trial function, shared by the block.
                  Real va[D] = {};
                  Real ta[D] = {};
                  Real sa[D] = {};
                  {
                    const Vector& alpha = aa[c * TrialCapacity + a];
                    const Vector& gradient = ga[c * TrialCapacity + a];
                    const Real volume = coefficients[3 * c + 1];
                    const Real twist = coefficients[3 * c + 2];
                    const Real shear = coefficients[3 * c];
                    RODIN_KOKKOS_UNROLL
                    for (Index s = 0; s < D; ++s)
                    {
                      if (s < vdim)
                      {
                        va[s] = volume * alpha(s);
                        ta[s] = twist * alpha(s);
                        sa[s] = shear * gradient(s);
                      }
                    }
                  }
                  RODIN_KOKKOS_UNROLL
                  for (Index j = 0; j < E; ++j)
                  {
                    const Index b = b0 + j;
                    if (b < nb)
                    {
                      const Vector& beta = ab[c * TestCapacity + b];
                      const Vector& gradient = gb[c * TestCapacity + b];
                      Real bt[D] = {};
                      Real shear = 0;
                      RODIN_KOKKOS_UNROLL
                      for (Index t = 0; t < D; ++t)
                      {
                        if (t < vdim)
                        {
                          bt[t] = beta(t);
                          shear += sa[t] * gradient(t);
                        }
                      }
                      RODIN_KOKKOS_UNROLL
                      for (Index t = 0; t < D; ++t)
                      {
                        RODIN_KOKKOS_UNROLL
                        for (Index s = 0; s < D; ++s)
                          local[j][t][s] += va[s] * bt[t] + ta[t] * bt[s];
                        local[j][t][t] += shear;
                      }
                    }
                  }
                }
              }
            }
            RODIN_KOKKOS_UNROLL
            for (Index j = 0; j < E; ++j)
            {
              const Index b = b0 + j;
              if (slot >= slots || b >= nb)
                continue;
              RODIN_KOKKOS_UNROLL
              for (Index t = 0; t < D; ++t)
              {
                RODIN_KOKKOS_UNROLL
                for (Index s = 0; s < D; ++s)
                {
                  if (t >= vdim || s >= vdim)
                    continue;
                  const size_t e =
                    base + size_t((b * vdim + t) * ntr + a * vdim + s) * spacing;
                  const Real v = sign * local[j][t][s];
                  if (accumulate)
                    values(e) += v;
                  else
                    values(e) = v;
                  if (index)
                  {
                    rows(e) = testIndices[b] + t * testDOFs.getStride();
                    cols(e) = trialIndices[a] + s * trialDOFs.getStride();
                  }
                }
              }
            }
          }
        }
      };
  };

  /**
   * @brief Device residual of the internal virtual work of the Neo-Hookean
   * law.
   *
   * Represents
   * @f[
   *   \int_\Omega P(F(u)) : \nabla v \ dx, \quad
   *   P = \mu F + (\lambda \ln J - \mu) F^{-T},
   * @f]
   * for vector device spaces, as NeoHookean::getFirstPiolaKirchhoffStress().
   * The quadrature order is that of the host, @f$ 2 \max(k_v, k_u) @f$
   * unless set by setQuadratureOrder().
   */
  template <Kokkos::Variational::VectorDeviceSpace TestFES,
    Kokkos::Variational::VectorDeviceSpace StateFES, class Data>
  class InternalVirtualWorkResidual<NeoHookean, Variational::TestFunction<TestFES>,
    Variational::GridFunction<StateFES, Data>>
    final : public Variational::LinearFormIntegratorBase<Real>,
            public Kokkos::Variational::DeviceLinearFormIntegrator<
              Kokkos::Variational::Internal::BackendOf<TestFES>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Constitutive law type.
      using LawType = NeoHookean;

      /// @brief Test function type.
      using TestType = Variational::TestFunction<TestFES>;

      /// @brief Current displacement state type.
      using StateType = Variational::GridFunction<StateFES, Data>;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Displacement state finite element space type.
      using StateFESType = StateFES;

      /// @brief Host integrator interface.
      using Parent = Variational::LinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = Kokkos::Variational::DeviceLinearFormIntegrator<
        Kokkos::Variational::Internal::BackendOf<TestFES>>;

      /// @brief Mesh type.
      using MeshType = typename TestFESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device degrees of freedom of the test space.
      using DOFsType = typename TestFESType::DOFsType;

      /// @brief Device quadrature tables.
      using QuadratureType = Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Kokkos::Variational::Internal::DevicePolytopes<MeshType>;

      /// @brief Device coefficients of the displacement.
      using StateViewType = Kokkos::Variational::GridFunctionView<StateFESType>;

      /**
       * @brief Constructs the residual at a displacement.
       * @param law Neo-Hookean law.
       * @param v Test displacement.
       * @param displacement Current displacement.
       */
      InternalVirtualWorkResidual(
        const LawType& law, const TestType& v, const StateType& displacement)
        : Parent(v),
          m_law(law),
          m_test(v),
          m_displacement(displacement),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      /// @brief Copy constructor. @param other Residual to copy.
      InternalVirtualWorkResidual(const InternalVirtualWorkResidual& other)
        : Parent(other),
          DeviceParent(other),
          m_law(other.m_law),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_quadOrder(other.m_quadOrder)
      {}

      /// @brief Sets the displacement. @param gf Current displacement. @return Reference to this.
      InternalVirtualWorkResidual& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        return *this;
      }

      /// @brief Sets the quadrature order. @param order Order. @return Reference to this.
      InternalVirtualWorkResidual& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /// @brief Unavailable on the device. @throws Alert::Exception Always.
      InternalVirtualWorkResidual& setInput([[maybe_unused]] InputFunction input)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Constitutive inputs are not available on Kokkos devices." << Alert::Raise;
        return *this;
      }

      /// @brief Gets the law. @return Neo-Hookean law.
      const LawType& getLaw() const
      {
        return m_law;
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      InternalVirtualWorkResidual& setPolytope(
        [[maybe_unused]] const Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t local) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&fes.getMesh(), 0, static_cast<Index>(fes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(), 0,
          fes.getDeviceDOFs().getStride(), 0,
          DOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /// @brief Locates the entries of every cell. @return Number of entries.
      size_t prepare() const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        return m_polytopes.prepare(fes.getMesh(), this->getRegion(),
          this->getAttributes(), Kokkos::Variational::ConstantDOFs{}, 1,
          fes.getDeviceDOFs(), fes.getVectorDimension(),
          [&](const Geometry::Polytope& polytope) { return getOrder(polytope); });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /// @brief Integrates every cell on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        const auto& displacement = m_displacement.get();
        const StateViewType state(
          displacement.getFiniteElementSpace(), displacement.getData());
        const Index vertices = m_polytopes.getQuadrature().getMaximumVertexCount();
        Kokkos::Variational::visitVertexCapacity(vertices, [&](auto n) {
          integratePolytopes<decltype(n)::value>(fes.getMesh(), m_polytopes,
            fes.getDeviceDOFs(), fes.getVectorDimension(), state,
            m_law.getLameFirstParameter(), m_law.getShearModulus(), values, offset, sign,
            accumulate, indices);
          return true;
        });
      }

      /// @brief Gets the integration region. @return Cells.
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      InternalVirtualWorkResidual* copy() const noexcept override
      {
        return new InternalVirtualWorkResidual(*this);
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * Local index a * vdim + c is component c of the scalar degree of
       * freedom a.
       */
      template <Index Vertices>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const DOFsType& dofs, Index vdim, const StateViewType& state, Real lambda,
        Real mu, const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& indices)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const PolytopeKernel<Vertices> kernel{mesh.getDeviceView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(), dofs,
          vdim, state, lambda, mu, values, offset, sign, accumulate,
          indices.extent(0) > 0, indices};
        ::Kokkos::parallel_for(
          "Rodin hyperelastic residual integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

    private:
      /// @brief Selects the quadrature order of a cell, as on the host.
      size_t getOrder(const Geometry::Polytope& polytope) const
      {
        if (m_quadOrder > 0)
          return m_quadOrder;
        const auto testOrder = m_test.get().getOrder(polytope);
        const auto stateOrder = m_displacement.get().getOrder(polytope);
        if (!testOrder || !stateOrder)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic constitutive quadrature order is established. "
            << "Specify it with setQuadratureOrder()." << Alert::Raise;
        }
        return 2 * std::max(*testOrder, *stateOrder);
      }

      void checkCompatibility(const StateType& displacement) const
      {
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& statefes = displacement.getFiniteElementSpace();
        const size_t d = testfes.getMesh().getSpaceDimension();
        if (&statefes.getMesh() != &testfes.getMesh() ||
          testfes.getVectorDimension() != d || statefes.getVectorDimension() != d ||
          testfes.getMesh().getDimension() != d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos hyperelastic residuals require displacement spaces of the "
            << "space dimension on one mesh of full dimension." << Alert::Raise;
        }
      }

      LawType m_law;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      size_t m_quadOrder;
      mutable PolytopesType m_polytopes;

    public:
      /**
       * @brief Kernel of integratePolytopes(), the local vector of one cell
       * per thread.
       */
      template <Index Vertices>
      struct PolytopeKernel
      {
        static constexpr Index Capacity =
          DOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION;

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        DOFsType dofs;
        Index vdim;
        StateViewType state;
        Real lambda;
        Real mu;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> indices;

        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto testIndices = dofs.bind(polytope);
          const Index nb = testIndices.size();
          const Index nte = nb * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          ::Kokkos::Array<Real, Capacity> local{};
          for (Index q = quadrature.begin(g); q < quadrature.end(g); ++q)
          {
            const Rodin::Kokkos::Variational::CellIntegrationPointView ip(polytope,
              quadrature.getReference(polytope, q),
              Rodin::Kokkos::Variational::DynamicGeometry{});
            const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
            const Rodin::Kokkos::Variational::NeoHookeanPoint point(
              state.interpolateJacobian(ip));
            // P = mu F + (lambda ln J - mu) F^{-T}
            const Real scale = lambda * point.logJ - mu;
            for (Index b = 0; b < nb; ++b)
            {
              const auto& gradient = dofs.getGradient(ip, b);
              for (Index c = 0; c < vdim; ++c)
              {
                Real v = 0;
                for (Index k = 0; k < vdim; ++k)
                  v += (mu * point.F(c, k) + scale * point.Finv(k, c)) * gradient(k);
                local[b * vdim + c] += wdet * v;
              }
            }
          }
          const size_t base = offset + offsets(p);
          for (Index te = 0; te < nte; ++te)
          {
            const size_t e = base + te * spacing;
            const Real v = sign * local[te];
            if (accumulate)
              values(e) += v;
            else
              values(e) = v;
            if (index)
              indices(e) = testIndices[te / vdim] + (te % vdim) * dofs.getStride();
          }
        }
      };
  };

  /**
   * @brief Device tangent of the internal virtual work of a hyperelastic
   * law of Kokkos::Variational::DeviceLaw.
   *
   * Represents
   * @f[
   *   \int_\Omega \mathbb{A}(F(u)) : \nabla \delta u : \nabla v \ dx,
   *   \quad \mathbb{A}_{iJkL} = \frac{\partial P_{iJ}}{\partial F_{kL}},
   * @f]
   * for vector device spaces. For the vector Lagrange bases
   * @f$ \phi_a e_s @f$ and @f$ \phi_b e_t @f$, the entry is
   * @f$ \sum_{J, L} \partial_J \phi_b \ \mathbb{A}_{tJsL} \ \partial_L \phi_a @f$,
   * the contraction of the host dense tangent. The elasticity tensor is
   * evaluated once per quadrature point.
   *
   * The quadrature order is that of the host,
   * @f$ 2 \max(k_{\delta u}, k_v, k_u) @f$ unless set by
   * setQuadratureOrder(). The displacement is copied to the device at
   * every assembly. Constitutive inputs are not available on the device.
   */
  template <Kokkos::Variational::DeviceHyperElasticLaw Law, class Solution,
    Kokkos::Variational::VectorDeviceSpace TrialFES,
    Kokkos::Variational::VectorDeviceSpace TestFES,
    Kokkos::Variational::VectorDeviceSpace StateFES, class Data>
  class InternalVirtualWorkTangent<Law, Variational::TrialFunction<Solution, TrialFES>,
    Variational::TestFunction<TestFES>, Variational::GridFunction<StateFES, Data>>
    final : public Variational::LocalBilinearFormIntegratorBase<Real>,
            public Kokkos::Variational::DeviceBilinearFormIntegrator<
              Kokkos::Variational::Internal::BackendOf<TrialFES>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Constitutive law type.
      using LawType = Law;

      /// @brief Device law type.
      using DeviceLawType = typename Kokkos::Variational::DeviceLaw<Law>::Type;

      /// @brief Trial function type.
      using TrialType = Variational::TrialFunction<Solution, TrialFES>;

      /// @brief Test function type.
      using TestType = Variational::TestFunction<TestFES>;

      /// @brief Current displacement state type.
      using StateType = Variational::GridFunction<StateFES, Data>;

      /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Displacement state finite element space type.
      using StateFESType = StateFES;

      /// @brief Host integrator interface.
      using Parent = Variational::LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = Kokkos::Variational::DeviceBilinearFormIntegrator<
        Kokkos::Variational::Internal::BackendOf<TrialFES>>;

      /// @brief Mesh type.
      using MeshType = typename TrialFESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device degrees of freedom of the trial space.
      using TrialDOFsType = typename TrialFESType::DOFsType;

      /// @brief Device degrees of freedom of the test space.
      using TestDOFsType = typename TestFESType::DOFsType;

      /// @brief Device quadrature tables.
      using QuadratureType = Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Kokkos::Variational::Internal::DevicePolytopes<MeshType>;

      /// @brief Device coefficients of the displacement.
      using StateViewType = Kokkos::Variational::GridFunctionView<StateFESType>;

      static_assert(std::is_trivially_copyable_v<DeviceLawType>);

      /**
       * @brief Constructs the tangent at a displacement.
       * @param law Hyperelastic law.
       * @param u Trial displacement increment.
       * @param v Test displacement.
       * @param displacement Current displacement.
       */
      InternalVirtualWorkTangent(const LawType& law, const TrialType& u,
        const TestType& v, const StateType& displacement)
        : Parent(u, v),
          m_law(law),
          m_trial(u),
          m_test(v),
          m_displacement(displacement),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      /// @brief Copy constructor. @param other Tangent to copy.
      InternalVirtualWorkTangent(const InternalVirtualWorkTangent& other)
        : Parent(other),
          DeviceParent(other),
          m_law(other.m_law),
          m_trial(other.m_trial),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_quadOrder(other.m_quadOrder)
      {}

      /// @brief Sets the displacement. @param gf Current displacement. @return Reference to this.
      InternalVirtualWorkTangent& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        return *this;
      }

      /// @brief Sets the quadrature order. @param order Order. @return Reference to this.
      InternalVirtualWorkTangent& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /// @brief Unavailable on the device. @throws Alert::Exception Always.
      InternalVirtualWorkTangent& setInput([[maybe_unused]] InputFunction input)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Constitutive inputs are not available on Kokkos devices." << Alert::Raise;
        return *this;
      }

      /// @brief Gets the law. @return Hyperelastic law.
      const LawType& getLaw() const
      {
        return m_law;
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      InternalVirtualWorkTangent& setPolytope(
        [[maybe_unused]] const Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&trialfes.getMesh(), static_cast<Index>(trialfes.getVectorDimension()),
          static_cast<Index>(testfes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(),
          trialfes.getDeviceDOFs().getStride(), testfes.getDeviceDOFs().getStride(),
          TrialDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES),
          TestDOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /**
       * @brief Locates the entries of every cell.
       * @return Number of entries, laid out as by the device integrals.
       */
      size_t prepare(bool) const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        return m_polytopes.prepare(trialfes.getMesh(), this->getRegion(),
          this->getAttributes(), trialfes.getDeviceDOFs(), trialfes.getVectorDimension(),
          testfes.getDeviceDOFs(), testfes.getVectorDimension(),
          [&](const Geometry::Polytope& polytope) { return getOrder(polytope); });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /// @brief Integrates every cell on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& displacement = m_displacement.get();
        const StateViewType state(
          displacement.getFiniteElementSpace(), displacement.getData());
        const DeviceLawType law = Kokkos::Variational::DeviceLaw<Law>::get(m_law);
        const Index vertices = m_polytopes.getQuadrature().getMaximumVertexCount();
        Kokkos::Variational::visitVertexCapacity(vertices, [&](auto n) {
          constexpr Index Vertices = decltype(n)::value;
          // Larger local matrices take one team per cell.
          if constexpr (getLocalSize<Vertices>() <= LocalCapacity)
            integratePolytopes<Vertices>(trialfes.getMesh(), m_polytopes,
              trialfes.getDeviceDOFs(), testfes.getDeviceDOFs(),
              testfes.getVectorDimension(), state, law, values, offset, sign, accumulate,
              rows, cols);
          else
            integrateTeams<Vertices>(trialfes.getMesh(), m_polytopes,
              trialfes.getDeviceDOFs(), testfes.getDeviceDOFs(),
              testfes.getVectorDimension(), state, law, values, offset, sign, accumulate,
              rows, cols);
          return true;
        });
      }

      /// @brief Gets the integration region. @return Cells.
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      InternalVirtualWorkTangent* copy() const noexcept override
      {
        return new InternalVirtualWorkTangent(*this);
      }

      /**
       * @brief Gets the number of entries of the local matrix of a cell.
       * @tparam Vertices Maximum vertex count of the cells.
       * @return Entry count of the local matrix in the kernel.
       */
      template <Index Vertices>
      static constexpr Index getLocalSize()
      {
        return TrialDOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION *
          TestDOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION;
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * The entry of local test index te and local trial index tr has rank
       * te * ntr + tr among the entries of the cell, the local index
       * a * vdim + r being component r of the scalar degree of freedom a.
       * At each quadrature point, the thread evaluates the elasticity
       * tensor, contracts it with @f$ \nabla \phi_a @f$ for each trial
       * function, then with @f$ \nabla \phi_b @f$ for each test function.
       */
      template <Index Vertices>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const TrialDOFsType& trialDOFs, const TestDOFsType& testDOFs, Index vdim,
        const StateViewType& state, const DeviceLawType& law,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const PolytopeKernel<Vertices> kernel{mesh.getDeviceView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(),
          trialDOFs, testDOFs, vdim, state, law, values, offset, sign, accumulate,
          rows.extent(0) > 0, rows, cols};
        ::Kokkos::parallel_for(
          "Rodin hyperelastic law tangent integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices with one
       * team per cell.
       *
       * The quadrature points are taken by groups of P points. The threads
       * of the team build the integration points of a group in team scratch
       * memory, one point per thread, then the gradients
       * @f$ \nabla \phi_a @f$, one pair of a point and a scalar basis
       * function per thread, which the test and state spaces share when
       * they have the type of the trial space. They then sum the
       * displacement gradients, one pair of a point and a basis function
       * of the state per thread, and evaluate the elasticity tensors, one
       * pair of a point and a row of the tensor per thread. Each thread
       * then takes a slot, a trial function and a block of TeamPairs test
       * functions, and accumulates in registers the vdim by vdim entries of
       * each of its pairs over the group, contracting each row of the
       * tensor with the trial gradient once for the whole block. The
       * entries are laid out as by integratePolytopes().
       */
      template <Index Vertices>
      static void integrateTeams(const MeshType& mesh, const PolytopesType& polytopes,
        const TrialDOFsType& trialDOFs, const TestDOFsType& testDOFs, Index vdim,
        const StateViewType& state, const DeviceLawType& law,
        const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset, Real sign,
        bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols)
      {
        using Kernel = TeamKernel<Vertices>;
        using Policy = typename Kernel::Policy;
        static_assert(std::is_trivially_destructible_v<typename Kernel::Point> &&
          std::is_trivially_copyable_v<typename Kernel::Vector>);
        constexpr Index TrialCapacity = Kernel::TrialCapacity;
        constexpr Index TestCapacity = Kernel::TestCapacity;
        constexpr Index E = Kernel::E;
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t scratch = Kernel::P * Kernel::PointBytes + 8 * Kernel::Alignment;
        const int level = scratch <= size_t(Policy::scratch_size_max(0)) ? 0 : 1;
        const Kernel kernel{mesh.getDeviceView(), polytopes.getQuadrature(),
          polytopes.getAttributes(), polytopes.getOffsets(), polytopes.getSpacing(),
          polytopes.getRegion(), polytopes.getDimension(), trialDOFs, testDOFs, vdim,
          state, law, values, offset, sign, accumulate, rows.extent(0) > 0, rows, cols,
          level};
        Policy policy(execution, polytopes.getPolytopeCount(), ::Kokkos::AUTO);
        policy.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        // One slot per thread, the fewest passes over the quadrature.
        constexpr Index Slots = TrialCapacity * ((TestCapacity + E - 1) / E);
        const Index wanted = (Slots + 31) / 32 * 32;
        const int maximum = policy.team_size_max(kernel, ::Kokkos::ParallelForTag());
        const int threads =
          std::max(1, std::min(maximum, int(std::min<Index>(wanted, TeamThreads))));
        Policy launch(execution, polytopes.getPolytopeCount(), threads);
        launch.set_scratch_size(level, ::Kokkos::PerTeam(scratch));
        ::Kokkos::parallel_for(
          "Rodin hyperelastic law tangent team integration", launch, kernel);
      }

    private:
      /// @brief Maximum number of entries of the local matrix of a thread,
      /// so larger local matrices take one team per cell.
      static constexpr Index LocalCapacity = 1024;

      /// @brief Largest group of quadrature points of a team.
      static constexpr Index TeamPoints = 32;

      /// @brief Bytes of team scratch memory of a group in the fast level.
      static constexpr size_t TeamScratch = 47104;

      /// @brief Fewest points of a group in the fast scratch level.
      static constexpr Index TeamMinimumPoints = 4;

      /// @brief Number of points of a group whose values exceed TeamScratch,
      /// in the slower scratch level.
      static constexpr Index TeamSlowPoints = 8;

      /// @brief Largest team, whose register bound the kernel is compiled for.
      static constexpr unsigned int TeamThreads = 512;

      /// @brief Test functions of the block of a thread, whose entries
      /// with one trial function it accumulates per pass over the quadrature.
      static constexpr Index TeamPairs = 4;

      /// @brief Selects the quadrature order of a cell, as on the host.
      size_t getOrder(const Geometry::Polytope& polytope) const
      {
        if (m_quadOrder > 0)
          return m_quadOrder;
        const auto trialOrder = m_trial.get().getOrder(polytope);
        const auto testOrder = m_test.get().getOrder(polytope);
        const auto stateOrder = m_displacement.get().getOrder(polytope);
        if (!trialOrder || !testOrder || !stateOrder)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic constitutive quadrature order is established. "
            << "Specify it with setQuadratureOrder()." << Alert::Raise;
        }
        return 2 * std::max({*trialOrder, *testOrder, *stateOrder});
      }

      void checkCompatibility(const StateType& displacement) const
      {
        const auto& trialfes = m_trial.get().getFiniteElementSpace();
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& statefes = displacement.getFiniteElementSpace();
        const size_t d = trialfes.getMesh().getSpaceDimension();
        if (&trialfes.getMesh() != &testfes.getMesh() ||
          &statefes.getMesh() != &testfes.getMesh() ||
          trialfes.getVectorDimension() != d || testfes.getVectorDimension() != d ||
          statefes.getVectorDimension() != d || trialfes.getMesh().getDimension() != d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos hyperelastic tangents require displacement spaces of the "
            << "space dimension on one mesh of full dimension." << Alert::Raise;
        }
      }

      LawType m_law;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      size_t m_quadOrder;
      mutable PolytopesType m_polytopes;

    public:
      /**
       * @brief Kernel of integratePolytopes(), the local matrix of one cell
       * per thread.
       */
      template <Index Vertices>
      struct PolytopeKernel
      {
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices);
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices);
        static constexpr Index D = Kokkos::Variational::ElasticityStride;
        static constexpr Index Row = Kokkos::Variational::ElasticityRow;

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        TrialDOFsType trialDOFs;
        TestDOFsType testDOFs;
        Index vdim;
        StateViewType state;
        DeviceLawType law;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> rows;
        ::Kokkos::View<Index*, MemorySpace> cols;

        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const Index na = trialIndices.size();
          const Index nb = testIndices.size();
          const Index ntr = na * vdim;
          const Index nte = nb * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          ::Kokkos::Array<Real,
            TrialCapacity * RODIN_MAXIMAL_SPACE_DIMENSION * TestCapacity *
              RODIN_MAXIMAL_SPACE_DIMENSION>
            local{};
          ::Kokkos::Array<Rodin::Kokkos::Geometry::Vector, TestCapacity> gb;
          Real elasticity[D * Row];
          for (Index q = quadrature.begin(g); q < quadrature.end(g); ++q)
          {
            const Rodin::Kokkos::Variational::CellIntegrationPointView ip(polytope,
              quadrature.getReference(polytope, q),
              Rodin::Kokkos::Variational::DynamicGeometry{});
            const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
            {
              const auto point = law.getPoint(state.interpolateJacobian(ip));
              for (Index t = 0; t < vdim; ++t)
                law.getElasticity(point, t, wdet, elasticity + t * Row);
            }
            for (Index b = 0; b < nb; ++b)
              gb[b] = testDOFs.getGradient(ip, b);
            for (Index a = 0; a < na; ++a)
            {
              const auto ga = trialDOFs.getGradient(ip, a);
              // m(t, J, s) = sum_L A(t, J, s, L) d_L phi_a
              Real m[D * D * D];
              for (Index t = 0; t < vdim; ++t)
              {
                for (Index J = 0; J < vdim; ++J)
                {
                  for (Index s = 0; s < vdim; ++s)
                  {
                    const Real* const A = elasticity + t * Row +
                      Rodin::Kokkos::Variational::getElasticityIndex(J, s, 0);
                    Real v = 0;
                    for (Index L = 0; L < vdim; ++L)
                      v += A[L] * ga(L);
                    m[(t * D + J) * D + s] = v;
                  }
                }
              }
              for (Index b = 0; b < nb; ++b)
              {
                for (Index t = 0; t < vdim; ++t)
                {
                  for (Index s = 0; s < vdim; ++s)
                  {
                    Real v = 0;
                    for (Index J = 0; J < vdim; ++J)
                      v += gb[b](J) * m[(t * D + J) * D + s];
                    local[(b * vdim + t) * ntr + a * vdim + s] += v;
                  }
                }
              }
            }
          }
          const size_t base = offset + offsets(p);
          for (Index te = 0; te < nte; ++te)
          {
            const Index row =
              testIndices[te / vdim] + (te % vdim) * testDOFs.getStride();
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
                cols(e) = trialIndices[tr / vdim] + (tr % vdim) * trialDOFs.getStride();
              }
            }
          }
        }
      };

      /**
       * @brief Kernel of integrateTeams(), the local matrix of one cell per
       * team.
       */
      template <Index Vertices>
      struct TeamKernel
      {
        using Point = Kokkos::Variational::CellIntegrationPointView;
        using Vector = Kokkos::Geometry::Vector;
        using ExecutionSpace = typename MeshType::ExecutionSpace;
        // The register bound keeps enough threads for one pass over the
        // quadrature, at the price of a few spills.
        using Policy = ::Kokkos::TeamPolicy<ExecutionSpace, ::Kokkos::IndexType<Index>,
          ::Kokkos::LaunchBounds<TeamThreads, 1>>;
        using Member = typename Policy::member_type;
        using StateDOFsType =
          std::decay_t<decltype(std::declval<const StateViewType&>().getDOFs())>;
        static constexpr Index TrialCapacity = TrialDOFsType::getCapacity(Vertices);
        static constexpr Index TestCapacity = TestDOFsType::getCapacity(Vertices);
        static constexpr Index D = Kokkos::Variational::ElasticityStride;
        static constexpr Index Row = Kokkos::Variational::ElasticityRow;
        static constexpr Index Square = D * D;
        static constexpr Index Tensor = D * Row;
        // Spaces of one type share their basis gradients.
        static constexpr bool SharedTest = std::is_same_v<TrialDOFsType, TestDOFsType>;
        static constexpr bool SharedState = std::is_same_v<TrialDOFsType, StateDOFsType>;
        static constexpr Index TestStorage = SharedTest ? 0 : TestCapacity;
        static constexpr size_t PointBytes = sizeof(Point) +
          (Square + Tensor) * sizeof(Real) +
          (TrialCapacity + TestStorage) * sizeof(Vector);
        // Groups that do not fit TeamScratch go to the slower scratch level.
        static constexpr Index Fit = TeamScratch / PointBytes;
        static constexpr Index P = Fit < TeamMinimumPoints ? TeamSlowPoints
          : Fit > TeamPoints                               ? TeamPoints
                                                           : Fit;
        static constexpr Index E = TeamPairs;
        static constexpr size_t Alignment =
          std::max({alignof(Point), alignof(Vector), alignof(Real)});

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        TrialDOFsType trialDOFs;
        TestDOFsType testDOFs;
        Index vdim;
        StateViewType state;
        DeviceLawType law;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> rows;
        ::Kokkos::View<Index*, MemorySpace> cols;
        int level;

        KOKKOS_INLINE_FUNCTION void operator()(const Member& team) const
        {
          const Index p = team.league_rank();
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto trialIndices = trialDOFs.bind(polytope);
          const auto testIndices = testDOFs.bind(polytope);
          const auto& stateDOFs = state.getDOFs();
          const auto& data = state.getData();
          const auto stateIndices = stateDOFs.bind(polytope);
          const Index ns = stateIndices.size();
          const Index na = trialIndices.size();
          const Index nb = testIndices.size();
          const Index ntr = na * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          const Index rank = team.team_rank();
          const Index threads = team.team_size();
          auto& memory = team.team_scratch(level);
          auto* const points =
            static_cast<Point*>(memory.get_shmem_aligned(P * sizeof(Point), Alignment));
          auto* const jacobians = static_cast<Real*>(
            memory.get_shmem_aligned(Square * P * sizeof(Real), Alignment));
          auto* const tensors = static_cast<Real*>(
            memory.get_shmem_aligned(Tensor * P * sizeof(Real), Alignment));
          auto* const ga = static_cast<Vector*>(
            memory.get_shmem_aligned(P * TrialCapacity * sizeof(Vector), Alignment));
          Vector* gb = ga;
          if (!SharedTest)
            gb = static_cast<Vector*>(
              memory.get_shmem_aligned(P * TestCapacity * sizeof(Vector), Alignment));
          // A slot is a trial function and a block of E test functions.
          const Index slots = na * ((nb + E - 1) / E);
          const Index begin = quadrature.begin(g);
          const Index end = quadrature.end(g);
          const size_t base = offset + offsets(p);
          // Groups of balanced sizes, at most P points each.
          const Index groups = (end - begin + P - 1) / P;
          const Index size = groups > 0 ? (end - begin + groups - 1) / groups : 0;
          for (Index first = 0; first < slots; first += threads)
          {
            const Index slot = first + rank;
            const Index a = slot % na;
            const Index b0 = (slot / na) * E;
            Real local[E][D][D] = {};
            for (Index q0 = begin; q0 < end; q0 += size)
            {
              const Index count = end - q0 < size ? end - q0 : size;
              // The values of the previous group are no longer read.
              team.team_barrier();
              for (Index c = rank; c < count; c += threads)
              {
                new (points + c)
                  Point(polytope, quadrature.getReference(polytope, q0 + c),
                    Rodin::Kokkos::Variational::DynamicGeometry{});
                for (Index k = 0; k < Square; ++k)
                  jacobians[c * Square + k] = 0;
              }
              team.team_barrier();
              for (Index i = rank; i < count * na; i += threads)
              {
                const Index c = i / na;
                const Index a = i - c * na;
                ga[c * TrialCapacity + a] = trialDOFs.getGradient(points[c], a);
              }
              if (!SharedTest)
              {
                for (Index i = rank; i < count * nb; i += threads)
                {
                  const Index c = i / nb;
                  const Index b = i - c * nb;
                  gb[c * TestCapacity + b] = testDOFs.getGradient(points[c], b);
                }
              }
              // The displacement gradient sums one term per pair of a point
              // and a basis function of the state, after the trial gradients
              // when the state shares them.
              const auto sum = [&](Index i) {
                const Index c = i / ns;
                const Index q = i - c * ns;
                Vector gradient;
                if (SharedState)
                  gradient = ga[c * TrialCapacity + q];
                else
                  gradient = stateDOFs.getGradient(points[c], q);
                for (Index r = 0; r < vdim; ++r)
                {
                  const Real value = data(stateIndices[q] + r * stateDOFs.getStride());
                  for (Index k = 0; k < vdim; ++k)
                    ::Kokkos::atomic_add(
                      &jacobians[c * Square + r * vdim + k], gradient(k) * value);
                }
              };
              if (!SharedState)
              {
                for (Index i = rank; i < count * ns; i += threads)
                  sum(i);
              }
              team.team_barrier();
              if (SharedState)
              {
                for (Index i = rank; i < count * ns; i += threads)
                  sum(i);
                team.team_barrier();
              }
              // Each row of an elasticity tensor rebuilds the kinematics of
              // its point.
              for (Index i = rank; i < count * vdim; i += threads)
              {
                const Index c = i / vdim;
                const Index t = i - c * vdim;
                Kokkos::Geometry::Matrix h(vdim, vdim);
                for (Index r = 0; r < vdim; ++r)
                  for (Index k = 0; k < vdim; ++k)
                    h(r, k) = jacobians[c * Square + r * vdim + k];
                const Real wdet =
                  quadrature.getWeight(q0 + c) * points[c].getDistortion();
                law.getElasticity(
                  law.getPoint(h), t, wdet, tensors + c * Tensor + t * Row);
              }
              team.team_barrier();
              if (slot < slots)
              {
                for (Index c = 0; c < count; ++c)
                {
                  Real alpha[D] = {};
                  {
                    const Vector& gradient = ga[c * TrialCapacity + a];
                    RODIN_KOKKOS_UNROLL
                    for (Index L = 0; L < D; ++L)
                      if (L < vdim)
                        alpha[L] = gradient(L);
                  }
                  Real beta[E][D] = {};
                  RODIN_KOKKOS_UNROLL
                  for (Index j = 0; j < E; ++j)
                  {
                    if (b0 + j < nb)
                    {
                      const Vector& gradient = gb[c * TestCapacity + b0 + j];
                      RODIN_KOKKOS_UNROLL
                      for (Index J = 0; J < D; ++J)
                        if (J < vdim)
                          beta[j][J] = gradient(J);
                    }
                  }
                  const Real* const A = tensors + c * Tensor;
                  // Row t of the tensor contracted with the trial gradient,
                  // shared by the block. Unused components are zero.
                  RODIN_KOKKOS_UNROLL
                  for (Index t = 0; t < D; ++t)
                  {
                    if (t < vdim)
                    {
                      Real m[D][D] = {};
                      RODIN_KOKKOS_UNROLL
                      for (Index J = 0; J < D; ++J)
                      {
                        RODIN_KOKKOS_UNROLL
                        for (Index s = 0; s < D; ++s)
                        {
                          const Real* const AJ = A + t * Row +
                            Rodin::Kokkos::Variational::getElasticityIndex(J, s, 0);
                          RODIN_KOKKOS_UNROLL
                          for (Index L = 0; L < D; ++L)
                            if (J < vdim && s < vdim && L < vdim)
                              m[J][s] += AJ[L] * alpha[L];
                        }
                      }
                      RODIN_KOKKOS_UNROLL
                      for (Index j = 0; j < E; ++j)
                      {
                        RODIN_KOKKOS_UNROLL
                        for (Index J = 0; J < D; ++J)
                        {
                          RODIN_KOKKOS_UNROLL
                          for (Index s = 0; s < D; ++s)
                            local[j][t][s] += beta[j][J] * m[J][s];
                        }
                      }
                    }
                  }
                }
              }
            }
            RODIN_KOKKOS_UNROLL
            for (Index j = 0; j < E; ++j)
            {
              const Index b = b0 + j;
              if (slot >= slots || b >= nb)
                continue;
              RODIN_KOKKOS_UNROLL
              for (Index t = 0; t < D; ++t)
              {
                RODIN_KOKKOS_UNROLL
                for (Index s = 0; s < D; ++s)
                {
                  if (t >= vdim || s >= vdim)
                    continue;
                  const size_t e =
                    base + size_t((b * vdim + t) * ntr + a * vdim + s) * spacing;
                  const Real v = sign * local[j][t][s];
                  if (accumulate)
                    values(e) += v;
                  else
                    values(e) = v;
                  if (index)
                  {
                    rows(e) = testIndices[b] + t * testDOFs.getStride();
                    cols(e) = trialIndices[a] + s * trialDOFs.getStride();
                  }
                }
              }
            }
          }
        }
      };
  };

  /**
   * @brief Device residual of the internal virtual work of a hyperelastic
   * law of Kokkos::Variational::DeviceLaw.
   *
   * Represents
   * @f[
   *   \int_\Omega P(F(u)) : \nabla v \ dx
   * @f]
   * for vector device spaces, with the first Piola-Kirchhoff stress of the
   * device law. The quadrature order is that of the host,
   * @f$ 2 \max(k_v, k_u) @f$ unless set by setQuadratureOrder().
   */
  template <Kokkos::Variational::DeviceHyperElasticLaw Law,
    Kokkos::Variational::VectorDeviceSpace TestFES,
    Kokkos::Variational::VectorDeviceSpace StateFES, class Data>
  class InternalVirtualWorkResidual<Law, Variational::TestFunction<TestFES>,
    Variational::GridFunction<StateFES, Data>>
    final : public Variational::LinearFormIntegratorBase<Real>,
            public Kokkos::Variational::DeviceLinearFormIntegrator<
              Kokkos::Variational::Internal::BackendOf<TestFES>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Constitutive law type.
      using LawType = Law;

      /// @brief Device law type.
      using DeviceLawType = typename Kokkos::Variational::DeviceLaw<Law>::Type;

      /// @brief Test function type.
      using TestType = Variational::TestFunction<TestFES>;

      /// @brief Current displacement state type.
      using StateType = Variational::GridFunction<StateFES, Data>;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Displacement state finite element space type.
      using StateFESType = StateFES;

      /// @brief Host integrator interface.
      using Parent = Variational::LinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = Kokkos::Variational::DeviceLinearFormIntegrator<
        Kokkos::Variational::Internal::BackendOf<TestFES>>;

      /// @brief Mesh type.
      using MeshType = typename TestFESType::MeshType;

      /// @brief Kokkos memory space.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Device degrees of freedom of the test space.
      using DOFsType = typename TestFESType::DOFsType;

      /// @brief Device quadrature tables.
      using QuadratureType = Kokkos::Variational::Internal::DeviceQuadrature<MeshType>;

      /// @brief Integrated polytopes.
      using PolytopesType = Kokkos::Variational::Internal::DevicePolytopes<MeshType>;

      /// @brief Device coefficients of the displacement.
      using StateViewType = Kokkos::Variational::GridFunctionView<StateFESType>;

      static_assert(std::is_trivially_copyable_v<DeviceLawType>);

      /**
       * @brief Constructs the residual at a displacement.
       * @param law Hyperelastic law.
       * @param v Test displacement.
       * @param displacement Current displacement.
       */
      InternalVirtualWorkResidual(
        const LawType& law, const TestType& v, const StateType& displacement)
        : Parent(v),
          m_law(law),
          m_test(v),
          m_displacement(displacement),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      /// @brief Copy constructor. @param other Residual to copy.
      InternalVirtualWorkResidual(const InternalVirtualWorkResidual& other)
        : Parent(other),
          DeviceParent(other),
          m_law(other.m_law),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_quadOrder(other.m_quadOrder)
      {}

      /// @brief Sets the displacement. @param gf Current displacement. @return Reference to this.
      InternalVirtualWorkResidual& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        return *this;
      }

      /// @brief Sets the quadrature order. @param order Order. @return Reference to this.
      InternalVirtualWorkResidual& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /// @brief Unavailable on the device. @throws Alert::Exception Always.
      InternalVirtualWorkResidual& setInput([[maybe_unused]] InputFunction input)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Constitutive inputs are not available on Kokkos devices." << Alert::Raise;
        return *this;
      }

      /// @brief Gets the law. @return Hyperelastic law.
      const LawType& getLaw() const
      {
        return m_law;
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      InternalVirtualWorkResidual& setPolytope(
        [[maybe_unused]] const Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t local) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Kokkos::Variational::DeviceLayout getLayout() const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        const auto& attributes = this->getAttributes();
        return {&fes.getMesh(), 0, static_cast<Index>(fes.getVectorDimension()),
          {attributes.begin(), attributes.end()}, this->getRegion(), 0,
          fes.getDeviceDOFs().getStride(), 0,
          DOFsType::getCapacity(RODIN_MAXIMUM_POLYTOPE_VERTICES)};
      }

      /// @brief Locates the entries of every cell. @return Number of entries.
      size_t prepare() const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        return m_polytopes.prepare(fes.getMesh(), this->getRegion(),
          this->getAttributes(), Kokkos::Variational::ConstantDOFs{}, 1,
          fes.getDeviceDOFs(), fes.getVectorDimension(),
          [&](const Geometry::Polytope& polytope) { return getOrder(polytope); });
      }

      /// @brief Gets the spacing of the entries. @return Spacing of the last integration.
      size_t getSpacing() const override
      {
        return m_polytopes.getSpacing();
      }

      /// @brief Integrates every cell on the device.
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices) const override
      {
        const auto& fes = m_test.get().getFiniteElementSpace();
        const auto& displacement = m_displacement.get();
        const StateViewType state(
          displacement.getFiniteElementSpace(), displacement.getData());
        const DeviceLawType law = Kokkos::Variational::DeviceLaw<Law>::get(m_law);
        const Index vertices = m_polytopes.getQuadrature().getMaximumVertexCount();
        Kokkos::Variational::visitVertexCapacity(vertices, [&](auto n) {
          integratePolytopes<decltype(n)::value>(fes.getMesh(), m_polytopes,
            fes.getDeviceDOFs(), fes.getVectorDimension(), state, law, values, offset,
            sign, accumulate, indices);
          return true;
        });
      }

      /// @brief Gets the integration region. @return Cells.
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      InternalVirtualWorkResidual* copy() const noexcept override
      {
        return new InternalVirtualWorkResidual(*this);
      }

      /**
       * @brief Integrates the cells of at most Vertices vertices.
       *
       * Local index a * vdim + c is component c of the scalar degree of
       * freedom a.
       */
      template <Index Vertices>
      static void integratePolytopes(const MeshType& mesh, const PolytopesType& polytopes,
        const DOFsType& dofs, Index vdim, const StateViewType& state,
        const DeviceLawType& law, const ::Kokkos::View<Real*, MemorySpace>& values,
        size_t offset, Real sign, bool accumulate,
        const ::Kokkos::View<Index*, MemorySpace>& indices)
      {
        const auto& execution = mesh.getContext().getExecutionSpace();
        const PolytopeKernel<Vertices> kernel{mesh.getDeviceView(),
          polytopes.getQuadrature(), polytopes.getAttributes(), polytopes.getOffsets(),
          polytopes.getSpacing(), polytopes.getRegion(), polytopes.getDimension(), dofs,
          vdim, state, law, values, offset, sign, accumulate, indices.extent(0) > 0,
          indices};
        ::Kokkos::parallel_for(
          "Rodin hyperelastic law residual integration",
          QuadratureType::policy(execution, polytopes.getPolytopeCount()), kernel);
      }

    private:
      /// @brief Selects the quadrature order of a cell, as on the host.
      size_t getOrder(const Geometry::Polytope& polytope) const
      {
        if (m_quadOrder > 0)
          return m_quadOrder;
        const auto testOrder = m_test.get().getOrder(polytope);
        const auto stateOrder = m_displacement.get().getOrder(polytope);
        if (!testOrder || !stateOrder)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No automatic constitutive quadrature order is established. "
            << "Specify it with setQuadratureOrder()." << Alert::Raise;
        }
        return 2 * std::max(*testOrder, *stateOrder);
      }

      void checkCompatibility(const StateType& displacement) const
      {
        const auto& testfes = m_test.get().getFiniteElementSpace();
        const auto& statefes = displacement.getFiniteElementSpace();
        const size_t d = testfes.getMesh().getSpaceDimension();
        if (&statefes.getMesh() != &testfes.getMesh() ||
          testfes.getVectorDimension() != d || statefes.getVectorDimension() != d ||
          testfes.getMesh().getDimension() != d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos hyperelastic residuals require displacement spaces of the "
            << "space dimension on one mesh of full dimension." << Alert::Raise;
        }
      }

      LawType m_law;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      size_t m_quadOrder;
      mutable PolytopesType m_polytopes;

    public:
      /**
       * @brief Kernel of integratePolytopes(), the local vector of one cell
       * per thread.
       */
      template <Index Vertices>
      struct PolytopeKernel
      {
        static constexpr Index Capacity =
          DOFsType::getCapacity(Vertices) * RODIN_MAXIMAL_SPACE_DIMENSION;

        Rodin::Kokkos::Geometry::MeshView view;
        QuadratureType quadrature;
        typename PolytopesType::AttributesType attributes;
        ::Kokkos::View<Index*, MemorySpace> offsets;
        size_t spacing;
        Geometry::Region region;
        Index dimension;
        DOFsType dofs;
        Index vdim;
        StateViewType state;
        DeviceLawType law;
        ::Kokkos::View<Real*, MemorySpace> values;
        size_t offset;
        Real sign;
        bool accumulate;
        bool index;
        ::Kokkos::View<Index*, MemorySpace> indices;

        KOKKOS_INLINE_FUNCTION void operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Variational::isIntegrated(
                view, region, polytope, attributes))
            return;
          KOKKOS_ASSERT(polytope.getVertices().size() <= Vertices);
          const auto testIndices = dofs.bind(polytope);
          const Index nb = testIndices.size();
          const Index nte = nb * vdim;
          const Index g = static_cast<Index>(polytope.getGeometry());
          ::Kokkos::Array<Real, Capacity> local{};
          for (Index q = quadrature.begin(g); q < quadrature.end(g); ++q)
          {
            const Rodin::Kokkos::Variational::CellIntegrationPointView ip(polytope,
              quadrature.getReference(polytope, q),
              Rodin::Kokkos::Variational::DynamicGeometry{});
            const Real wdet = quadrature.getWeight(q) * ip.getDistortion();
            Rodin::Kokkos::Geometry::Matrix stress;
            law.getStress(law.getPoint(state.interpolateJacobian(ip)), stress);
            for (Index b = 0; b < nb; ++b)
            {
              const auto& gradient = dofs.getGradient(ip, b);
              for (Index c = 0; c < vdim; ++c)
              {
                Real v = 0;
                for (Index k = 0; k < vdim; ++k)
                  v += stress(c, k) * gradient(k);
                local[b * vdim + c] += wdet * v;
              }
            }
          }
          const size_t base = offset + offsets(p);
          for (Index te = 0; te < nte; ++te)
          {
            const size_t e = base + te * spacing;
            const Real v = sign * local[te];
            if (accumulate)
              values(e) += v;
            else
              values(e) = v;
            if (index)
              indices(e) = testIndices[te / vdim] + (te % vdim) * dofs.getStride();
          }
        }
      };
  };
}

#endif
