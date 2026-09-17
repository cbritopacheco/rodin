/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Finite-element integrator kernel and assembly benchmarks.
 *
 * The kernel cases isolate element binding, quadrature, basis evaluation, and
 * local entry generation. The assembly cases include mesh traversal,
 * local-to-global scatter, and sparse-operator construction. Keeping both
 * levels prevents assembly costs from hiding a regression in a specialized
 * QuadratureRule and prevents a faster kernel from being mistaken for a
 * faster assembled form.
 */

#include <benchmark/benchmark.h>

#include <cmath>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/H1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{

  namespace
  {
    template <size_t Order>
    struct H1Family
    {
        template <class Range>
        using FES = H1<Order, Range, LocalMesh>;

        static FES<Real> scalar(const LocalMesh& mesh)
        {
          return FES<Real>(std::integral_constant<size_t, Order>{}, mesh);
        }

        static FES<Math::SpatialVector<Real>> vector(
          const LocalMesh& mesh, size_t dimension)
        {
          return FES<Math::SpatialVector<Real>>(
            std::integral_constant<size_t, Order>{}, mesh, dimension);
        }
    };

    struct P1Family
    {
        template <class Range>
        using FES = P1<Range, LocalMesh>;

        static FES<Real> scalar(const LocalMesh& mesh)
        {
          return FES<Real>(mesh);
        }

        static FES<Math::SpatialVector<Real>> vector(
          const LocalMesh& mesh, size_t dimension)
        {
          return FES<Math::SpatialVector<Real>>(mesh, dimension);
        }
    };

    template <size_t Dimension>
    LocalMesh makeMesh()
    {
      if constexpr (Dimension == 2)
      {
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {10, 10});
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
      else
      {
        auto mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {4, 4, 4});
        mesh.getConnectivity().compute(3, 2);
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        mesh.getConnectivity().compute(2, 3);
        return mesh;
      }
    }

    LocalMesh makeGeometryMesh(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Segment:
          return LocalMesh::UniformGrid(geometry, {10});
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          return LocalMesh::UniformGrid(geometry, {10, 10});
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Wedge:
        case Polytope::Type::Pyramid:
          return LocalMesh::UniformGrid(geometry, {4, 4, 4});
        default:
          assert(false);
          return {};
      }
    }

    template <size_t Dimension>
    void installP2Transformations(LocalMesh& mesh)
    {
      const Polytope::Type geometry =
        Dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron;
      RealH1Element<2> geometryElement(geometry);
      Math::SpatialPoint physical;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        PointCloud nodes(Dimension, geometryElement.getCount());
        for (size_t a = 0; a < geometryElement.getCount(); ++a)
        {
          cell->getTransformation().transform(physical, geometryElement.getNode(a));
          Real warp = 0.05 * std::sin(0.3 * physical(0));
          for (size_t d = 1; d < Dimension; ++d)
            warp *= std::sin(0.3 * physical(d));
          physical(Dimension - 1) += warp;
          for (size_t d = 0; d < Dimension; ++d)
            nodes(d, a) = physical(d);
        }
        mesh.setPolytopeTransformation({Dimension, cell->getIndex()},
          new ParametricTransformation<RealH1Element<2>>(
            std::move(nodes), geometryElement));
      }
    }

    template <class IntegralType, class TrialFES, class TestFES>
    void runBilinearKernel(benchmark::State& state, IntegralType& integral,
      const TrialFES& trialFES, const TestFES& testFES)
    {
      const auto& mesh = trialFES.getMesh();
      const size_t d = mesh.getDimension();
      const Index cells = mesh.getPolytopeCount(d);
      int64_t entries = 0;

      for (Index i = 0; i < cells; ++i)
      {
        entries += static_cast<int64_t>(trialFES.getDOFs(d, i).size()) *
          static_cast<int64_t>(testFES.getDOFs(d, i).size());
      }

      for (auto _ : state)
      {
        Real checksum = 0;
        for (Index i = 0; i < cells; ++i)
        {
          const auto polytope = mesh.getPolytope(d, i);
          integral.setPolytope(*polytope);
          const auto& trialDOFs = trialFES.getDOFs(d, i);
          const auto& testDOFs = testFES.getDOFs(d, i);
          for (size_t l = 0; l < static_cast<size_t>(testDOFs.size()); ++l)
            for (size_t m = 0; m < static_cast<size_t>(trialDOFs.size()); ++m)
              checksum += integral.integrate(m, l);
        }
        benchmark::DoNotOptimize(checksum);
      }

      state.SetItemsProcessed(state.iterations() * entries);
      state.counters["cells"] = static_cast<double>(cells);
      state.counters["entries"] = static_cast<double>(entries);
    }

    template <class IntegralType, class FES>
    void runLinearKernel(benchmark::State& state, IntegralType& integral, const FES& fes)
    {
      const auto& mesh = fes.getMesh();
      const size_t d = mesh.getDimension();
      const Index cells = mesh.getPolytopeCount(d);
      int64_t entries = 0;
      for (Index i = 0; i < cells; ++i)
        entries += static_cast<int64_t>(fes.getDOFs(d, i).size());

      for (auto _ : state)
      {
        Real checksum = 0;
        for (Index i = 0; i < cells; ++i)
        {
          const auto polytope = mesh.getPolytope(d, i);
          integral.setPolytope(*polytope);
          const auto& dofs = fes.getDOFs(d, i);
          for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
            checksum += integral.integrate(l);
        }
        benchmark::DoNotOptimize(checksum);
      }

      state.SetItemsProcessed(state.iterations() * entries);
      state.counters["cells"] = static_cast<double>(cells);
      state.counters["entries"] = static_cast<double>(entries);
    }

    template <class IntegralType, class TrialFES, class TestFES>
    void runGlobalBilinearKernel(benchmark::State& state, IntegralType& integral,
      const TrialFES& trialFES, const TestFES& testFES)
    {
      const auto& mesh = trialFES.getMesh();
      const size_t d = mesh.getDimension();
      const Index cells = mesh.getPolytopeCount(d);
      int64_t entries = 0;
      for (Index te = 0; te < cells; ++te)
        for (Index tr = 0; tr < cells; ++tr)
          entries += static_cast<int64_t>(trialFES.getDOFs(d, tr).size()) *
            static_cast<int64_t>(testFES.getDOFs(d, te).size());

      for (auto _ : state)
      {
        Real checksum = 0;
        for (Index te = 0; te < cells; ++te)
        {
          const auto testPolytope = mesh.getPolytope(d, te);
          const auto& testDOFs = testFES.getDOFs(d, te);
          for (Index tr = 0; tr < cells; ++tr)
          {
            const auto trialPolytope = mesh.getPolytope(d, tr);
            integral.setPolytope(*trialPolytope, *testPolytope);
            const auto& trialDOFs = trialFES.getDOFs(d, tr);
            for (size_t l = 0; l < static_cast<size_t>(testDOFs.size()); ++l)
              for (size_t m = 0; m < static_cast<size_t>(trialDOFs.size()); ++m)
                checksum += integral.integrate(m, l);
          }
        }
        benchmark::DoNotOptimize(checksum);
      }

      state.SetItemsProcessed(state.iterations() * entries);
      state.counters["cell_pairs"] = static_cast<double>(cells * cells);
      state.counters["entries"] = static_cast<double>(entries);
    }

    template <class Form>
    void runAssembly(benchmark::State& state, Form& form, Index cells)
    {
      form.assemble();
      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
      state.SetItemsProcessed(state.iterations() * cells);
      state.counters["cells"] = static_cast<double>(cells);
      state.counters["nonzeros"] = static_cast<double>(form.getOperator().nonZeros());
    }

    template <class Form>
    void runLinearAssembly(benchmark::State& state, Form& form, Index cells)
    {
      form.assemble();
      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getVector().norm());
      }
      state.SetItemsProcessed(state.iterations() * cells);
      state.counters["cells"] = static_cast<double>(cells);
      state.counters["dofs"] = static_cast<double>(form.getVector().size());
    }

    template <class Family, size_t Dimension>
    void BM_BasisLoadKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TestFunction v(fes);
      auto integral = Integral(v);
      runLinearKernel(state, integral, fes);
    }

    template <class Family, size_t Dimension>
    void BM_SourceLoadKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TestFunction v(fes);
      RealFunction f([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(f, v);
      runLinearKernel(state, integral, fes);
    }

    template <class Family, size_t Dimension>
    void BM_FluxLoadGenericKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TestFunction v(fes);
      if constexpr (Dimension == 2)
      {
        VectorFunction flux{1.0, 0.5};
        auto integral = Integral(flux, Grad(v));
        runLinearKernel(state, integral, fes);
      }
      else
      {
        VectorFunction flux{1.0, 0.5, 0.25};
        auto integral = Integral(flux, Grad(v));
        runLinearKernel(state, integral, fes);
      }
    }

    template <class Family, size_t Dimension>
    void BM_MassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_WeightedMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      RealFunction c([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(c * u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_OuterWeightedMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      RealFunction c([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(c * Dot(u, v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_GridFunctionWeightedMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      GridFunction coefficient(fes);
      coefficient = RealFunction([](const Point& p) { return 1.0 + p.x(); });
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(coefficient * u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_GradGradKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Grad(u), Grad(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_WeightedGradGradKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      RealFunction c([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(c * Grad(u), Grad(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_VectorMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_VectorSourceLoadKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TestFunction v(fes);
      if constexpr (Dimension == 2)
      {
        VectorFunction f{1.0, 0.5};
        auto integral = Integral(f, v);
        runLinearKernel(state, integral, fes);
      }
      else
      {
        VectorFunction f{1.0, 0.5, 0.25};
        auto integral = Integral(f, v);
        runLinearKernel(state, integral, fes);
      }
    }

    template <class Family, size_t Dimension>
    void BM_AnisotropicMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      Math::Matrix<Real> matrix = Math::Matrix<Real>::Identity(Dimension, Dimension);
      matrix(0, Dimension - 1) = 0.125;
      MatrixFunction coefficient(matrix);
      auto integral = Integral(coefficient * u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_JacobianJacobianKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Jacobian(u), Jacobian(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_WeightedJacobianJacobianKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      RealFunction c([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(c * Jacobian(u), Jacobian(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_AdvectionKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      if constexpr (Dimension == 2)
      {
        VectorFunction beta{1.0, 0.5};
        auto integral = Integral(Dot(Jacobian(u) * beta, v));
        runBilinearKernel(state, integral, fes, fes);
      }
      else
      {
        VectorFunction beta{1.0, 0.5, 0.25};
        auto integral = Integral(Dot(Jacobian(u) * beta, v));
        runBilinearKernel(state, integral, fes, fes);
      }
    }

    template <class Family, size_t Dimension>
    void BM_GridFunctionAdvectionKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      GridFunction beta(fes);
      if constexpr (Dimension == 2)
        beta = VectorFunction{1.0, 0.5};
      else
        beta = VectorFunction{1.0, 0.5, 0.25};
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Dot(Jacobian(u) * beta, v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_DivPressureKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto velocity = Family::vector(mesh, Dimension);
      auto pressure = Family::scalar(mesh);
      TrialFunction u(velocity);
      TestFunction q(pressure);
      auto integral = Integral(Div(u), q);
      runBilinearKernel(state, integral, velocity, pressure);
    }

    template <class Family, size_t Dimension>
    void BM_PressureDivKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto pressure = Family::scalar(mesh);
      auto velocity = Family::vector(mesh, Dimension);
      TrialFunction p(pressure);
      TestFunction v(velocity);
      auto integral = Integral(p, Div(v));
      runBilinearKernel(state, integral, pressure, velocity);
    }

    template <class Family, size_t Dimension>
    void BM_DivDivGenericKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Div(u), Div(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_ElasticityGenericKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral =
        Integral(Jacobian(u) + Jacobian(u).T(), 0.5 * (Jacobian(v) + Jacobian(v).T()));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <size_t Dimension>
    void BM_P1PotentialKernel(benchmark::State& state)
    {
      auto mesh = Dimension == 2
        ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5})
        : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
      P1<Real, LocalMesh> fes(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      const auto kernel = [](const Point&, const Point&) { return 1.0; };
      auto integral = Integral(Potential(kernel, u), v);
      runGlobalBilinearKernel(state, integral, fes, fes);
    }

    template <size_t Dimension>
    void BM_P1VectorPotentialKernel(benchmark::State& state)
    {
      auto mesh = Dimension == 2
        ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5})
        : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
      P1<Math::SpatialVector<Real>, LocalMesh> fes(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      const auto kernel = [](Math::SpatialMatrix<Real>& value, const Point&,
                            const Point&) { value.setIdentity(); };
      auto integral = Integral(Potential(kernel, u), v);
      runGlobalBilinearKernel(state, integral, fes, fes);
    }

    template <size_t TrialOrder, size_t TestOrder, size_t Dimension>
    void BM_MixedOrderMassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto trialFES = H1Family<TrialOrder>::scalar(mesh);
      auto testFES = H1Family<TestOrder>::scalar(mesh);
      TrialFunction u(trialFES);
      TestFunction v(testFES);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, trialFES, testFES);
    }

    template <size_t TrialOrder, size_t TestOrder, size_t Dimension>
    void BM_MixedOrderGradGradKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto trialFES = H1Family<TrialOrder>::scalar(mesh);
      auto testFES = H1Family<TestOrder>::scalar(mesh);
      TrialFunction u(trialFES);
      TestFunction v(testFES);
      auto integral = Integral(Grad(u), Grad(v));
      runBilinearKernel(state, integral, trialFES, testFES);
    }

    void BM_P1GeometryMassKernel(benchmark::State& state)
    {
      const auto geometry = static_cast<Polytope::Type>(state.range(0));
      auto mesh = makeGeometryMesh(geometry);
      P1<Real, LocalMesh> fes(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, fes, fes);
      state.counters["geometry"] = static_cast<double>(state.range(0));
    }

    void BM_P1GeometryGradGradKernel(benchmark::State& state)
    {
      const auto geometry = static_cast<Polytope::Type>(state.range(0));
      auto mesh = makeGeometryMesh(geometry);
      P1<Real, LocalMesh> fes(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Grad(u), Grad(v));
      runBilinearKernel(state, integral, fes, fes);
      state.counters["geometry"] = static_cast<double>(state.range(0));
    }

    template <size_t Dimension>
    void BM_CurvedH1P2MassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      installP2Transformations<Dimension>(mesh);
      auto fes = H1Family<2>::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <size_t Dimension>
    void BM_CurvedH1P2GradGradKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      installP2Transformations<Dimension>(mesh);
      auto fes = H1Family<2>::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(Grad(u), Grad(v));
      runBilinearKernel(state, integral, fes, fes);
    }

    template <class Family, size_t Dimension>
    void BM_MassAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(u, v);
      runAssembly(state, form, mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_GradGradAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(Grad(u), Grad(v));
      runAssembly(state, form, mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_ReactionDiffusionAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(Grad(u), Grad(v)) + Integral(u, v);
      runAssembly(state, form, mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_ColdMassAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      for (auto _ : state)
      {
        BilinearForm form(u, v);
        form = Integral(u, v);
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
      state.SetItemsProcessed(state.iterations() * mesh.getCellCount());
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_ColdGradGradAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      for (auto _ : state)
      {
        BilinearForm form(u, v);
        form = Integral(Grad(u), Grad(v));
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
      state.SetItemsProcessed(state.iterations() * mesh.getCellCount());
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_ElasticityAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::vector(mesh, Dimension);
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(Div(u), Div(v)) +
        Integral(Jacobian(u) + Jacobian(u).T(), 0.5 * (Jacobian(v) + Jacobian(v).T()));
      runAssembly(state, form, mesh.getCellCount());
    }

    template <class Family, size_t Dimension>
    void BM_BoundaryMassAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = BoundaryIntegral(u, v);
      runAssembly(state, form, mesh.getFaceCount());
    }

    template <class Family, size_t Dimension>
    void BM_BoundaryLoadAssembly(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      auto fes = Family::scalar(mesh);
      TestFunction v(fes);
      RealFunction f([](const Point& p) { return 1.0 + p.x(); });
      LinearForm form(v);
      form = BoundaryIntegral(f, v);
      runLinearAssembly(state, form, mesh.getFaceCount());
    }

    template <size_t Dimension>
    void BM_P0MassKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      P0<Real, LocalMesh> fes(mesh);
      TrialFunction u(fes);
      TestFunction v(fes);
      auto integral = Integral(u, v);
      runBilinearKernel(state, integral, fes, fes);
    }

    template <size_t Dimension>
    void BM_P0SourceLoadKernel(benchmark::State& state)
    {
      auto mesh = makeMesh<Dimension>();
      P0<Real, LocalMesh> fes(mesh);
      TestFunction v(fes);
      RealFunction f([](const Point& p) { return 1.0 + p.x(); });
      auto integral = Integral(f, v);
      runLinearKernel(state, integral, fes);
    }

/**
 * @brief Registers the local-kernel cases of one family and dimension.
 * @param FAMILY Finite element family providing scalar() and vector()
 * @param LABEL Name the family appears under in the benchmark output
 * @param DIMENSION Spatial dimension of the mesh
 *
 * Covers every expression a QuadratureRule specializes, and the generic
 * controls alongside them, so a regression cannot hide behind the optimized
 * cases alone.
 */
#define RODIN_REGISTER_KERNELS(FAMILY, LABEL, DIMENSION)                                 \
  BENCHMARK_TEMPLATE(BM_BasisLoadKernel, FAMILY, DIMENSION)                              \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/BasisLoad");                     \
  BENCHMARK_TEMPLATE(BM_SourceLoadKernel, FAMILY, DIMENSION)                             \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/SourceLoad");                    \
  BENCHMARK_TEMPLATE(BM_FluxLoadGenericKernel, FAMILY, DIMENSION)                        \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GenericFluxLoad");               \
  BENCHMARK_TEMPLATE(BM_MassKernel, FAMILY, DIMENSION)                                   \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/Mass");                          \
  BENCHMARK_TEMPLATE(BM_WeightedMassKernel, FAMILY, DIMENSION)                           \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/WeightedMass");                  \
  BENCHMARK_TEMPLATE(BM_OuterWeightedMassKernel, FAMILY, DIMENSION)                      \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/OuterWeightedMass");             \
  BENCHMARK_TEMPLATE(BM_GridFunctionWeightedMassKernel, FAMILY, DIMENSION)               \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GridFunctionWeightedMass");      \
  BENCHMARK_TEMPLATE(BM_GradGradKernel, FAMILY, DIMENSION)                               \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GradGrad");                      \
  BENCHMARK_TEMPLATE(BM_WeightedGradGradKernel, FAMILY, DIMENSION)                       \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/WeightedGradGrad");              \
  BENCHMARK_TEMPLATE(BM_VectorMassKernel, FAMILY, DIMENSION)                             \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/VectorMass");                    \
  BENCHMARK_TEMPLATE(BM_VectorSourceLoadKernel, FAMILY, DIMENSION)                       \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/VectorSourceLoad");              \
  BENCHMARK_TEMPLATE(BM_AnisotropicMassKernel, FAMILY, DIMENSION)                        \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/AnisotropicMass");               \
  BENCHMARK_TEMPLATE(BM_JacobianJacobianKernel, FAMILY, DIMENSION)                       \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/JacobianJacobian");              \
  BENCHMARK_TEMPLATE(BM_WeightedJacobianJacobianKernel, FAMILY, DIMENSION)               \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/WeightedJacobianJacobian");      \
  BENCHMARK_TEMPLATE(BM_AdvectionKernel, FAMILY, DIMENSION)                              \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/Advection");                     \
  BENCHMARK_TEMPLATE(BM_GridFunctionAdvectionKernel, FAMILY, DIMENSION)                  \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GridFunctionAdvection");         \
  BENCHMARK_TEMPLATE(BM_DivPressureKernel, FAMILY, DIMENSION)                            \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/DivPressure");                   \
  BENCHMARK_TEMPLATE(BM_PressureDivKernel, FAMILY, DIMENSION)                            \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/PressureDiv");                   \
  BENCHMARK_TEMPLATE(BM_DivDivGenericKernel, FAMILY, DIMENSION)                          \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GenericDivDiv");                 \
  BENCHMARK_TEMPLATE(BM_ElasticityGenericKernel, FAMILY, DIMENSION)                      \
    ->Name("Integrator/Kernel/" LABEL "/D" #DIMENSION "/GenericElasticity")

/**
 * @brief Registers the full-assembly cases of one family and dimension.
 * @param FAMILY Finite element family providing scalar() and vector()
 * @param LABEL Name the family appears under in the benchmark output
 * @param DIMENSION Spatial dimension of the mesh
 *
 * A subset of the expressions, since assembly adds mesh traversal and sparse
 * construction to every case and the matrix would otherwise grow past what a
 * runner can afford.
 */
#define RODIN_REGISTER_ASSEMBLY(FAMILY, LABEL, DIMENSION)                                \
  BENCHMARK_TEMPLATE(BM_MassAssembly, FAMILY, DIMENSION)                                 \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/ReassemblyMass");              \
  BENCHMARK_TEMPLATE(BM_GradGradAssembly, FAMILY, DIMENSION)                             \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/ReassemblyGradGrad");          \
  BENCHMARK_TEMPLATE(BM_ReactionDiffusionAssembly, FAMILY, DIMENSION)                    \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/ReactionDiffusion");           \
  BENCHMARK_TEMPLATE(BM_ColdMassAssembly, FAMILY, DIMENSION)                             \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/ColdMass");                    \
  BENCHMARK_TEMPLATE(BM_ColdGradGradAssembly, FAMILY, DIMENSION)                         \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/ColdGradGrad");                \
  BENCHMARK_TEMPLATE(BM_ElasticityAssembly, FAMILY, DIMENSION)                           \
    ->Name("Integrator/Assembly/" LABEL "/D" #DIMENSION "/Elasticity")

    // These run on CI, so the matrix is kept to what distinguishes the
    // integrators rather than to what the machine can be made to do. Meshes are
    // small deliberately: every case reports its cell and entry counts, so cost
    // per element is what is being compared, and a larger mesh multiplies the
    // bill without sharpening that comparison. Cubic elements in three
    // dimensions are registered for assembly only; as a local kernel one such
    // case took over a second an iteration against milliseconds for the rest,
    // which buys a number nobody reads at a price paid on every run.
    RODIN_REGISTER_KERNELS(P1Family, "P1", 2);
    RODIN_REGISTER_KERNELS(H1Family<1>, "H1P1", 2);
    RODIN_REGISTER_KERNELS(H1Family<2>, "H1P2", 2);
    RODIN_REGISTER_KERNELS(H1Family<3>, "H1P3", 2);
    RODIN_REGISTER_KERNELS(P1Family, "P1", 3);
    RODIN_REGISTER_KERNELS(H1Family<2>, "H1P2", 3);

    RODIN_REGISTER_ASSEMBLY(P1Family, "P1", 2);
    RODIN_REGISTER_ASSEMBLY(H1Family<2>, "H1P2", 2);
    RODIN_REGISTER_ASSEMBLY(H1Family<3>, "H1P3", 2);
    RODIN_REGISTER_ASSEMBLY(P1Family, "P1", 3);
    RODIN_REGISTER_ASSEMBLY(H1Family<2>, "H1P2", 3);

    BENCHMARK_TEMPLATE(BM_BoundaryMassAssembly, P1Family, 2)
      ->Name("Integrator/Assembly/P1/D2/BoundaryMass");
    BENCHMARK_TEMPLATE(BM_BoundaryLoadAssembly, P1Family, 2)
      ->Name("Integrator/Assembly/P1/D2/BoundaryLoad");
    BENCHMARK_TEMPLATE(BM_BoundaryMassAssembly, H1Family<2>, 2)
      ->Name("Integrator/Assembly/H1P2/D2/BoundaryMass");
    BENCHMARK_TEMPLATE(BM_BoundaryLoadAssembly, H1Family<2>, 2)
      ->Name("Integrator/Assembly/H1P2/D2/BoundaryLoad");
    BENCHMARK_TEMPLATE(BM_BoundaryMassAssembly, P1Family, 3)
      ->Name("Integrator/Assembly/P1/D3/BoundaryMass");
    BENCHMARK_TEMPLATE(BM_BoundaryLoadAssembly, P1Family, 3)
      ->Name("Integrator/Assembly/P1/D3/BoundaryLoad");

    BENCHMARK_TEMPLATE(BM_P0MassKernel, 2)->Name("Integrator/Kernel/P0/D2/Mass");
    BENCHMARK_TEMPLATE(BM_P0SourceLoadKernel, 2)
      ->Name("Integrator/Kernel/P0/D2/SourceLoad");
    BENCHMARK_TEMPLATE(BM_P0MassKernel, 3)->Name("Integrator/Kernel/P0/D3/Mass");
    BENCHMARK_TEMPLATE(BM_P0SourceLoadKernel, 3)
      ->Name("Integrator/Kernel/P0/D3/SourceLoad");

    BENCHMARK_TEMPLATE(BM_P1PotentialKernel, 2)
      ->Name("Integrator/Kernel/P1/D2/GlobalPotential");
    BENCHMARK_TEMPLATE(BM_P1PotentialKernel, 3)
      ->Name("Integrator/Kernel/P1/D3/GlobalPotential");
    BENCHMARK_TEMPLATE(BM_P1VectorPotentialKernel, 2)
      ->Name("Integrator/Kernel/P1/D2/VectorGlobalPotential");
    BENCHMARK_TEMPLATE(BM_P1VectorPotentialKernel, 3)
      ->Name("Integrator/Kernel/P1/D3/VectorGlobalPotential");

    BENCHMARK_TEMPLATE(BM_MixedOrderMassKernel, 2, 1, 2)
      ->Name("Integrator/Kernel/H1P2-P1/D2/Mass");
    BENCHMARK_TEMPLATE(BM_MixedOrderGradGradKernel, 2, 1, 2)
      ->Name("Integrator/Kernel/H1P2-P1/D2/GradGrad");
    BENCHMARK_TEMPLATE(BM_MixedOrderMassKernel, 3, 2, 2)
      ->Name("Integrator/Kernel/H1P3-P2/D2/Mass");
    BENCHMARK_TEMPLATE(BM_MixedOrderGradGradKernel, 3, 2, 2)
      ->Name("Integrator/Kernel/H1P3-P2/D2/GradGrad");
    BENCHMARK_TEMPLATE(BM_MixedOrderMassKernel, 2, 1, 3)
      ->Name("Integrator/Kernel/H1P2-P1/D3/Mass");
    BENCHMARK_TEMPLATE(BM_MixedOrderGradGradKernel, 2, 1, 3)
      ->Name("Integrator/Kernel/H1P2-P1/D3/GradGrad");

    BENCHMARK(BM_P1GeometryMassKernel)
      ->Name("Integrator/Kernel/P1/Geometry/Mass")
      ->Arg(static_cast<int64_t>(Polytope::Type::Segment))
      ->Arg(static_cast<int64_t>(Polytope::Type::Triangle))
      ->Arg(static_cast<int64_t>(Polytope::Type::Quadrilateral))
      ->Arg(static_cast<int64_t>(Polytope::Type::Tetrahedron))
      ->Arg(static_cast<int64_t>(Polytope::Type::Hexahedron))
      ->Arg(static_cast<int64_t>(Polytope::Type::Wedge))
      ->Arg(static_cast<int64_t>(Polytope::Type::Pyramid));
    BENCHMARK(BM_P1GeometryGradGradKernel)
      ->Name("Integrator/Kernel/P1/Geometry/GradGrad")
      ->Arg(static_cast<int64_t>(Polytope::Type::Segment))
      ->Arg(static_cast<int64_t>(Polytope::Type::Triangle))
      ->Arg(static_cast<int64_t>(Polytope::Type::Quadrilateral))
      ->Arg(static_cast<int64_t>(Polytope::Type::Tetrahedron))
      ->Arg(static_cast<int64_t>(Polytope::Type::Hexahedron))
      ->Arg(static_cast<int64_t>(Polytope::Type::Wedge))
      ->Arg(static_cast<int64_t>(Polytope::Type::Pyramid));

    BENCHMARK_TEMPLATE(BM_CurvedH1P2MassKernel, 2)
      ->Name("Integrator/Kernel/H1P2/D2/CurvedMass");
    BENCHMARK_TEMPLATE(BM_CurvedH1P2GradGradKernel, 2)
      ->Name("Integrator/Kernel/H1P2/D2/CurvedGradGrad");
    BENCHMARK_TEMPLATE(BM_CurvedH1P2MassKernel, 3)
      ->Name("Integrator/Kernel/H1P2/D3/CurvedMass");
    BENCHMARK_TEMPLATE(BM_CurvedH1P2GradGradKernel, 3)
      ->Name("Integrator/Kernel/H1P2/D3/CurvedGradGrad");

#undef RODIN_REGISTER_ASSEMBLY
#undef RODIN_REGISTER_KERNELS
  }
}
