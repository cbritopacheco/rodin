/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>

#include <boost/filesystem.hpp>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace RodinBenchmark
{
  struct Poisson_UniformGrid_16x16 : public benchmark::Fixture
  {
    public:
      static constexpr const Geometry::Attribute dirichletAttr = 1;

      void SetUp(const benchmark::State&)
      {
        mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
        mesh.getConnectivity().compute(1, 2);
        vhPtr.reset(new P1<Scalar, Context::Serial>(mesh));
      }

      void TearDown(const benchmark::State&)
      {}

      boost::filesystem::path meshfile;
      Mesh<Context::Serial> mesh;
      std::unique_ptr<P1<Scalar, Context::Serial>> vhPtr;
  };

  BENCHMARK_F(Poisson_UniformGrid_16x16, Assembly_ConstantCoefficient_ConstantSource)
  (benchmark::State& st)
  {
    assert(vhPtr);
    const auto& vh = *vhPtr;
    TrialFunction u(vh);
    TestFunction  v(vh);
    ScalarFunction f(1.0);
    ScalarFunction zero(0.0);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, zero).on(dirichletAttr);

    for (auto _ : st)
      poisson.assemble();
  }
}