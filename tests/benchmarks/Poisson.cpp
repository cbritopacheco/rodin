/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <Corrade/Utility/Debug.h>
#include <Corrade/Utility/Resource.h>
#include <Corrade/TestSuite/Tester.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace RodinTest
{
  struct Poisson : Corrade::TestSuite::Tester
  {
    using FES = H1<Scalar, Context::Serial>;

    static constexpr const size_t batchSize = 100;
    static constexpr const size_t batchCount = 100;
    static constexpr const Geometry::Attribute dirichletAttr = 1;

    explicit Poisson();

    void Assembly_ConstantCoefficient_ConstantSource();

    Corrade::Utility::Resource      rs;
    std::string                     meshfile;
    const Mesh<Context::Serial>     mesh;
    FES                             vh;
    TrialFunction<FES>              u;
    TestFunction<FES>               v;
    ScalarFunction<Scalar>          f;
    ScalarFunction<Scalar>          g;
  };

  Poisson::Poisson()
    : rs("rodin"),
      mesh(meshfile),
      vh(mesh),
      u(vh),
      v(vh),
      f(1.0),
      g(0.0)

  {
    addBenchmarks({ &Poisson::Assembly_ConstantCoefficient_ConstantSource }, batchCount);
  }

  void Poisson::Assembly_ConstantCoefficient_ConstantSource()
  {
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, g).on(dirichletAttr);
    poisson.assemble();
  }
}

CORRADE_TEST_MAIN(RodinTest::Poisson)

