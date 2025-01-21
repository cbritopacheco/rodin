#include <gtest/gtest.h>

#include "Rodin/Variational/RelativeError.h"
#include "Rodin/Test/Random.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;

TEST(Rodin_Manufactured_P1_Triangle, Poisson)
{
  auto pi = Rodin::Math::Constants::pi();

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 64, 64 });
  mesh.scale(1.0 / 63);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);

  RandomFloat gen(0.0, 100.0);
  for (size_t i = 0; i < 25; i++)
  {
    RealFunction c = gen();

    auto f = 2 * c * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto error = RelativeError::lInf(u.getSolution(), c * sin(pi * F::x) * sin(pi * F::y));
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}


TEST(Rodin_Manufactured_P1_Quadrilateral, Poisson)
{
  auto pi = Rodin::Math::Constants::pi();

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });
  mesh.scale(1.0 / 31);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);

  RandomFloat gen(0.0, 100.0);
  for (size_t i = 0; i < 25; i++)
  {
    RealFunction c = gen();

    auto f = 2 * c * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto error = RelativeError::lInf(u.getSolution(), c * sin(pi * F::x) * sin(pi * F::y));
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}

