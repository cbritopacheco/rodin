#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char** argv)
{
  Mesh mesh =
    Mesh<Rodin::Context::Sequential>::Builder()
    .initialize(2)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Type::Triangle, {0, 1, 2})
    .polytope(Polytope::Type::Triangle, {1, 3, 2})
    .finalize();

  P1 fes(mesh);
  GridFunction gf(fes);

  ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
  gf.project(c);

  mesh.save("miaow.mesh");
  gf.save("miaow.gf");
}

