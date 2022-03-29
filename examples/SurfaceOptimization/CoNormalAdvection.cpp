#include <Rodin/Mesh.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main(int, char**)
{
  const char* meshFile = "rodin.mesh";
  const double pi = std::atan(1) * 4;

  Mesh Omega = Mesh::load(meshFile);

  // Sphere radius
  double r = 1;

  // Hole radius
  double rr = 0.2;

  // Hole centers on sphere
  std::vector<std::array<double, 3>> cs = {
    { r * sin(0) * cos(0), r * sin(0) * sin(0), r * cos(0) },
    { r * sin(pi / 2) * cos(pi), r * sin(pi / 2) * sin(pi), r * cos(pi / 2) },
    { r * sin(pi / 2) * cos(pi / 2), r * sin(pi / 2) * sin(pi / 2), r * cos(pi / 2) }
  };

  // Geodesic distance
  auto gd = [&](const double* x, const double* c, int)
            {
              return std::acos((x[0] * c[0] + x[1] * c[1] + x[2] * c[2]));
            };

  // Function for generating holes
  auto f = ScalarCoefficient(
      [&](const double* x, int dim)
      {
        double d = std::numeric_limits<double>::max();
        for (const auto& c : cs)
        {
          double dd = gd(x, c.data(), dim) - rr;
          d = std::min(d, dd);
        }
        return -d;
      });

  H1 Vh(Omega);
  GridFunction dist(Vh);
  dist = f;

  Omega.save("Omega.mesh");
  dist.save("dist.gf");

  auto mmgMesh = Cast(Omega).to<MMG::MeshS>();
  auto mmgLs = Cast(dist).to<MMG::IncompleteScalarSolutionS>().setMesh(mmgMesh);

  mmgMesh.save("mmg.mesh");
  mmgLs.save("mmg.sol");

  return 0;
}
