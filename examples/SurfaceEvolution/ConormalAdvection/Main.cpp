#include <Rodin/Math.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace std;
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main()
{
  const char* meshFile = "miaow.mesh";

  Mesh Omega;
  Omega.load(meshFile);

  H1 Vh(Omega);

  GridFunction dist(Vh);
  // dist.projectOnBoundary(f);
  std::ofstream fout("obj.txt");

  size_t N = 100;
  for (size_t i = 0; i < N; i++)
  {
    // fout << err << "\n";
    fout.flush();
  }

  Omega = MMG::ImplicitDomainMesher()
                             .split(6, {3, 6})
                             .noSplit(2)
                             .setHMax(0.05)
                             .surface()
                             .discretize(dist);

  return 0;
}
