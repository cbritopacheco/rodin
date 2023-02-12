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

constexpr double dt = 0.01;
constexpr size_t maxIt = 100;
constexpr Geometry::Attribute sphereCap = 3;
constexpr char meshFile[] =
  "../resources/examples/SurfaceEvolution/ConormalAdvection/SphereCap.medit.mesh";

int main()
{
  Math::Tensor<3> t(2, 4, 7);
  for (int i = 0; i < t.size(); i++)
    t(i) = i;
  for (int i = 0; i < t.size(); i++)
    std::cout << *(t.data() + i) << " ";
  std::cout << std::endl;

  std::cout << "cout\n";
  std::cout << t << std::endl;

  std::cout << "mat" << std::endl;

  Math::Matrix c = Math::slice(t, 1, 1);
  std::cout << c << std::endl;
  std::exit(1);

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  // Evolution
  for (size_t i = 0; i < maxIt; i++)
  {
    // Build finite element space on the mesh
    H1 vh(th);
    H1 uh(th, th.getSpaceDimension());

    // Distance the subdomain
    GridFunction dist(vh);
    dist = MMG::Distancer(vh).setInteriorDomain(sphereCap)
                             .distance(th);

    // Compute gradient of signed distance function
    Grad gd(dist);
    GridFunction conormal(uh);
    conormal = gd / Frobenius(gd);

    // Advect
    MMG::Advect(dist, conormal).step(dt);

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(0.05)
                                    .setHausdorff(0.01)
                                    .discretize(dist);

    // Save results
    th.save("out/SphereCap." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
  }
  return 0;
}
