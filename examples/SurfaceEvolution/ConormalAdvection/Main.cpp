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

constexpr Real T = M_PI / 2 - 0.1;
constexpr Real dt = 1;
constexpr Real hmax = 1;
constexpr Real hmin = 0.1 * hmax;
constexpr size_t maxIt = 100;
constexpr Geometry::Attribute sphereCap = 3;
constexpr char meshFile[] =
  "../resources/examples/SurfaceEvolution/ConormalAdvection/SphereCap.medit.mesh";

int main()
{
  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  // Evolution
  size_t i = 0;
  Real t = 0;
  while (true)
  {
    if (i == 0)
      th.save("out/SphereCap.medit." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "t: " << t << '\n'
                  << "i: " << i
                  << Alert::Raise;

    // Build finite element space on the mesh
    P1 vh(th);
    P1 uh(th, th.getSpaceDimension());

    // Distance the subdomain
    GridFunction dist(vh);
    dist = MMG::Distancer(vh).setInteriorDomain(sphereCap)
                             .distance(th);

    // Compute gradient of signed distance function
    Grad gd(dist);
    GridFunction conormal(uh);
    conormal = gd / Frobenius(gd);

    // Advect
    MMG::Advect(dist, conormal).step(std::min(dt, T - t));
    t += std::min(dt, T - t);
    i++;

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMin(hmin)
                                    .setHMax(hmax)
                                    .setHausdorff(0.01)
                                    .discretize(dist);

    MMG::Optimizer().setAngleDetection(false)
                        .setHMin(hmin)
                        .setHMax(hmax)
                        .setHausdorff(0.01)
                        .optimize(th);

    // Save results
    th.save("out/SphereCap.medit." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    th.save("out/SphereCap.mfem." + std::to_string(i) + ".mesh", IO::FileFormat::MFEM);

    if (t + std::numeric_limits<double>::epsilon() > T)
      break;
  }
  return 0;
}
