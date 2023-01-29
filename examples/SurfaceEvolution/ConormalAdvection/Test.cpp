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

constexpr double T = M_PI;
constexpr double dt = 0.01;
constexpr size_t maxIt = T / dt;
constexpr double azimuth = 0.1;
constexpr Geometry::Attribute sphereCap = 3;
constexpr char meshFile[] =
  "../resources/examples/SurfaceEvolution/ConormalAdvection/SphereCap.medit.mesh";

double phi(double t, const Point& p);

int main()
{
  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  std::ofstream fout("err.txt");
  double t = 0;
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "--------------\n"
                  << "i: " << i << '\n'
                  << "t: " << t << '\n'
                  << Alert::Raise;

    // Build finite element space on the mesh
    H1 vh(th);
    H1 uh(th, th.getSpaceDimension());

    auto phit = ScalarFunction([&](const Point& p) -> double { return phi(t, p); });

    // Distance the subdomain
    GridFunction dist(vh);

    if (i == 0)
    {
      dist = phit;
    }
    else
    {
      dist = MMG::Distancer(vh).setInteriorDomain(sphereCap)
                               .distance(th);
    }

    // Compute L2 error
    GridFunction diff(vh);
    diff = Pow(dist - phit, 2);
    double error = Integral(diff).compute();
    diff.save("diff.gf");

    // Compute gradient of signed distance function
    Grad gd(dist);
    GridFunction conormal(uh);
    conormal = gd / Frobenius(gd);
    conormal.save("conormal.gf");

    // Advect
    MMG::Advect(dist, conormal).step(dt);
    dist.getFiniteElementSpace().getMesh().save("miaow.mesh");
    dist.save("miaow.gf");

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(0.02)
                                    .setHausdorff(0.01)
                                    .discretize(dist);

    // Save results
    th.save("out/SphereCap." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    Alert::Info() << "l2: " << error << '\n' << Alert::Raise;
    fout << error << '\n' << std::flush;
    t += dt;
  }
  return 0;
}

double phi(double t, const Point& p)
{
  const double x = p.x();
  const double y = p.y();
  const double z = p.z();
  assert(x != 0 && y != 0);
  assert(x != 0 && y != 0 && z != 0);
  const double alpha = std::acos(z);
  assert(0 <= alpha && alpha <= M_PI);
  const double beta = std::atan2(y, x) + M_PI;
  assert(0 <= beta && beta <= 2 * M_PI);
  return alpha - t - azimuth;
}

