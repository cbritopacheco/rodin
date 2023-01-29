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
constexpr double azimuth = 0.3;
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
    dist.getFiniteElementSpace().getMesh().save("miaow.mesh");
    dist.save("miaow.gf");

    GridFunction diff(vh);
    auto phit = ScalarFunction([&](const Point& p) -> double { return phi(t, p); });
    diff = phit;
    diff.save("diff.gf");
    std::exit(1);

    diff = Pow(dist - phit, 2);
    double error = Integral(diff);

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(0.05)
                                    .setHausdorff(0.01)
                                    .discretize(dist);

    // Save results
    th.save("out/SphereCap." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    Alert::Info() << "l2: " << error << '\n' << Alert::Raise;
    fout << error << '\n' << std::flush;
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
  const double dom =
    std::sin(alpha) * std::sin(azimuth + t) + std::cos(alpha) * std::cos(azimuth + t);
  assert(-1.0 <= dom && dom <= 1.0);
  const double res = std::acos(dom);
  assert(std::isfinite(res));
  return -Math::sgn(azimuth + t - alpha) * res;
}

