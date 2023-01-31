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

constexpr Geometry::Attribute sphereCap = 3;
constexpr char meshFile[] =
  "../resources/examples/SurfaceEvolution/ConormalAdvection/SphereCap.medit.mesh";

struct Experiment
{
  const double azimuth;
  const double hmax;
  const double c;
  const double T;
  const double hausd;
};

std::vector<Experiment> experiments =
{
  {0.1, 1, 0.1, M_PI, 0.01},
  {0.1, 1, 0.01, M_PI, 0.01},
  {0.1, 0.1, 0.1, M_PI, 0.01},
  {0.1, 0.1, 0.01, M_PI, 0.01},
  {0.1, 0.02, 0.1, M_PI, 0.01},
  {0.1, 0.02, 0.01, M_PI, 0.01}
};

double phi(double t, const Point& p, const double azimuth);

int main(int argc, char** argv)
{
  const size_t experimentId = std::atoi(argv[1]);
  Experiment experiment = experiments[experimentId - 1];

  const double dt = experiment.c * experiment.hmax;
  const size_t maxIt = (experiment.T - experiment.azimuth) / dt;

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  std::ofstream fout("L2Error_" + std::to_string(experimentId) + ".csv");
  fout << "t,$\\mathcal{E}(t)$\n" << std::flush;
  double t = 0;
  for (size_t i = 0; i < maxIt; i++)
  {
    // Alert::Info() << "--------------\n"
    //               << "i: " << i << " / " << maxIt << '\n'
    //               << "t: " << t << '\n'
    //               << Alert::Raise;

    // Build finite element space on the mesh
    H1 vh(th);
    H1 uh(th, th.getSpaceDimension());

    auto phit = ScalarFunction(
        [&](const Point& p) -> double
        { return phi(t, p, experiment.azimuth); });

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
    fout << t << "," << error << '\n' << std::flush;

    // Compute gradient of signed distance function
    Grad gd(dist);
    GridFunction conormal(uh);
    conormal = gd / Frobenius(gd);

    // Advect
    MMG::Advect(dist, conormal).step(dt);

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(experiment.hmax)
                                    .setHausdorff(experiment.hausd)
                                    .discretize(dist);

    // th.save("out/SphereCap." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    // Alert::Info() << "l2: " << error << '\n' << Alert::Raise;
    t += dt;
  }
  return 0;
}

double phi(double t, const Point& p, const double azimuth)
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

