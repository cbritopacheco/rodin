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

std::vector<Experiment> experiments = {};

double phi(double t, const Point& p, const double azimuth);

int main(int argc, char** argv)
{
  size_t N = 32;
  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.1, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.2, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.3, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.4, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.5, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.6, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.7, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.8, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 0.9, M_PI / 2 - 0.1, 0.01});

  for (size_t i = 0; i < N; i++)
    experiments.push_back(
        {0.1, 0.02 + (1 - 0.02) * float(i) / float(N), 1.0, M_PI / 2 - 0.1, 0.01});

  if (argc == 1)
    Alert::Exception() << "Experiment id not specified." << Alert::Raise;

  const size_t experimentId = std::atoi(argv[1]);

  if (experimentId <= 0 || experimentId > experiments.size())
    std::exit(EXIT_FAILURE);

  Experiment experiment = experiments[experimentId - 1];

  const bool physical = false;
  double dt = NAN;
  if (physical)
    dt = experiment.c * experiment.hmax;
  else
    dt = experiment.c;

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  std::string filename;
  if (physical)
    filename = "L2ErrorPhysical_" + std::to_string(experimentId) + ".csv";
  else
    filename = "L2Error_" + std::to_string(experimentId) + ".csv";

  std::ofstream fout(filename);
  fout << "t,$\\mathcal{E}(t)$\n" << std::flush;
  double t = 0;
  size_t i = 0;
  while (true)
  {
    // Build finite element space on the mesh
    P1 vh(th);
    P1 uh(th, th.getSpaceDimension());

    // Distance the subdomain
    GridFunction dist(vh);

    if (i == 0)
    {
      dist = [&](const Point& p) { return phi(0, p, experiment.azimuth); };
    }
    else
    {
      dist = MMG::Distance(vh).setInteriorDomain(sphereCap)
                               .distance(th);

      // Compute gradient of signed distance function
      Grad gd(dist);
      GridFunction conormal(uh);
      conormal = gd / Frobenius(gd);

      // Advect
      MMG::Advect(dist, conormal).step(std::min(dt, experiment.T - t));
      t += std::min(dt, experiment.T - t);
    }

    auto phit = ScalarFunction(
        [&](const Point& p) { return phi(t, p, experiment.azimuth); });

    // Compute L2 error
    GridFunction diff(vh);
    diff = Pow(dist - phit, 2);
    double error = Integral(diff).compute();
    fout << t << "," << error << '\n' << std::flush;

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(experiment.hmax)
                                    .setHausdorff(experiment.hausd)
                                    .discretize(dist);

    th.save("out/SphereCap.mfem." + std::to_string(i) + ".mesh");

    if (t + std::numeric_limits<double>::epsilon() > experiment.T)
      break;

    // Alert::Info() << "l2: " << error << '\n' << Alert::Raise;
    i++;
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

