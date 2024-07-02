#include <thread>
#include <queue>
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

class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        condition.wait(lock, [this] { return stop || !tasks.empty(); });

                        if (stop && tasks.empty()) {
                            return;
                        }

                        task = std::move(tasks.front());
                        tasks.pop();
                    }

                    task();
                }
            });
        }
    }

    template <typename Func>
    void enqueue(Func func) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.emplace(func);
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers) {
            worker.join();
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
};

struct Experiment
{
  const double azimuth;
  const double hmax;
  const double c;
  const double T;
  const double hausd;
};

std::vector<Experiment> experiments = {};

void run(size_t expId, const std::vector<Experiment>& experiments);

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

  const size_t hwc = std::thread::hardware_concurrency();
  const size_t n = hwc - 2;
  ThreadPool threadPool(n);
  for (size_t i = 0; i < experiments.size(); i++)
    threadPool.enqueue([i]{ run(i, experiments); });

  return 0;
}

void run(size_t experimentId, const std::vector<Experiment>& experiments)
{
  Experiment experiment = experiments[experimentId];

  const bool physical = true;
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
    filename = "physical/L2ErrorPhysical_" + std::to_string(experimentId) + ".csv";
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
      dist = MMG::Distancer(vh).setInteriorDomain(sphereCap)
                               .distance(th);

      // Compute gradient of signed distance function
      GridFunction conormal(uh);
      auto gd = Grad(dist);
      conormal = gd / Frobenius(gd);

      // Advect
      try
      {
        MMG::Advect(dist, conormal).step(std::min(dt, experiment.T - t));
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning()
          << "Advection failed for experiment " << experimentId
          << Alert::NewLine
          << "Retrying..."
          << Alert::Raise;
        continue;
      }

      t += std::min(dt, experiment.T - t);
    }

    auto phit = RealFunction(
        [&](const Point& p) { return phi(t, p, experiment.azimuth); });

    // Compute L2 error
    GridFunction diff(vh);
    diff = Pow(dist - phit, 2);
    diff.setWeights();
    double error = Integral(diff).compute();
    fout << t << "," << error << '\n' << std::flush;

    // Generate mesh to subdomain
    th = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                    .setHMax(experiment.hmax)
                                    .setHausdorff(experiment.hausd)
                                    .discretize(dist);

    MMG::Optimizer().setAngleDetection(false)
                   .setHMax(experiment.hmax)
                   .setHausdorff(experiment.hausd)
                   .optimize(th);

    // th.save("out/SphereCap.mfem." + std::to_string(i) + ".mesh");

    if (t + std::numeric_limits<double>::epsilon() > experiment.T)
      break;

    // Alert::Info() << "l2: " << error << '\n' << Alert::Raise;
    i++;
  }
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

