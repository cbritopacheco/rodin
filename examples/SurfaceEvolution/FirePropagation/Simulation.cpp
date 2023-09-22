#include <cmath>
#include <chrono>

#include <Rodin/Math.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>
#include <Rodin/Variational/LazyEvaluator.h>

using namespace std;
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

static const char* meshfile =
  "../resources/examples/SurfaceEvolution/FirePropagation/Topography.mfem.mesh";
static const Scalar hmax = 400;
static Scalar hmin = 0.01 * hmax;
static const Scalar hausdorff = 20;
// static const Scalar hmax = 600;
// static const Scalar hmin = 100;

class Environment
{
  public:
    using Context = Context::Serial;
    using ScalarFES = P1<Scalar, Context>;
    using VectorFES = P1<Math::Vector, Context>;

    struct Plane
    {
      Scalar a, b, c, d;
    };

    enum Terrain
    {
      WorldBorder = 1,
      Vegetation = 2,
      Burnt = 3,
      Fire = 5
    };

    struct VegetalStratum
    {
      // Energy ratio between incident radiant energy and ignition energy of
      // dry fuel
      Scalar A0;

      // Vertical velocity
      Scalar u00;

      // Energy ratio between incident radiant energy emitted from the flame
      // base and ignition energy of dry fuel
      Scalar R00;

      // Moisture factor
      Scalar a;

      // Rate of spread for no slope and no wind (m/s)
      std::function<Scalar(const Point&)> R0;

      // Buoyancy velocity component for a zero slope (m/s)
      std::function<Scalar(const Point&)> u0;

      // Energy ratio between incident radiant energy and ignition energy of wet fuel
      std::function<Scalar(const Point&)> A;

      // ROS coefficient (m/s)
      std::function<Scalar(const Point&)> v0;

      // Moisture content (%)
      std::function<Scalar(const Point&)> m;

      // Thickness of the vegetal stratum (m)
      std::function<Scalar(const Point&)> e;

      // Surface density of vegetal fuel, (kg/m^2)
      std::function<Scalar(const Point&)> sigma;

      // Residence time (s)
      std::function<Scalar(const Point&)> tau;

      // Vegetal fuel surface to volume ratio (1/m)
      std::function<Scalar(const Point&)> sv;

      // Gas flame density (kg/m^3)
      std::function<Scalar(const Point&)> pv;

      // Absorption coefficient
      std::function<Scalar(const Point&)> v;

      // Flame gas temperature (K)
      std::function<Scalar(const Point&)> T;

      // Air temperature (K)
      std::function<Scalar(const Point&)> Ta;
    };

    class Flame
    {
      public:
        Flame(Environment& env)
          : m_env(env), m_direction(env.m_vfes)
        {}

        Flame& step(Scalar dt)
        {
          assert(dt > 0);
          const auto& env = m_env.get();

          // Slope vector
          auto p = Grad(env.m_terrainHeight);
          p.traceOf(Terrain::Burnt);

          // Angle between slope and ground
          auto alpha =
            [&](const Point& v)
            {
              Scalar angle = std::acos(p.z()(v) / Frobenius(p)(v));
              assert(angle >= 0);
              // assert(angle <= Math::Constants::pi<Scalar>() / 2.0);
              return Math::Constants::pi() / 2.0 - angle;
            };

          Grad gdist(env.m_fireDist);
          gdist.traceOf(Terrain::Vegetation);

          const auto conormal = gdist / Frobenius(gdist);

          // Angle between slope and conormal
          auto phi =
            [&](const Point& v) -> Scalar
            {
              Scalar fv = p(v).dot(conormal(v)) / p(v).norm();
              if (std::isfinite(fv))
                return std::acos(fv);
              else
                return 0.0;
            };

          // Angle between wind and conormal
          auto psi =
            [&](const Point& v) -> Scalar
            {
              Scalar fv =
                env.m_wind(v).dot(conormal(v)) / env.m_wind(v).norm();
              if (std::isfinite(fv))
                return std::acos(fv);
              else
                return 0.0;
            };

          // Compute the tilt angle
          auto gamma =
            [&](const Point& v) -> Scalar
            {
              Scalar rhs =
                std::tan(alpha(v)) * std::cos(phi(v)
                    ) + env.m_wind(v).norm() * std::cos(psi(v));
              return std::atan(rhs);
            };

          // Compute rate of spread
          auto R =
            [&](const Point& v) -> Scalar
            {
              const Scalar g = gamma(v);
              if (g > 0)
              {
                const Scalar R0 = env.m_vegetalStratum.R0(v);
                const Scalar v0 = env.m_vegetalStratum.v0(v);
                const Scalar A = env.m_vegetalStratum.A(v);
                const Scalar Ra = R0 + A * v0 * (
                    1 + std::sin(g) - std::cos(g)) / std::cos(g) - v0 / std::cos(g);
                return 0.5 * (Ra + std::sqrt(Ra * Ra + (4 * v0 * R0) / std::cos(g)));
              }
              else
              {
                return env.m_vegetalStratum.R0(v);
              }
            };

          m_direction = GridFunction(env.m_vfes);
          m_direction = dt * ScalarFunction(R) * conormal;

          return *this;
        }

        const GridFunction<VectorFES>& getDirection() const
        {
          return m_direction;
        }

      private:
        std::reference_wrapper<Environment> m_env;
        GridFunction<VectorFES> m_direction;
    };

    Environment(MMG::Mesh& topography, const VegetalStratum& vegetalStratum)
      : m_topography(topography),
        m_sfes(m_topography),
        m_vfes(m_topography, m_topography.getSpaceDimension()),
        m_wind(m_vfes),
        m_terrainHeight(m_sfes),
        m_vegetalStratum(vegetalStratum),
        m_flame(*this),
        m_gravity(-9.8),
        m_fireDist(m_sfes),
        m_elapsedTime(0.0)
      {
        m_terrainHeight = [](const Point& v) { return v.z(); };
      }

    Environment& step(Scalar dt)
    {
      std::cout << "wind\n";
      // auto wind = VectorFunction{0, 0, 0};
      // auto wind = VectorFunction{
      //     [](const Geometry::Point& p){ return p.y() - 25000; },
      //     [](const Geometry::Point& p){ return -p.x() + 25000; },
      //     0};

      // m_wind = wind / Frobenius(wind);

      std::cout << "dist\n";
      m_fireDist =
        MMG::Distancer(m_sfes).setInteriorDomain(Terrain::Burnt)
                              .distance(m_topography);

      std::cout << "flame\n";
      m_flame.step(1);

      std::cout << "advect: " << dt << "\n";
      MMG::Advect(m_fireDist, m_flame.getDirection()).step(dt);

      std::cout << "mmgls\n";
      m_topography = MMG::ImplicitDomainMesher().setHMax(hmax)
                                                .setHMin(hmin)
                                                .setHausdorff(hausdorff)
                                                .setAngleDetection(false)
                                                // .split(Terrain::Burnt,
                                                //     {Terrain::Burnt, Terrain::Vegetation})
                                                // .split(Terrain::Vegetation,
                                                //     {Terrain::Burnt, Terrain::Vegetation})
                                                .setBoundaryReference(Terrain::Fire)
                                                .discretize(m_fireDist);

      std::cout << "mmgoptim\n";
      MMG::Optimize().setAngleDetection(false)
                          .setHausdorff(hausdorff)
                          .setHMin(hmin)
                          .setHMax(hmax)
                          .optimize(m_topography);


      // Rebuild finite element spaces with new topography
      m_sfes = ScalarFES(m_topography);
      m_vfes = VectorFES(m_topography, m_topography.getSpaceDimension());

      std::cout << "terrain\n";
      m_terrainHeight = GridFunction(m_sfes);
      m_terrainHeight = [](const Point& v) { return v.z(); };
      m_wind = GridFunction(m_vfes);
      m_elapsedTime += dt;
      std::cout << "end\n";
      return *this;
    }

    const Flame& getFlame() const
    {
      return m_flame;
    }

    template <class FunctionDerived>
    Environment& setWind(const FunctionBase<FunctionDerived>& fn)
    {
      m_wind.project(fn);
      return *this;
    }

    const Mesh<Context>& getTopography() const
    {
      return m_topography;
    }

    Scalar getGravity() const
    {
      return m_gravity;
    }

  private:
    MMG::Mesh& m_topography;

    ScalarFES m_sfes;
    VectorFES m_vfes;

    GridFunction<VectorFES> m_wind;

    GridFunction<ScalarFES> m_terrainHeight;

    VegetalStratum m_vegetalStratum;

    Flame m_flame;

    const Scalar m_gravity;

    GridFunction<ScalarFES> m_fireDist;

    Scalar m_elapsedTime;
};


int main()
{
  MMG::Mesh topography;
  // topography.load("/Users/carlos/Projects/rodin/build/Topography.mesh");
  topography.load(meshfile);

  topography.getConnectivity().compute(1, 0);

  for (auto it = topography.getVertex(); !it.end(); ++it)
  {
    // Scalar avgHeight = 0;
    // for (auto vit = it->getVertex(); !vit.end(); ++vit)
    //   avgHeight += vit->z();
    // avgHeight /= 2;

    if (it->z() > 200)
      topography.setCorner(it->getIndex());
  }

  MMG::Optimize().setAngleDetection(false)
                      .setHausdorff(hausdorff)
                      .setHMin(hmin)
                      .setHMax(hmax)
                      .optimize(topography);

  // Make a fire somewhere
  Alert::Info() << "Initializing fire..." << Alert::Raise;
  {
    P1 fes(topography);

    // Compute elevation
    GridFunction elevation(fes);
    elevation = [](const Point& p) { return p.z(); };
    // topography.save("Elevation.mesh");
    // elevation.save("Elevation.gf");

    GridFunction phi(fes);
    phi = [](const Point& p)
    {
      const Scalar r = 1000;
      Scalar rr = 0;

      Math::FixedSizeVector<3> c0;
      c0 << 19500, 19500, p.z();
      rr = (p.getCoordinates() - c0).norm() - r;

      // Math::FixedSizeVector<3> c1;
      // c1 << 37500, 37500, p.z();
      // rr = std::min(rr, (p.getCoordinates() - c1).norm() - r);

      // Math::FixedSizeVector<3> c2;
      // c2 << 37500, 12500, p.z();
      // rr = std::min(rr, (p.getCoordinates() - c2).norm() - r);

      // Math::FixedSizeVector<3> c3;
      // c3 << 12500, 37500, p.z();
      // rr = std::min(rr, (p.getCoordinates() - c3).norm() - r);

      return rr;
    };

    topography = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                            .setHMax(hmax)
                                            .setHMin(hmin)
                                            .setHausdorff(hausdorff)
                                            .discretize(phi);

  }

  // Define vegetal stratum
  Environment::VegetalStratum stratum;

  stratum.a = 0.05;
  stratum.A0 = 2.25;
  stratum.R00 = 0.05;
  stratum.u00 = 80;

  stratum.pv = [](const Point&) { return 680; };
  stratum.m = [](const Point&) { return 0.1; };
  stratum.sv = [](const Point&) { return 4550; };
  stratum.sigma = [](const Point&) { return 0.5; };
  stratum.tau = [](const Point&) { return 20; };
  stratum.e = [](const Point&) { return 4.0; };
  stratum.v =
    [&](const Point& v) -> Scalar
    {
      Scalar pv = stratum.pv(v);
      assert(pv > 0);
      return std::min(1.0, stratum.sv(v) * stratum.sigma(v) / (4.0 * pv));
    };
  stratum.A =
    [&](const Point& v) -> Scalar
    {
      Scalar absorption = stratum.v(v);
      return absorption * (stratum.A0 / (1 + stratum.a * stratum.m(v)));
    };
  stratum.R0 =
    [&](const Point& v) -> Scalar
    {
      return stratum.e(v) / stratum.sigma(v) * stratum.R00 / (1 + stratum.a * stratum.m(v));
    };
  stratum.u0 =
    [&](const Point& v) -> Scalar
    {
      return stratum.u00 * stratum.sigma(v) / stratum.tau(v);
    };
  stratum.v0 =
    [&](const Point& v) -> Scalar
    {
      return 12 * stratum.R0(v);
    };

  // Define environment and step through it
  Alert::Info() << "Starting simulation..." << Alert::Raise;
  Environment environment(topography, stratum);
  const Scalar dt = 1200;
  Scalar t = 0;
  for (int i = 0; i < std::numeric_limits<int>::max(); i++)
  {
    topography.save("FirePropagation.mfem.mesh");
    Alert::Info info;
    info << "i: " << i << " | "
         << "t: " << (t / 60) << "m";

    topography.save("out/FirePropagation.mfem."  + std::to_string(i) + ".mesh", IO::FileFormat::MFEM);
    topography.save("out/FirePropagation.medit." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);

    auto start = std::chrono::steady_clock::now();
    environment.step(dt);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    info << " | cpu: " << std::setw(4) << diff.count() << "s";
    info << Alert::NewLine << topography.getElementCount();
    info.raise();
    t += dt;
  }

  return 0;
}

