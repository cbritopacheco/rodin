#include <cmath>
#include <boost/bimap.hpp>

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

static constexpr double tgv = std::numeric_limits<double>::max();

class Environment
{
  public:
    using Context = Context::Serial;
    using FES = H1<Context>;

    struct Plane
    {
      double a, b, c, d;
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
      double A0;

      // Vertical velocity
      double u00;

      // Energy ratio between incident radiant energy emitted from the flame
      // base and ignition energy of dry fuel
      double R00;

      // Moisture factor
      double a;

      // Rate of spread for no slope and no wind (m/s)
      std::function<double(const Point&)> R0;

      // Buoyancy velocity component for a zero slope (m/s)
      std::function<double(const Point&)> u0;

      // Energy ratio between incident radiant energy and ignition energy of wet fuel
      std::function<double(const Point&)> A;

      // ROS coefficient (m/s)
      std::function<double(const Point&)> v0;

      // Moisture content (%)
      std::function<double(const Point&)> m;

      // Thickness of the vegetal stratum (m)
      std::function<double(const Point&)> e;

      // Surface density of vegetal fuel, (kg/m^2)
      std::function<double(const Point&)> sigma;

      // Residence time (s)
      std::function<double(const Point&)> tau;

      // Vegetal fuel surface to volume ratio (1/m)
      std::function<double(const Point&)> sv;

      // Gas flame density (kg/m^3)
      std::function<double(const Point&)> pv;

      // Absorption coefficient
      std::function<double(const Point&)> v;

      // Flame gas temperature (K)
      std::function<double(const Point&)> T;

      // Air temperature (K)
      std::function<double(const Point&)> Ta;
    };

    class Flame
    {
      public:
        Flame(Environment& env)
          : m_env(env),
            m_direction(new GridFunction(*env.m_vfes))
        {}

        Flame& step(double dt)
        {
          assert(dt > 0);

          // Slope vector
          auto p = Grad(*m_env.m_terrainHeight);
          p.traceOf(Terrain::Burnt);

          // Angle between slope and ground
          auto alpha =
            [&](const Point& v)
            {
              double angle = std::acos(p.z()(v) / Frobenius(p)(v));
              assert(angle >= 0);
              // assert(angle <= Math::Constants::pi<double>() / 2.0);
              return Math::Constants::pi<double>() / 2.0 - angle;
            };

          Solver::UMFPack solver;

          TrialFunction d(*m_env.m_vfes);
          TestFunction  v(*m_env.m_vfes);

          const auto& dist = *m_env.m_fireDist;
          auto gdist = Grad(dist);
          gdist.traceOf(Terrain::Burnt);

          double lambda = 1.0;
          Problem cnd(d, v);
          cnd = Integral(lambda * Jacobian(d), Jacobian(v))
              + Integral(d, v)
              + DirichletBC(d, gdist).on(Terrain::Fire);
          solver.solve(cnd);

          const auto conormal = d.getGridFunction() / Frobenius(d.getGridFunction());

          // Angle between slope and conormal
          auto phi =
            [&](const Point& v) -> double
            {
              auto fn = Dot(p, conormal) / (Frobenius(p) * Frobenius(conormal));
              double fv = fn(v);
              if (std::isfinite(fv))
                return std::acos(fv);
              else
                return 0.0;
            };

          // Angle between wind and conormal
          const auto& wind = m_env.getWind();
          auto psi =
            [&](const Point& v) -> double
            {
              auto fn = Dot(wind, conormal) / (Frobenius(wind) * Frobenius(conormal));
              double fv = fn(v);
              if (std::isfinite(fv))
                return std::acos(fv);
              else
                return 0.0;
            };

          // Compute the tilt angle
          auto gamma =
            [&](const Point& v) -> double
            {
              auto w = Frobenius(wind)(v);
              double rhs =
                std::tan(alpha(v)) * std::cos(phi(v)) + w * std::cos(psi(v));
              return std::atan(rhs);
            };

          // Compute rate of spread
          auto R =
            [&](const Point& v) -> double
            {
              const double g = gamma(v);
              if (g > 0)
              {
                const double R0 = m_env.m_vegetalStratum.R0(v);
                const double v0 = m_env.m_vegetalStratum.v0(v);
                const double A = m_env.m_vegetalStratum.A(v);
                const double Ra = R0 + A * v0 * (
                    1 + std::sin(g) - std::cos(g)) / std::cos(g) - v0 / std::cos(g);
                return 0.5 * (Ra + std::sqrt(Ra * Ra + (4 * v0 * R0) / std::cos(g)));
              }
              else
              {
                return m_env.m_vegetalStratum.R0(v);
              }
            };

          TrialFunction g(*m_env.m_vfes);
          TestFunction  w(*m_env.m_vfes);
          Problem hilbert(g, w);
          hilbert = Integral(lambda * Jacobian(g), Jacobian(w))
                  + Integral(g, w)
                  + DirichletBC(g, VectorFunction{0, 0, 0}).on(Terrain::WorldBorder)
                  + DirichletBC(g, ScalarFunction(R) * conormal).on(Terrain::Fire)
                  ;
          solver.solve(hilbert);

          m_direction.reset(new GridFunction(g.getGridFunction()));
          return *this;
        }

        const GridFunction<FES>& getDirection() const
        {
          return *m_direction;
        }

      private:
        Environment& m_env;
        std::unique_ptr<GridFunction<FES>> m_direction;
    };

    Environment(MMG::Mesh& topography, const VegetalStratum& vegetalStratum)
      : m_topography(topography),
        m_sfes(new FES(m_topography)),
        m_vfes(new FES(m_topography, m_topography.getSpaceDimension())),
        m_wind(new VectorFunction{0.0, 0.0, 0.0}),
        m_terrainHeight(new GridFunction(*m_sfes)),
        m_vegetalStratum(vegetalStratum),
        m_flame(*this),
        m_gravity(-9.8),
        m_fireDist(new GridFunction(*m_sfes)),
        m_elapsedTime(0.0)
    {
      *m_terrainHeight = [](const Point& v) { return v.z(); };
    }

    Environment& step(double dt)
    {
      m_fireDist.reset(
          new GridFunction(
            MMG::Distancer(*m_sfes)
              .setInteriorDomain(Terrain::Burnt)
              .distance(m_topography)));

      // m_topography.save("dist.mesh", IO::FileFormat::MEDIT);
      // m_fireDist->save("dist.sol", IO::FileFormat::MEDIT);

      m_flame.step(dt);
      MMG::Advect(*m_fireDist, m_flame.getDirection()).step(dt);

      m_topography = MMG::ImplicitDomainMesher().setHMax(0.05)
                                                .setHMin(0.002)
                                                .setHausdorff(0.001)
                                                .setAngleDetection(false)
                                                // .split(Terrain::Burnt,
                                                //     {Terrain::Burnt, Terrain::Vegetation})
                                                // .split(Terrain::Vegetation,
                                                //     {Terrain::Burnt, Terrain::Vegetation})
                                                .setBoundaryReference(Terrain::Fire)
                                                .discretize(*m_fireDist);

      // Rebuild finite element spaces with new topography
      m_sfes.reset(new FES(m_topography));
      m_vfes.reset(new FES(m_topography, m_topography.getSpaceDimension()));
      m_terrainHeight.reset(new GridFunction(*m_sfes));
      *m_terrainHeight = [](const Point& v) { return v.z(); };
      m_elapsedTime += dt;
      return *this;
    }

    const Flame& getFlame() const
    {
      return m_flame;
    }

    const VectorFunctionBase& getWind() const
    {
      return *m_wind;
    }

    const Mesh<Context>& getTopography() const
    {
      return m_topography;
    }

    double getGravity() const
    {
      return m_gravity;
    }

  private:
    MMG::Mesh& m_topography;

    std::unique_ptr<FES> m_sfes;
    std::unique_ptr<FES> m_vfes;

    std::unique_ptr<VectorFunctionBase> m_wind;

    std::unique_ptr<GridFunction<FES>> m_terrainHeight;

    VegetalStratum m_vegetalStratum;

    Flame m_flame;

    const double m_gravity;

    std::unique_ptr<GridFunction<FES>> m_fireDist;

    double m_elapsedTime;
};


int main()
{
  const char* meshfile = "topo.mesh";
  MMG::Mesh topography;
  topography.load(meshfile, IO::FileFormat::MEDIT);
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
    [&](const Point& v) -> double
    {
      double pv = stratum.pv(v);
      assert(pv > 0);
      return std::min(1.0, stratum.sv(v) * stratum.sigma(v) / (4.0 * pv));
    };
  stratum.A =
    [&](const Point& v) -> double
    {
      double absorption = stratum.v(v);
      return absorption * (stratum.A0 / (1 + stratum.a * stratum.m(v)));
    };
  stratum.R0 =
    [&](const Point& v) -> double
    {
      return stratum.e(v) / stratum.sigma(v) * stratum.R00 / (1 + stratum.a * stratum.m(v));
    };
  stratum.u0 =
    [&](const Point& v) -> double
    {
      return stratum.u00 * stratum.sigma(v) / stratum.tau(v);
    };
  stratum.v0 =
    [&](const Point& v) -> double
    {
      return 12 * stratum.R0(v);
    };

  Environment environment(topography, stratum);
  for (int i = 0; i < 1000; i++)
  {
    std::cout << "i: " << i << std::endl;

    environment.step(0.01);
    topography.save("out/evolution." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
    topography.save("woof.mesh", IO::FileFormat::MEDIT);
    // environment.getFlame().getDirection().save("direction.gf");
  }

  return 0;
}

