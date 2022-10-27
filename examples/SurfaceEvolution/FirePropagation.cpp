#include <cmath>
#include <boost/bimap.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace std;
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

class Environment
{
  public:
    using Context = Context::Serial;
    using FES = H1<Context>;

    enum class Terrain
    {
      Vegetation,
      Fire,
      Burnt
    };

    struct Plane
    {
      double a, b, c, d;
    };

    struct Parameters
    {
      double R0, u0, A, v0;
    };

    class Flame
    {
      public:
        Flame(Environment& env)
          : m_env(env),
            m_height(env.m_sfes),
            m_tiltAngle(env.m_sfes),
            m_rateOfSpread(env.m_sfes)
        {}

        Flame& step(double dt)
        {
          double alpha = 0.0;

          auto p = Grad(m_env.m_terrainHeight);

          double u = m_env.getParameters().u0 / std::cos(alpha);


          // Angle between wind and conormal
          auto phi =
            ScalarFunction(
                [&](const Vertex& v) -> double
                {
                  return 1.0;
                });

          // Angle between gradient and conormal
          auto psi =
            ScalarFunction(
                [](const Vertex& v) -> double
                {
                  return 1.0;
                });

          // Compute the tilt angle
          auto gamma =
            ScalarFunction(
              [&](const Vertex& v) -> double
              {
                double rhs =
                  std::tan(alpha) * std::cos(phi(v)) + m_env.getWind()(v) * std::cos(psi(v));
                return std::atan(rhs);
              });
          return *this;
        }

        const GridFunction<FES>& getHeight() const
        {
          return m_height;
        }

        const GridFunction<FES>& getRateOfSpread() const
        {
          return m_rateOfSpread;
        }

      private:
        Environment& m_env;
        GridFunction<FES> m_height;
        GridFunction<FES> m_tiltAngle;
        GridFunction<FES> m_rateOfSpread;
    };

    Environment(const Parameters& params, Mesh<Context>& topography)
      : m_params(params),
        m_topography(topography),
        m_sfes(m_topography),
        m_vfes(m_topography, m_topography.getSpaceDimension()),
        m_wind(m_sfes),
        m_terrainHeight(m_sfes),
        m_ground({0, 0, 0, 0}),
        m_groundSlope(m_vfes),
        m_flame(*this)
    {
      m_terrainMap.insert({1, Terrain::Vegetation});
      m_terrainMap.insert({2, Terrain::Fire});
      m_terrainMap.insert({3, Terrain::Burnt});

      m_terrainHeight =
        ScalarFunction(
          [](const Vertex& v) -> double
          {
            return v.z();
          });

      m_groundSlope = Grad(m_terrainHeight);

      m_terrainHeight.save("terrainHeight.gf");
      topography.save("groundSlope.mesh");
      m_groundSlope.save("groundSlope.gf");
    }

    Environment& step(double dt)
    {
      m_flame.step(dt);
      return *this;
    }

    Environment& setGround(const Plane& ground)
    {
      m_ground = ground;
      return *this;
    }

    Environment& setTerrainMap(const boost::bimap<int, Terrain>& terrainMap)
    {
      m_terrainMap = terrainMap;
      return *this;
    }

    const Flame& getFlame() const
    {
      return m_flame;
    }

    const GridFunction<FES>& getWind() const
    {
      return m_wind;
    }

    const Parameters& getParameters() const
    {
      return m_params;
    }

  private:
    Parameters m_params;
    Mesh<Context>& m_topography;

    FES m_sfes;
    FES m_vfes;

    GridFunction<FES> m_wind;

    GridFunction<FES> m_terrainHeight;
    boost::bimap<int, Terrain> m_terrainMap;

    Plane m_ground;
    GridFunction<FES> m_groundSlope;

    Flame m_flame;
};


int main()
{
  const char* meshfile = "topo.mesh";
  Mesh topography;
  topography.load(meshfile);
  Environment environment({0, 0, 0, 0}, topography);
  return 0;
}
