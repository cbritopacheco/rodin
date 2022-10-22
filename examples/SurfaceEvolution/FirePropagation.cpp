#include <cmath>
#include <boost/bimap.hpp>
#include <boost/math/tools/roots.hpp>

#include "Rodin/Mesh.h"
#include "Rodin/Variational.h"
#include "RodinExternal/MMG.h"

using namespace std;
using namespace Rodin;
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

    struct Flame
    {
      Flame(FES& sfes)
        : m_height(sfes),
          m_rateOfSpread(sfes)
      {}

      const GridFunction<FES>& getHeight() const
      {
        return m_height;
      }

      const GridFunction<FES>& getRateOfSpread() const
      {
        return m_rateOfSpread;
      }

      private:
        GridFunction<FES> m_height;
        GridFunction<FES> m_rateOfSpread;
    };

    Environment(Mesh<Context>& topography, const Parameters& params)
      : m_topography(topography),
        m_ground({0, 0, 0, 0}),
        m_params(params),
        m_sfes(m_topography),
        m_vfes(m_topography, m_topography.getSpaceDimension()),
        m_wind(m_sfes),
        m_flame(m_sfes)
    {
      m_terrainMap.insert({1, Terrain::Vegetation});
      m_terrainMap.insert({2, Terrain::Fire});
      m_terrainMap.insert({3, Terrain::Burnt});
    }

    void step(double dt);

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

    const GridFunction<FES>& getWind() const
    {
      return m_wind;
    }

  private:
    Mesh<Context>& m_topography;
    Plane m_ground;
    Parameters m_params;
    FES m_sfes;
    FES m_vfes;
    GridFunction<FES> m_wind;
    Flame m_flame;
    boost::bimap<int, Terrain> m_terrainMap;
};


int main()
{
  const char* meshfile = "topo.mesh";
  Mesh topography;
  topography.load(meshfile);
  // topography.initialize(2, 3)
  //           .vertex({-0.5, -0.5, 0})
  //           .vertex({-0.5,  0.5, 0})
  //           .vertex({0.5,  -0.5, 0})
  //           .vertex({0.5,   0.5, 0})
  //           .element(Geometry::Triangle, {0, 1, 2})
  //           .element(Geometry::Triangle, {1, 2, 3})
  //           .finalize()
  //           .refine()
  //           .refine()
  //           .refine()
  //           .refine()
  //           .refine()
  //           .refine()
  //           ;

  // H1 fes(topography, 3);
  // GridFunction disp(fes);
  // disp = VectorFunction{
  //   0,
  //   0,
  //   [](const Vertex& v) -> double
  //   {
  //     return std::exp(-40 * (v.x() * v.x() + v.y() * v.y()));
  //   }
  // };

  // disp.save("disp.gf");
  // topography.displace(disp);
  // topography.save("topo.mesh");
  Environment environment(topography, {0, 0, 0, 0});

  return 0;
}

