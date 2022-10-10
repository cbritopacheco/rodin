#include <cmath>
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

    struct Ground
    {
      double a, b, c, d;
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

    Environment(Mesh<Context>& topography)
      : m_topography(topography),
        m_sfes(m_topography),
        m_vfes(m_topography, m_topography.getSpaceDimension()),
        m_wind(m_sfes),
        m_flame(m_sfes)
    {}

    void step(double dt);

    const GridFunction<FES>& getWind() const
    {
      return m_wind;
    }

  private:
    Mesh<Context>& m_topography;

    FES m_sfes;
    FES m_vfes;
    GridFunction<FES> m_wind;

    Flame m_flame;
};


int main()
{
  const char* meshfile = "../resources/mfem/fire-propagation.mesh";
  Mesh topography;
  topography.load(meshfile);

  Environment environment(topography);
  // Environment::Flame flame(environment);

  return 0;
}

