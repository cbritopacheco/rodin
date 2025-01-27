#include "GeometryIndexed.h"
#include "RippleSimplexification.h"

namespace Rodin::Geometry
{
  GeometryIndexed<std::vector<IndexArray>>
  RippleSimplexification<Mesh<Context::Local>>::s_simplices =
  {
    {
      Polytope::Type::Triangle,
      {
        { {0, 1, 2} },
      }
    },
    {
      Polytope::Type::Quadrilateral,
      {
        { {0, 1, 2}, {1, 2, 3} },
        { {0, 1, 2}, {0, 2, 3} }
      }
    },
    {
      Polytope::Type::TriangularPrism,
      {
        { {0, 0, 0, 1, 0} },
        { {0, 0, 1, 1, 0} },
        { {0, 1, 0, 0, 0} },
        { {0, 1, 0, 1, 0} },
        { {0, 1, 1, 0, 0} },
      }
    }
  };
}
