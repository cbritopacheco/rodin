#include "QFGG.h"


namespace Rodin::Variational
{
  std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_VARIATIONAL_QFGG_MAX_ORDER> QFGG::s_points = {};
  std::array<Geometry::GeometryIndexed<Math::Vector>, RODIN_VARIATIONAL_QFGG_MAX_ORDER> QFGG::s_weights = {};
}

