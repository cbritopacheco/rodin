#include "QFGM.h"

namespace Rodin::QF
{
  std::vector<std::array<Math::Vector, RODIN_QFGM_MAX_ORDER>>
  QFGM::initializeWeights()
  {
    return {};
  }

  std::vector<std::array<Math::Matrix, RODIN_QFGM_MAX_ORDER>>
  QFGM::initializePoints()
  {
     return {};
  }

  const std::vector<std::array<Math::Vector, RODIN_QFGM_MAX_ORDER>> QFGM::s_weights = initializeWeights();
  const std::vector<std::array<Math::Matrix, RODIN_QFGM_MAX_ORDER>> QFGM::s_points = initializePoints();
}

