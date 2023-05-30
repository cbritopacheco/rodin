#include "P1Element.h"

namespace Rodin::Variational
{
  const std::array<Math::Matrix, 3> P1Element<Scalar>::s_dofs =
  {
    Math::Matrix{{0}},
    Math::Matrix{{0}, {1}},
    Math::Matrix{{0, 1, 0},
                 {0, 0, 1}}
  };
}

