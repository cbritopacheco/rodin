#ifndef RODIN_EXTERNAL_MMG_UTILITY_H
#define RODIN_EXTERNAL_MMG_UTILITY_H

#include <mfem.hpp>

#include "Common.h"

namespace Rodin::External::MMG
{
  /**
   * @internal
   * @brief Performs a deep copy of the MMG5_Mesh objects.
   *
   * This method expects that the destination has already been allocated.
   *
   * @param[in] src Source mesh
   * @oaram[in] dst Destination mesh
   */
  void MMG5_Mesh_Copy(const MMG5_pMesh src, MMG5_pMesh dst);
}

#endif
