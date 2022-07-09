#ifndef RODIN_EXTERNAL_MMG_UTILITY_H
#define RODIN_EXTERNAL_MMG_UTILITY_H

#include <mfem.hpp>

#include "Rodin/Variational/GridFunction.h"

#include "Common.h"

namespace Rodin::External::MMG
{
  MMG5_pMesh MMG5_Create_Mesh(int version, int spaceDim);

  MMG5_pMesh MMG5_Load_Mesh(const boost::filesystem::path& filename);

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

  MMG5_pSol MMG5_Create_Sol(MMG5_pMesh mesh, int vdim);

  void MMG5_Sol_Copy(const MMG5_pSol src, MMG5_pSol dst);

  void MMG5_Sol_Swap(MMG5_pSol a, MMG5_pSol b);

  void MMG5_Sol_Free(MMG5_pSol sol);

  void MMG5_Sol_To_Rodin_GridFunction(
      const MMG5_pSol src, Variational::GridFunctionBase& dst);

  /**
   * @internal
   *
   * This method will allocate the memory required for the MMG5_pSol->m field.
   */
  void Rodin_GridFunction_To_MMG5_Sol(
      const Variational::GridFunctionBase& gf, MMG5_pSol sol);

  template <class Trait>
  void Rodin_Mesh_To_MMG5_Mesh(const Rodin::Mesh<Trait>& src, MMG5_pMesh dst);

  template <>
  void Rodin_Mesh_To_MMG5_Mesh<Traits::Serial>(
      const Rodin::Mesh<Traits::Serial>& src, MMG5_pMesh dst);

  /**
   * @internal
   */
  void Load_MMG5_Sol(const boost::filesystem::path&, int meshDim, MMG5_pSol sol);
}

#endif
