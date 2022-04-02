#ifndef RODIN_EXTERNAL_MMG_UTILITY_H
#define RODIN_EXTERNAL_MMG_UTILITY_H

#include <mfem.hpp>

#include "Rodin/Variational/GridFunction.h"

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

  void MMG5_Sol_Copy(const MMG5_pSol src, MMG5_pSol dst);

  void MMG5_Sol_Swap(MMG5_pSol a, MMG5_pSol b);

  void MMG5_Sol_Free(MMG5_pSol sol);

  void MMG5_Sol_To_Rodin_IncompleteGridFunction(
      const MMG5_pSol src, Variational::IncompleteGridFunction& dst);

  /**
   * @internal
   *
   * This method will allocate the memory required for the MMG5_pSol->m field.
   */
  template <class FES>
  void Rodin_GridFunction_To_MMG5_Sol(const Variational::GridFunction<FES>& gf, MMG5_pSol sol)
  {
    auto [data, size] = gf.getData();
    if (size)
    {
      int vdim = gf.getFiniteElementSpace().getVectorDimension();
      assert(size % vdim == 0);
      size_t n = size / vdim;
      sol->np  = n;
      sol->npi = n;
      sol->npmax = std::max({static_cast<int>(1.5 * sol->np), MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
      MMG5_SAFE_CALLOC(sol->m, (sol->size * (sol->npmax + 1)), double,
          Alert::Exception("Failed to allocate memory for MMG5_pSol->m").raise());
      if (vdim == 1)
      {
        std::copy(data, data + size, sol->m + 1);
      }
      else
      {
        switch (gf.getFiniteElementSpace().getFES().GetOrdering())
        {
          case mfem::Ordering::byNODES:
          {
            for (size_t i = 0; i < n; i++)
              for (size_t j = 0; j < vdim; j++)
                sol->m[(i + 1) * sol->size + j] = data[i + j * n];
            break;
          }
          case mfem::Ordering::byVDIM:
          {
            std::copy(data, data + size, sol->m + sol->size);
            break;
          }
        }
      }
    }
  }

  /**
   * @internal
   */
  void Load_MMG5_Sol(const boost::filesystem::path&, int meshDim, MMG5_pSol sol);
}

#endif
