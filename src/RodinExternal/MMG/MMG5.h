/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODINEXTERNAL_MMG_MMG5_H
#define RODIN_RODINEXTERNAL_MMG_MMG5_H

#include <boost/filesystem.hpp>
#include <mfem.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"

namespace Rodin::External::MMG
{
  class MMG5
  {
   public:
    MMG5()
      : m_ridgeDetection(true)
    {}

    static constexpr int s_meshVersionFormatted = 2;

    // ---- Mesh methods ---------------------------------------------------
    static MMG5_pMesh createMesh(int version, int dim, std::optional<int> spaceDim = {});

    static void copyMesh(const MMG5_pMesh src, MMG5_pMesh dst);

    static bool isSurfaceMesh(MMG5_pMesh mesh);

    static void destroyMesh(MMG5_pMesh);

    static MMG5_pMesh rodinToMesh(const MMG::Mesh& src);

    static MMG::Mesh meshToRodin(const MMG5_pMesh src);

    // ---- Solution methods -----------------------------------------------
    static MMG5_pSol createSolution(MMG5_pMesh mesh, int vdim);

    static void copySolution(const MMG5_pSol src, MMG5_pSol dst);

    template <class FES>
    static void copySolution(const MMG5_pSol src, Variational::GridFunction<FES>& dst)
    {
      assert(src->type == MMG5_Scalar);
      double* data = new double[src->np];
      // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
      // + np + 1.
      std::copy(src->m + 1, src->m + src->np + 1, data);
      dst.getHandle().SetDataAndSize(data, src->np);
      dst.getHandle().MakeDataOwner();
    }

    template <class FES>
    static void copySolution(const Variational::GridFunction<FES>& src, MMG5_pSol dst)
    {
      assert(dst);
      auto [data, size] = src.getData();
      if (size)
      {
       int vdim = src.getFiniteElementSpace().getVectorDimension();
       assert(size % vdim == 0);
       assert(dst->size == vdim);
       size_t n = size / vdim;
       assert(n > 0);
       dst->np  = n;
       dst->npi = n;
       dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
       assert(dst->np < dst->npmax);
       if (!dst->m)
       {
        // So (dst->size + 1) * (dst->np + 1) seems to work for most
        // applications
        MMG5_SAFE_CALLOC(dst->m, (dst->size + 1) * (dst->np + 1), double,
            Alert::Exception("Failed to allocate memory for MMG5_pSol->m").raise());
       }
       if (vdim == 1)
       {
         std::copy(data, data + size, dst->m + 1);
       }
       else
       {
        switch (src.getFiniteElementSpace().getHandle().GetOrdering())
        {
          case mfem::Ordering::byNODES:
          {
            for (size_t i = 0; i < n; i++)
              for (size_t j = 0; j < vdim; j++)
                dst->m[(i + 1) * dst->size + j] = data[i + j * n];
           break;
          }
          case mfem::Ordering::byVDIM:
          {
           std::copy(data, data + size, dst->m + dst->size);
           break;
          }
        }
       }
      }
    }

    static void swapSolution(MMG5_pSol a, MMG5_pSol b);

    static void destroySolution(MMG5_pSol sol);

   protected:
    MMG5& setAngleDetection(bool b = true)
    {
      m_ridgeDetection = b;
      return *this;
    }

    /**
     * @brief Sets the minimal edge size.
     *
     * @param[in] hmin Minimal edge size.
     *
     * @returns Reference to self (for method chaining)
     *
     * The `hmin` option allows to truncate the edge sizes to be lower than the
     * `hmax` one.
     *
     * The default values for this parameters are computed from the mesh
     * bounding box or, if provided, from the given metric.
     *
     * - Without metric, the minimal edge size is set to 0.01 of the bounding
     * box size.
     *
     * - With metric, the minimal edge size is set to 0.1 of the
     * smallest prescribed size.
     *
     * @see setHMax(double)
     */
    MMG5& setHMin(double hmin)
    {
      m_hmin = hmin;
      return *this;
    }

    /**
     * @brief Sets the maximal edge size parameter.
     *
     * @param[in] hmax Maximal edge size.
     *
     * @returns Reference to self (for method chaining)
     *
     * The `hmax` option allows to truncate the edge sizes to be greater than
     * the `hmin` parameter.
     *
     * The default values for this parameters are computed from the mesh
     * bounding box or, if provided, from the given metric.
     *
     * - Without metric, the maximal edge size is set to two times the
     *  bounding box size.
     *
     * - With metric, the maximal one is set to 10 times the maximal
     *  prescribed size.
     *
     * @see setHMin(double)
     */
    MMG5& setHMax(double hmax)
    {
      m_hmax = hmax;
      return *this;
    }

    /**
     * @brief Sets the Hausdorff parameter.
     *
     * @param[in] hausd Hausdorff parameter.
     *
     * @returns Reference to self (for method chaining)
     *
     * The Hausdorff parameter controls the boundary approximation.  It
     * imposes the maximal distance between the piecewise linear
     * representation of the boundary and the reconstructed ideal boundary.
     * Thus, a low Hausdorff parameter leads to the refinement of high
     * curvature areas.
     *
     * By default, the Hausdorff value is set to 0.01, which is a suitable
     * value for an object of size 1 in each direction. For smaller (resp.
     * larger) objects, you may need to decrease (resp. increase) the
     * Hausdorff parameter.
     *
     */
    MMG5& setHausdorff(double hausd)
    {
      m_hausd = hausd;
      return *this;
    }

    /**
     * @brief Sets the gradation parameter
     *
     * @param[in] hgrad Gradation parameter
     *
     * @returns Reference to self (for method chaining)
     *
     * The `hgrad` option allows to set the gradation value. It controls the
     * ratio between two adjacent edges. With a gradation of @f$ h @f$, two
     * adjacent edges @f$ e_1 @f$ and @f$ e_2 @f$ must respect the following
     * constraint:
     *
     * @f[
     *  \dfrac{1}{h} \leq \dfrac{ |e_1| }{ |e_2| } \leq h
     * @f]
     *
     * By default, the gradation value is 1.3.
     *
     */
    MMG5& setGradation(double hgrad)
    {
      m_hgrad = hgrad;
      return *this;
    }

   protected:
    MMG5& setParameters(MMG5_pMesh mesh);

   private:
    std::optional<double> m_hmin,
                   m_hmax,
                   m_hausd,
                   m_hgrad;
    bool m_ridgeDetection;
  };
}
#endif
