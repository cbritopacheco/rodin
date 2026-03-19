/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MMG5.h
 * @brief Core MMG wrapper utilities for mesh/solution conversion and parameters.
 */
#ifndef RODIN_RODINEXTERNAL_MMG_MMG5_H
#define RODIN_RODINEXTERNAL_MMG_MMG5_H

#include <boost/filesystem.hpp>

#include "Rodin/Math.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Utility/IsSpecialization.h"

#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"
#include "GridFunction.h"

namespace Rodin::MMG
{
  /**
   * @brief Exception type for MMG wrapper failures.
   * @tparam FuncName Type of the function-name object passed by caller.
   *
   * This exception prepends a standardized "class/function" prefix to MMG
   * diagnostic messages, then can be streamed with additional details before
   * raising through Rodin's alert system.
   */
  template <class FuncName>
  class MMG5Exception : public Alert::Exception
  {
    public:
      /**
       * @brief Builds an MMG wrapper exception message prefix.
       * @param[in] funcName Function identifier used in the diagnostic message.
       */
      MMG5Exception(const FuncName& funcName)
      {
        const auto& className = boost::typeindex::type_id_with_cvr<MMG5>().pretty_name();
        *this << "In " << Alert::Identifier::Function(funcName)
              << " of class " << Alert::Identifier::Class(className) << ": ";
      }
  };

  /**
   * @brief Return-code type used by MMG C API calls.
   *
   * MMG kernels conventionally return `MMG5_SUCCESS`, `MMG5_LOWFAILURE`, or
   * `MMG5_STRONGFAILURE`.
   */
  using ReturnCode = int;

  /**
   * @brief Base class for MMG-backed operations in Rodin.
   *
   * Provides:
   * - allocation/deallocation helpers for MMG mesh and solution structures,
   * - conversion between Rodin meshes/grid functions and MMG native data,
   * - shared remeshing parameter configuration (hmin/hmax/hausd/hgrad/angle).
   *
   * High-level operators such as @ref Adapt, @ref Optimizer, and
   * @ref LevelSetDiscretizer build on this class.
   */
  class MMG5
  {
    public:
      /// MMG mesh version tag used in generated MMG mesh objects.
      static constexpr int s_meshVersionFormatted = 2;

      // ---- Mesh methods ---------------------------------------------------
      /**
       * @internal
       * @brief Allocates and initializes an MMG mesh object.
       * @param[in] version MMG mesh version.
       * @param[in] dim Topological mesh dimension.
       * @param[in] spaceDim Optional embedding space dimension.
       * @returns Newly allocated MMG mesh pointer.
       */
      static MMG5_pMesh createMesh(size_t version, size_t dim, Optional<size_t> spaceDim = {});

      /**
       * @internal
       * @brief Copies source mesh to a destination mesh.
       * @param[in] src Source MMG mesh.
       * @param[out] dst Destination MMG mesh.
       *
       * Performs deep-copy allocation for dynamically owned MMG arrays.
       */
      static void copyMesh(const MMG5_pMesh src, MMG5_pMesh dst);

      /**
       * @internal
       * @brief Determines if a mesh is surface or not.
       * @param[in] mesh Mesh handle.
       * @returns `true` if the mesh is a surface manifold.
       */
      static bool isSurfaceMesh(const MMG5_pMesh mesh);

      /**
       * @internal
       * @brief Destroys the mesh object and frees the allocated memory.
       * @param[in] mesh Pointer to mesh.
       */
      static void destroyMesh(MMG5_pMesh);

      /**
       * @brief Converts a Rodin local mesh to a native MMG mesh.
       * @param[in] src Source Rodin mesh.
       * @returns Newly allocated MMG mesh.
       *
       * The conversion preserves MMG-specific tags when `src` is an
       * @ref MMG::Mesh (corners, ridges, required entities).
       */
      static MMG5_pMesh rodinToMesh(const Rodin::Geometry::LocalMesh& src);

      /**
       * @brief Converts a native MMG mesh to @ref MMG::Mesh.
       * @param[in] src Source MMG mesh handle.
       * @returns Converted Rodin mesh.
       */
      static MMG::Mesh meshToRodin(const MMG5_pMesh src);

      // ---- Solution methods -----------------------------------------------

      /**
       * @internal
       * @brief Constructs a solution and allocates space for it.
       * @param[in] mesh Owning MMG mesh.
       * @param[in] vdim Value dimension (`1` for scalar, `>1` for vector).
       * @returns Newly allocated MMG solution.
       */
      static MMG5_pSol createSolution(MMG5_pMesh mesh, size_t vdim);

      /**
       * @internal
       * @brief Deep-copies an MMG solution.
       * @param[in] src Source solution.
       * @param[out] dst Destination solution.
       */
      static void copySolution(const MMG5_pSol src, MMG5_pSol dst);

      /**
       * @internal
       * @brief Copies values from MMG solution to Rodin MMG grid function.
       * @tparam Range Grid function range type (`Real` or `Math::Vector<Real>`).
       * @param[in] src Source MMG solution.
       * @param[out] dst Destination grid function.
       *
       * The transfer preserves MMG node ordering and writes into the matrix
       * layout used by Rodin's @ref Variational::GridFunction data container.
       */
      template <class Range>
      static void copySolution(const MMG5_pSol src, MMG::GridFunction<Range>& dst)
      {
        assert(src);
        if constexpr (std::is_same_v<Real, Range>)
        {
          assert(src->type == MMG5_Scalar);
          assert(dst.getFiniteElementSpace().getVectorDimension() == 1);
          Math::Matrix<Real>& data = dst.getData();
          assert(data.rows() == 1);
          data.resize(1, src->np);
          // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
          // + np + 1.
          std::copy(src->m + 1, src->m + src->np + 1, data.data());
        }
        else if constexpr (std::is_same_v<Math::Vector<Real>, Range>)
        {
          const size_t vdim = src->size;
          assert(src->type == MMG5_Vector);
          assert(vdim == dst.getFiniteElementSpace().getVectorDimension());
          Math::Matrix<Real>& data = dst.getData();
          assert(data.rows() >= 0);
          assert(static_cast<size_t>(data.rows()) == vdim);
          data.resize(vdim, src->np);
          // MMG5_pSol->m is 1 indexed. We must start at m + vdim and finish at
          // m + vdim * (src->np + 1).
          std::copy(src->m + vdim, src->m + vdim * (src->np + 1), data.data());
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @internal
       * @brief Copies values from Rodin MMG grid function to MMG solution.
       * @tparam Range Grid function range type (`Real` or `Math::Vector<Real>`).
       * @param[in] src Source grid function.
       * @param[out] dst Destination MMG solution.
       *
       * Allocates MMG storage buffers when needed and copies nodal values in the
       * memory layout expected by the selected MMG kernel.
       */
      template <class Range>
      static void copySolution(const MMG::GridFunction<Range>& src, MMG5_pSol dst)
      {
        assert(dst);
        if constexpr (std::is_same_v<Real, Range>)
        {
          assert(dst->type == MMG5_Scalar);
          assert(src.getFiniteElementSpace().getVectorDimension() == 1);
          const Math::Matrix<Real>& data = src.getData();
          assert(data.cols() == 1);
          assert(dst->size == 1);
          const size_t n = data.size();
          if (n)
          {
            dst->np  = n;
            dst->npi = n;
            dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
            assert(dst->size == 1);
            assert(dst->np < dst->npmax);
            if (!dst->m)
            {
              // So 2 * (dst->np + 1) seems to work for most applications
              MMG5_SAFE_CALLOC(dst->m, 2 * (dst->npmax + 1), Real,
                Alert::Exception() << "Failed to allocate memory for MMG5_pSol->m." << Alert::Raise);
            }
            std::copy(data.data(), data.data() + n, dst->m + 1);
          }
          else
          {
            dst->np  = 0;
            dst->npi = 0;
            dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
          }
        }
        else if constexpr (std::is_same_v<Math::Vector<Real>, Range>)
        {
          assert(dst->type == MMG5_Vector);
          const size_t vdim = src.getFiniteElementSpace().getVectorDimension();
          assert(dst->size >= 0);
          assert(vdim == static_cast<size_t>(dst->size));
          const Math::Matrix<Real>& data = src.getData();
          assert(dst->size == data.rows());
          const size_t n = data.cols();
          assert(n > 0);
          dst->np  = n;
          dst->npi = n;
          dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
          assert(dst->np < dst->npmax);
          if (n)
          {
            if (!dst->m)
            {
              // So (dst->size + 1) * (dst->np + 1) seems to work for most
              // applications
              MMG5_SAFE_CALLOC(dst->m, (dst->size + 1) * (dst->npmax + 1), Real,
                Alert::Exception() << "Failed to allocate memory for MMG5_pSol->m" << Alert::Raise);
            }
            std::copy(data.data(), data.data() + data.size(), dst->m + dst->size);
          }
          else
          {
            dst->np  = 0;
            dst->npi = 0;
            dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
          }
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @internal
       * @brief Swaps the data between two instances of type MMG5_pSol.
       * @param[in, out] a First solution.
       * @param[in, out] b Second solution.
       */
      static void swapSolution(MMG5_pSol a, MMG5_pSol b);

      /**
       * @internal
       * @brief Destroys and frees the allocated memory for a MMG5_pSol object.
       * @param[in] sol MMG solution pointer.
       */
      static void destroySolution(MMG5_pSol sol);

      /**
       * @brief Default constructor.
       */
      MMG5();

      /**
       * @brief Enables/disables ridge angle detection.
       * @param[in] b If `true`, enable angle detection.
       * @returns Reference to this object.
       */
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
       * @see setHMax(Real)
       */
      MMG5& setHMin(Real hmin)
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
       * @see setHMin(Real)
       */
      MMG5& setHMax(Real hmax)
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
      MMG5& setHausdorff(Real hausd)
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
      MMG5& setGradation(Real hgrad)
      {
        m_hgrad = hgrad;
        return *this;
      }

    protected:
      /**
       * @brief Applies configured remeshing parameters to an MMG mesh.
       * @param[in, out] mesh Target MMG mesh.
       * @returns Reference to this object.
       */
      MMG5& setParameters(MMG5_pMesh mesh);

    private:
      Optional<Real> m_hmin, m_hmax, m_hausd, m_hgrad; ///< Optional MMG scalar parameters.
      bool m_ridgeDetection; ///< Whether ridge angle detection is enabled.
  };
}
#endif
