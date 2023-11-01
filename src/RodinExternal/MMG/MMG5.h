/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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

namespace Rodin::External::MMG
{
  template <class FuncName>
  class MMG5Exception : public Alert::Exception
  {
    public:
      MMG5Exception(const FuncName& funcName)
      {
        const auto& className = boost::typeindex::type_id_with_cvr<MMG5>().pretty_name();
        *this << "In " << Alert::Identifier::Function(funcName)
              << " of class " << Alert::Identifier::Class(className) << ": ";
      }
  };

  /// Type of return code used by the MMG functions.
  using ReturnCode = int;

  /**
   * @brief Class representing MMG5 objects utilized in the MMG framework.
   *
   * This class is used for wrappping the functionality of the MMG library.
   */
  class MMG5
  {
    public:
      static constexpr int s_meshVersionFormatted = 2;

      // ---- Mesh methods ---------------------------------------------------
      /**
       * @internal
       */
      static MMG5_pMesh createMesh(size_t version, size_t dim, std::optional<size_t> spaceDim = {});

      /**
       * @internal
       * @brief Copies source mesh to a destination mesh.
       *
       * This method performs the necessary memory allocations when copying the
       * data.
       */
      static void copyMesh(const MMG5_pMesh src, MMG5_pMesh dst);

      /**
       * @internal
       * @brief Determines if a mesh is surface or not.
       * @param[in] mesh Mesh.
       */
      static bool isSurfaceMesh(const MMG5_pMesh mesh);

      /**
       * @internal
       * @brief Destroys the mesh object and frees the allocated memory.
       */
      static void destroyMesh(MMG5_pMesh);

      /**
       * @brief Converts an MMG::Mesh object to a MMG5_pMesh object.
       */
      static MMG5_pMesh rodinToMesh(const Rodin::Geometry::SerialMesh& src);

      /**
       * @brief Converts an MMG5_pMesh object to a MMG::Mesh object.
       */
      static MMG::Mesh meshToRodin(const MMG5_pMesh src);

      // ---- Solution methods -----------------------------------------------

      /**
       * @internal
       * @brief Constructs a solution and allocates space for it.
       */
      static MMG5_pSol createSolution(MMG5_pMesh mesh, size_t vdim);

      /**
       * @internal
       */
      static void copySolution(const MMG5_pSol src, MMG5_pSol dst);

      /**
       * @internal
       * @brief Copies the solution from MMG5_pSol object to an MMG::GridFunction object.
       * @tparam Range type of value
       */
      template <class Range>
      static void copySolution(const MMG5_pSol src, MMG::GridFunction<Range>& dst)
      {
        assert(src);
        if constexpr (std::is_same_v<Scalar, Range>)
        {
          assert(src->type == MMG5_Scalar);
          assert(dst.getFiniteElementSpace().getVectorDimension() == 1);
          Math::Matrix& data = dst.getData();
          assert(data.rows() == 1);
          data.resize(1, src->np);
          // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
          // + np + 1.
          std::copy(src->m + 1, src->m + src->np + 1, data.data());
        }
        else if constexpr (std::is_same_v<Math::Vector, Range>)
        {
          const size_t vdim = src->size;
          assert(src->type == MMG5_Vector);
          assert(vdim == dst.getFiniteElementSpace().getVectorDimension());
          Math::Matrix& data = dst.getData();
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
       * @brief Copies the solution from MMG::GridFunction object to an MMG5_pSol object.
       * @tparam Range type of value
       */
      template <class Range>
      static void copySolution(const MMG::GridFunction<Range>& src, MMG5_pSol dst)
      {
        assert(dst);
        if constexpr (std::is_same_v<Scalar, Range>)
        {
          assert(dst->type == MMG5_Scalar);
          assert(src.getFiniteElementSpace().getVectorDimension() == 1);
          const Math::Matrix& data = src.getData();
          assert(data.rows() == 1);
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
              MMG5_SAFE_CALLOC(dst->m, 2 * (dst->npmax + 1), double,
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
        else if constexpr (std::is_same_v<Math::Vector, Range>)
        {
          assert(dst->type == MMG5_Vector);
          const size_t vdim = src.getFiniteElementSpace().getVectorDimension();
          assert(dst->size >= 0);
          assert(vdim == static_cast<size_t>(dst->size));
          const Math::Matrix& data = src.getData();
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
              MMG5_SAFE_CALLOC(dst->m, (dst->size + 1) * (dst->npmax + 1), double,
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
       */
      static void swapSolution(MMG5_pSol a, MMG5_pSol b);

      /**
       * @internal
       * @brief Destroys and frees the allocated memory for a MMG5_pSol object.
       */
      static void destroySolution(MMG5_pSol sol);

      /**
       * @brief Default constructor.
       */
      MMG5();

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
      std::optional<double> m_hmin, m_hmax, m_hausd, m_hgrad;
      bool m_ridgeDetection;
  };
}
#endif
