/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file V5.h
 * @brief Native MMG version 5 mesh and solution operations.
 */
#ifndef RODIN_MMG_V5_H
#define RODIN_MMG_V5_H
#include "Base.h"
#include "GridFunction.h"

/**
 * @namespace Rodin::MMG::V5::Mesh
 * @brief Native MMG version 5 mesh handle operations and conversions.
 */
namespace Rodin::MMG::V5::Mesh
{
  /// MMG mesh version tag used for generated native mesh objects.
  constexpr int sMeshVersionFormatted = 2;

  /**
   * @brief Allocates and initializes an MMG mesh object.
   * @param[in] version MMG mesh version.
   * @param[in] dim Topological dimension of a planar, surface, or volume mesh.
   * @param[in] spaceDim Embedding dimension, defaulting to @p dim. Supported
   * topology/embedding pairs are @f$(2,2)@f$, @f$(2,3)@f$ and @f$(3,3)@f$.
   * @returns Newly allocated mesh, owned by the caller and released with
   * @ref destroy.
   * @throws Alert::Exception if the supplied dimensions are unsupported.
   */
  MMG5_pMesh create(size_t version, size_t dim, Optional<size_t> spaceDim = {});

  /**
   * @brief Copies source mesh to a destination mesh.
   * @param[in] src Source MMG mesh.
   * @param[out] dst Destination MMG mesh.
   *
   * Performs deep-copy allocation for dynamically owned MMG arrays.
   * Cached cell adjacencies use MMG's three entries per triangle or four
   * entries per tetrahedron, with five additional array entries.
   */
  void copy(const MMG5_pMesh src, MMG5_pMesh dst);

  /**
   * @brief Determines if a mesh is surface or not.
   * @param[in] mesh Mesh handle.
   * @returns `true` if the mesh is a surface manifold.
   */
  bool isSurface(const MMG5_pMesh mesh);

  /**
   * @brief Destroys the mesh object and frees the allocated memory.
   * @param[in] mesh Pointer to mesh.
   */
  void destroy(MMG5_pMesh mesh);

  /**
   * @brief Converts a Rodin local mesh to a native MMG mesh.
   * @param[in] src Source Rodin mesh.
   * @returns Newly allocated MMG mesh.
   *
   * The conversion preserves MMG-specific tags when `src` is an
   * @ref MMG::Mesh (corners, ridges, required entities).
   */
  MMG5_pMesh toV5(const Rodin::Geometry::LocalMesh& src);

  /**
   * @brief Converts a native MMG mesh to @ref MMG::Mesh.
   * @param[in] src Source MMG mesh handle.
   * @returns Converted Rodin mesh.
   */
  MMG::Mesh fromV5(const MMG5_pMesh src);

}

/**
 * @namespace Rodin::MMG::V5::Solution
 * @brief Native MMG version 5 solution handle operations and transfers.
 */
namespace Rodin::MMG::V5::Solution
{

  /**
   * @brief Constructs a solution and allocates space for it.
   * @param[in] mesh Owning MMG mesh.
   * @param[in] vdim Value dimension (`1` for scalar, `>1` for vector).
   * @returns Newly allocated MMG solution.
   */
  MMG5_pSol create(MMG5_pMesh mesh, size_t vdim);

  /**
   * @brief Deep-copies an MMG solution.
   * @param[in] src Source solution.
   * @param[out] dst Destination solution.
   */
  void copy(const MMG5_pSol src, MMG5_pSol dst);

  /**
   * @brief Copies values from MMG solution to Rodin MMG grid function.
   * @tparam Range Grid function range type (`Real` or `Math::SpatialVector<Real>`).
   * @param[in] src Source MMG solution.
   * @param[out] dst Destination grid function.
   *
   * The transfer preserves MMG node ordering and writes into the matrix
   * layout used by Rodin's @ref Variational::GridFunction data container.
   */
  template <class Range>
  void copy(const MMG5_pSol src, MMG::GridFunction<Range>& dst)
  {
    assert(src);
    if constexpr (std::is_same_v<Real, Range>)
    {
      assert(src->type == MMG5_Scalar);
      assert(dst.getFiniteElementSpace().getVectorDimension() == 1);
      auto& data = dst.getData();
      data.resize(src->np);
      // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
      // + np + 1.
      std::copy(src->m + 1, src->m + src->np + 1, data.data());
    }
    else if constexpr (std::is_same_v<Math::SpatialVector<Real>, Range>)
    {
      const size_t vdim = src->size;
      assert(src->type == MMG5_Vector);
      assert(vdim == dst.getFiniteElementSpace().getVectorDimension());
      auto& data = dst.getData();
      data.resize(vdim * src->np);
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
   * @brief Copies values from Rodin MMG grid function to MMG solution.
   * @tparam Range Grid function range type (`Real` or `Math::SpatialVector<Real>`).
   * @param[in] src Source grid function.
   * @param[out] dst Destination MMG solution.
   *
   * Allocates MMG storage buffers when needed and copies nodal values in the
   * memory layout expected by the selected MMG kernel.
   */
  template <class Range>
  void copy(const MMG::GridFunction<Range>& src, MMG5_pSol dst)
  {
    assert(dst);
    if constexpr (std::is_same_v<Real, Range>)
    {
      assert(dst->type == MMG5_Scalar);
      assert(src.getFiniteElementSpace().getVectorDimension() == 1);
      const auto& data = src.getData();
      assert(data.cols() == 1);
      assert(dst->size == 1);
      const size_t n = data.size();
      if (n)
      {
        dst->np = n;
        dst->npi = n;
        dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
        assert(dst->size == 1);
        assert(dst->np < dst->npmax);
        if (!dst->m)
        {
          // So 2 * (dst->np + 1) seems to work for most applications
          MMG5_SAFE_CALLOC(dst->m, 2 * (dst->npmax + 1), Real,
            Alert::Exception() << "Failed to allocate memory for MMG5_pSol->m."
                               << Alert::Raise);
        }
        std::copy(data.data(), data.data() + n, dst->m + 1);
      }
      else
      {
        dst->np = 0;
        dst->npi = 0;
        dst->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
      }
    }
    else if constexpr (std::is_same_v<Math::SpatialVector<Real>, Range>)
    {
      assert(dst->type == MMG5_Vector);
      [[maybe_unused]] const size_t vdim =
        src.getFiniteElementSpace().getVectorDimension();
      assert(dst->size >= 0);
      assert(vdim == static_cast<size_t>(dst->size));
      const auto& data = src.getData();
      assert(data.size() % vdim == 0);
      const size_t n = data.size() / vdim;
      dst->np = n;
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
            Alert::Exception() << "Failed to allocate memory for MMG5_pSol->m"
                               << Alert::Raise);
        }
        std::copy(data.data(), data.data() + data.size(), dst->m + dst->size);
      }
      else
      {
        dst->np = 0;
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
   * @brief Swaps the data between two instances of type MMG5_pSol.
   * @param[in, out] a First solution.
   * @param[in, out] b Second solution.
   */
  void swap(MMG5_pSol a, MMG5_pSol b);

  /**
   * @brief Destroys and frees the allocated memory for a MMG5_pSol object.
   * @param[in] sol MMG solution pointer.
   */
  void destroy(MMG5_pSol sol);

}
#endif
